"""
Tumour normal SV calling pipeline using manta
"""

import random
from pathlib import Path
import numpy as np
import pandas as pd

### Helper functions
def get_normal_file_id(tumour_file_id, donor_table):
    """
    Given the file id of a tumour bam, get the corresponding normal bam ID
    """
    donor = donor_table[donor_table['File ID'] == tumour_file_id]['ICGC Donor'].values[0]
    donor_df = donor_table[(donor_table['ICGC Donor'] == donor)]
    normal_file_id = donor_df[
        donor_df['Specimen Type'] \
        .str.contains('normal', case=False)]['File ID'].values[0]
    return normal_file_id

def get_tumour_normal_pairs(manifest_table, donor_table):
    """
    returns a dictionary of tumour/normal bam file ids.
    """
    # exclude donors in the donor_table that arent in the manifest
    donors = set(manifest_table['donor_id/donor_count'])
    donor_table = donor_table[donor_table['ICGC Donor'].isin(donors)]

    # get all files with tumour specimen type
    is_tumour = donor_table.apply(
        lambda x: 'tumour' in x['Specimen Type'].lower(), axis=1)
    tumour_file_ids = donor_table[is_tumour]['File ID']

    # for each tumour sample, pair it with the corresponding normal sample
    tumour_normal_pairs = {
        tumour_file_id: get_normal_file_id(tumour_file_id, donor_table)
        for tumour_file_id in tumour_file_ids
    }

    # do some quick testing
    for tid, nid in tumour_normal_pairs.items():
        assert(len(donor_table[donor_table['File ID'] == tid]) == 1)
        assert(len(donor_table[donor_table['File ID'] == nid]) == 1)
        assert('tumour' in donor_table[
            donor_table['File ID'] == tid]['Specimen Type'].values[0].lower())
        assert('normal' in donor_table[
            donor_table['File ID'] == nid]['Specimen Type'].values[0].lower())

    return tumour_normal_pairs


### Setup
###################################################
outdir = config['outdir']
manifest_table = pd.read_csv(config['manifest'], sep='\t')
file_ids = manifest_table['file_id'].tolist()
tumour_normal_pairs = get_tumour_normal_pairs(manifest_table, donor_table)
### TODO pick a single tumour/normal pair to test with 
file_ids = ['FI10014', 'FI10013']
tumour_normal_pairs = {'FI10014', tumour_normal_pairs['FI10014']}

### Rules
################################################################################
rule all:
    input:
        expand(f'{outdir}/{{tumour_file_id}}/{{tumour_file_id}}.diploidSV.vcf.gz',
               tumour_file_id=list(tumour_normal_pairs.keys())),
        expand(f'{outdir}/{{tumour_file_id}}/{{tumour_file_id}}.somaticSV.vcf.gz',
               tumour_file_id=list(tumour_normal_pairs.keys())),
        expand(f'{outdir}/{{tumour_file_id}}/{{tumour_file_id}}.candidateSV.vcf.gz',
               tumour_file_id=list(tumour_normal_pairs.keys())),
        expand(f'{outdir}/{{tumour_file_id}}/{{tumour_file_id}}.candidateSmallIndels.vcf.gz',
               tumour_file_id=list(tumour_normal_pairs.keys())),

rule GetReference:
    output:
        fasta = temp(f'{outdir}/ref/hs37d5.fa'),
        fai = temp(f'{outdir}/ref/hs37d5.fa.fai')
    shell:
        """
        aws s3 cp s3://layerlabcu/ref/genomes/hs37d5/hs37d5.fa {output.fasta}
        aws s3 cp s3://layerlabcu/ref/genomes/hs37d5/hs37d5.fa.fai {output.fai}
        """

rule GetManifest:
    output:
        f'{outdir}/manifests/{{file_id}}-manifest.tsv'
    run:
        Path(f'{outdir}/manifests').mkdir(exist_ok=True)
        m = manifest_table[manifest_table['file_id'] == wildcards.file_id]
        m.to_csv(output[0], index=False, sep='\t')


rule GetBam:
    input:
        manifest = f'{outdir}/manifests/{{file_id}}-manifest.tsv'
    output:
        bam = temp(f"{outdir}/bam/{{file_id}}.bam"),
        bai = temp(f"{outdir}/bam/{{file_id}}.bai")
    run:
        Path(f'{outdir}/bam').mkdir(exist_ok=True)
        shell(f"""score-client --quiet download \\
        --validate false \\
        --output-dir {outdir}/bam \\
        --manifest {input.manifest}""")

        # remove bam and rename index with file_id
        bam = f'{outdir}/bam/{manifest_table[manifest_table.file_id == wildcards.file_id].file_name.values[0]}'
        bai = f'{bam}.bai'
        Path(bam).rename(output.bam)
        Path(bai).rename(output.bai)


rule DownloadsComplete:
    input:
        bams = expand(f'{outdir}/bam/{{file_id}}.bam'),
        bais = expand(f'{outdir}/bam/{{file_id}}.bai')
    output:
        f'{outdir}/downloads_complete.txt'
    shell:
        'touch {output}'


# change inputs/outputs
rule RunManta:
    input:
        downloads_complete = rules.DownloadsComplete.output,
        tumour_bam = f'{outdir}/bam/{{tumour_file_id}}.bam',
        tumour_bai = f'{outdir}/bam/{{tumour_file_id}}.bai',
        normal_bam = lambda w: f'{outdir}/bam/{tumour_normal_pairs[w.tumour_file_id]}.bam',
        normal_bai = lambda w: f'{outdir}/bam/{tumour_normal_pairs[w.tumour_file_id]}.bai',
        fasta = rules.GetReference.output.fasta,
        fai = rules.GetReference.output.fai,
    output:
        germline = f'{outdir}/{{tumour_file_id}}/{{tumour_file_id}}.diploidSV.vcf.gz',
        somatic = f'{outdir}/{{tumour_file_id}}/{{tumour_file_id}}.somaticSV.vcf.gz',
        candidate = f'{outdir}/{{tumour_file_id}}/{{tumour_file_id}}.candidateSV.vcf.gz',
        candidate_smallindels = f'{outdir}/{{tumour_file_id}}/{{tumour_file_id}}.candidateSmallIndels.vcf.gz',
    threads:
        workflow.cores
    shell:
        f"""
        mkdir -p {outdir}/{{wildcards.tumour_file_id}}
        rm -f {outdir}/{wildcards.tumour_file_id}}/runWorkflow.py.pickle
        /mnt/local/bin/configManta.py \
            --normalBam {input.normal_bam} \
            --tumorBam {input.tumour_bam} \
            --referenceFasta {input.fasta} \
            --runDir {outdir}/{{wildcards.tumour_file_id}}
        {outdir}/{{wildcards.tumour_file_id}}/runWorkflow.py -j {{threads}}
        
        cp {outdir}/{{wildcards.tumour_file_id}}/results/variants/diploidSV.vcf.gz {{output.germline}}
        cp {outdir}/{{wildcards.tumour_file_id}}/results/variants/somaticSV.vcf.gz {{output.somatic}}
        cp {outdir}/{{wildcards.tumour_file_id}}/results/variants/candidateSV.vcf.gz {{output.candidate}}
        cp {outdir}/{{wildcards.tumour_file_id}}/results/variants/candidateSmallIndels.vcf.gz {{output.candidate_smallindels}}
        find {outdir}/{{wildcards.tumour_file_id}}/ -type f -not -name '{{wildcards.tumour_file_id}}.*.vcf.gz'
        """
