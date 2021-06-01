import random
import numpy as np
import pandas as pd
from pathlib import Path

### Setup
###################################################
outdir = config['outdir']
manifest_table = pd.read_csv(config['manifest'], sep='\t')
file_ids = manifest_table['file_id'].tolist()

if (i := config['partition']) > 0:
    random.seed(42)
    random.shuffle(file_ids)
    n = len(file_ids)
    file_ids = file_ids[(i-1)*n//4 : i*n//4]

### Rules
###################################################
rule all:
    input:
        f'{outdir}/sites.smoove.square.vcf.gz'

rule GetVcf:
    """
    for each file_id, get the smoove/manta vcfs
    """
    output:
        smoove = f'{outdir}/vcf/{{file_id}}-smoove.genotyped.vcf',
        manta = f'{outdir}/vcf/{{file_id}}-manta.vcf'
    params:
        smoove_fname = lambda w: f'{w.file_id}-smoove.genotyped.vcf.gz',
        manta_fname = lambda w: f'{w.file_id}-manta.vcf.gz',
    shell:
        f"""
        aws s3 cp s3://layerlabcu/icgc/smoove.2/{{params.smoove_fname}} {outdir}/vcf/
        aws s3 cp s3://layerlabcu/icgc/manta.2/{{params.manta_fname}} {outdir}/vcf/
        
        # SURVIVOR does not work right with vcf.gz input
        bcftools view {outdir}/vcf/{{params.smoove_fname}} > {{output.smoove}}
        bcftools view {outdir}/vcf/{{params.manta_fname}} > {{output.manta}}
        """


rule SurvivorMergeVCFs:
    threads:
        1
    input:
        smoove = expand(f'{outdir}/vcf/{{file_id}}-smoove.genotyped.vcf', file_id=file_ids),
        manta = expand(f'{outdir}/vcf/{{file_id}}-manta.vcf', file_id=file_ids),
    output:
        vcf = f'{outdir}/survivor-merged.vcf'
    shell:
        f"""
        cat <(echo {{input.smoove_normal}}) <(echo {{input.manta_normal}}) |
            tr ' ' '\n' > {outdir}/vcf-list.txt
        
        max_dist_between_breakpoints=0.2 # fraction of SVLEN
        min_support=1
        take_type_into_account=1
        take_strand_into_account=0
        estimate_dist_from_SV_size=1
        min_size=50

        SURVIVOR merge {outdir}/vcf-list.txt \\
            $max_dist_between_breakpoints \\
            $min_support \\
            $take_type_into_account \\
            $take_strand_into_account \\
            $estimate_dist_from_SV_size \\
            $min_size \\
            {{output.vcf}}
        """

rule GetReference:
    output:
        fasta = temp(f'{outdir}/ref/hs37d5.fa'),
        fai = temp(f'{outdir}/ref/hs37d5.fa.fai'),
    shell:
        f"""
        aws s3 cp s3://layerlabcu/ref/genomes/hs37d5/hs37d5.fa {{output.fasta}}
        aws s3 cp s3://layerlabcu/ref/genomes/hs37d5/hs37d5.fa.fai {{output.fai}}
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
        bam = temp(f'{outdir}/bam/{{file_id}}.bam'),
        bai = temp(f'{outdir}/bam/{{file_id}}.bai'),
    resources:
        num_downloads = 1
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


### TODO
rule SmooveGenotype:
    input:
        fasta = rules.GetReference.output.fasta,
        fai = rules.GetReference.output.fai,
        bam = rules.GetBam.output.bam,
        bai = rules.GetBam.output.bai,
        vcf = rules.SurvivorMergeVCFs.output.vcf
    output:
        f'{outdir}/svtyper-vcf/{{file_id}}-smoove.genotyped.vcf.gz'
    conda:
        'envs/smoove.yaml'
    shell:
        f"""
        smoove genotype -p 1 -f {{input.fasta}} -v {{input.vcf}}\\
            -n {{wildcards.file_id}} \\
            -o {outdir}/svtyper-vcf/ \\
            {{input.bam}}
        """

rule SmoovePasteVCFs:
    conda:
        'envs/smoove.yaml'
    input:
        vcfs = expand(f'{outdir}/svtyper-vcf/{{file_id}}-smoove.genotyped.vcf.gz',
                      file_id=file_ids)
    output:
        f'{outdir}/sites.smoove.square.vcf.gz'
    shell:
        f"""
        smoove paste -o {outdir} --name sites {{input.vcfs}}
        """
