import subprocess
import numpy as np
import pandas as pd
from pathlib import Path


### Setup
###############################################################################
### TODO add config to replace hardcoded paths
outdir = "/home/much8161/data/stix/1kg/manta-tumour-normal-analysis"

# TODO put the base location in a config
# or create a rule to download the base data
icgc_bedpe_dir = '/home/much8161/data/stix/1kg/pcawg/icgc/open'

# table containing donor ids, file ids, specimen types, etc.
donor_table = pd.read_csv(
    '/home/much8161/Repositories/icgc-data/data_listings/donor_table.tsv',
    sep='\t')

# file id => bam sample name
fileid2sample = {
    line.split()[0] : line.split()[1]
    for line in open('/home/much8161/Repositories/icgc-data/data_listings/sample_names.txt').readlines()}   

# list the files in manta-tumour-normal somaticSV vcfs
file_listing = subprocess.check_output(                                                                                                                                   
    ['aws s3 ls s3://layerlabcu/icgc/manta-tumour-normal/'], shell=True,)
tumour_file_ids = list(set(
    line.split()[3].split('.')[0]                                                                                                                      
    for line in file_listing.decode().rstrip().split('\n')))

### Rules
###############################################################################
rule All:
    input:
        expand(f'{outdir}/bed/{{fid}}.stix.single_sample.bed',
               fid=tumour_file_ids)

rule GetSingleSampleVCFs:
    """
    get the single sample mode manta vcfs from s3
    """
    output:
        f'{outdir}/single_sample_vcf/{{fid}}-manta.vcf.gz'
    shell:
        """
        aws s3 cp s3://layerlabcu/icgc/manta.2/{wildcards.fid}-manta.vcf.gz {output}
        """

rule GetSomaticVCFs:
    """
    get the somaticSV vcfs from s3
    """
    output:
        f'{outdir}/somaticSV_vcf/{{fid}}.somaticSV.vcf.gz'
    shell:
        """
        aws s3 cp s3://layerlabcu/icgc/manta-tumour-normal/{wildcards.fid}.somaticSV.vcf.gz {output}
        """

rule GetSingleSampleDels:
    input:
        rules.GetSingleSampleVCFs.output
    output:
        f'{outdir}/single_sample_vcf/{{fid}}.del.vcf.gz'
    conda:
        'envs/pysam.yaml'
    shell:
        'python3 scripts/get_dels.py {input} | bgzip -c > {output}'

rule StixQuerySingleSample:
    threads:
        workflow.cores
    input:
        rules.GetSingleSampleDels.output
    output:
        f'{outdir}/bed/{{fid}}.stix.single_sample.bed'
    conda:
        'envs/pysam.yaml'
    shell:
        """
        bash scripts/stix_cmd.sh {input} {output} {threads}
        """

rule GetIcgcSampleDels:
    """
    get del regions for a given sample from the icgc bedpe
    svclass = column 11
    """
    input:
        lambda w: f'{icgc_bedpe_dir}/{fileid2sample[w.fid]}.pcawg_consensus_1.6.161116.somatic.sv.bedpe.gz'
    output:
        f'{outdir}/icgc_bed/{{fid}}.del.bed'
    shell:
        """
        bash scripts/icgc_bedpe2bed.sh {input} {output} DEL
        """

rule ThresholdCalledRegions:
    """
    filter out called del regions with:
    - gt 0 stix hits
    - gt 1 stix hits
    output to separate bed files
    """
    input:
        rules.StixQuerySingleSample.output
    output:
        gt0_bed = f'{outdir}/bed/{{fid}}.gt0.stix.bed',
        gt1_bed = f'{outdir}/bed/{{fid}}.gt1.stix.bed'
    shell:
        """
        bash scripts/threshold_called_regions.sh {input} {output.gt0_bed} 0
        bash scripts/threshold_called_regions.sh {input} {output.gt1_bed} 1
        """
