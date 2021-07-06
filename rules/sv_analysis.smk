import subprocess
import pandas as pd
import numpy as np
from pathlib import Path

### Setup
################################################################################
outdir = config['outdir']
donor_table = pd.read_csv(config['donor_table'], sep='\t')
genes_bed = config['genes_bed']
regions_18p16q = config['regions_18p16q']
somaticSV_bucket = config['somaticSV_bucket']
file_listing = subprocess.check_output([f'aws s3 ls {somaticSV_bucket}/'], shell=True)
tumour_file_ids = list(set(
    line.split()[3].split('.')[0]
    for line in file_listing.decode().rstrip().split('\n')))

### Rules
################################################################################
rule All:
    input:
        ## TODO
        expand(f'{outdir}/intersect_8p16q/{{fid}}.8p16q.vcf.gz',
               fid=tumour_file_ids)

rule GetSomaticVCFs:
    """
    get the manta somaticSV vcfs from s3
    """
    output:
        f'{outdir}/somaticSV/{{fid}}.somaticSV.vcf.gz'
    shell:
        f"""
        aws s3 cp {somaticSV_bucket}/{{wildcards.fid}}.somaticSV.vcf.gz {{output}}
        """

rule Intersect8p16q:
    """
    Intersect the somatic variants with the regions in 8p/16q
    """
    input:
        vcf = rules.GetSomaticVCFs.output
    output:
        f'{outdir}/intersect_8p16q/{{fid}}.8p16q.vcf.gz'
    conda:
        'envs/bedtools.yaml'
    shell:
        f"""
        bedtools intersect -header -a {{input.vcf}} -b {regions_18p16q} |
        bgzip -c > {{output}}
        """

    
