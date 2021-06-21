import subprocess
import numpy as np
import pandas as pd
from pathlib import Path


### Setup
###############################################################################
### TODO add config later
outdir = "/data/stix/1kg/manta-tumour-normal-analysis"

# list the files in manta-tumour-normal somaticSV vcfs
file_listing = subprocess.check_output(                                                                                                                                   
    ['aws s3 ls s3://layerlabcu/icgc/manta-tumour-normal/'], shell=True,)
tumour_file_ids = list(set(
    line.split()[3].split('.')[0]                                                                                                                      
    for line in file_listing.decode().rstrip().split('\n')))

### Rules
###############################################################################
### TODO
rule All:
    input:
        single_sample = expand(
            f'{outdir}/single_sample_vcf/{{fid}}-manta.vcf.gz',
            fid=tumour_file_ids),
        somaticSV = expand(
            f'{outdir}/somaticSV_vcf/{{fid}}.somaticSV.vcf.gz',
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

### TODO
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
    
