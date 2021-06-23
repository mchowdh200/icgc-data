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
        gt0_icgc = expand(f'{outdir}/intersections/{{fid}}.gt0.icgc.del.bed',
                          fid=tumour_file_ids),
        gt1_icgc = expand(f'{outdir}/intersections/{{fid}}.gt1.icgc.del.bed',
                          fid=tumour_file_ids),
        gnomad_icgc = expand(f'{outdir}/intersections/{{fid}}.gnomad-sub.icgc.del.bed',
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
        f"""
        mkdir -p {outdir}/icgc_bed
        bash scripts/icgc_bedpe2bed.sh {{input}} {{output}} DEL
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
        gt0_bed = f'{outdir}/thresholded/{{fid}}.gt0.stix.bed',
        gt1_bed = f'{outdir}/thresholded/{{fid}}.gt1.stix.bed'
    shell:
        f"""
        mkdir -p {outdir}/thresholded
        bash scripts/threshold_called_regions.sh {{input}} {{output.gt0_bed}} 0
        bash scripts/threshold_called_regions.sh {{input}} {{output.gt1_bed}} 1
        """

### TODO get gnomad subtraction regions for comparison
rule SubtractGnomadRegions:
    """
    remove regions in the sv callset that overlap with gnomad
    """
    input:
        stix_bed = rules.StixQuerySingleSample.output,
        gnomad_bed = "/home/much8161/data/stix/1kg/gnomad.DEL.bed" # TODO don't hardcode
    output:
        f'{outdir}/gnomad_subtracted/{{fid}}.gnomad_subtracted.del.bed'
    conda:
        'envs/bedtools.yaml'
    shell:
        f"""
        mkdir -p {outdir}/gnomad_subtracted
        bedtools intersect -v -r -f 0.9 -a {{input.stix_bed}} -b {{input.gnomad_bed}} |
            grep -v hs | grep -v GL | grep -v X | grep -v Y > {{output}}
        """
    
rule IntersectICGC:
    """
    Intersect the thresholded SV calls with ICGC truth regions
    TODO Intersect the gnomad subtracted SV calls as well
    """
    input:
        gt0_bed = rules.ThresholdCalledRegions.output.gt0_bed,
        gt1_bed = rules.ThresholdCalledRegions.output.gt1_bed,
        icgc_bed = rules.GetIcgcSampleDels.output,
        gnomad_sub_bed = rules.SubtractGnomadRegions.output
    output:
        gt0_icgc = f'{outdir}/intersections/{{fid}}.gt0.icgc.del.bed',
        gt1_icgc = f'{outdir}/intersections/{{fid}}.gt1.icgc.del.bed',
        gnomad_icgc = f'{outdir}/intersections/{{fid}}.gnomad-sub.icgc.del.bed'
    conda:
        'envs/bedtools.yaml'
    shell:
        f"""
        mkdir -p {outdir}/intersections
        scripts/intersect_icgc.sh {{input.gt0_bed}} {{input.icgc_bed}} {{output.gt0_icgc}}
        scripts/intersect_icgc.sh {{input.gt1_bed}} {{input.icgc_bed}} {{output.gt1_icgc}}
        scripts/intersect_icgc.sh {{input.gnomad_sub_bed}} {{input.icgc_bed}} {{output.gnomad_icgc}}
        """
