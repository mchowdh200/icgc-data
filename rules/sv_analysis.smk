import subprocess
import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict
from types import SimpleNamespace

### Setup
################################################################################
# don't want to use dict notation to get stuff from config
conf = SimpleNamespace(**config)
donor_table = pd.read_csv(conf.donor_table, sep='\t')
# outdir = config['outdir']
# somaticSV_bucket = config['somaticSV_bucket']
# genes_bed = config['genes_bed']
# regions8p16q = config['regions8p16q']
# cytoband_bed = config['cytoband_bed']

file_listing = subprocess.check_output([f'aws s3 ls {conf.somaticSV_bucket}/'], shell=True)
tumour_file_ids = list(set(
    line.split()[3].split('.')[0]
    for line in file_listing.decode().rstrip().split('\n')))

# map the file id to other useful information
fid2donor = {
    fid : donor_table[donor_table['File ID'] == fid]['ICGC Donor'].values[0]
    for fid in tumour_file_ids
}
fid2sample = {
    (x:=line.rstrip().split())[0] : x[1]
    for line in open(config['sample_names'])
}

### Rules
################################################################################
rule All:
    input:
        expand(f'{conf.outdir}/del_fusions/{{fid}}.del_fusions.bedpe',
               fid=tumour_file_ids),
        f'{conf.outdir}/combined.bedpe'
        # expand(f'{conf.outdir}/intersect_8p16q/{{fid}}.8p16q.bedpe',
               # fid=tumour_file_ids)
        # expand(f'{conf.outdir}/intersect_8p16q/{{fid}}.8p16q.vcf.gz',
        #        fid=tumour_file_ids)

rule GetSomaticVCFs:
    """
    get the manta somaticSV vcfs from s3
    """
    output:
        f'{conf.outdir}/somaticSV/{{fid}}.somaticSV.vcf.gz'
    shell:
        f"""
        aws s3 cp {conf.somaticSV_bucket}/{{wildcards.fid}}.somaticSV.vcf.gz {{output}}
        """

rule VCF2Bedpe:
    input:
        rules.GetSomaticVCFs.output
    output:
        f'{conf.outdir}/bedpe/{{fid}}.bedpe'
    conda:
        'envs/bedtools.yaml'
    shell:
        f"""
        mkdir -p {conf.outdir}/bedpe
        bcftools view {{input}} | svtools vcftobedpe > {{output}}
        """

rule IntersectGenes:
    """
    Intersect the bedpe of variants with a genes bed.
    Gene will be on col 27.
    Strand of gene will be on col 28
    """
    input:
        rules.VCF2Bedpe.output
    output:
        f'{conf.outdir}/intersect_genes/{{fid}}.genes.bedpe'
    conda:
        'envs/bedtools.yaml'
    shell:
        f"""
        mkdir -p {conf.outdir}/intersect_genes
        pairToBed -a {{input}} -b {conf.genes_bed} > {{output}}
        """

rule IntersectCytoband:
    """
    Intersect variants with cytoband notated regions
    """
    input:
        rules.VCF2Bedpe.output
    output:
        f'{conf.outdir}/intersect_cytoband/{{fid}}.cytoband_regions.bed'
    conda:
        'envs/bedtools.yaml'
    shell:
        f"""
        mkdir -p {conf.outdir}/intersect_cytoband
        pairToBed -a {{input}} -b {{conf.cytoband_bed}} > {{output}}
        """

rule GetDelFusions:
    ## TODO check the deletion strand as well
    # as a final filtering step
    """
    Collapse del regions into a comma separated list of genes they intersect with.
    ideally, for bedpe del breakpoints, this should intersect with at most two genes per region.
    We will also ensure that the genes must be on the same strand to be considered
    a deletion based fusion.
    """
    input:
        rules.IntersectGenes.output
    output:
        f'{conf.outdir}/del_fusions/{{fid}}.del_fusions.bedpe'
    conda:
        'envs/bedtools.yaml'
    shell:
        f"""
        mkdir -p {conf.outdir}/del_fusions
        bash scripts/del_fusions.sh {{input}} {{output}}
        """

rule Intersect8p16q:
    """
    Intersect the somatic variants with the regions in 8p/16q.
    """
    input:
        rules.VCF2Bedpe.output
        # rules.GetSomaticVCFs.output
    output:
        f'{conf.outdir}/intersect_8p16q/{{fid}}.8p16q.bedpe'
        # f'{conf.outdir}/intersect_8p16q/{{fid}}.8p16q.vcf.gz'
    conda:
        'envs/bedtools.yaml'
    shell:
        f"""
        bedtools intersect -header -a {{input}} -b {conf.regions8p16q} > {{output}}
        # bedtools intersect -header -a {{input}} -b {conf.regions8p16q} | bgzip -c > {{output}}
        """

rule CombineIntersections:
    """
    cat all the intersected bedpe files.
    In the process, get (or add) only the following fields:
    - genomic regions: cols [0:6)
    - strand A and B: cols 9, 10
    - SVTYPE: col 11
    - Donor ID: add
    - File ID: add
    - Specimen Type: add
    """
    input:
        expand(f'{conf.outdir}/intersect_8p16q/{{fid}}.8p16q.bedpe',
               fid=tumour_file_ids)
    output:
        f'{conf.outdir}/combined.bedpe'
    run:
        header = '#chromA\tstartA\tstartB\tchromB\tstartB\tendB\tstrandA\tstrandB\tSVType\tdonorID\tfileID\tspecimenType\n'
        with open(output[0], 'w') as o:
            o.write(header)
            for bedpe_fname in input:
                with open(bedpe_fname) as f:
                    for line in f:
                        if line[0] == '#':
                            continue
                        A = line.rstrip().split()
                        genomic_region = A[:6]
                        strands = A[8:10]
                        svtype = A[10]
                        print(bedpe_fname)
                        fid = Path(bedpe_fname).name.split('.')[0]
                        donor_id = donor_table[
                            donor_table['File ID'] == fid]['ICGC Donor'].values[0]
                        specimen_type = donor_table[
                            donor_table['File ID'] == fid]['Specimen Type'].values[0]
                        o.write('\t'.join(genomic_region + strands +
                                          [svtype, donor_id, fid, specimen_type]))
                        o.write('\n')

        
    
