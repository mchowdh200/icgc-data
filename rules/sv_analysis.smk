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
        f'{conf.outdir}/fusion_support.bedpe',
        # f'{conf.outdir}/fusion_list.bedpe',
        # expand(f'{conf.outdir}/filtered_fusions/{{fid}}.fusions.filtered.bedpe',
        #        fid=tumour_file_ids),
        # expand(f'{conf.outdir}/fusion_candidates/{{fid}}.fusion_candidates.bedpe',
        #        fid=tumour_file_ids),
        # f'{conf.outdir}/wordcloud.png',
        # f'{conf.outdir}/gene_cooccurrence.npz',
        # f'{conf.outdir}/gene_cooccurrence.graphml'
        # expand(f'{conf.outdir}/intersect_cytoband/{{fid}}.cytoband.bed',
        #        fid=tumour_file_ids),
        # expand(f'{conf.outdir}/intersect_windows/{{fid}}.bins.bed',
        #        fid=tumour_file_ids),
        # expand(f'{conf.outdir}/del_fusions/{{fid}}.del_fusions.bedpe',
        #        fid=tumour_file_ids),
        # f'{conf.outdir}/combined.bedpe',
        # expand(f'{conf.outdir}/intersect_8p16q/{{fid}}.8p16q.bedpe',
        #        fid=tumour_file_ids)
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

rule GetGeneOccurrence:
    """
    get # of variant occurrences for each gene
    """
    input:
        genes=conf.genes_bed,
        vcf = rules.GetSomaticVCFs.output
    output:
        f'{conf.outdir}/gene_occurrence/{{fid}}.gene_occurrence.bed'
    shell:
        'bedtools intersect -a {input.genes} -b {input.vcf} -c > {output}'

rule GeneTopicModel:
    input:
        expand(f'{conf.outdir}/intersect_cytoband/{{fid}}.cytoband.bed',
               fid=tumour_file_ids)
        # expand(f'{conf.outdir}/gene_occurrence/{{fid}}.gene_occurrence.bed',
        #        fid=tumour_file_ids)
    output:
        f'{conf.outdir}/wordcloud.png'
    params:
        # columns to get bed info from (0-based)
        feature_column = 3,
        count_column = 5
    shell:
        """
        python scripts/hdp_model.py {params.feature_column} \\
                                    {params.count_column} \\
                                    {output} \\
                                    {input}
        """

rule GetGeneCooccurrenceMatrix:
    """
    Compute (sparse) cooccurrence matrix and write to binary format on disk
    """
    input:
        expand(f'{conf.outdir}/gene_occurrence/{{fid}}.gene_occurrence.bed',
               fid=tumour_file_ids)
    output:
        matrix = f'{conf.outdir}/gene_cooccurrence.npz',
        gene_features = f'{conf.outdir}/gene_features.txt'
    params:
        # columns to get bed info from (0-based)
        feature_column = 3,
        count_column = 5
    threads:
        workflow.cores
    shell:
        """
        NUMBA_NUM_THREADS={{threads}}
        python scripts/cooccurrence.py {params.feature_column} \\
                                       {params.count_column} \\
                                       {output.matrix} \\
                                       {output.gene_features} \\
                                       {input}
        """
rule CoocurrenceMatrix2Graph:
    """
    Convert cooccurrence matrix to graphml format for visualization
    """
    input:
        gene_matrix = rules.GetGeneCooccurrenceMatrix.output.matrix,
        gene_features = rules.GetGeneCooccurrenceMatrix.output.gene_features
    output:
        f'{conf.outdir}/gene_cooccurrence.graphml'
    shell:
        """
        python scripts/co_occ2graph.py {input.gene_matrix} {input.gene_features} {output}
        """
        



rule GenomeMakeWindows:
    input:
        conf.contig_lengths
    output:
        f'{conf.outdir}/hs37d5.windows.bed'
    conda:
        'envs/bedtools.yaml'
    shell:
        'bedtools makewindows -g {input} -w 1000000 > {output}'

rule IntersectWindows:
    """
    Take the vcfs and intersect with the genome windows to get
    # of overlaps within each window. A = windows, B = vcf
    """
    input:
        windows = rules.GenomeMakeWindows.output,
        vcf = rules.GetSomaticVCFs.output
    output:
        f'{conf.outdir}/intersect_windows/{{fid}}.bins.bed'
    shell:
        'bedtools intersect -a {input.windows} -b {input.vcf} -c > {output}'


rule IntersectCytobandWindows:
    """
    take the vcfs and intersect with the cytoband regions to get #
    of overlaps within each window.  A = windows, B = vcf
    """
    input:
        windows = conf.cytoband_bed,
        vcf = rules.GetSomaticVCFs.output
    output:
        f'{conf.outdir}/intersect_cytoband/{{fid}}.cytoband.bed'
    shell:
        f"""
        mkdir -p {conf.outdir}/intersect_cytoband
        bedtools intersect -a {{input.windows}} -b {{input.vcf}} -c > {{output}}
        """


rule PlotCooccurrenceMatrix:
    input:
        expand(f'{conf.outdir}/intersect_cytoband/{{fid}}.cytoband.bed',
               fid=tumour_file_ids)
    output:
        f'{conf.outdir}/cooccurrence_matrix.png'
    shell:
        'python3 scripts/plot_cooccurrence.py {output} {input}'



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
    Breakpoints = cols 0-5 inclusive
    Breakend orientation = cols 8, 9
    Gene = col 26.
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


rule GetFusionCandidates:
    """
    Collapse del regions into a comma separated list of genes they intersect with.
    ideally, for bedpe del breakpoints, this should intersect with at most two genes per region.
    We will also ensure that the genes must be on the same strand to be considered
    a deletion based fusion.
    """
    input:
        rules.IntersectGenes.output
    output:
        f'{conf.outdir}/fusion_candidates/{{fid}}.fusion_candidates.bedpe'
    conda:
        'envs/bedtools.yaml'
    shell:
        f"""
        mkdir -p {conf.outdir}/del_fusions
        bash scripts/fusion_candidates.sh {{input}} {{output}}
        """

rule FilterFusionCandidates:
    """
    Using genes strand, breakpoint orientation, and SVTYPE, filter
    the fusion candidates to ensure that only valid fusions are present.
    """
    input:
        candidates = rules.GetFusionCandidates.output,
        genes = conf.genes_bed
    output:
        f'{conf.outdir}/filtered_fusions/{{fid}}.fusions.filtered.bedpe'
    shell:
        """
        python scripts/filter_fusions.py {input.candidates} {input.genes} > {output}
        """

rule GetFusionsList:
    """
    From our set of fusions, compile a list of unique 
    * fusion genes,
    * their loci 
    * and the SVTYPE that gave rise to them

    
    * genes bed input columns
    0,1,2 -- region
    3 -- gene
    * fusion bedpe relevant columns
    col 6 -- SVTYPE
    col 9 -- GENE_A,GENE_B

    * output format
    CHR_A START_A END_A  CHR_B START_B END_B GENE_A,GENE_B SVTYPE
    """
    input:
        genes_bed = conf.genes_bed,
        fusion_bedpe = expand(
            f'{conf.outdir}/filtered_fusions/{{fid}}.fusions.filtered.bedpe',
            fid=tumour_file_ids)
    output:
        f'{conf.outdir}/fusion_list.bedpe'
    shell:
        """
        python scripts/get_fusions.py {input} |
        bedtools sort | uniq > {output}"""

rule StixGetFusionSupport:
    """
    get total evidence in support for all fusions
    """
    input:
        rules.GetFusionsList.output
    output:
        f'{conf.outdir}/fusion_support.bedpe'
    threads:
        workflow.cores
    shell:
        'bash scripts/stix_fusion_query.sh {input} {output} {threads}'




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
        bedtools intersect -header -a {{input}} -b {conf.regions_8p16q} > {{output}}
        # bedtools intersect -header -a {{input}} -b {conf.regions_8p16q} | bgzip -c > {{output}}
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

        
    
