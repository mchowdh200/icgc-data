
outdir = config['outdir']
manifest_dir = config['manifest_dir']
with open(config['donor_list']) as f:
    donors = [x.rstrip() for x in f.readlines()]


rule all:
    ## this should be the final SVTyper sites vcf.
    # Take care when mixing global variables and Snakemake
    # wildcards in a formatted string. In general, surround
    # wildcards with double curly braces.
    # EXAMPLE: BOOK_FILE = f'{INPUT_DIR}{{book}}.txt'
    input:
        f'{outdir}/sites-SVTyper.vcf.gz'
    shell:
        """
        aws s3 cp {input} s3://layerlabcu/icgc/
        """

### TODO
rule GetSmooveVCFs:

### TODO
rule RenameSmooveSamples:

### TODO
rule GetMantaVCFs:

### TODO
rule RenameMantaVCFs:

### TODO
rule SurvivorMergeVCFs:

### TODO
rule GetBams:

### TODO
rule SmooveGenotype:

### TODO
rule SmoovePasteVCFs:


