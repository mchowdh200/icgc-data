
outdir = config['outdir']
manifest_dir = config['manifest_dir']
with open(config['donor_list']) as f:
    donors = [x.rstrip() for x in f.readlines()]


rule all:
    ## TODO this should be the final SVTyper sites vcf.
    # input:
    #     f'{outdir}/sites-SVTyper.vcf.gz'
    # shell:
    #     """
    #     aws s3 cp {input} s3://layerlabcu/icgc/
    #     """
    input:
        directory(expand(f'{outdir}/smoove-vcf/{{donor}}', donor=donors))

### TODO
rule GetSmooveVCFs:
    output:
        directory(f'{outdir}/smoove-vcf/{{donor}}')
    shell:
        f"""
        aws s3 cp --recursive \
            s3://layerlabcu/icgc/smoove/{{wildcards.donor}}/ \
            {outdir}/smoove-vcf/{{wildcards.donor}}
        """


### TODO
rule RenameSmooveSamples:

### TODO
rule GetMantaVCFs:

### TODO
rule RenameMantaVCFs:

### TODO
rule SurvivorMergeVCFs:

### TODO
# see manta snakefile for similar rule
# and copy the resource management method
rule GetBams:

### TODO
# can use as many cores as we need.
rule SmooveGenotype:

### TODO
rule SmoovePasteVCFs:


