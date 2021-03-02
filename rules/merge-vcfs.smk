
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
        normal = expand(f'{outdir}/smoove-vcf/{{donor}}/{{donor}}.normal.vcf.gz',
                        donor=donors),
        tumour = expand(f'{outdir}/smoove-vcf/{{donor}}/{{donor}}.tumour.vcf.gz',
                        donor=donors)

rule GetSmooveVCFs:
    output:
        directory(f'{outdir}/smoove-vcf/{{donor}}')
    shell:
        f"""
        aws s3 cp --recursive \
            s3://layerlabcu/icgc/smoove/{{wildcards.donor}}/ \
            {outdir}/smoove-vcf/{{wildcards.donor}}
        """

rule RenameSmooveSamples:
    input:
        dir = f'{outdir}/smoove-vcf/{{donor}}'
    output:
        # normal = f'{outdir}/smoove-vcf/{{donor}}/{{donor}}.normal.vcf.gz',
        # tumour = f'{outdir}/smoove-vcf/{{donor}}/{{donor}}.tumour.vcf.gz'
        normal = f'{input.dir}/{{donor}}.normal.vcf.gz',
        tumour = f'{input.dir}/{{donor}}/{{donor}}.tumour.vcf.gz'
    conda:
        'envs/gatk.yaml'
    shell:
        # grep for the normal/tumour vcfs
        # run gatk to rename samples
        f"""
        normal_in=$(find {outdir}/smoove-vcf/{{wildcards.donor}} -name '*.vcf.gz' |
                    grep -i normal)
        tumour_in=$(find {outdir}/smoove-vcf/{{wildcards.donor}} -name '*.vcf.gz' |
                    grep -i tumour)
        echo $normal_in
        echo $tumour_in
        
        gatk RenameSampleInVcf \
            --INPUT $normal_in \
            --OUTPUT {{output.normal}} \\
            --NEW_SAMPLE_NAME {{wildcards.donor}}-normal
        gatk RenameSampleInVcf \
            --INPUT $tumour_in \
            --OUTPUT {{output.tumour}} \\
            --NEW_SAMPLE_NAME {{wildcards.donor}}-tumour
        """
    
### TODO
# not for this workflow of course, but looks like there is a bioconda
# version of manta.  Look into using that for the manta workflow
# instead of installing with the setup script.
rule GetMantaVCFs:

### TODO
rule RenameMantaVCFs:

### TODO
# try using the bioconda version of survivor.
# it uses the same version that we currently have installed
# then we can remove that from installation
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


