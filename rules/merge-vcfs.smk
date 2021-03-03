
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

### TODO handle if there is a missing vcf
# this rule also checks if there are missing vcfs, and produces a dummy output
# using touch.  This will have to be handled down the line before SURVIVOR
# is run.
rule RenameSmooveSamples:
    output:
        normal = temp(f'{outdir}/smoove-vcf/{{donor}}/{{donor}}.normal.vcf.gz'),
        tumour = temp(f'{outdir}/smoove-vcf/{{donor}}/{{donor}}.tumour.vcf.gz')
    conda:
        '../envs/samtools.yaml'
    shell:
        # get stored smoove vcf
        # grep for the normal/tumour vcfs
        # run gatk to rename samples
        f"""
        aws s3 cp --recursive \
            s3://layerlabcu/icgc/smoove/{{wildcards.donor}}/ \
            {outdir}/smoove-vcf/{{wildcards.donor}}/
        
        set +euo pipefail
        normal_in=$(find {outdir}/smoove-vcf/{{wildcards.donor}} -name '*.vcf.gz' |
                    grep -i normal)
        tumour_in=$(find {outdir}/smoove-vcf/{{wildcards.donor}} -name '*.vcf.gz' |
                    grep -i tumour | head -1)

        if [[ ! -z $normal_in ]]; then
            echo {{wildcards.donor}}-normal \
                > {outdir}/smoove-vcf/{{wildcards.donor}}/normal-sample.txt
            bcftools reheader \
                -s {outdir}/smoove-vcf/{{wildcards.donor}}/normal-samples.txt \
                -o {{output.normal}} \
                $normal_in
        else
            touch {{output.normal}}
        fi
        if [[ ! -z $tumour_in ]]; then
            echo {{wildcards.donor}}-tumour \
                > {outdir}/smoove-vcf/{{wildcards.donor}}/tumour-sample.txt
            bcftools reheader \
                -s {outdir}/smoove-vcf/{{wildcards.donor}}/tumour-sample.txt \
                -o {{output.normal}} \
                $tumour_in
        else
            touch {{output.tumour}}
        fi
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


