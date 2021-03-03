
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
# right now the grep command is returning a non-zero if there isn't a file
# present. Do an if statement to check if there is a missing vcf.
rule RenameSmooveSamples:
    output:
        normal = f'{outdir}/smoove-vcf/{{donor}}/{{donor}}.normal.vcf.gz',
        tumour = f'{outdir}/smoove-vcf/{{donor}}/{{donor}}.tumour.vcf.gz'
    conda:
        'envs/gatk.yaml'
    shell:
        # get stored smoove vcf
        # grep for the normal/tumour vcfs
        # run gatk to rename samples
        f"""
        aws s3 cp --recursive \
            s3://layerlabcu/icgc/smoove/{{wildcards.donor}}/ \
            {outdir}/smoove-vcf/{{wildcards.donor}}/
        
        set +e
        normal_in=$(find {outdir}/smoove-vcf/{{wildcards.donor}} -name '*.vcf.gz' |
                    grep -i normal)
        tumour_in=$(find {outdir}/smoove-vcf/{{wildcards.donor}} -name '*.vcf.gz' |
                    grep -i tumour)
        set -e

        if [[ -z $normal_in ]]; then
            gatk RenameSampleInVcf \
                --INPUT $normal_in \
                --OUTPUT {{output.normal}} \\
                --NEW_SAMPLE_NAME {{wildcards.donor}}-normal
        fi
        if [[ -z $tumour_in ]]; then
            gatk RenameSampleInVcf \
                --INPUT $tumour_in \
                --OUTPUT {{output.tumour}} \\
                --NEW_SAMPLE_NAME {{wildcards.donor}}-tumour
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


