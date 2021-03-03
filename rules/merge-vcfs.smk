
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
        f'{outdir}/survivor-merged.vcf.gz'

### TODO handle if there is a missing vcf
# this rule also checks if there are missing vcfs, and produces a dummy output
# using touch.  This will have to be handled down the line before SURVIVOR
# is run (use the "file" command).
rule RenameSmooveSamples:
    output:
        normal = f'{outdir}/smoove-vcf/{{donor}}/{{donor}}.normal.vcf.gz',
        tumour = f'{outdir}/smoove-vcf/{{donor}}/{{donor}}.tumour.vcf.gz'
    shell:
        f"""
        aws s3 cp --recursive \\
            s3://layerlabcu/icgc/smoove/{{wildcards.donor}}/ \\
            {outdir}/smoove-vcf/{{wildcards.donor}}/
        
        set +e # grep returns non zero if search is empty
        normal_in=$(find {outdir}/smoove-vcf/{{wildcards.donor}} -name '*.vcf.gz' |
                    grep -i normal)
        tumour_in=$(find {outdir}/smoove-vcf/{{wildcards.donor}} -name '*.vcf.gz' |
                    grep -i tumour | head -1)
        set e

        if [[ ! -z $normal_in ]]; then
            bcftools reheader \\
                -s <(echo {{wildcards.donor}}-normal) \\
                -o {{output.normal}} \\
                $normal_in
        else
            # dummy output if search was empty
            touch {{output.normal}} 
        fi
        if [[ ! -z $tumour_in ]]; then
            bcftools reheader \\
                -s <(echo {{wildcards.donor}}-tumour) \\
                -o {{output.tumour}} \\
                $tumour_in
        else
            touch {{output.tumour}}
        fi
        """
    

### TODO
rule RenameMantaSamples:
    output:
        normal = f'{outdir}/manta-vcf/{{donor}}/{{donor}}.normal.vcf.gz',
        tumour = f'{outdir}/manta-vcf/{{donor}}/{{donor}}.tumour.vcf.gz'
    shell:
        f"""
        aws s3 cp s3://layerlabcu/icgc/manta/{{wildcards.donor}}/diploidSV.vcf.gz \\
            {outdir}/manta-vcf/{{wildcards.donor}}/diploidSV.vcf.gz
        aws s3 cp s3://layerlabcu/icgc/manta/{{wildcards.donor}}/somaticSV.vcf.gz \\
            {outdir}/manta-vcf/{{wildcards.donor}}/somaticSV.vcf.gz
        
        bcftools reheader \\
            -s <(echo {{wildcards.donor}}-normal) \\
            -o {{output.normal}} \\
            {outdir}/manta-vcf/{{wildcards.donor}}/diploidSV.vcf.gz
        bcftools reheader \\
            -s <(echo {{wildcards.donor}}-tumour) \\
            -o {{output.tumour}} \\
            {outdir}/manta-vcf/{{wildcards.donor}}/somaticSV.vcf.gz
        """
        

### TODO
# try using the bioconda version of survivor.
# it uses the same version that we currently have installed
# then we can remove that from installation
rule SurvivorMergeVCFs:
    input:
        smoove_normal = expand(f'{outdir}/smoove-vcf/{{donor}}/{{donor}}.normal.vcf.gz',
                               donor=donors),
        smoove_tumour = expand(f'{outdir}/smoove-vcf/{{donor}}/{{donor}}.tumour.vcf.gz',
                               donor=donors),
        manta_normal = expand(f'{outdir}/manta-vcf/{{donor}}/{{donor}}.normal.vcf.gz',
                              donor=donors),
        manta_tumour = expand(f'{outdir}/manta-vcf/{{donor}}/{{donor}}.tumour.vcf.gz',
                              donor=donors)
    output:
        f'{outdir}/survivor-merged.vcf.gz'
    shell:
        f"""
        cat <(echo {{input.smoove_normal}}) <(echo {{input.smoove_tumour}}) \\
            <(echo {{input.manta_normal}})  <(echo {{input.manta_tumour}}) |
            tr ' ' '\n' | xargs file | grep gzip | cut -d':' -f1 \\
            > {outdir}/vcf-list.txt

        touch {output} # TODO dummy temp file
        """
### TODO
# see manta snakefile for similar rule
# and copy the resource management method
rule GetBams:

### TODO
# can use as many cores as we need.
rule SmooveGenotype:

### TODO
rule SmoovePasteVCFs:


