
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
        f'{outdir}/survivor-merged.vcf'

### TODO handle if there is a missing vcf
# this rule also checks if there are missing vcfs, and produces a dummy output
# using touch.  This will have to be handled down the line before SURVIVOR
# is run (use the "file" command).
rule RenameSmooveSamples:
    output:
        normal = f'{outdir}/smoove-vcf/{{donor}}/{{donor}}.normal.vcf',
        tumour = f'{outdir}/smoove-vcf/{{donor}}/{{donor}}.tumour.vcf'
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
                $normal_in |
                zcat > {{output.normal}}
        else
            # dummy output if search was empty
            touch {{output.normal}} 
        fi
        if [[ ! -z $tumour_in ]]; then
            bcftools reheader \\
                -s <(echo {{wildcards.donor}}-tumour) \\
                $tumour_in |
                zcat > {{output.tumour}}
        else
            touch {{output.tumour}}
        fi
        """
    

### TODO
rule RenameMantaSamples:
    output:
        normal = f'{outdir}/manta-vcf/{{donor}}/{{donor}}.normal.vcf',
        tumour = f'{outdir}/manta-vcf/{{donor}}/{{donor}}.tumour.vcf'
    shell:
        f"""
        aws s3 cp s3://layerlabcu/icgc/manta/{{wildcards.donor}}/diploidSV.vcf.gz \\
            {outdir}/manta-vcf/{{wildcards.donor}}/diploidSV.vcf.gz
        aws s3 cp s3://layerlabcu/icgc/manta/{{wildcards.donor}}/somaticSV.vcf.gz \\
            {outdir}/manta-vcf/{{wildcards.donor}}/somaticSV.vcf.gz
        
        bcftools reheader \\
            -s <(echo {{wildcards.donor}}-normal) \\
            {outdir}/manta-vcf/{{wildcards.donor}}/diploidSV.vcf.gz |
            zcat > {{output.normal}}
        # somaticSV doesn't genotype, I just want the region.
        # it contains two samples, but Ill just extract the first
        somaticSV={outdir}/manta-vcf/{{wildcards.donor}}/somaticSV.vcf.gz
        bcftools view -s $(bcftools query -l $somaticSV | head -1) $somaticSV |
            bgzip -c | 
            bcftools reheader -s <(echo {{wildcards.donor}}-tumour) |
            zcat > {{output.tumour}}
        """
        

### TODO
# try using the bioconda version of survivor.
# it uses the same version that we currently have installed
# then we can remove that from installation
rule SurvivorMergeVCFs:
    input:
        smoove_normal = expand(f'{outdir}/smoove-vcf/{{donor}}/{{donor}}.normal.vcf',
                               donor=donors),
        smoove_tumour = expand(f'{outdir}/smoove-vcf/{{donor}}/{{donor}}.tumour.vcf',
                               donor=donors),
        manta_normal = expand(f'{outdir}/manta-vcf/{{donor}}/{{donor}}.normal.vcf',
                              donor=donors),
        manta_tumour = expand(f'{outdir}/manta-vcf/{{donor}}/{{donor}}.tumour.vcf',
                              donor=donors)
    output:
        f'{outdir}/survivor-merged.vcf'
    shell:
        f"""
        cat <(echo {{input.smoove_normal}}) <(echo {{input.smoove_tumour}}) \\
            <(echo {{input.manta_normal}})  <(echo {{input.manta_tumour}}) |
            tr ' ' '\n' | xargs file | grep -v empty | cut -d':' -f1 \\
            > {outdir}/vcf-list.txt
        
        max_dist_between_breakpoints=0.1 # fraction of SVLEN
        min_support=1
        take_type_into_account=1
        take_strand_into_account=0
        estimate_dist_from_SV_size=1
        min_size=50

        SURVIVOR merge {outdir}/vcf-list.txt \\
            $max_dist_between_breakpoints \\
            $min_support \\
            $take_type_into_account \\
            $take_strand_into_account \\
            $estimate_dist_from_SV_size \\
            $min_size \\
            {{output}}
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


