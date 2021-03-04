
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

rule SurvivorMergeByDonor:
    ## Get merged regions by caller tumour-noral,
    ## then merge into a single donor vcf.
    input:
        smoove_normal = f'{outdir}/smoove-vcf/{{donor}}/{{donor}}.normal.vcf',
        smoove_tumour = f'{outdir}/smoove-vcf/{{donor}}/{{donor}}.tumour.vcf',
        manta_normal = f'{outdir}/manta-vcf/{{donor}}/{{donor}}.normal.vcf',
        manta_tumour = f'{outdir}/manta-vcf/{{donor}}/{{donor}}.tumour.vcf'
    output:
        f'{outdir}/{{donor}}/{{donor}}-merged.vcf'
    shell:
        f"""
        ## SURVIVOR Params
        max_dist_between_breakpoints=0.3 # fraction of SVLEN
        min_support=1
        take_type_into_account=1
        take_strand_into_account=0
        estimate_dist_from_SV_size=1
        min_size=50
        
        ## Merge the smoove tumour/normal vcfs
        # ugh. we need to handle potential empty vcfs
        printf "{{input.smoove_normal}}\n{{input.smoove_tumour}}\n" |
            xargs file | grep -v empty | cut -d':' -f1 \\
            > {outdir}/{{wildcards.donor}}/smoove-vcf-list.txt
        if [[ $(cat {outdir}/{{wildcards.donor}}/smoove-vcf-list.txt | wc -l) -eq 2 ]]; then
            SURVIVOR merge {outdir}/{{wildcards.donor}}/smoove-vcf-list.txt \\
                $max_dist_between_breakpoints \\
                $min_support \\
                $take_type_into_account \\
                $take_strand_into_account \\
                $estimate_dist_from_SV_size \\
                $min_size \\
                {outdir}/{{wildcards.donor}}/smoove-merged.vcf
        else
            cp $(head -1 {outdir}/{{wildcards.donor}}/smoove-vcf-list.txt) \\
               {outdir}/{{wildcards.donor}}/smoove-merged.vcf
        fi
        
        ## Merge manta tumour/normal vcfs
        printf "{{input.manta_normal}}\n{{input.manta_tumour}}\n" \\
            > {outdir}/{{wildcards.donor}}/manta-vcf-list.txt
        SURVIVOR merge {outdir}/{{wildcards.donor}}/manta-vcf-list.txt \\
            $max_dist_between_breakpoints \\
            $min_support \\
            $take_type_into_account \\
            $take_strand_into_account \\
            $estimate_dist_from_SV_size \\
            $min_size \\
            {outdir}/{{wildcards.donor}}/manta-merged.vcf
        
        ## Merge caller-merged vcfs
        printf "{outdir}/{{wildcards.donor}}/smoove-merged.vcf\n{outdir}/{{wildcards.donor}}/manta-merged.vcf\n" \\
            > {outdir}/{{wildcards.donor}}/caller-merged-vcf-list.txt
        SURVIVOR merge {outdir}/{{wildcards.donor}}/caller-merged-vcf-list.txt \\
            $max_dist_between_breakpoints \\
            $min_support \\
            $take_type_into_account \\
            $take_strand_into_account \\
            $estimate_dist_from_SV_size \\
            $min_size \\
            {{output}}
        """

rule SurvivorMergeDonors:
    ## merge all donor vcfs
    input:
        expand(f'{outdir}/{{donor}}/{{donor}}-merged.vcf', donor=donors)
    output:
        f'{outdir}/survivor-merged.vcf'
    shell:
        f"""
        ## SURVIVOR params
        max_dist_between_breakpoints=0.3 # fraction of SVLEN
        min_support=4
        take_type_into_account=1
        take_strand_into_account=0
        estimate_dist_from_SV_size=1
        min_size=50

        ## Merge all donor-vcfs
        echo {{input}} > {outdir}/donor-vcf-list.txt
        SURVIVOR merge {outdir}/donor-vcf-list.txt \\
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


