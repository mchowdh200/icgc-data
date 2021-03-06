outdir = config['outdir']
manifest_dir = config['manifest_dir']
max_bams = config['max_bams']
ref_s3_path = config['ref_s3_path']

### Full donor list for SURVIVOR
with open(config['donor_list']) as f:
    donors_all = [x.rstrip() for x in f.readlines()]

### TODO Donor list part for SVTyper
part = config['donor_list_part'][-1]
with open(config['donor_list_part']) as f:
    donors_part = [x.rstrip() for x in f.readlines()]


rule all:
    ## TODO this should be the final SVTyper sites vcf.
    # input:
    #     f'{outdir}/sites.smoove.square.vcf.gz'
    # shell:
    #     """
    #     aws s3 cp {input} s3://layerlabcu/icgc/
    #     """
    input:
        svtyper = expand(f'{outdir}/svtyper-vcf/{{donor_part}}-{{specimen_type}}-smoove.genotyped.vcf.gz',
                        donor_part=donors_part, specimen_type=['normal', 'tumour']),
        survivor = f'{outdir}/survivor-merged.vcf'
   

rule RenameSmooveSamples:
    threads:
        1
    output:
        normal = f'{outdir}/smoove-vcf/{{donor_all}}/{{donor_all}}.normal.vcf',
        tumour = f'{outdir}/smoove-vcf/{{donor_all}}/{{donor_all}}.tumour.vcf'
    shell:
        f"""
        aws s3 cp --recursive \\
            s3://layerlabcu/icgc/smoove/{{wildcards.donor_all}}/ \\
            {outdir}/smoove-vcf/{{wildcards.donor_all}}/
        
        set +e # grep returns non zero if search is empty
        normal_in=$(find {outdir}/smoove-vcf/{{wildcards.donor_all}} -name '*.vcf.gz' |
                    grep -i normal)
        tumour_in=$(find {outdir}/smoove-vcf/{{wildcards.donor_all}} -name '*.vcf.gz' |
                    grep -i tumour | head -1)
        set e

        if [[ ! -z $normal_in ]]; then
            bcftools reheader \\
                -s <(echo {{wildcards.donor_all}}-normal) \\
                $normal_in |
                zcat > {{output.normal}}
        else
            # dummy output if search was empty
            touch {{output.normal}} 
        fi
        if [[ ! -z $tumour_in ]]; then
            bcftools reheader \\
                -s <(echo {{wildcards.donor_all}}-tumour) \\
                $tumour_in |
                zcat > {{output.tumour}}
        else
            touch {{output.tumour}}
        fi
        """

rule RenameMantaSamples:
    threads:
        1
    output:
        normal = f'{outdir}/manta-vcf/{{donor_all}}/{{donor_all}}.normal.vcf',
        tumour = f'{outdir}/manta-vcf/{{donor_all}}/{{donor_all}}.tumour.vcf'
    shell:
        f"""
        aws s3 cp s3://layerlabcu/icgc/manta/{{wildcards.donor_all}}/diploidSV.vcf.gz \\
            {outdir}/manta-vcf/{{wildcards.donor_all}}/diploidSV.vcf.gz
        aws s3 cp s3://layerlabcu/icgc/manta/{{wildcards.donor_all}}/somaticSV.vcf.gz \\
            {outdir}/manta-vcf/{{wildcards.donor_all}}/somaticSV.vcf.gz
        
        bcftools reheader \\
            -s <(echo {{wildcards.donor_all}}-normal) \\
            {outdir}/manta-vcf/{{wildcards.donor_all}}/diploidSV.vcf.gz |
            zcat > {{output.normal}}
        # somaticSV doesn't genotype, I just want the region.
        # it contains two samples, but Ill just extract the first
        somaticSV={outdir}/manta-vcf/{{wildcards.donor_all}}/somaticSV.vcf.gz
        bcftools view -s $(bcftools query -l $somaticSV | head -1) $somaticSV |
            bgzip -c | 
            bcftools reheader -s <(echo {{wildcards.donor_all}}-tumour) |
            zcat > {{output.tumour}}
        """

rule SurvivorMergeVCFs:
    threads:
        1
    input:
        smoove_normal = expand(f'{outdir}/smoove-vcf/{{donor_all}}/{{donor_all}}.normal.vcf',
                               donor_all=donors_all),
        smoove_tumour = expand(f'{outdir}/smoove-vcf/{{donor_all}}/{{donor_all}}.tumour.vcf',
                               donor_all=donors_all),
        manta_normal = expand(f'{outdir}/manta-vcf/{{donor_all}}/{{donor_all}}.normal.vcf',
                              donor_all=donors_all),
        manta_tumour = expand(f'{outdir}/manta-vcf/{{donor_all}}/{{donor_all}}.tumour.vcf',
                              donor_all=donors_all)
    output:
        f'{outdir}/survivor-merged.vcf'
    shell:
        f"""
        cat <(echo {{input.smoove_normal}}) <(echo {{input.manta_normal}}) \\
            <(echo {{input.smoove_tumour}})  <(echo {{input.manta_tumour}}) |
            tr ' ' '\n' | xargs file | grep -v empty | cut -d':' -f1 \\
            > {outdir}/vcf-list.txt
        
        max_dist_between_breakpoints=0.3 # fraction of SVLEN
        min_support=4
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
        aws s3 cp {{output}} s3://layerlabcu/icgc/SURVIVOR/
        """

rule GetReference:
    output:
        fasta = temp(f'{outdir}/ref/hs37d5.fa'),
        fai = temp(f'{outdir}/ref/hs37d5.fa.fai'),
    shell:
        f"""
        aws s3 cp {ref_s3_path} {{output.fasta}}
        aws s3 cp {ref_s3_path}.fai {{output.fai}}
        """

rule GetBam:
    threads:
        max(workflow.cores//6, 1)
    resources:
        num_downloads = 1
    input:
        manifest = f'{manifest_dir}/{{donor_part}}-tumour-normal.tsv'
    output:
        bam = temp(f'{outdir}/{{donor_part}}/{{donor_part}}-{{specimen_type}}.bam'),
        bai = temp(f'{outdir}/{{donor_part}}/{{donor_part}}-{{specimen_type}}.bam.bai')
    shell:
        f"""
        if [[ ! -d {outdir}/{{wildcards.donor_part}} ]]; then
            mkdir {outdir}/{{wildcards.donor_part}}
        fi

        # wait for resources to be available
        while [[ $(find {outdir} -name '*.bam' -or -name '*.cram' | wc -l) -ge {max_bams} ]]; do
            sleep 5
        done
        
        # get manifest containing just specimen type we are after
        # after the header, normal is the first entry and tumour is the second
        # the bam filename is on the fifth column of the manifest.
        specimen_type="{{wildcards.specimen_type}}"
        if [[ $specimen_type == "normal" ]]; then
           sed '3d' {{input.manifest}} > {outdir}/{{wildcards.donor_part}}/{{wildcards.donor_part}}-{{wildcards.specimen_type}}.tsv
        else
           sed '2d' {{input.manifest}} > {outdir}/{{wildcards.donor_part}}/{{wildcards.donor_part}}-{{wildcards.specimen_type}}.tsv
        fi
        bam_fname={outdir}/{{wildcards.donor_part}}/$(tail -1 {outdir}/{{wildcards.donor_part}}/{{wildcards.donor_part}}-{{wildcards.specimen_type}}.tsv | cut -f5)
        
        score-client download \\
            --validate false \\
            --output-dir {outdir}/{{wildcards.donor_part}} \\
            --manifest {outdir}/{{wildcards.donor_part}}/{{wildcards.donor_part}}-{{wildcards.specimen_type}}.tsv
        
        mv $bam_fname {{output.bam}}
        mv $bam_fname.bai {{output.bai}}
        """

# we need to do this in order for the samples to match when running SVTyper :(
# also, I'm making the output crams to save space.
rule ReplaceReadGroups:
    threads:
        max(workflow.cores//4, 1)
    input:
        fasta = f'{outdir}/ref/hs37d5.fa',
        fai = f'{outdir}/ref/hs37d5.fa.fai',
        bam = f'{outdir}/{{donor_part}}/{{donor_part}}-{{specimen_type}}.bam',
        bai = f'{outdir}/{{donor_part}}/{{donor_part}}-{{specimen_type}}.bam.bai'
    output:
        cram = temp(f'{outdir}/{{donor_part}}/{{donor_part}}-{{specimen_type}}.cram'),
        crai = temp(f'{outdir}/{{donor_part}}/{{donor_part}}-{{specimen_type}}.cram.crai')
    shell:
        f"""
        sample_name="{{wildcards.donor_part}}-{{wildcards.specimen_type}}"
        samtools addreplacerg -r "ID:$sample_name" \\
                              -r "SM:$sample_name" \\
                              -r "LB:$sample_name" \\
                              -@ {{threads}} \\
                              --reference {{input.fasta}} \\
                              -O CRAM -o {{output.cram}} \\
                              {{input.bam}}
        samtools index -b -@ {{threads}} {{output.cram}}
        """

### TODO
rule SmooveGenotype:
    threads:
        max(workflow.cores//3, 1)
    resources:
        SVTyper_instances=1
    conda:
        'envs/smoove.yaml'
    input:
        fasta = f'{outdir}/ref/hs37d5.fa',
        fai = f'{outdir}/ref/hs37d5.fa.fai',
        cram = f'{outdir}/{{donor_part}}/{{donor_part}}-{{specimen_type}}.cram',
        crai = f'{outdir}/{{donor_part}}/{{donor_part}}-{{specimen_type}}.cram.crai',
        vcf = f'{outdir}/survivor-merged.vcf'
    output:
        f'{outdir}/svtyper-vcf/{{donor_part}}-{{specimen_type}}-smoove.genotyped.vcf.gz'
    shell:
        f"""
        smoove genotype -p {{threads}} -f {{input.fasta}} -v {{input.vcf}}\\
            -n {{wildcards.donor_part}}-{{wildcards.specimen_type}} \\
            -o {outdir}/svtyper-vcf/ \\
            {{input.cram}}
        aws s3 cp {{output}} s3://layerlabcu/icgc/svtyper/$(basename {{output}})
        """

### This is so fast I'll just do it manually later.
# rule SmoovePasteVCFs:
#     conda:
#         'envs/smoove.yaml'
#     input:
#         expand(f'{outdir}/svtyper-vcf/{{donor_all}}-{{specimen_type}}-smoove.genotyped.vcf.gz',
#                donor_all=donors_all, specimen_type=['normal', 'tumour'])
#     output:
#         f'{outdir}/sites.smoove.square.vcf.gz'
#     shell:
#         f"""
#         smoove paste -o {outdir} --name sites {outdir}/svtyper-vcf/*.vcf.gz
#         """
