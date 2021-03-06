outdir = config['outdir']
manifest_dir = config['manifest_dir']
max_bams = config['max_bams']
ref_s3_path = config['ref_s3_path']
with open(config['donor_list']) as f:
    donors = [x.rstrip() for x in f.readlines()]
### TODO testing purposes
donors=donors[:2]
rule all:
    ## TODO this should be the final SVTyper sites vcf.
    input:
        f'{outdir}/sites.smoove.square.vcf.gz'
    # shell:
    #     """
    #     aws s3 cp {input} s3://layerlabcu/icgc/
    #     """


rule RenameSmooveSamples:
    threads:
        1
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
    threads:
        1
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

rule SurvivorMergeVCFs:
    threads:
        1
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
        max(workflow.cores//12, 1)
    resources:
        num_downloads = 1
    input:
        manifest = f'{manifest_dir}/{{donor}}-tumour-normal.tsv'
    output:
        bam = temp(f'{outdir}/{{donor}}/{{donor}}-{{specimen_type}}.bam'),
        bai = temp(f'{outdir}/{{donor}}/{{donor}}-{{specimen_type}}.bam.bai')
    shell:
        f"""
        if [[ ! -d {outdir}/{{wildcards.donor}} ]]; then
            mkdir {outdir}/{{wildcards.donor}}
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
           sed '3d' {{input.manifest}} > {outdir}/{{wildcards.donor}}/{{wildcards.donor}}-{{wildcards.specimen_type}}.tsv
        else
           sed '2d' {{input.manifest}} > {outdir}/{{wildcards.donor}}/{{wildcards.donor}}-{{wildcards.specimen_type}}.tsv
        fi
        bam_fname={outdir}/{{wildcards.donor}}/$(tail -1 {outdir}/{{wildcards.donor}}/{{wildcards.donor}}-{{wildcards.specimen_type}}.tsv | cut -f5)
        
        score-client download \\
            --validate false \\
            --output-dir {outdir}/{{wildcards.donor}} \\
            --manifest {outdir}/{{wildcards.donor}}/{{wildcards.donor}}-{{wildcards.specimen_type}}.tsv
        
        mv $bam_fname {{output.bam}}
        mv $bam_fname.bai {{output.bai}}
        """

# we need to do this in order for the samples to match when running SVTyper :(
# also, I'm making the output crams to save space.
rule ReplaceReadGroups:
    threads:
        max(workflow.cores//6, 1)
    input:
        fasta = f'{outdir}/ref/hs37d5.fa',
        fai = f'{outdir}/ref/hs37d5.fa.fai',
        bam = f'{outdir}/{{donor}}/{{donor}}-{{specimen_type}}.bam',
        bai = f'{outdir}/{{donor}}/{{donor}}-{{specimen_type}}.bam.bai'
    output:
        cram = temp(f'{outdir}/{{donor}}/{{donor}}-{{specimen_type}}.cram'),
        crai = temp(f'{outdir}/{{donor}}/{{donor}}-{{specimen_type}}.cram.crai')
    shell:
        f"""
        sample_name="{{wildcards.donor}}-{{wildcards.specimen_type}}"
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
        cram = f'{outdir}/{{donor}}/{{donor}}-{{specimen_type}}.cram',
        crai = f'{outdir}/{{donor}}/{{donor}}-{{specimen_type}}.cram.crai',
        vcf = f'{outdir}/survivor-merged.vcf'
    output:
        f'{outdir}/svtyper-vcf/{{donor}}-{{specimen_type}}-smoove.genotyped.vcf.gz'
    shell:
        f"""
        smoove genotype -p {{threads}} -f {{input.fasta}} -v {{input.vcf}}\\
            -n {{wildcards.donor}}-{{wildcards.specimen_type}} \\
            -o {outdir}/svtyper-vcf/ \\
            {{input.cram}}
        """

rule SmoovePasteVCFs:
    conda:
        'envs/smoove.yaml'
    input:
        expand(f'{outdir}/svtyper-vcf/{{donor}}-{{specimen_type}}-smoove.genotyped.vcf.gz',
               donor=donors, specimen_type=['normal', 'tumour'])
    output:
        f'{outdir}/sites.smoove.square.vcf.gz'
    shell:
        f"""
        smoove paste -o {outdir} --name sites {outdir}/svtyper-vcf/*.vcf.gz
        """
