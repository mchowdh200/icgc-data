
outdir = config['outdir']
manifest_dir = config['manifest_dir']
with open(config['donor_list']) as f:
    donors = [x.rstrip() for x in f.readlines()]
get_from_s3 = config['get_from_s3']

rule all:
    input:
        f'{outdir}/covviz_report.html',
        expand(f'{outdir}/{{donor}}/covviz_report.html', donor=donors)

rule RunCovviz:
    threads:
        workflow.cores
    input:
        bai = expand(f'{outdir}/indices/{{donor}}-{{specimen_type}}.bai',
                     specimen_type=['normal', 'tumour'], donor=donors),
        fasta = f'{outdir}/ref/hs37d5.fa',
        fai = f'{outdir}/ref/hs37d5.fa.fai',
        svtyper_variants=f'{outdir}/annotations/squared.sites.vcf.gz',
        LNCaP_variants=f'{outdir}/annotations/LNCAPEXP_REFINEFINAL1.vcf',
        genes_hg19=f'{output}/annotations/genes_hg19.bed'
    output:
        f'{outdir}/covviz_report.html'
    shell:
        f"""
        bcftools annotate -x '^INFO/SVTYPE,INFO/SVLEN,INFO/SUPP' \\
            {{input.svtyper_variants}} > {outdir}/annotations/svtyper_filtered_info.vcf

        goleft indexcov \\
            --directory {outdir}/covviz-all \\
            --fai {{input.fai}} \\
            {outdir}/indices/*.bai
        covviz --output {outdir}/covviz_report.html \\
               --vcf {outdir}/annotations/svtyper_filtered_info.vcf \\
               --vcf {{input.LNCaP_variants}} \\
               --bed {{input.genes_hg19}} \\
               --ped {outdir}/covviz-all/covviz-all-indexcov.ped \\
               {outdir}/covviz-all/covviz-all-indexcov.bed.gz
        aws s3 cp {{output}} s3://layerlabcu/icgc/covviz/
        aws s3 cp {{output}} s3://icgc-vis/
        """

rule RunCovvizPairwise:
    input:
        bai = expand(f'{outdir}/indices/{{donor}}-{{specimen_type}}.bai',
                     specimen_type=['normal', 'tumour'], donor='{donor}'),
        fasta = f'{outdir}/ref/hs37d5.fa',
        fai = f'{outdir}/ref/hs37d5.fa.fai',
        svtyper_variants=f'{outdir}/annotations/squared.sites.vcf.gz',
        LNCaP_variants=f'{outdir}/annotations/LNCAPEXP_REFINEFINAL1.vcf',
        genes_hg19=f'{output}/annotations/genes_hg19.bed'
    output:
        f'{outdir}/{{donor}}/covviz_report.html'
    shell:
        f"""
        bcftools view -s {{wildcards.donor}}-tumour,{{wildcards.donor}}-normal \\
            {{input.svtyper_variants}} |
        bcftools annotate -x '^INFO/SVTYPE,INFO/SVLEN,INFO/SUPP' > {outdir}/{{wildcards.donor}}/donor-variants.vcf
        goleft indexcov \\
            --directory {outdir}/{{wildcards.donor}} \\
            --fai {{input.fai}} \\
            {outdir}/indices/{{wildcards.donor}}-*.bai
        covviz --output {outdir}/{{wildcards.donor}}/covviz_report.html \\
               --vcf {outdir}/{{wildcards.donor}}/donor-variants.vcf \\
               --vcf {{input.LNCaP_variants}} \\
               --bed {{input.genes_hg19}} \\
               --ped {outdir}/{{wildcards.donor}}/{{wildcards.donor}}-indexcov.ped \\
               {outdir}/{{wildcards.donor}}/{{wildcards.donor}}-indexcov.bed.gz
        aws s3 cp {{output}} s3://layerlabcu/icgc/covviz/{{wildcards.donor}}/
        aws s3 cp {{output}} s3://icgc-vis/{{wildcards.donor}}/
        """

rule GetAnnotationRegions:
    output:
        svtyper_variants=f'{outdir}/annotations/squared.sites.vcf.gz',
        LNCaP_variants=f'{outdir}/annotations/LNCAPEXP_REFINEFINAL1.vcf',
        genes_hg19=f'{output}/annotations/genes_hg19.bed'
    shell:
        f"""
        aws s3 cp --recursive s3://layerlabcu/icgc/misc_annotations/ {outdir}/annotations/
        aws s3 cp s3://layerlabcu/icgc/svtyper/squared.sites.vcf.gz {outdir}/annotations/
        """
    
rule GetReference:
    output:
        fasta = outdir+'/ref/hs37d5.fa',
        fai = outdir+'/ref/hs37d5.fa.fai'
    shell:
        """
        aws s3 cp s3://layerlabcu/ref/genomes/hs37d5/hs37d5.fa {output.fasta}
        aws s3 cp s3://layerlabcu/ref/genomes/hs37d5/hs37d5.fa.fai {output.fai}
        """

if get_from_s3:
    rule S3GetBamIndex:
        output:
            expand(f'{outdir}/indices/{{donor}}-{{specimen_type}}.bai',
                   donor=donors, specimen_type=['normal', 'tumour']),
        shell:
            f"""
            aws s3 cp --recursive s3://layerlabcu/icgc/bam_indices/ {outdir}/indices/
            """
else:
    rule CombineManifests:
        input:
            expand(f'{manifest_dir}/{{donor}}-tumour-normal.tsv', donor=donors)
        output:
            combined_manifest = f'{outdir}/combined-manifest.tsv'
        run:
            # get top line of a manifest
            with open(input[0]) as manifest:
                header = manifest.readline()

            # write the rest 
            with open(output.combined_manifest, 'w') as combined_manifest:
                combined_manifest.write(header)
                for fname in input:
                    with open(fname, 'r') as f:
                        lines = f.readlines()
                        for line in lines[1:]: #skip header
                            combined_manifest.write(line)

    rule MountDirectory:
        input:
            f'{outdir}/combined-manifest.tsv'
        output:
            temp(f'{outdir}/mounted-successfully.out')
        shell:
            f"""
            [[ ! -d {outdir}/temp ]] && mkdir {outdir}/temp
            if mount | grep -q {outdir}/temp ; then
                touch {{output}}
                exit 0
            fi

            score-client mount --daemonize \\
                --mount-point {outdir}/temp \\
                --manifest {{input}} 
            sleep 30s # give time for the directory to be mounted
            touch {{output}}
            """

            
    rule ScoreClientGetBamIndex:
        threads:
            workflow.cores # i don't want to overwhelm the FUSE filesystem.
        input:
            manifest = f'{manifest_dir}/{{donor}}-tumour-normal.tsv',
            receipt = f'{outdir}/mounted-successfully.out',
        output:
            normal = f'{outdir}/indices/{{donor}}-normal.bai',
            tumour = f'{outdir}/indices/{{donor}}-tumour.bai'
        shell:
            f"""
            if [[ ! -d {outdir}/indices ]]; then
                mkdir {outdir}/indices
            fi
            normal_bam=$(sed '2q;d' {{input.manifest}} | cut -f5)
            tumour_bam=$(sed '3q;d' {{input.manifest}} | cut -f5)
            normal_bai=$(find {outdir}/temp -name '*.bai' | grep $normal_bam)
            tumour_bai=$(find {outdir}/temp -name '*.bai' | grep $tumour_bam)
            sleep 1s
            cp --no-preserve=mode $normal_bai {{output.normal}}
            sleep 1s
            cp --no-preserve=mode $tumour_bai {{output.tumour}}

            aws s3 cp {{output.normal}} s3://layerlabcu/icgc/bam_indices/
            aws s3 cp {{output.tumour}} s3://layerlabcu/icgc/bam_indices/
            """
