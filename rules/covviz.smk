
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
    ## TODO just combine into one report
    threads:
        workflow.cores
    input:
        bai = expand(outdir+'/indices/{donor}-{specimen_type}.bam.bai',
                     specimen_type=['normal', 'tumour'], donor=donors),
        fasta = outdir+'/ref/hs37d5.fa',
        fai = outdir+'/ref/hs37d5.fa.fai'
    params:
        outdir = outdir,
        baidir = outdir+'/*/*.bai',
    output:
        outdir+'/covviz_report.html'
    shell:
        """
        nextflow run brwnj/covviz -latest \
            -w /mnt/local \\
            --indexes '{params.baidir}' \
            --fai {input.fai} \
            --outdir {params.outdir}
        aws s3 cp {output} s3://layerlabcu/icgc/covviz/
        """

rule RunCovvizPairwise:
    input:
        bai = expand(f'{outdir}/indices/{{donor}}-{{specimen_type}}.bam.bai',
                     specimen_type=['normal', 'tumour'], donor='{donor}'),
        fasta = f'{outdir}/ref/hs37d5.fa',
        fai = f'{outdir}/ref/hs37d5.fa.fai'
    output:
        f'{outdir}/{{donor}}/covviz_report.html'
    shell:
        f"""
        nextflow run brwnj/covviz -latest \\
            -w /mnt/local \\
            --indexes '{outdir}/indices/{{wildcards.donor}}-*.bam.bai' \\
            --fai {{input.fai}} \\
            --outdir $(dirname {{output}})
        aws s3 cp {{output}} s3://layerlabcu/icgc/covviz/{{wildcards.donor}}/
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
            normal = f'{outdir}/indices/{{donor}}-normal.bam.bai',
            tumour = f'{outdir}/indices/{{donor}}-tumour.bam.bai'
        shell:
            """
            aws s3 cp s3://layerlabcu/icgc/bam_indices/$(basename {output.normal}) \
                {output.normal}
            aws s3 cp s3://layerlabcu/icgc/bam_indices/$(basename {output.tumour}) \
                {output.tumour}
            """
else:
    ## TODO update directory structure if I need to redownload
    rule CombineManifests:
        input:
            expand(manifest_dir+'/{donor}-tumour-normal.tsv', donor=donors)
        output:
            combined_manifest = temp(outdir+'/combined-manifest.tsv')
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
            outdir+'/combined-manifest.tsv'
        output:
            temp(outdir+'/mounted-successfully.out')
        params:
            mountdir = outdir+'/temp'
        shell:
            """
            [[ ! -d {params.mountdir} ]] && mkdir {params.mountdir}
            if mount | grep -q {params.mountdir} ; then
                touch {output}
                exit 0
            fi

            score-client mount --daemonize \
                --mount-point {params.mountdir} \
                --manifest {input} 
            sleep 30s # give time for the directory to be mounted
            touch {output}
            """

    rule CreateOutdirs:
        output:
            normal = directory(outdir+'/normal'),
            tumour = directory(outdir+'/tumour')
        shell:
            """
            [[ ! -d {output.normal} ]] && mkdir {output.normal}
            [[ ! -d {output.tumour} ]] && mkdir {output.tumour}
            """
            
    rule ScoreClientGetBamIndex:
        input:
            manifest = manifest_dir+'/{donor}-tumour-normal.tsv',
            receipt = outdir+'/mounted-successfully.out',
            normal_dir = directory(outdir+'/normal'),
            tumour_dir = directory(outdir+'/tumour')
        output:
            normal = outdir+'/normal/{donor}-normal.bam.bai',
            tumour = outdir+'/tumour/{donor}-tumour.bam.bai'
        params:
            mountdir = outdir+'/temp'
        shell:
            """
            normal_bam=$(sed '2q;d' {input.manifest} | cut -f5)
            tumour_bam=$(sed '3q;d' {input.manifest} | cut -f5)
            normal_bai=$(find {params.mountdir} -name '*.bai' | grep $normal_bam)
            tumour_bai=$(find {params.mountdir} -name '*.bai' | grep $tumour_bam)
            cp --no-preserve=mode $normal_bai {output.normal}
            cp --no-preserve=mode $tumour_bai {output.tumour}

            if ! aws s3 ls s3://layerlabcu/icgc/bam_indices | 
                grep -q $(basename {output.normal}); then
                aws s3 cp {output.normal} s3://layerlabcu/icgc/bam_indices/
            fi

            if ! aws s3 ls s3://layerlabcu/icgc/bam_indices | 
                grep -q $(basename {output.tumour}); then
                aws s3 cp {output.tumour} s3://layerlabcu/icgc/bam_indices/
            fi
            """
