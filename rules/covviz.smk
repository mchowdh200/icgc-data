
outdir = config['outdir']
manifest_dir = config['manifest_dir']
with open(config['donor_list']) as f:
    donors = [x.rstrip() for x in f.readlines()]
get_from_s3 = config['get_from_s3']

rule all:
    input:
        expand(outdir+'/{specimen_type}/results/covviz_report.html',
               specimen_type=['normal', 'tumour'])

rule RunCovviz:
    threads:
        workflow.cores
    input:
        bai = expand(outdir+'/{specimen_type}/{donor}-normal.bam.bai',
                     specimen_type='{specimen_type}', donor=donors),
        fasta = outdir+'/ref/hs37d5.fa',
        fai = outdir+'/ref/hs37d5.fa.fai'
    params:
        baidir = outdir+'/{specimen_type}',
    output:
        outdir+'/{specimen_type}/results/covviz_report.html'
    shell:
        """
        nexflow run brwnj/covviz -latest \
            --indexes '{params.baidir}/*.bai'
            --fai {input.fasta} \
            --outdir {params.baidir}
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
            normal = outdir+'/normal/{donor}-normal.bam.bai',
            tumour = outdir+'/tumour/{donor}-tumour.bam.bai'
        shell:
            """
            aws s3 cp s3://layerlabcu/icgc/bam_indices/$(basename {output.normal}) \
                {output.normal}
            aws s3 cp s3://layerlabcu/icgc/bam_indices/$(basename {output.tumour}) \
                {output.tumour}
            """
else:
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
