
outdir = config['outdir']
manifest_dir = config['manifest_dir']
with open(config['donor_list']) as f:
    donors = [x.rstrip() for x in f.readlines()]

rule all:
    input:
        expand(outdir+'/{donor}/{donor}-normal.bai', donor=donors),
        expand(outdir+'/{donor}/{donor}-tumour.bai', donor=donors)

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
        [[ ! -d {params.mountdir}]] && mkdir {params.mountdir}
        mount | grep -q {params.mountdir} || unmount {params.mountdir}

        score-client mount --daemonize \
            --mount-point {params.mountdir} \
            --manifest {input} 
        touch {output}
        """

rule get_bam_index:
    input:
        manifest = manifest_dir+'/{donor}-tumour-normal.tsv',
        receipt = outdir+'/mounted-successfully.out'
    output:
        normal = outdir+'/{donor}/{donor}-normal.bai',
        tumour = outdir+'/{donor}/{donor}-tumour.bai'
    params:
        mountdir = outdir+'/temp'
    shell:
        """
        normal_bam=$(sed '2q;d' {input.manifest} | cut -f5)
        tumour_bam=$(sed '3q;d' {input.manifest} | cut -f5)
        normal_bai=$(find {params.mountdir} -name *.bai | grep $normal_bam)
        tumour_bai=$(find {params.mountdir} -name *.bai | grep $tumor_bam)
        cp $normal_bai {output.normal}
        cp $tumour_bai {output.tumour}
        """
