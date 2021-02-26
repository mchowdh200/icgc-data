
outdir = config['outdir']
manifest_dir = config['manifest_dir']
with open(config['donor_list']) as f:
    donors = [x.rstrip() for x in f.readlines()]

rule all:
    input:
        expand(outdir+'/{donor}/{donor}-normal.bai', donor=donors),
        expand(outdir+'/{donor}/{donor}-tumour.bai', donor=donors)

rule get_bam_index:
    input:
        manifest = manifest_dir+'/{donor}-tumour-normal.tsv'
    output:
        normal = outdir+'/{donor}/{donor}-normal.bai',
        tumour = outdir+'/{donor}/{donor}-normal.bai'
    shell:
        """
        [[ ! -d /mnt/{donor} ]] && sudo mkdir /mnt/{donor}
        sudo chown ubuntu /mnt/{donor}

        score-client mount --daemonize \
            --mount-point /mnt/{donor} \
            --manifest {input.manifest} 
        
        normal_bam=$(sed '2q;d {input.manifest} | cut -f5)
        tumour_bam=$(sed '3q;d {input.manifest} | cut -f5)
        
        normal_bai=$(find /mnt/{donor} -name *.bai | grep $normal_bam)
        tumour_bai=$(find /mnt/{donor} -name *.bai | grep $tumor_bam)
        
        cp $normal_bai {output.normal}
        cp $tumour_bai {output.tumour}
        """
