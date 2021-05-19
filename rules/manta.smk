from pathlib import Path
import numpy as np
import pandas as pd

### Setup
###################################################
outdir = config['outdir']
manifest_table = pd.read_csv(config['manifest'], sep='\t')
file_ids = manifest_table['file_id'].tolist()

### Helper Functions
###################################################
def bam_disk_mb(wildcards):
    return int(manifest_table[
        manifest_table.file_id == wildcards.file_id].file_size * 1e-6)

### Rules
###################################################

rule All:
    input:
        expand(f'{outdir}/{{file_id}}/{{file_id}}-manta.vcf.gz',
               file_id=file_ids)

rule GetReference:
    output:
        fasta = temp(f'{outdir}/ref/hs37d5.fa'),
        fai = temp(f'{outdir}/ref/hs37d5.fa.fai')
    shell:
        """
        aws s3 cp s3://layerlabcu/ref/genomes/hs37d5/hs37d5.fa {output.fasta}
        aws s3 cp s3://layerlabcu/ref/genomes/hs37d5/hs37d5.fa.fai {output.fai}
        """

rule GetManifest:
    output:
        f'{outdir}/manifests/{{file_id}}-manifest.tsv'
    run:
        Path(f'{outdir}/manifests').mkdir(exist_ok=True)
        m = manifest_table[manifest_table['file_id'] == wildcards.file_id]
        m.to_csv(output[0], index=False, sep='\t')

rule GetBam:
    input:
        manifest = f'{outdir}/manifests/{{file_id}}-manifest.tsv'
    output:
        bam = temp(f"{outdir}/bam/{{file_id}}.bam"),
        bai = temp(f"{outdir}/bam/{{file_id}}.bai")
    resources:
        disk_mb = bam_disk_mb,
        num_downloads = 1
    run:
        Path(f'{outdir}/bam').mkdir(exist_ok=True)
        shell(f"""score-client download \\
        --validate false \\
        --output-dir {outdir}/bam \\
        --manifest {input.manifest}""")

        # remove bam and rename index with file_id
        bam = f'{outdir}/bam/{manifest_table[manifest_table.file_id == wildcards.file_id].file_name.values[0]}'
        bai = f'{bam}.bai'
        Path(bam).rename(output.bam)
        Path(bai).rename(output.bai)

rule RunManta:
    threads:
        np.ceil(workflow.cores/2)
    input:
        fasta = rules.GetReference.output.fasta,
        fai = rules.GetReference.output.fai,
        bam = rules.GetBam.output.bam,
        bai = rules.GetBam.output.bai
    output:
        f'{outdir}/{{file_id}}/{{file_id}}-manta.vcf.gz'
    shell:
        f"""
        mkdir -p {outdir}/{{wildcards.file_id}}
        /mnt/local/manta/bin/configManta.py \\
            --bam {{input.bam}} \\
            --referenceFasta {{input.fasta}} \\
            --runDir {outdir}/{{wildcards.file_id}}
        python2 {outdir}/{{wildcards.file_id}}/runWorkflow.py -j {{threads}}

        cp {outdir}/{{wildcards.file_id}}/results/variants/diploidSV.vcf.gz {{output}}
        find {outdir}/{{wildcards.file_id}}/ -type f -not -name '{{wildcards.file_id}}-manta.vcf.gz' -delete
        """
        
        
