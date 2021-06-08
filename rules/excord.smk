import random
from pathlib import Path
import numpy as np
import pandas as pd

### Setup
###################################################
outdir = config['outdir']
manifest_table = pd.read_csv(config['manifest'], sep='\t')
file_ids = manifest_table['file_id'].tolist()

if (i := config['partition']) > 0:
    random.seed(42)
    random.shuffle(file_ids)
    n = len(file_ids)
    file_ids = file_ids[(i-1)*n//4 : i*n//4]

### Rules
###################################################
rule All:
    input:
        expand(f'{outdir}/bed/{{file_id}}.excord.bed.gz', file_id=file_ids)
    shell:
        f'aws s3 cp {outdir}/bed/ s3://layerlabcu/icgc/excord.2/'

    
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
    run:
        Path(f'{outdir}/bam').mkdir(exist_ok=True)
        shell(f"""score-client --quiet download \\
        --validate false \\
        --output-dir {outdir}/bam \\
        --manifest {input.manifest}""")

        # remove bam and rename index with file_id
        bam = f'{outdir}/bam/{manifest_table[manifest_table.file_id == wildcards.file_id].file_name.values[0]}'
        bai = f'{bam}.bai'
        Path(bam).rename(output.bam)
        Path(bai).rename(output.bai)

rule GetExcord:
    output:
        excord = f'{outdir}/excord'
    shell:
        """
        wget -O {output.excord} https://github.com/brentp/excord/releases/download/v0.2.2/excord_linux64
        chmod +x {output.excord}
        """
    
rule RunExcord:
    input:
        bam = rules.GetBam.output.bam,
        bai = rules.GetBam.output.bai,
        fasta = rules.GetReference.output.fasta,
        fai = rules.GetReference.output.fai,
        excord = rules.GetExcord.output.excord,
    output:
        f'{outdir}/bed/{{file_id}}.excord.bed.gz'
    log:
        'logs/{file_id}-excord.log'
    shell:
        'bash scripts/excord_cmd.sh {input.bam} {input.fasta} {output} {input.excord} &> {log}'
        

