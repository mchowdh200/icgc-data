"""
Call SVs using smoove.
Label output vcfs with the icgc data portal file ID.
"""
from pathlib import Path
import pandas as pd

### Setup
###################################################
# TODO put in a separate python file as a function and import?
outdir = config['outdir']
manifest_table = pd.read_csv(config['manifest'], sep='\t')
file_ids = manifest_table['file_id'].tolist()
donor_table = pd.read_csv(config['donor_table'], sep='\t')
donors = donor_table['ICGC Donor'].tolist()

### Helper Functions
###################################################
# TODO put this in a utils python file
def bam_disk_mb(wildcards):
    return int(manifest_table[
        manifest_table.file_id == wildcards.file_id].file_size * 1e-6)

### Rules
###################################################

rule All:
    input:
        expand(f'{outdir}/{{file_id}}/{{file_id}}-smoove.genotyped.vcf.gz',
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

rule GetExcludeRegions:
    output:
        f'{outdir}/bed/exclude.bed'
    shell:
        'aws s3 cp s3://layerlabcu/BED/ceph18.b37.lumpy.exclude.2014-01-15.bed {output}'


# TODO put this in a separate snakefile and import?
# - how can i do this with a config.
# - maybe one config file for the whole project?
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
        disk_mb = bam_disk_mb
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


rule SmooveCall:
    input:
        fasta = rules.GetReference.output.fasta,
        fai = rules.GetReference.output.fai,
        exclude = rules.GetExcludeRegions.output,
        bam = rules.GetBam.output.bam,
        bai = rules.GetBam.output.bai
    output:
        f'{outdir}/{{file_id}}/{{file_id}}-smoove.genotyped.vcf.gz'
    resources:
        disk_mb = bam_disk_mb
    conda:
        'envs/smoove.yaml'
    shell:
        f"""
        mkdir -p {outdir}/{{wildcards.file_id}}
        smoove call -p 1 \\
            --name {{wildcards.file_id}} \\
            --fasta {{input.fasta}} \\
            --exclude {{input.exclude}} \\
            --outdir {outdir}/{{wildcards.file_id}} \\
            --genotype --duphold \\
            {{input.bam}}
        find {outdir}/{{wildcards.file_id}}/ -type f -not -name '*.vcf.gz*' -delete
        """
