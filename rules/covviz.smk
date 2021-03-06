import re
import shutil
from pathlib import Path
import pandas as pd

### Setup
######################################################################

# get variables from config file
outdir = config['outdir']
manifest_table = pd.read_csv(config['manifest'], sep='\t')
file_ids = manifest_table['file_id'].tolist()
donor_table = pd.read_csv(config['donor_table'], sep='\t')
donors = donor_table['ICGC Donor'].tolist()
get_from_s3 = config['get_from_s3']


### helper functions
######################################################################
def bam_disk_mb(wildcards):
    return int(manifest_table[
        manifest_table.file_id == wildcards.file_id].file_size * 1e-6)

### Rules
######################################################################
rule all:
    input:
        f'{outdir}/covviz_report.html'

if get_from_s3:
    rule S3GetIndex:
        output:
            f'{outdir}/indices/{{file_id}}.bai'
        shell:
            f"""
            mkdir -p {outdir}/indices
            aws s3 cp s3://layerlabcu/icgc/indices.05.09.21/{{wildcards.file_id}}.bai {{output}}
            """
else:
    rule GetManifest:
        """For each donor, create manifest to download files"""
        output:
            f'{outdir}/manifests/{{file_id}}-manifest.tsv'
        run:
            manifest = pd.read_csv(config['manifest'], sep='\t')
            manifest = manifest[manifest['file_id'] == wildcards.file_id]
            Path(f'{outdir}/manifests').mkdir(exist_ok=True)
            manifest.to_csv(output[0], index=False, sep='\t')


    rule ScoreClientGetIndex:
        input:
            manifest = f'{outdir}/manifests/{{file_id}}-manifest.tsv'
        output:
            bai = f"{outdir}/indices/{{file_id}}.bai"
        resources:
            disk_mb = bam_disk_mb
        run:
            Path(f'{outdir}/indices').mkdir(exist_ok=True)
            shell(f"""score-client download \\
            --validate false \\
            --output-dir {outdir}/indices \\
            --manifest {input.manifest}""")

            # remove bam and rename index with file_id
            bam = f'{outdir}/indices/{manifest_table[manifest_table.file_id == wildcards.file_id].file_name.values[0]}'
            bai = f'{bam}.bai'
            Path(bam).unlink()
            Path(bai).rename(f'{outdir}/indices/{wildcards.file_id}.bai')

rule GetRefIndex:
    output:
        f'{outdir}/ref/hs37d5.fa.fai'
    shell:
        """
        aws s3 cp s3://layerlabcu/ref/genomes/hs37d5/hs37d5.fa.fai {output}
        """

rule GoleftIndexcov:
    """
    Get coverage data from indices.  First rename indices to have metadata in them.
    """
    input:
        indices = expand(f'{outdir}/indices/{{file_id}}.bai', file_id=file_ids),
        fai = f'{outdir}/ref/hs37d5.fa.fai'
    output:
        temp(f'{outdir}/goleft/goleft-indexcov.bed.gz'),
        temp(f'{outdir}/goleft/goleft-indexcov.ped')
    run:
        # get dict of to-be-renamed bai files indexed by file id
        new_filenames = dict()
        for _, row in donor_table.iterrows():
            donor = row['ICGC Donor']
            fileid = row['File ID']
            specimen_type = re.compile('(\s|-)+').sub('.', row['Specimen Type'])
            new_filenames[fileid] = '_'.join([donor, fileid, specimen_type]) + '.bai'

        Path(f'{outdir}/renamed_indices').mkdir(exist_ok=True)
        for index in input.indices:
            # rename using the dict
            fileid = Path(index).stem
            shutil.copy(
                index, f'{outdir}/renamed_indices/{new_filenames[fileid]}')

        shell(f"""goleft indexcov --directory {outdir}/goleft \\
        --fai {{input.fai}} {outdir}/renamed_indices/*.bai""")

rule Covviz:
    input:
        bed = f'{outdir}/goleft/goleft-indexcov.bed.gz',
        ped = f'{outdir}/goleft/goleft-indexcov.ped'
        # svtyper_variants=f'{outdir}/annotations/squared.sites.vcf.gz',
        # LNCaP_variants=f'{outdir}/annotations/LNCAPEXP_REFINEFINAL1.vcf',
        # genes_hg19=f'{outdir}/annotations/genes_hg19.bed'
    output:
        f'{outdir}/covviz_report.html'
    shell:
        f"""
        covviz --output {{output}} --ped {{input.ped}} {{input.bed}}
        aws s3 cp {{output}} s3://icgc-vis/
        """

# rule GetAnnotationRegions:
#     output:
#         svtyper_variants=f'{outdir}/annotations/squared.sites.vcf.gz',
#         LNCaP_variants=f'{outdir}/annotations/LNCAPEXP_REFINEFINAL1.vcf',
#         genes_hg19=f'{outdir}/annotations/genes_hg19.bed'
#     shell:
#         f"""
#         aws s3 cp --recursive s3://layerlabcu/icgc/misc_annotations/ {outdir}/annotations/
#         aws s3 cp s3://layerlabcu/icgc/svtyper/squared.sites.vcf.gz {outdir}/annotations/
#         """
    
