from pathlib import Path
import pandas as pd

### Setup
################################################################################
outdir = '/mnt/local/data'
manifest_table = pd.read_csv(
    '/home/ubuntu/icgc-data/data_listings/donor_manifest.tsv', sep='\t')
file_ids = manifest_table.file_id.tolist()


### Rules
################################################################################
rule All:
    input:
        expand(f'{outdir}/bam_headers/{{file_id}}_bam_header.txt',
               file_id=file_ids)

rule GetBamHeader:
    # for each file id, use score-client view command
    output:
        f'{outdir}/bam_headers/{{file_id}}_bam_header.txt'
    run:
        Path(f'{outdir}/bam_headers').mkdir(exist_ok=True)
        object_id = str(manifest_table[
            manifest_table.file_id == wildcards.file_id].object_id)
        shell(f'score-client view --object-id {object_id} --header-only > {output[0]}')

