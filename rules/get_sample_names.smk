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
        f'{outdir}/sample_names.txt'

rule GetBamHeader:
    # for each file id, use score-client view command
    output:
        f'{outdir}/bam_headers/{{file_id}}_bam_header.txt'
    run:
        Path(f'{outdir}/bam_headers').mkdir(exist_ok=True)
        object_id = manifest_table[
            manifest_table.file_id == wildcards.file_id].object_id.values[0]
        shell(f'score-client view --object-id {object_id} --header-only > {output[0]}')

rule GetSampleName:
    # parse bam header and get the sample name from the RG info
    input:
        f'{outdir}/bam_headers/{{file_id}}_bam_header.txt'
    output:
        f'{outdir}/{{file_id}}_sample_name.txt'
    shell:
        f"grep RG {{input}} | head -1 | tr '\t' '\n' | grep SM | cut -d':' -f2 > {{output}}"

rule Combine:
    """
    Write a file with tab delimited file id, sample name mappings.
    """
    input:
        expand(f'{outdir}/{{file_id}}_sample_name.txt', file_id=file_ids)
    output:
        f'{outdir}/sample_names.txt'
    run:
        with open(output[0], 'w') as out:
            for file in input:
                file_id = file.split('_')[0]
                with open(file, 'r') as f:
                    sample_name = f.readline().rstrip()
                    out.write(f'{file_id}\t{sample_name}\n')
        
