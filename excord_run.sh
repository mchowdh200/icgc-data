#!/bin/bash
set -eu

function setup() {
    local data_dir=$1

    # get the ref genome
    aws s3 cp s3://layerlabcu/ref/genomes/hs37d5/hs37d5.fa $data_dir/ref/
    aws s3 cp s3://layerlabcu/ref/genomes/hs37d5/hs37d5.fa.fai $data_dir/ref/

    # setup excord
    wget -O /mnt/local/excord https://github.com/brentp/excord/releases/download/v0.2.4/excord
    chmod +x /mnt/local/excord

}

function download_and_run() {
    
    local data_dir=$1
    local manifest_line=$2
    local donor_table=$3

    ## get file_name from manifest
    local file_name=$(cut -f5 <<<"$manifest_line")

    ## match to donor_id/specimen_type
    local donor_id=$(grep "$file_name" "$donor_table" | grep -v mini | cut -f4)

    # remove hyphen and spaces, sep words with dots
    local specimen_type=$(grep "$file_name" "$donor_table" | grep -v mini |
                          cut -f5 | sed 's/-//g' | sed -E 's/\s+/\./g')

    ## create mini-manifest
    # first create working directory based on the donor_id/specimen_type
    # then place the mini manifest in the working directory
    [[ ! -d $data_dir/$donor_id ]] &&
        mkdir $data_dir/$donor_id
    [[ ! -d $data_dir/$donor_id/$specimen_type ]] &&
        mkdir $data_dir/$donor_id/$specimen_type
    cat manifest_header.tsv <(echo "$manifest_line") \
        > $data_dir/$donor_id/$specimen_type/manifest.tsv

    ## Download BAM
    score-client/bin/score-client download \
        --verify-connection false \
        --validate false \
        --output-dir $data_dir/$donor_id/$specimen_type \
        --manifest $data_dir/$donor_id/$specimen_type/manifest.tsv

    ## Run excord on bam
    samtools view -b $data_dir/$donor_id/$specimen_type/$file_name |
    /mnt/local/excord  \
        --discordantdistance 500 \
        --fasta $data_dir/ref/hs37d5.fa \
        /dev/stdin |
    bgzip -c > $data_dir/$donor_id/$specimen_type/$donor_id.$specimen_type.excord.bed.gz

    ## Upload results
    aws s3 cp $data_dir/$donor_id/$specimen_type/$donor_id.$specimen_type.excord.bed.gz \
              s3://layerlabcu/icgc/excord/$donor_id/$specimen_type/

    ## Cleanup working directory
    rm -r $data_dir/$donor_id/$specimen_type
}
export -f download_and_run

### Run everything
manifest=$1
setup /mnt/local/data
cat $manifest | gargs  -p 5 'download_and_run /mnt/local/data "{}" object_donor_specimen.tsv'

