#!/bin/bash
set -eu

function setup() {
    local data_dir=$1
    local manifest=$2

    # get the header of the manifset
    head -1 $manifest > $data_dir/manifest_header.tsv

    # get the ref genome
    aws s3 cp s3://layerlabcu/ref/genomes/hs37d5/hs37d5.fa $data_dir/ref/
    aws s3 cp s3://layerlabcu/ref/genomes/hs37d5/hs37d5.fa.fai $data_dir/ref/

    # get exclude bed
    aws s3 cp s3://layerlabcu/BED/ceph18.b37.lumpy.exclude.2014-01-15.bed $data_dir/BED
}

function download_and_call() {
    ## Inputs:
    # 0. base $data_dir
    # 1. line of manifest to download with
    # 2. path to the object_donor_specimen.tsv
    #
    ## Procedure:
    # 1. DONE get the $file_name from the manifest line
    # 2. DONE match the $file_name to the $donor_id/$specimen_type
    # 3. DONE cat the manifest line to the header to create a "mini-manifest"
    # 4. DONE create output directory using $donor_id/$specimen_type
    # 5. DONE use the score client to download file using the mini manifest
    # 6. DONE call svs with smoove
    # 7. TODO upload results vcf to s3
    # 8. TODO clean up working directory
    
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
    cat $data_dir/manifest_header.tsv <(echo "$manifest_line") \
        > $data_dir/$donor_id/$specimen_type/manifest.tsv

    ## Download BAM
    score-client/bin/scoreclient download \
        --verify-connection false \
        --validate false \
        --output-dir $data_dir/$donor_id/$specimen_type \
        --manifest $data_dir/$donor_id/$specimen_type/manifest.tsv \

    ## Call SVs
        smoove call \
            --name $donor_id.$specimen_type \
            --fasta $data_dir/ref/hs37d5.fa \
            --exclude $data_dir/BED/ceph18.b37.lumpy.exclude.2014-01-15.bed \
            --outdir $donor_id/$specimen_type \
            --genotype \
            --duphold

    ## Upload results VCF TODO
    ## Cleanup working directory TODO

}

function testing() {
    local data_dir=~/Repositories/icgc-data
    local manifest=$data_dir/prostate_cancer_manifest.tsv
    local donor_table=$data_dir/object_donor_specimen.tsv

    # get the second line from the manifest for testing
    local manifest_line=$(head -2 $manifest | tail -1)

    # setup $data_dir $manifest

    download_and_call "$data_dir" "$manifest_line" "$donor_table"
}

### Run everything
testing
