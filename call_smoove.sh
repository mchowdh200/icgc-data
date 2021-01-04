#!/bin/bash
set -u

function setup() {
    local data_dir=$1

    # get the ref genome
    aws s3 cp s3://layerlabcu/ref/genomes/hs37d5/hs37d5.fa $data_dir/ref/
    aws s3 cp s3://layerlabcu/ref/genomes/hs37d5/hs37d5.fa.fai $data_dir/ref/

    # get exclude bed
    aws s3 cp s3://layerlabcu/BED/ceph18.b37.lumpy.exclude.2014-01-15.bed $data_dir/BED/
}

function download_and_call() {
    
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
    score-client/bin/score-client download \
        --verify-connection false \
        --validate false \
        --output-dir $data_dir/$donor_id/$specimen_type \
        --manifest $data_dir/$donor_id/$specimen_type/manifest.tsv \

    ## Call SVs
    smoove call \
        --name $donor_id.$specimen_type \
        --fasta $data_dir/ref/hs37d5.fa \
        --exclude $data_dir/BED/ceph18.b37.lumpy.exclude.2014-01-15.bed \
        --outdir $data_dir/$donor_id/$specimen_type \
        --genotype \
        --duphold \
        $data_dir/$donor_id/$specimen_type/$file_name

    ## Upload results VCF
    aws s3 cp $data_dir/$donor_id/$specimen_type/$donor_id.$specimen_type-smoove.genotyped.vcf.gz \
        s3://layerlabcu/icgc/smoove/$donor_id/$specimen_type/

    ## Cleanup working directory
    rm -r $data_dir/$donor_id/$specimen_type

}
export -f download_and_call

### Run everything
local manifest=$1

# get the second line from the manifest for testing
# local manifest_line=$(head -2 $manifest | tail -1)

setup $data_dir

cat $manifest | gargs --dry-run -p 4 'download_and_call /mnt/local/data "{}" object_donor_specimen.tsv'
