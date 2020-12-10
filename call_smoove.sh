#!/bin/bash
set -eu

function call_svs() {
    object_id=$(sed 's/-//g' <<<"$1")
    file_name=$2
    donor_id=$3
    specimen_type=$(sed 's/-//g' <<<"$4" | sed 's/\s+/\./g')
    fasta=$5
    exclude_bed=$6
    bam=/mnt/icgc/$object_id/$file_name

    smoove call \
        --name $donor_id.$specimen_type \
        --fasta $fasta \
        --exclude $exclude_bed \
        --outdir $donor_id.$specimen_type \
        --genotype \
        --duphold
}
export -f call_svs

# get ref genome
aws s3 cp s3://layerlabcu/ref/genomes/hs37d5/hs37d5.fa /mnt/local/data/ref/
aws s3 cp s3://layerlabcu/ref/genomes/hs37d5/hs37d5.fa.fai /mnt/local/data/ref/
fasta=/mnt/local/data/ref/hs37d5.fa

# get exclude bed
aws s3 cp s3://layerlabcu/BED/ceph18.b37.lumpy.exclude.2014-01-15.bed /mnt/local/data/BED
exclude_bed=/mnt/local/data/BED/ceph18.b37.lumpy.exclude.2014-01-15.bed

# mount icgc bams using manifest
sudo mkdir /mnt/icgc
sudo chmod 777 /mnt/icgc
score-client mount \
    --mount-point /mnt/icgc \
    --manifest prostate_cancer_manifest.tsv \
    --cache-metadata \
    --daemonize

# call svs with smoove -----
# format of data list:
# file_id,object_id,file_name,donor_id,specimen_type
# {0}    ,{1}      ,{2}      ,{3}     ,{4}
tail -n+2 object_donor_speciment.csv |
    gargs -s ',' -p 1 "call_svs {1} {2} {3} {4} $fasta $exclude_bed"



