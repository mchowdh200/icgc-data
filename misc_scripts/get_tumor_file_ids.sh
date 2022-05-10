#!/bin/env bash

set -eu

# pipe list of files that contain icgc file ids and separate into
# lists by if they are tumor or normal associated

indir=$1
outdir=$2

donor_table="../data_listings/donor_table.tsv"
fid2type=$(tail -n+2 $donor_table | cut -f2,7)
normal_ids=$(echo "$fid2type" | grep -i normal | cut -f1 | tr '\n' '|')
tumor_ids=$(echo "$fid2type" | grep -i tumour | cut -f1 | tr '\n' '|')

ls $indir | grep -i -E \"$normal_ids\" > $outdir/normals.txt
ls $indir | grep -i -E \"$tumor_ids\" > $outdir/tumors.txt

