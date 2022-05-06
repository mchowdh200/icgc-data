#!/bin/env bash

out_normal=$1
out_tumor=$2

donor_table="../data_listings/donor_table.tsv"
fid2type=$(tail -n+2 $donor_table | cut -f2,7)
echo "$fid2type" | wc -l
echo "$fid2type" | grep -i normal > $out_normal
echo "$fid2type" | grep -i tumour > $out_tumor
wc -l $out_normal
wc -l $out_tumor

