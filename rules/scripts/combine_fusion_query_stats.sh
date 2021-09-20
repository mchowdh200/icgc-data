#!/bin/env bash

pca_support=$1
onekg_support=$2

# sort/get unique lines form both inputs
# sort by genomic position afterwards
pca_unique=$(sort $pca_support | uniq | bedtools sort)
onekg_unique=$(sort $onekg_support | uniq | bedtools sort |
             awk '{print $10}') # we only need the support number

printf "#CHROM_A\tSTART_A\tEND_A\tCHROM_B\tSTART_B\tEND_B\tGENE_A\tGENE_B\tSVTYPE\tPCA_index_num_reads\tonekg_index_num_samples\t"
paste <(echo "$pca_unique") <(echo "$onekg_unique")
