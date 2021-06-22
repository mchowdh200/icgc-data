#!/bin/env bash
set euo pipefail
# genomic interval
c=$1
s=$2
e=$3
data_dir=$4
d="1kg.ped.db" # ped database
i="alt_sort_b" # index directory

# TODO seems to be a bug in stix where you need to call it
# from the directory with the index in it.
hit=$(cd $data_dir && stix -d 1kg.ped.db -t DEL -s 500 -i alt_sort_b -l $c:$s-$s -r $c:$e-$e |
      tail -n+2 | awk '{print $7+$8}' | paste -sd " " - )
nz=$(echo "$hit" | tr ' ' '\n' | awk '$1>0' | wc -l)
max=$(echo "$hit" | tr ' ' '\n' | sort -n | tail -n 1)
total=$( echo $hit | python scripts/sum.py )
echo -e "$c\t$s\t$e\t$nz\t$max\t$total"
