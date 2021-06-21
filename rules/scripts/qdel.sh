#!/bin/bash

# genomic interval
c=$1
s=$2
e=$3

d="/data/stix/1kg/1kg.ped.db" # ped database
i="/data/stix/1kg/alt_sort_b" # index directory

hit=$(stix -d $d -t DEL -s 500 -i $i -l $c:$s-$s -r $c:$e-$e |
      tail -n+2 | awk '{print $7+$8}' | paste -sd " " - )
nz=$(echo "$hit" | tr ' ' '\n' | awk '$1>0' | wc -l)
max=$(echo "$hit" | tr ' ' '\n' | sort -n | tail -n 1)
total=$( echo $hit | python scripts/sum.py )
#echo -e "$nz\t$max\t$total\t$hit"
echo -e "$nz\t$max\t$total"
