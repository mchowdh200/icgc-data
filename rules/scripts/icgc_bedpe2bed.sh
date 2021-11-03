#!/bin/env bash
set -exuo pipefail

input=$1
output=$2
svtype=$3

if [[ $svtype=='INV' ]]; then
    zcat $input | tail -n+2 |
        awk '{OFS="\t"; if ($11=="h2hINV" || $11=="t2tINV") print $1,$2,$6,$11;}' > $output
else
    zcat $input | tail -n+2 |
        awk -v svtype=$svtype '{OFS="\t"; if ($11==svtype) print $1,$2,$6,$11;}' > $output
fi
