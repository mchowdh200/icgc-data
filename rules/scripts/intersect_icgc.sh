#!/bin/env bash

set -exuo pipefail

input_bed=$1
icgc_bed=$2
output=$3

bedtools intersect -a $input_bed -b $icgc_bed -f 0.9 -r -wb | sort | uniq > $output
