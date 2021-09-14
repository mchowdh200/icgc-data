#!/bin/env bash

input=$1
output=$2

# BEDTOOLS STATEMENT
# group by the breakpoints
# operate on the gene column
# generate a comma separated list of distinct genes
# this will not filter by breakend orientation or gene strand
# 1-6 = breakpoints
# 11 = SVTYPE
# 9,10 = Breakend orientation
# 27 = gene
# AWK STATEMENT:
# $10 will be the column containing the gene comma separatedlist
bedtools groupby -i $input -g 1-6,11,9,10 -c 27 -o distinct |
    awk '{ split($10,a,","); if (length(a) == 2) print $0; }'> $output
