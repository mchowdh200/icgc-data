#!/bin/env bash

input=$1
output=$2

bedtools groupby -i $input -g 1-7,28 -c 27 -o distinct | grep DEL |
    awk '{ split($9,a,","); if (length(a) == 2) print $0; }'> $output
