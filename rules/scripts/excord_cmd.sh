#!/bin/env bash
set -euo pipefail
bam=$1
fasta=$2
output=$3
excord=$4

# get insert-sizes of first 100k reads, to compute discordant distance
samtools view -f66 $bam | head -100000 |
    awk '{print sqrt($9^2)}' \
    > $bam.insert-sizes.txt
disc_dist=$(python3 scripts/get_disc_distance.py $bam.insert-sizes.txt)

echo "running excord with discordant distance = $disc_dist"
samtools view -b $bam |
    $excord --discordantdistance $disc_dist --fasta $fasta /dev/stdin |
    bgzip -c > $output
