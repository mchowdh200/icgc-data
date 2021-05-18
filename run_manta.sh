#!/bin/env bash

cores=$(grep -c ^processor /proc/cpuinfo)
num_downloads=$(python3 -c "print(int(0.75*$cores))")

snakemake -s rules/manta.smk \
          --configfile conf/manta-tumour-normal.yaml \
          --resources num_downloads=$num_downloads \
          -j $cores
