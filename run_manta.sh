#!/bin/env bash

cores=$(grep -c ^processor /proc/cpuinfo)
snakemake -s rules/manta.smk \
          --configfile rules/conf/config.yaml \
          --resources num_downloads=$(python3 -c "print($cores//2)") \
          -j $cores
