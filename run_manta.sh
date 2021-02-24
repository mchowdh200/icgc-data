#!/bin/env bash

snakemake -s rules/manta-tumour-normal.smk \
          --configfile conf/manta-tumour-normal.yaml \
          --resources num_downloads=2 bams_on_disk=6 manta_running=1 \
          -j $(grep -c ^processor /proc/cpuinfo)
