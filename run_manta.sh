#!/bin/env bash

snakemake -s rules/manta-tumour-normal.smk \
          --configfile conf/manta-tumour-normal.yaml \
          --resources bams_on_disk=2 manta_running=1 \
          -j 4
