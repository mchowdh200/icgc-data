#!/bin/env bash

snakemake -s manta-tumour-normal.smk \
          --configfile conf/config.yaml \
          -j $(grep -c ^processor /proc/cpuinfo) \
          --use-conda --conda-frontend mamba
