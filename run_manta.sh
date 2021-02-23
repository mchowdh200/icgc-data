#!/bin/env bash

snakemake -s rules/manta-tumour-normal.smk \
          --configfile conf/manta-tumour-normal.yaml \
          -j 4 -n
