#!/bin/env bash

snakemake -s merge-vcfs.smk \
          --configfile conf/config.yaml \
          -j $(grep -c ^processor /proc/cpuinfo) \
          --use-conda --conda-frontend mamba
          
