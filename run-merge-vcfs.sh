#!/bin/env bash

snakemake -s rules/merge-vcfs.smk \
          --configfile rules/conf/config.yaml \
          -j $(grep -c ^processor /proc/cpuinfo) \
          --use-conda \
          --resources num_downloads=4
          
