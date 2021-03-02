#!/bin/env bash

snakemake -n -s rules/merge-vcfs.smk \
          --configfile conf/merge-vcfs.yaml \
          -j $(grep -c ^processor /proc/cpuinfo)
          
