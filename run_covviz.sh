#!/bin/env bash

snakemake -s rules/covviz.smk \
          --configfile conf/covviz.yaml \
          -j $(grep -c ^processor /proc/cpuinfo)
