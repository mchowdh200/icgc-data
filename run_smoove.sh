#!/bin/env bash

snakemake -s rules/smoove.smk \
          --configfile rules/conf/config.yaml \
          --use-conda --conda-frontend mamba \
          -j $(grep -c ^processor /proc/cpuinfo)
