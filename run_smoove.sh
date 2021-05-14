#!/bin/env bash

snakemake -s rules/smoove.smk \
          --configfile rules/conf/config.yaml \
          --resources disk_mb=$(df -m | grep /mnt/local | awk '{print $4*0.95}') \
          --scheduler greedy \
          --use-conda --conda-frontend mamba \
          -j $(grep -c ^processor /proc/cpuinfo)
