#!/bin/env bash

snakemake -s sv_analysis.smk -j 8 \
          --use-conda --conda-frontend mamba \
          --config outdir=~/data/icgc/sv_analysis/ \
          donor_table=~/Repositories/icgc-data/data_listings/donor_table.tsv \
          genes_bed='~/Repositories/icgc-data/data_listings/genes_hg19.bed' \
          regions_18p16q='~/Repositories/icgc-data/data_listings/8p16q_regions.bed' \
          somaticSV_bucket='s3://layerlabcu/icgc/manta-tumour-normal'
          
