#!/bin/env bash

snakemake -s sv_analysis.smk -j 8 \
          --use-conda --conda-frontend mamba \
          --config outdir='/home/murad/data/icgc/sv_analysis' \
          donor_table='/home/murad/Repositories/icgc-data/data_listings/donor_table.tsv' \
          sample_names='/home/murad/Repositories/icgc-data/data_listings/sample_names.txt' \
          genes_bed='/home/murad/Repositories/icgc-data/data_listings/genes_hg19.bed' \
          regions_8p16q='/home/murad/Repositories/icgc-data/data_listings/8p16q_regions.bed' \
          cytoband_bed='/home/murad/Repositories/icgc-data/data_listings/hg37_cytoband.bed' \
          somaticSV_bucket='s3://layerlabcu/icgc/manta-tumour-normal'
          
