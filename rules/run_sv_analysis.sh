#!/bin/env bash

snakemake -s ~/Repositories/icgc-data/rules/sv_analysis.smk -j 8 \
          --use-conda --conda-frontend mamba \
          --directory /home/murad/Repositories/icgc-data/rules \
          --config outdir='/home/murad/data/icgc/sv_analysis' \
          donor_table='/home/murad/Repositories/icgc-data/data_listings/donor_table.tsv' \
          sample_names='/home/murad/Repositories/icgc-data/data_listings/sample_names.txt' \
          genes_bed='/home/murad/Repositories/icgc-data/data_listings/genes.hg19.collapsed.bed' \
          regions_8p16q='/home/murad/Repositories/icgc-data/data_listings/8p16q_regions.bed' \
          cytoband_bed='/home/murad/Repositories/icgc-data/data_listings/hg37_cytoband.bed' \
          somaticSV_bucket='s3://layerlabcu/icgc/manta-tumour-normal' \
          contig_lengths='/home/murad/Repositories/icgc-data/data_listings/hs37d5.contig_lengths.txt'
          
