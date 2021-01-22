#!/bin/bash
set -eu

manifest=$1

# get reference
aws s3 cp s3://layerlabcu/ref/genomes/hs37d5/hs37d5.fa /mnt/local/data/ref/
aws s3 cp s3://layerlabcu/ref/genomes/hs37d5/hs37d5.fa.fai /mnt/local/data/ref/

# setup excord
wget -O /mnt/local/excord https://github.com/brentp/excord/releases/download/v0.2.4/excord
chmod +x /mnt/local/excord

