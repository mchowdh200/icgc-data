"""
After merging sv callsets from multiple sv callers,
Some samples have ref mismatches, this is a hacky way of
getting around that before squaring off the vcfs
"""
import sys
import pysam

vcf = pysam.VariantFile(sys.argv[1], 'r')
outdir = sys.argv[2]

print(vcf.header, end='')

for variant in vcf:
    variant.ref = '.'
    print(variant, end='')
