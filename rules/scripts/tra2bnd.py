"""
SURVIVOR infers svtype of TRA from some BND regions.
SVTyper doesn't recognize TRA, so this script changes them back to BND.
"""
import sys
import pysam

vcf = pysam.VariantFile(sys.argv[1], 'r')
outdir = sys.argv[2]

print(vcf.header, end='')

for variant in vcf:
    if variant.info['SVTYPE'] == 'TRA':
        variant.info['SVTYPE'] == 'BND'
        print(variant, end='')
