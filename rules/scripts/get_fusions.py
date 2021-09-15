import sys

genes_bed = sys.argv[1]
fusion_bedpe_list = sys.argv[2:]

### get the loci of all genes
genes = dict()
with open(genes_bed) as f:
    for line in f:
        A = line.rstrip().split()
        g = A[3]
        c, s, e = A[0], A[1], A[2]
        genes[g] = [c, s, e]

### print gene A, gene B in bedpe format along with svtype
for bedpe in fusion_bedpe_list:
    with open(bedpe) as f:
        for line in f:
            A = line.rstrip().split()
            ga, gb = A[9].split(',')
            ca, sa, ea = genes[ga]
            cb, sb, eb = genes[gb]
            sv = A[6]

            print('\t'.join([
                ca, sa, ea,
                cb, sb, eb,
                sv, ga, gb
            ]))
