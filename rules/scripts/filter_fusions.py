import sys

candidate_bedpe = sys.argv[1]
genes_bed = sys.argv[2]

### Get the gene -> strand mapping
gene2strand = {
    (x:=line.rstrip().split())[3] : x [4]
    for line in open(genes_bed).readlines()
}


with open(candidate_bedpe) as f:
    for line in f:
        A = line.rstrip().split()
        svtype = A[6]
        orientation = A[7], A[8]
        genes = A[9].split(',')
        # print(f'{svtype=} {orientation=} {genes=}')
        gene_strands = [gene2strand[x] for x in genes]

        if svtype == 'DUP':
            continue
        if svtype == 'DEL':
            if gene_strands[0] != gene_strands[1]: # must be same strand
                continue
        elif svtype == 'BND':
            if gene_strands[0] != gene_strands[1]:
                if orientation[0] != orientation[1]: # needs to be ++ or --
                    continue
            else:
                if orientation[0] == orientation[1]: # needs to be +- or -+
                    continue

        print(line, end='')
