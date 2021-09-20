#!/bin/env/bash

function stix_get_support
{
    while (( "$#" )); do
        case "$1" in
            -c|--chr)
                chr=$2
                shift 2;;
            -s|--start)
                start=$2
                shift 2;;
            -e|--end)
                end=$2
                shift 2;;
            -g|--gene)
                gene=$2
                shift 2;;
            -t|--svtype)
                svtype=$2
                shift 2;;
            --) # end argument parsing
                shift
                break;;
            -*|--*=) # unsupported flags
                echo "Error: Unsupported flag $1" >&2
                exit 1;;
        esac
    done
    # this part hurts me
    cd /mnt/local

    ## TODO separate summation from the stix query
    # we want to count the number of samples that support a query
    # ie count the number of lines in the result first. Then sum it.
    # yep its all hard coded...
    samples=$(/mnt/local/bin/stix \
         -i tumour_index \
         -d pca.ped.db \
         -t $svtype -s 500 \
         -l "$chr:$start-$end" \
         -r "$chr:$start-$end" |
         tail -n+3 | # skip first two lines
         cut -f4,5)
    n_samples=$(echo "$samples" | awk '{if ($1+$2 > 0) print}' | wc -l)
    support=$(echo "$samples" |
              python3 -c 'import sys; print(sum(int(i) for line in sys.stdin
                                            for i in line.rstrip().split()))')
    printf "$chr\t$start\t$end\t$gene\t$svtype\t$n_samples\t$support\n"
}
export -f stix_get_support

set -eu

gene_list=$1
output=$2
threads=$3

cat $gene_list |
    gargs -p $threads \
          "stix_get_support --chr {0} --start {1} --end {2} \\
                            --gene {3} --svtype {4}" > $output

