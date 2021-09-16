#!/bin/bash

function stix_get_support
{
    while (( "$#" )); do
        case "$1" in
            -ca|--chrA)
                chrA=$2
                shift 2;;
            -sa|--startA)
                startA=$2
                shift 2;;
            -ea|--endA)
                endA=$2
                shift 2;;
            -cb|--chrB)
                chrB=$2
                shift 2;;
            -sb|--startB)
                startB=$2
                shift 2;;
            -eb|--endB)
                endB=$2
                shift 2;;
            -sv|--svtype)
                svtype=$2
                shift 2;;
            -ga|--geneA)
                geneA=$2
                shift 2;;
            -gb|--geneB)
                geneB=$2
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

    # yep its all hard coded...
    support=$(/mnt/local/bin/stix \
         -i tumour_index \
         -d pca.ped.db \
         -t $svtype -s 500 \
         -l "$chrA:$startA-$endA" \
         -r "$chrB:$startB-$endB" |
    tail -n+3 | # skip first two lines
    cut -f4,5 | # get the last two columns (integers)
    python3 -c 'import sys; print(sum(int(i) for line in sys.stdin
                                for i in line.rstrip().split()))')

    printf "$chrA\t$startA\t$endA\t$chrB\t$startB\t$endB\t$geneA\t$geneB\t$svtype\t$support\n"
}
export -f stix_get_support

set -eu

fusion_list=$1
output=$2
threads=$3

cat $fusion_list |
    gargs -p $threads \
          "stix_get_support --chrA {0} --startA {1} --endA {2} \\
                            --chrB {3} --startB {4} --endB {5} \\
                            --svtype {6} --geneA {7} --geneB {8}" > $output

