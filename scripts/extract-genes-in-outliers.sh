#!/usr/bin/env bash

# Simple script to parse gene identifiers within curated regions

# $1 = TSV file with coordinates to extract genes between
# $2 = GFF file to extract genes from
# $3 = Output file

while IFS=$'\t' read -r -a lineArray; do
    # Arrage: 0 = chromosome, 1 = start, 2 = end, 3 = classifcation
    ID=(
        "$(
            gawk \
                -v chr="${lineArray[0]}" -v start="${lineArray[1]}" -v end="${lineArray[2]}" \
                '$1 == chr && $3 == "gene" && $4 >= start && $5 <= end {print $9}' "$2" |
                grep -Eo "FUN_.*;" | sed 's/;.*//'
        )"
    )

    # Loop over IDs for each region
    for i in ${ID[@]}; do
        echo -e "${lineArray[3]}\t${lineArray[0]}\t${lineArray[1]}\t${lineArray[2]}\t${i}" >>"$3"
    done
done <"$1"
