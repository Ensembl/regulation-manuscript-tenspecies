#!/bin/bash
if [[ $# -ne 3 ]]; then
  echo "Usage: $0 <species> <input.bed> <output.bed>"
  exit 1
fi
SPECIES=$1
INPUT_BED=$2
OUTPUT_BED=$3

patterns=$(cut -f1 supporting/"$SPECIES".chrom.sizes | sed 's/^/^/; s/$/[[:space:]]/' | paste -sd '|' -)
grep -E "$patterns" "$INPUT_BED" > processing/ctcf_canonical.bed

bedtools slop -b 3 -i processing/ctcf_canonical.bed -g supporting/$SPECIES.chrom.sizes > processing/slop.bed
bedtools merge -s -d 15 -c 4,5,6 -o collapse,sum,distinct -i processing/slop.bed > processing/merged_1.bed 
bedtools merge -c 4,5,6 -o collapse,sum,distinct -i processing/merged_1.bed > processing/merged_2.bed
sed 's/\+,-/\./' processing/merged_2.bed > "$OUTPUT_BED"