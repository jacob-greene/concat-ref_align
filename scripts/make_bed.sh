#!/bin/bash

ml BEDTools

samtools view -B -h "${1}" \
| bedtools bamtobed -bedpe -i stdin \
| cut -f1,2,6,7 | sort -k1,1 -k2n,2n -k3n,3n \
| awk -v OFS='\t' '{len = $3 - $2 ; print $0, len }' > ${2}