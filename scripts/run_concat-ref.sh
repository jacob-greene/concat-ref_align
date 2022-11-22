#!/bin/sh
export mm10=/shared/biodata/ngs/Reference/iGenomes/Mus_musculus/UCSC/mm10/Sequence/BWAIndex/genome.fa
export hg38=/shared/biodata/ngs/Reference/iGenomes/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/genome.fa
export OUT=/fh/fast/henikoff_s/user/jgreene/refs/hg38_plus_mm10_BWA/GRCh38_plus_UCSC_mm10.fa


#concat and index
./concatenate_reference.py --humanRef $hg38 --mouseRef $mm10 --concatRef $OUT --tag UCSC_mm10 --human_tag UCSC_hg38
bwa index $OUT
