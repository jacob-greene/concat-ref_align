#!/bin/bash

samtools view "${1}" \
    | perl -nle '@reads=split(/\t/,$_); {$start=$reads[3]-1; $end=$reads[3]+length($reads[9])-1; print "$reads[2]\t","$start\t","$end\t",$reads[0];}' \
    | awk '{gsub(":","\t",$4);gsub("_","\t",$4); print $1"\t"$2"\t"$3"\t"$4"\t"$5}' \
    | awk '{print $1"\t"$2"\t"$3"\t"$9"_"$10"_"$11"_"$12}' | sort -k1,1 -k2,2n -k4,4 \
    | uniq -c  | awk '{print $2"\t"$3"\t"$4"\t"$5"\t"$1}' | bgzip > "${2}"