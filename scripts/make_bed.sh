#!/bin/bash

samtools view "${1}" \
    | perl -nle '@reads=split(/\t/,$_); {$start=$reads[3]-1; $end=$reads[3]+length($reads[9])-1; print "$reads[2]\t","$start\t","$end\t",$reads[0];}' \
    | bgzip > "${2}"