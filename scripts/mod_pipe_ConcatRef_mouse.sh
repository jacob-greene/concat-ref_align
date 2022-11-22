#! /bin/bash
#$ -S /bin/bash
#$ -cwd


AlignedPath=$1 #path/to/bam_file/

INPUT=$2 #bam_file.bam #aligned to ConcatRef

human_tag=$3 #human_tag="UCSC_hg38"

ID=${INPUT%%.bam*}


# Extract mm10 aligned reads 

#this command gets the input chromosomes from idxstats, cuts out the other columns, removes anything with _NCBI_GRCm38 or * 
chr_list=$(samtools idxstats  $AlignedPath$INPUT | cut -f 1 | grep -v $human_tag | grep -vw \* )

echo $chr_list

#AL note: this next line just takes all reads that aligned with the MOUSE genome and puts them in a new file
#however it caused problems because there are pairs where one read is aligned to a normal chromosome and one is aligned to a 'junk chromosome' this lead to only a single read from the pair getting retained
#to fix this, I altered the chrom list
samtools view -b $AlignedPath$INPUT $chr_list > $AlignedPath$ID.mouse.bam #AL modified filenames to remove .hg19 and add .human_unfiltered


echo 'human reads removed'

#AL note: next step is to bai index the file
samtools index $AlignedPath$ID.mouse.bam $AlignedPath$ID.mouse.bam.bai

echo 'index created'


# Remove remained RNEXT mouse reads using tag information

#AL notes:
#index($2,"_NCBI_GRCm38")==0 is intended to remove the mouse contigs from the header #i guess awk is 1 indexed
#index($7,"_NCBI_GRCm38")==0 removes reads where the read pair is aligned to the mouse genome. This column (RNEXT) has an = for pairs on the same contig
samtools view -h $AlignedPath$ID.mouse.bam | awk -F '\t' -vTAG=${human_tag}  '{ if (index($2,TAG)==0 && index($7,TAG)==0) print $0 }' | samtools view -bS - > $AlignedPath$ID.mouse.cleaned.bam

echo 'reads with human-aligned mates removed'
samtools index $AlignedPath$ID.mouse.cleaned.bam $AlignedPath$ID.mouse.cleaned.bam.bai

echo 'index created'

