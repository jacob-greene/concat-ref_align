#!/usr/bin/env python3
# coding: utf-8

import argparse

# Initiate the parser
parser = argparse.ArgumentParser()
parser.add_argument('--humanRef', action="store", dest="human_ref_path", help="Path to human reference genome")
parser.add_argument('--mouseRef', action="store", dest="mouse_ref_path", help="Path to mouse reference genome")
parser.add_argument('--concatRef', action="store", dest="out_file", help="Output concated reference file name")
parser.add_argument('--tag', action="store", dest="tag", help="tag info to rename mouse contig")
parser.add_argument('--human_tag', action="store", dest="human_tag", help="tag info to rename human contig") #AL mod

args = parser.parse_args()

#define references
####################
# ##change the path to reference genome (human  & mouse), out_file and tag below
# human_ref_path = '/fh/fast/ha_g/grp/reference/GRCh38/10XG/refdata-GRCh38-2.1.0/fasta/genome.fa'
# mouse_ref_path = '/fh/fast/ha_g/grp/reference/Mus_musculus_UCSC_mm10/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa'
# out_file = '../fasta/GRCh38_plus_UCSC_mm10.fa'
# tag ='UCSC_mm10' #rename mouse contig
####################

#read in the files
with open(args.human_ref_path,'r') as f:
    human_ref = f.readlines()
    
with open(args.mouse_ref_path,'r') as f:
    mouse_ref = f.readlines()

#view the first few lines
print('# of lines in human:',len(human_ref))
print(human_ref[0:10])
print('\n')
print('# of lines in mouse:',len(mouse_ref))
print(mouse_ref[0:10])


#rename the mouse contigs
for i in range(len(mouse_ref)):
    if mouse_ref[i].startswith('>'):
        mouse_ref[i] = mouse_ref[i].strip('\n')+(f'_{args.tag}\n')

#rename the human contigs #AL added this section
for i in range(len(human_ref)):
    if human_ref[i].startswith('>'):
        human_ref[i] = human_ref[i].strip('\n')+(f'_{args.human_tag}\n')

#export to new file
with open (args.out_file, 'w+') as f:
    f.write(''.join(human_ref))
    f.write(''.join(mouse_ref))






