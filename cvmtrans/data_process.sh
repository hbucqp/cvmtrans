#!/usr/bin/bash

# concatenate paried fastq files

for file in *_1.clean.fq
do
    base=$(basename $file _1.clean.fq)
    cat ${base}_1.clean.fq ${base}_2.clean.fq > ${base}.fq
done

# change fastq file header

for file in *.fq
do
    base=$(basename $file .fq)
    awk '{print (NR%4 == 1) ? "@""'"${base}_"'"++i : $0}' ${base}.fq > ${base}.fastq
done

# compress fasta files

gzip *.fastq



# extract reads contains the transposon tags

python3 cvmtrans_extract_reads.py -i T3.fastq.gz -o test.fq.gz -tag "CCGGGGACTTATCAGCCAACCTGT"


