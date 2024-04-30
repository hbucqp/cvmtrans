#!/usr/bin/python3

from Bio import SeqIO

count = SeqIO.convert("05ZYH33.gbk", "genbank", "05ZYH33.fa", "fasta")
print("Converted %i records" % count)
