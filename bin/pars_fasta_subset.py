#!/usr/bin/env python

import pandas as pd
from Bio import SeqIO
import numpy as np
from collections import Counter
import argparse

parser = argparse.ArgumentParser(description="Parse a fasta file into a table understandable by pastml and select sites to simulate")
parser.add_argument("align", help="alignment in fasta")
parser.add_argument("length", help = "nb sites" )
parser.add_argument("min_seq", help = "nb sequences a mutation should be" )
parser.add_argument("output", help="name of the output")
args = parser.parse_args()


align_path = str(args.align)
len_seq = int(args.length)
minseq = int(args.min_seq)
output = str(args.output)

align=dict()
for rec in SeqIO.parse(align_path, 'fasta'):
    align[rec.name] = rec.seq[0:len_seq]

align_df = pd.DataFrame([list(i) for i in align.values()], index = align.keys())
align_df.replace('X', np.nan, inplace = True)
align_df.replace('-', np.nan, inplace = True)
align_df.replace('*', np.nan, inplace = True)


AA = ['A','R','N','D','C','Q', 'E', 'G' ,'H' ,'I' ,'L' ,'K' ,'M' ,'F', 'P', 'S', 'T', 'W', 'Y', 'V']

positions_to_test = []
for pos in range(len_seq):
    if len([Counter(align_df[pos])[i] for i in AA if Counter(align_df[pos])[i] >= minseq]) >= 2:
        positions_to_test.append(pos) #numerotation starts at 0


subset_align = align_df[positions_to_test]
B = pd.DataFrame([[i]*len(positions_to_test) for i in AA], index =  ["".join(['fake_', i]) for i in AA])
pd.concat([subset_align, B]).to_csv(output+"pastml_input.tsv.gz", sep='\t', compression = "gzip")


with open("positions_to_test.txt", "w") as wf:
    wf.write("\n".join([str(i+1) for i in positions_to_test])+"\n") ## numerotation starts at 1
