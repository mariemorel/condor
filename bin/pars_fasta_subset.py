#!/usr/bin/env python

import pandas as pd
from Bio import SeqIO
import numpy as np
from collections import Counter
import argparse
import os.path

parser = argparse.ArgumentParser(description="Parse a fasta file into a table understandable by pastml and select sites to simulate")
parser.add_argument("align", help="alignment in fasta")
parser.add_argument("length", help = "nb sites" )
parser.add_argument("min_seq", help = "nb sequences a mutation should be" )
parser.add_argument("output", help="name of the output")
parser.add_argument("positions", help="file with specific positions to test")
args = parser.parse_args()


align_path = str(args.align)
len_seq = int(args.length)
minseq = int(args.min_seq)
output = str(args.output)

#import alignment and store it in a dataframe

align=dict()
for rec in SeqIO.parse(align_path, 'fasta'):
    align[rec.name] = rec.seq[0:len_seq]

align_df = pd.DataFrame([list(i) for i in align.values()], index = align.keys())

all_content = set()
for i in align_df.values:
    all_content.update(i)

to_replace = [a for a in Counter(all_content).keys() if not a in list("ARNDCQEGHILKMFPSTWYV")]

for r in to_replace:
    align_df.replace(r, np.nan, inplace = True)

#replace unknown amino-acids by empty for pastml

AA = ['A','R','N','D','C','Q', 'E', 'G' ,'H' ,'I' ,'L' ,'K' ,'M' ,'F', 'P', 'S', 'T', 'W', 'Y', 'V']

#test positions with mutations present in at least minseq and it is not a constant position.
positions_to_test = []
if os.path.isfile(args.positions):
    with open(args.positions) as f:
        for line in f:
            positions_to_test.append(int(line.strip())) #numerotation at 0
else:
    for pos in range(len_seq):
        if len([Counter(align_df[pos])[i] for i in AA if Counter(align_df[pos])[i] >= minseq]) >= 2:
            positions_to_test.append(pos) #numerotation starts at 0

    
#create subset of alignment
subset_align = align_df[positions_to_test]

#pastml won't run if not all amino acids are present in the column. So we create fake lines in the table to have a representation of all amino acids.
B = pd.DataFrame([[i]*len(positions_to_test) for i in AA], index =  ["".join(['fake_', i]) for i in AA])
pd.concat([subset_align, B]).to_csv(output+"pastml_input.tsv.gz", sep='\t', compression = "gzip")

#write positions to test with numerotation starting at one (easier for user)
with open("positions_to_test.txt", "w") as wf:
    wf.write("\n".join([str(i+1) for i in positions_to_test])+"\n") ## numerotation starts at 1
