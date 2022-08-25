#!/usr/bin/env python

import pandas as pd

# alignment
from Bio import AlignIO

# handle arguments
import argparse
import os

parser = argparse.ArgumentParser(description="BayesTraits script")

parser.add_argument("tips_file", help="path towards fasta alignment", type=str)

parser.add_argument(
    "pos_mut", help="path towards positions and mutations to test", type=str
)

parser.add_argument(
    "pheno_file",
    help="path towards list of species with convergent phenotype",
    type=str,
)

args = parser.parse_args()

tips_file = args.tips_file
positions = args.pos_mut
pheno_file = args.pheno_file


path_name = os.path.dirname(os.path.abspath(positions))

Align_tips = AlignIO.read(open(tips_file), "fasta")
Tips_df = pd.DataFrame(
    [list(rec.seq) for rec in Align_tips], index=[rec.id for rec in Align_tips]
)

pos_mut = []
with open(positions) as f:
    for line in f:
        pos_mut.append(line.strip())

binary_df = pd.DataFrame(index=Tips_df.index, columns=[i for i in pos_mut])
for i in pos_mut:
    pos = int(i[:-1])
    mut = i[-1]
    binary_df[i] = [1 if a == mut else 0 for a in Tips_df[pos - 1].values]


phenotype = []
with open(pheno_file) as f:
    for line in f:
        phenotype.append(line.strip())

binary_df["phenotype"] = [1 if i in phenotype else 0 for i in Tips_df.index]

binary_df.to_csv(path_name + "/binary_tested_sites.tsv", sep="\t")

