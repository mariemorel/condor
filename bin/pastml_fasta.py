#!/usr/bin/env python

import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from iteround import saferound
import numpy as np

import argparse

parser = argparse.ArgumentParser(description="pastmlfile into fasta")
parser.add_argument("pastml_acr", help="pastml output")
parser.add_argument("positions", help = "file with positions to test" )
parser.add_argument("length", help = "nb sites" )
parser.add_argument("marginal", help = "marginal probabilities of root")
parser.add_argument("nb_simu", help = "number of simulations to perform")
parser.add_argument("output", help="name of the output")
args = parser.parse_args()

pastml_acr = str(args.pastml_acr)
path_pos = str(args.positions)
nb_sites = int(args.length)
output = str(args.output)
nb_simu = int(args.nb_simu)

Pos_list = []
with open(path_pos) as f:
    for line in f:
        Pos_list.append(int(line.strip()))

#list with all positions to test

characters = []
marginal_list = []
with open(args.marginal) as f:
    for line in f:
        if line.strip().isdigit(): 
            characters.append(line.strip())
        else:
            marginal_list.append(line.strip().split("\t")[1:])

##list of marginal proba for root
##list of different amino acids given at root


AA = list("ARNDCQEGHILKMFPSTWYV")

root_dict = {}
for a,pos in enumerate(characters): 
    root= []

    probs = []
    for proba in marginal_list[a]:
        probs.append(nb_simu*float(proba))
    for i, count in enumerate(saferound(probs,0)):
        root.extend( int(count) * [AA[i]])

    root_dict[pos] = root

##dictionary of each pos with a vector of amino acids corresponding at the root amino acids mutliplied by the marginal proba * nb simulations

#A = pd.read_csv(pastml_acr, sep='\t', index_col = 0)
A = pd.read_csv(pastml_acr, sep='\t', index_col = 0, compression = "gzip")
B = A[[str(i-1) for i in Pos_list]]

## adjust position since different counting

#create fasta sequences of acr output.
record = []
for n in set(B.index):
    record.append(SeqRecord(Seq("".join(list(B.loc[n]))), id = n, name = n, description = n ))
    if n == "root":
        reconstructed_root = "".join(list(B.loc[n]))
    elif n == "N000000001":
        reconstructed_root = "".join(list(B.loc[n]))
#retrieve reconstructed root

SeqIO.write(record, output+"pastml_acr.fasta", "fasta")
with open("reconstructed_root", "w") as wf:
    wf.write(reconstructed_root)


with open(output+"marginal_posterior.txt", "w") as wf:
    for i in root_dict:
        wf.write("\t".join([i, "".join(root_dict[i])])+"\n")


