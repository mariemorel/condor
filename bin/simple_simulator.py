#!/usr/bin/env python

import numpy as np
import pandas as pd
from scipy import linalg

import argparse

from ete3 import PhyloTree

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from scipy import sparse
import time


def parseRates(rates_file):
    rates = []

    with open(rates_file) as f:
        lines = (line.strip() for line in f if line)
        rates = [float(line) for line in lines]

    return rates

parser = argparse.ArgumentParser(
    description="Create an amino acid alignment corresponding to a tree and an ancestral sequence and count the number of substitutions from the root only if transmitted to leaves"
)
parser.add_argument(
    "root", help="root sequence"
) 
parser.add_argument(
    "tree_file", help="rooted tree with names at the internal nodes")
parser.add_argument("output", help="output file name")
parser.add_argument(
    "rates_file", help="file containing one site rate per line") 
parser.add_argument(
    "freq", help="file with the freqs amino acid from the model")
parser.add_argument("matrix_model", help="matrix of the model")
args = parser.parse_args()

tree = PhyloTree(args.tree_file, format=1)
amino_acid = list("ARNDCQEGHILKMFPSTWYV")
rates = parseRates(args.rates_file)

with open(args.root) as f:
    root_list =  f.read().splitlines()

root = root_list[0] #probabilistic root at one position
nb_residues = int(len(root)) #should equal to nb simulations

frequencies = []
with open(args.freq) as f:
    for line in f:
        frequencies.append(float(line.split("\t")[1].strip()))

#list of frequencies by aa

S_df = pd.read_csv(args.matrix_model, sep="\t", index_col=0)
# PI = np.array([0.060490222, 0.066039665, 0.044127815, 0.042109048, 0.020075899, 0.053606488,
# 0.071567447, 0.072308239, 0.022293943, 0.069730629, 0.098851122, 0.056968211, 0.019768318,
# 0.028809447, 0.046025282, 0.05060433, 0.053636813, 0.033011601, 0.028350243, 0.061625237])
PI = np.array(frequencies)

S = np.array(S_df)

def get_normalised_generator(frequencies, rate_matrix=None):
    """
    Calculates the normalised generator from the rate matrix and character state frequencies.

    :param frequencies: character state frequencies.
    :type frequencies: numpy.array
    :param rate_matrix: (optional) rate matrix (by default an all-equal-rate matrix is used)
    :type rate_matrix: numpy.ndarray
    :return: normalised generator 1/mu Q
    :rtype: numpy.ndarray
    """
    if rate_matrix is None:
        n = len(frequencies)
        rate_matrix = np.ones(shape=(n, n), dtype=np.float64) - np.eye(n)
    generator = rate_matrix * frequencies
    generator -= np.diag(generator.sum(axis=1))
    mu = -generator.diagonal().dot(frequencies)
    generator /= mu
    return generator


Q_norm = get_normalised_generator(PI, S)


class Simulation:
    align = {}  # attribut de la classe,
    fasta_align = []

    def _random_choice_prob_index(
        self, a, axis=1
    ):  # choose the amino acid in the simulation
        r = np.random.rand(1)
        return (a.cumsum() > r).argmax()

    def Simulator_choice(self, T, root, rates):
        def process_node(node, parent_seq):
            target_seq = [0]*len(parent_seq)
            for siteidx, siteaa in enumerate(parent_seq):
                P = np.array(linalg.expm(Q_norm * node.dist * rates[siteidx]))  # 20,10 #matrix exp of normalized substitution matrix * length of branch *evolution rate
                target_seq[siteidx] = self._random_choice_prob_index(
                    P[siteaa, :], axis=1
                )  # create a list for the simulations: random sample of a vector of indexes corresponding to target amino acids knowing the parent ones
            sequence_aa = "".join([amino_acid[i] for i in target_seq]) #transform it in amino acids
            if node.is_leaf():
                self.align[node.name] = sequence_aa #attribute simulated sequence to node
                self.fasta_align.append(
                    SeqRecord(Seq(sequence_aa), node.name,
                              name=node.name, description="")
                ) #add simulated sequence in fasta
            for child in node.children:  # recursive
                process_node(child, target_seq)

        parent_seq = np.array(
            [amino_acid.index(i) for i in root]
        )  # fix the root as a array
        for child in T.children:
            process_node(child, parent_seq) #first iteration from root

start_time = time.time()
align_simu = Simulation()
align_simu.Simulator_choice(tree, root, rates)
#print("--- %s seconds ---" % (time.time() - start_time))

for rec in align_simu.fasta_align:
    print(f'>{rec.name}')
    print(f'{rec.seq}')
