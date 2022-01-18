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


parser = argparse.ArgumentParser(
    description="Create an amino acid alignment corresponding to a tree and an ancestral sequence and count the number of substitutions from the root only if transmitted to leaves"
)
parser.add_argument(
    "root", help="probabilistic root for one position, we proceed position per position"
)  # one position
parser.add_argument(
    "tree_file", help="rooted tree with names at the internal nodes")
parser.add_argument("output", help="output file name")
parser.add_argument(
    "rates", help="value of the rate for the position")  # one rate
parser.add_argument(
    "freq", help="file with the freqs amino acid from the model")
parser.add_argument("matrix_model", help="matrix of the model")
args = parser.parse_args()

tree = PhyloTree(args.tree_file, format=1)
amino_acid = list("ARNDCQEGHILKMFPSTWYV")

rate = float(args.rates)  # should be one rate only

with open(args.root) as f:
    root_list =  f.read().splitlines()


root = root_list[0] #probabilistic root at one position
nb_residues = int(len(root)) #should equal to nb simulations

#print(root)


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
        r = np.expand_dims(np.random.rand(a.shape[1 - axis]), axis=axis)
        return (a.cumsum(axis=axis) > r).argmax(axis=axis)

    def Simulator_choice(self, T, root):
        def process_node(node, parent_seq):
            P = np.array(linalg.expm(Q_norm * node.dist * rate))  # 20,10 #matrix exp of normalized substitution matrix * length of branch *evolution rate
            target_seq = self._random_choice_prob_index(
                P[parent_seq, :], axis=1
            )  # create a list for the simulations: random sample of a vector of indexes corresponding to target amino acids knowing the parent ones
            sequence_aa = "".join([amino_acid[i] for i in target_seq]) #transform it in amino acids
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
align_simu.Simulator_choice(tree, root)
#print("--- %s seconds ---" % (time.time() - start_time))

align_txt = {}
for rec in align_simu.fasta_align:
    align_txt[rec.name] = rec.seq


# we did not simulate the root at it is an input --> add it in the tree
root_sequence = SeqRecord(Seq(root),
                          tree.name, name=tree.name, description="")
align_txt[root_sequence.name] = root_sequence.seq

# function to count EEMs from tips
class Count_lineages:
    counter_emergence = np.zeros(
        shape=(len(amino_acid), nb_residues))  # attribut de la classe,

    def _annotate_lineages(self, node): ##test if ancestral sequences have same aa as tip.
        if node.is_leaf():
            test = dict()
            for aa in amino_acid:
                test[aa] = np.array([i == aa for i in align_txt[node.name]])
            node.add_feature("lineage_bool", test)
            return test
        else:
            Answers = dict()
            child1, child2 = node.children
            p = self._annotate_lineages(child1)  # recursive
            q = self._annotate_lineages(child2)  # recursive
            for aa in amino_acid:
                X = p[aa] | q[aa]
                Answers[aa] = np.array(
                    [i == aa for i in align_txt[node.name]]) & X
            node.add_feature("lineage_bool", Answers)
            return Answers

    def Traverse_count(self, tree):
        self._annotate_lineages(tree)
        for node in tree.get_descendants():
            seq = align_txt[node.name]
            bool_dict = node.lineage_bool
            diff = np.array(
                [i != j for i, j in zip(seq, align_txt[node.up.name])])
            for aa in amino_acid:
                test_aa = np.array([i == aa for i in seq])
                X = diff & bool_dict[aa]
                to_add = X & test_aa
                self.counter_emergence[amino_acid.index(
                    aa), :] += to_add.astype(int)


start_time = time.time()
Count_lineages().Traverse_count(tree)
#print("--- %s seconds ---" % (time.time() - start_time))

ONE_POS = pd.DataFrame(
    Count_lineages.counter_emergence.astype(int), index=amino_acid)
# should be always the same position but times 100000 #!!!!!!!!!!!!!!!!!!!!!!! changes the index !!!!!!!!!!!! ARN --> ACD

# create sparse matrix from One Pos
data_csr = sparse.csr_matrix(np.array(ONE_POS))
sparse.save_npz(
    "".join(["count_", args.output, "{:04d}".format(nb_residues), ".npz"]), data_csr, "csr")
