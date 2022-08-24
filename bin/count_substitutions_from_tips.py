#!/usr/bin/env python

from ete3 import PhyloTree
import argparse
import pandas as pd
from Bio.Align import AlignInfo
from Bio import AlignIO


parser = argparse.ArgumentParser(
    description="Create from alignment and tree a matrix with nb substitutions"
)
parser.add_argument(
    "alignment_file", help="full path of alignment in FASTA from pastml"
)
parser.add_argument("tree_file", help="rooted tree with names internal nodes")
parser.add_argument("positions", help="positions to test") #1250 for synthetic 347 ror HIV real

parser.add_argument("min_eem", help = "minimum emergence-events to test the position")

args = parser.parse_args()

t = PhyloTree(args.tree_file, format=1)
t.link_to_alignment(args.alignment_file, alg_format="fasta")

align_Name = args.alignment_file.split(".")[0]

alignment = AlignIO.read(args.alignment_file, "fasta")
residues = AlignInfo.SummaryInfo(alignment)._get_all_letters()

eem = int(args.min_eem)

Pos_list = []
with open(args.positions) as f:
    for line in f:
        Pos_list.append(int(line.strip())) #numerotation at 1

amino_acid = list("ARNDCQEGHILKMFPSTWYV")


class Count_apparitions:
    # CHANGES = {j:{i:0 for i in amino_acid} for j in range(nb_residues)} #attribut de la classe,

    def __init__(self):
        self.changes = {
            j: {i: 0 for i in residues} for j in Pos_list
        }  # attribut de la classe,
        self.substitutions = {j: {i: {k: 0 for k in residues}
                                  for i in residues} for j in Pos_list}

    def apparitions_leaves(self, t):
        for posnum,pos in enumerate(Pos_list):
            T = self._attribution(t)
            for leaf in T.get_leaves():  # for each leaf in the tree
                # if leaf.sequence[pos] != T.sequence[pos]: #if the leaf amino acid is different from root amino acid
                self._recursive_count(leaf, posnum, pos)  # recursive count
        return([self.changes, self.substitutions])

    def _attribution(self, t):  # need to do this for each position
        for node in t.traverse():
            node.add_feature("traversed", set())
        return t

    def _recursive_count(self, node, posnum, pos):  # need to do this for each position
        if node.is_root():
            return
        if node.name in node.up.traversed: #if node was already seen
            return
        if node.sequence[posnum] != node.up.sequence[posnum]: #if there is a change
            node.up.traversed.add(node.name)
            # add a change towards amino acid at node and pos
            self.changes[pos][node.sequence[posnum]] += 1
            self.substitutions[pos][node.up.sequence[posnum]
                                    ][node.sequence[posnum]] += 1
        else:
            node = (
                node.up
            )  # if the node was never seen and there is no change, we go up
            self._recursive_count(node, posnum, pos)


results = Count_apparitions().apparitions_leaves(t)


df = pd.DataFrame(results[0])

DF = pd.DataFrame(df.describe().loc['max'] > eem )
pos_to_test = DF[DF[max]].index.to_list()

A = []
for i in df.columns:
    if df[i].describe().loc['max'] > eem :
        A.append(["".join([str(i),j]) for j in df[df[i] > eem].index.values])


with open("pos_mut_to_test_eem.txt", "w") as wf:
    for i in A:
        for j in i: 
            wf.write("".join(str(j))+"\n")


with open("positions_to_test_eem.txt", "w") as wf:
    wf.write("\n".join(list(map(str, pos_to_test))))

df.to_csv(
    "".join([align_Name, "substitutions_even_root.tsv"]), sep="\t", encoding="utf-8"
) ##numerotation from 1

per_base_tips = []
for pos in results[1]:
    for anc in results[1][pos]:
        for mut in results[1][pos][anc]:
            if results[1][pos][anc][mut] > 0:
                per_base_tips.append(
                    [pos, anc, mut, results[1][pos][anc][mut]])

pd.DataFrame(per_base_tips).to_csv("".join(
    [align_Name, 'substitutions_aa_tips_per_base.tsv']), sep='\t', index=False, encoding="utf-8") #numerotation from 1
