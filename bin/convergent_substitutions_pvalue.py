#!/usr/bin/env python

import pandas as pd
from scipy import sparse

# mutliple correction
import statsmodels as sm
import statsmodels.api

# handle alignment
from Bio import AlignIO
from Bio.Align import AlignInfo
from collections import Counter
import numpy as np

# handle arguments
import argparse
import glob


parser = argparse.ArgumentParser(
    description="count the number of substitutions")
parser.add_argument("positions", help="positions to test")
parser.add_argument("root", help="root sequence")
parser.add_argument("rates", help="file with the rates per position")
parser.add_argument("align", help="file with the tips alignment in fasta")
parser.add_argument(
    "ref_matrix", help="matrix of the counts for the reference alignment")
parser.add_argument(
    "substitutions", help="list of every substitutions per position")
parser.add_argument(
    "nb_simu", help="Nb of simulations that have been performed")
parser.add_argument(
    "freq", help="frequencies of the model")
parser.add_argument(
    "matrix_model", help="model of evolution representating the data")
parser.add_argument(
    "min_seq", help="nb of sequences min to test convergence")
parser.add_argument(
    "min_eems", help="nb of emergence min to test convergence")
parser.add_argument(
    "alpha", help="threshold for multiple tests correction"
)
parser.add_argument(
    "correction", help="name of multiple tests correction: holm or fdr_bh"
)

args = parser.parse_args()
root_df = args.root
alignment_file = args.align
ref_counting = args.ref_matrix
nb_simu = int(args.nb_simu)
minseq = int(args.min_seq)
mineems = int(args.min_eems)
risk = float(args.alpha)
test_type = args.correction


list_pos = []
with open(args.positions) as f:
    for line in f:
        list_pos.append(int(line.strip())) 
#numerotation from 1

AA = list("ARNDCQEGHILKMFPSTWYV")

frequencies = []
with open(args.freq) as f:
    for line in f:
        frequencies.append(float(line.split("\t")[1].strip()))


root_seq = []
root_DF = pd.read_csv(root_df, sep="\t", header=None, names = ["pos", "root"])
for p in list_pos : 
    root_seq.append(root_DF[root_DF.pos == p].root.values[0])

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


# same order than frequencies ARNDC...
Q_norm = get_normalised_generator(PI, S)
Q_norm_df = pd.DataFrame(Q_norm, index=AA, columns=AA)


align = AlignIO.read(open(alignment_file), "fasta")
align_df = pd.DataFrame(align, index = [i.name for i in align])

summary_align = AlignInfo.SummaryInfo(align)
consensus = list(summary_align.dumb_consensus(threshold=0.1))

rate = [float(line.strip()) for line in open(args.rates)]

subset_rate = [rate[pos-1] for pos in list_pos]

substitutions_ref_df = pd.read_csv(ref_counting, sep="\t", index_col=0)
substitutions_ref_eem_df = substitutions_ref_df[[str(i) for i in list_pos]]


ref_df = np.array(substitutions_ref_eem_df.reindex(AA))

Substitutions_list_df = pd.read_csv(args.substitutions, sep="\t")
Substitutions_list_df.columns = ["position", "anc", "mut", "nb"]


Substitutions_list_df["pos_mut"] = ["".join([str(i), j]) for i, j in zip(
    Substitutions_list_df["position"], Substitutions_list_df["mut"])]
Substitutions_list_df["pos_anc"] = ["".join([str(i), j]) for i, j in zip(
    Substitutions_list_df["position"], Substitutions_list_df["anc"])]

# Proportion = pd.read_csv(args.proportion, sep="\t")
# Proportion.columns = ["position", "mut", "leaves",
#                       "transmitted", "proportion", "sum_proportion"]

# Proportion["pos_mut"] = ["".join([str(i), j]) for i, j in zip(
#     Proportion["position"], Proportion["mut"])]


DISTANCE_MATRIX = np.zeros(shape=(20, 20), dtype=np.int64)
DISTANCE_MATRIX[np.tril_indices(20, k=-1)] = \
    [2,
     2, 2,
     1, 2, 1,
     2, 1, 2, 2,
     2, 1, 2, 2, 3,
     1, 2, 2, 1, 3, 1,
     1, 1, 2, 1, 1, 2, 1,
     2, 1, 1, 1, 2, 1, 2, 2,
     2, 1, 1, 2, 2, 2, 2, 2, 2,
     2, 1, 2, 2, 2, 1, 2, 2, 1, 1,
     2, 1, 1, 2, 3, 1, 1, 2, 2, 1, 2,
     2, 1, 2, 3, 3, 2, 2, 2, 3, 1, 1, 1,
     2, 2, 2, 2, 1, 3, 3, 2, 2, 1, 1, 3, 2,
     1, 1, 2, 2, 2, 1, 2, 2, 1, 2, 1, 2, 2, 2,
     1, 1, 1, 2, 1, 2, 2, 1, 2, 1, 1, 2, 2, 1, 1,
     1, 1, 1, 2, 2, 2, 2, 2, 2, 1, 2, 1, 1, 2, 1, 1,
     2, 1, 3, 3, 1, 2, 2, 1, 3, 3, 1, 2, 2, 2, 2, 1, 2,
     2, 2, 1, 1, 1, 2, 2, 2, 1, 2, 2, 2, 3, 1, 2, 1, 2, 2,
     1, 2, 2, 1, 2, 2, 1, 1, 2, 1, 1, 2, 1, 1, 2, 2, 2, 2, 2]
DISTANCE_MATRIX = np.maximum(DISTANCE_MATRIX, DISTANCE_MATRIX.T)

distance_matrix_df = pd.DataFrame(DISTANCE_MATRIX, index=AA, columns=AA)


pvalues = []
all_results = []

def compute_pvalue(posnumber, position):
    anc = root_seq[posnumber] #min seq + eem
    cons_aa = consensus[position-1] #no filter on positions
    rate_pos = rate[position-1] #no filter on positions
    counts_simu = glob.glob("".join(["count_", str(position), "_named_tree_*.npz"]))[0] 
    load_csr = sparse.load_npz(counts_simu)
    dense = np.array(load_csr.todense()) 
    align_count = Counter(align_df[position-1]) # no filter on positions
    for aa, nb_seq in align_count.items():
        if (nb_seq >= minseq) and (aa in AA):
            aa_index = AA.index(aa)
            change = ref_df.T[posnumber][aa_index]
            if change > mineems: #strictly more than x EEMs
                variance = np.var(dense[aa_index])
                mean = np.mean(dense[aa_index])
                max_occur = max(load_csr[aa_index].data)
                higher_simu = len([i for i in load_csr[aa_index].data if i >= change])
                pval = higher_simu+0.5/(nb_simu+1)
                all_results.append(
                    [anc, cons_aa, position, aa, change, nb_seq,  max_occur, pval, variance, mean, rate_pos])


for posnum, pos in enumerate(list_pos):
    compute_pvalue(posnum, pos) #numerotation from 1


all_tests = pd.DataFrame(
    all_results,
    columns=[
        "pastml_root",
        "consensus_root",
        "position",
        "mut",
        "ref_EEM",
        "nbseq",
        "max_simu",
        "pvalue_raw",
        "variance",
        "mean",
        "evol_rate"
    ],
)


# all_tests = pd.merge(all_tests, Proportion, how='inner',
#                      on=[ "position",  'mut'])

# Substitutions behind
Type = []
Reversion_nb = []
details = []
Rev_details = []
genetic_distance = []
findability = []
exchangeability = []
max_anc = []

for acr, pos, mut in zip(all_tests.pastml_root,all_tests.position, all_tests.mut ):
    posmut = "".join([str(pos), mut])
    data = Substitutions_list_df[Substitutions_list_df["pos_mut"] == posmut]
    if len(data) ==1: #if only one aa leads to this mut. 
        if acr == mut :
            Type.append("reversion, parallel")
        else:
            Type.append("parallel")
    else:
        if acr == mut :
            Type.append("reversion, convergent")
        else:
            Type.append("convergent")
    Reversion_nb.append(
        sum(Substitutions_list_df[Substitutions_list_df["pos_anc"] == posmut].nb))
    Rev_details.append("; ".join([":".join([str(i), str(j)]) for i, j in dict(
        Substitutions_list_df[Substitutions_list_df["pos_anc"] == posmut][["mut", "nb"]].values).items()]))
    anc_dict = {anc:int(nb) for anc, nb in zip(data.anc, data.nb)}
    MaxKey = max(anc_dict, key=anc_dict.get)
    findability.append(np.log10(1/Q_norm_df[MaxKey][mut]))
    exchangeability.append(Q_norm_df[MaxKey][mut])
    genetic_distance.append(distance_matrix_df[MaxKey][mut])
    max_anc.append(MaxKey)
    details.append("; ".join([":".join([str(i), str(j)]) for i, j in dict(data[["anc", "nb"]].values).items()]))


all_tests["max_anc"] = max_anc
all_tests["type_substitution"] = Type
all_tests["details"] = details
all_tests["loss"] = Reversion_nb
all_tests["loss_details"] = Rev_details
all_tests["genetic_distance"] = genetic_distance
all_tests["substitution_rate"] = exchangeability
all_tests["findability"] = findability


#############################################
# at this step we have all the info on the substitutions
# but we don't known which ones pass the test
#############################################



pvals_holm = sm.stats.multitest.multipletests(all_tests.pvalue_raw, alpha=risk, method="holm", is_sorted=False, returnsorted=False)
pvals_fdr = sm.stats.multitest.multipletests(all_tests.pvalue_raw, alpha=risk, method="fdr_bh", is_sorted=False, returnsorted=False)


if test_type == "holm":
    all_tests["adjust_pvalue"] = list(pvals_holm[1])
    all_tests["adjust_pvalue_fdr"] = list(pvals_fdr[1])
    all_tests["detected_EEM"] = ["PASS" if i else "NOT PASS" for i in pvals_holm[0]]
    cols = ["pastml_root",
            "consensus_root",
            "position",
            "mut",
            "max_anc", 
            "ref_EEM",
            "nbseq",
            "evol_rate",
            "genetic_distance",
            "substitution_rate",
            "findability",
            "type_substitution",
            "details",
            "loss",
            "loss_details",
            "max_simu",
            "variance",
            "mean",
            "pvalue_raw",
            "adjust_pvalue",
            "adjust_pvalue_fdr",
            "detected_EEM"]
    all_tests = all_tests[cols]

else:
    all_tests["adjust_pvalue"] = list(pvals_fdr[1])
    all_tests["adjust_pvalue_holm"] = list(pvals_holm[1])
    all_tests["detected_EEM"] = ["PASS" if i else "NOT PASS" for i in pvals_fdr[0]]
    cols = ["pastml_root",
        "consensus_root",
        "position",
        "mut",
        "max_anc", 
        "ref_EEM",
        "nbseq",
        "evol_rate",
        "genetic_distance",
        "substitution_rate",
        "findability",
        "type_substitution",
        "details",
        "loss",
        "loss_details",
        "max_simu",
        "variance",
        "mean",
        "pvalue_raw",
        "adjust_pvalue",
        "adjust_pvalue_holm",
        "detected_EEM"]
    all_tests = all_tests[cols]

     
all_tests.sort_values('pvalue_raw', inplace=True)
all_tests.to_csv("all_results_metrics.tsv", sep="\t", index=False)
detected = all_tests[all_tests["detected_EEM"] == "PASS"]
detected.to_csv("detected_metrics.tsv", sep="\t", index=False)