#!/usr/bin/env python

import pandas as pd
from collections import Counter
import argparse

parser = argparse.ArgumentParser(description="Merge Correlation and Emergence results")

parser.add_argument("bayes_file", help="output BayesTraits", type=str)
parser.add_argument(
    "binary_file", help="table with phenotype and mutations per sequence", type=str
)
parser.add_argument(
    "lim_bayes", help="lim value for log Bayes factor", type=float, default=10
)
parser.add_argument(
    "all_results", help="path towards all_results_metrics.tsv", type=str
)

args = parser.parse_args()

bayes_file = args.bayes_file
binary = args.binary_file
lim_log = args.lim_bayes
tested = args.all_results


tested_df = pd.read_csv(tested, sep="\t")
binary_df = pd.read_csv(binary, sep="\t")
bayes_df = pd.read_csv(
    bayes_file,
    header=None,
    sep="\t",
    names=["x1", "log-dep", "x2", "log-indep", "BF", "posmut"],
)
bayes_df.drop(["x1", "x2"], axis=1, inplace=True)

treated = Counter(binary_df.phenotype)[1]
non_treated = Counter(binary_df.phenotype)[0]

corr = []
for pos_mut in binary_df.columns[1:-1]:
    counts = Counter(binary_df[(binary_df[pos_mut]) == 1].phenotype)
    if counts[1] / treated > counts[0] / non_treated:
        corr.append("positive")
    elif counts[1] / treated < counts[0] / non_treated:
        corr.append("negative")
    else:
        corr.append("equivalent")

bayes_df["correlation"] = corr
bayes_df.to_csv("BayesFactor.txt", sep="\t", index=None)

tested_df["posmut"] = [
    "".join([str(i), j]) for i, j in zip(tested_df.position, tested_df.mut)
]

All_bayes = pd.merge(tested_df, bayes_df, on="posmut", how="left")
ConDor_detected = All_bayes[
    (All_bayes.detected == "PASS")
    & (All_bayes.BF > lim_log)
    & (All_bayes.correlation == "positive")
]

All_bayes.to_csv("condor_tested_results.tsv", sep="\t", index=None)
ConDor_detected.to_csv("condor_detected_results.tsv", sep="\t", index=None)

