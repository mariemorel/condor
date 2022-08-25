#!/usr/bin/env python

import pandas as pd
from collections import Counter
import argparse

parser = argparse.ArgumentParser(description="Correlation results")

parser.add_argument("bayes_file", help="output BayesTraits", type=str)
parser.add_argument(
    "binary_file", help="table with phenotype and mutations per sequence", type=str
)
parser.add_argument(
    "lim_bayes", help="lim value for log Bayes factor", type=float, default=10
)

args = parser.parse_args()

bayes_file = args.bayes_file
binary = args.binary_file
lim_log = args.lim_bayes



binary_df = pd.read_csv(binary, sep="\t")
bayes_df = pd.read_csv(
    bayes_file,
    header=None,
    sep="\t",
    names=["x1", "log-dep", "x2", "log-indep", "BF", "posmut"],
)
bayes_df.drop(["x1", "x2"], axis=1, inplace=True)


cols = ["posmut", "log-dep", "log-indep", "BF"]

bayes_df = bayes_df[cols] 

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
bayes_df.to_csv("tested_results.tsv", sep="\t", index=None)


Bayes_detected = bayes_df[(bayes_df.BF > lim_log)
    & (bayes_df.correlation == "positive")
]

Bayes_detected.to_csv("detected_results.tsv", sep="\t", index=None)