#!/usr/bin/env python

import pandas as pd 
import numpy as np
import argparse


model_list = ["BLOSUM62","CPREV","DAYHOFF","FLU","HIVB","HIVW","JTT","DCMUT","JTTDCMUT","LG","MTART","MTMAM","MTREV","MTZOA","MTMET","MTVER","MTINV",
"POISSON","PMB","RTREV","VT","WAG","GTR20","Q_PFAM","Q_PFAM_GB","Q_LG","Q_BIRD","Q_INSECT","Q_MAMMAL","Q_PLANT","Q_YEAST","FLAVI"]

parser = argparse.ArgumentParser(description="find the model matrix and adapt to pastml format")
parser.add_argument("model", help="best model found by iqtree MFP, default = LG")

args = parser.parse_args()

model_name=str(args.model)

AA_1_LETTER_CODES = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G',
                     'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S',
                     'T', 'W', 'Y', 'V']

NUM_AA = 20

model_dict={}
list_rates = []
for line in open("/pasteur/zeus/projets/p01/Evolbioinfo/users/mamorel/Projet_Convergence/Code/detecting-convergent-evolution/protein_model.txt"):
    if line.startswith("model"):
        if len(list_rates) > 0:
            model_dict[model] = list_rates
            list_rates = []
            model = line.split(" ")[1].strip("\n")
        else:
            model = line.split(" ")[1].strip("\n")
    elif line.strip().endswith(";"):
        freq = [float(i) for i in line.strip().strip(";").split(" ")]
    elif len(line.strip()) == 0:
        pass
    else:
        list_rates.extend([float(i) for i in line.strip().split(" ")])
    model_dict[model] = list_rates



MODEL_RATE_MATRIX = np.zeros(shape=(NUM_AA, NUM_AA), dtype=np.float64)
MODEL_RATE_MATRIX[np.tril_indices(NUM_AA, k=-1)] = model_dict[model_name]
MODEL_RATE_MATRIX = np.maximum(MODEL_RATE_MATRIX, MODEL_RATE_MATRIX.T)



pd.DataFrame(MODEL_RATE_MATRIX, index = AA_1_LETTER_CODES, columns = AA_1_LETTER_CODES).to_csv("".join([model_name, "simulator_matrix.model"]), sep="\t")

with open("".join([model_name,"pastml_matrix"]), 'w') as wf:
    wf.write(" ".join(["#", " ".join(AA_1_LETTER_CODES)])+"\n")
    for i in MODEL_RATE_MATRIX:
        wf.write(" ".join([str(a) for a in list(i)])+"\n")

