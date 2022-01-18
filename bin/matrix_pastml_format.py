#!/usr/bin/env python3

import pandas as pd 
import numpy as np
import argparse


model_list = ["BLOSUM62","CPREV","DAYHOFF","FLU","HIVB","HIVW","JTT","DCMUT","JTTDCMUT","LG","MTART","MTMAM","MTREV","MTZOA","MTMET","MTVER","MTINV",
"POISSON","PMB","RTREV","VT","WAG","GTR20","Q_PFAM","Q_PFAM_GB","Q_LG","Q_BIRD","Q_INSECT","Q_MAMMAL","Q_PLANT","Q_YEAST","FLAVI"]

parser = argparse.ArgumentParser(description="find the model matrix and adapt to pastml format")
parser.add_argument("model", help="best model found by iqtree MFP, default = LG")
parser.add_argument("matrices", help="File containing all the matrices")

args = parser.parse_args()

model_name=str(args.model)
matrices=str(args.matrices)


AA_1_LETTER_CODES = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G',
                     'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S',
                     'T', 'W', 'Y', 'V']

NUM_AA = 20

model_dict={}
list_rates = []
for line in open(matrices): #open iqtree matrices file
    if line.startswith("model"):
        if len(list_rates) > 0:
            model_dict[model] = list_rates # incremente for each model a dictionnary containing the list of substitution rates. 
            list_rates = [] # empty list for this model and start with a new one
            model = line.split(" ")[1].strip("\n") #retrieve model name for the second model and after
        else:
            model = line.split(" ")[1].strip("\n") #retrieve first model name 
    elif line.strip().endswith(";"):
        freq = [float(i) for i in line.strip().strip(";").split(" ")] # once we have all the subsitution rates --> frequencies
    elif len(line.strip()) == 0: 
        pass #if line empty pass 
    else:
        list_rates.extend([float(i) for i in line.strip().split(" ")]) # fill the list of rates with rates
    model_dict[model] = list_rates #last model

MODEL_RATE_MATRIX = np.zeros(shape=(NUM_AA, NUM_AA), dtype=np.float64) #create empty matrix 20*20
MODEL_RATE_MATRIX[np.tril_indices(NUM_AA, k=-1)] = model_dict[model_name] # fill lower triangle matrix with substitution rates 
MODEL_RATE_MATRIX = np.maximum(MODEL_RATE_MATRIX, MODEL_RATE_MATRIX.T) #fill upper triangle matrix

pd.DataFrame(MODEL_RATE_MATRIX, index = AA_1_LETTER_CODES, columns = AA_1_LETTER_CODES).to_csv("".join([model_name, "simulator_matrix.model"]), sep="\t") #save substitution matrix in a format readable by simulator

with open("".join([model_name,"pastml_matrix"]), 'w') as wf:
    wf.write(" ".join(["#", " ".join(AA_1_LETTER_CODES)])+"\n")
    for i in MODEL_RATE_MATRIX:
        wf.write(" ".join([str(a) for a in list(i)])+"\n")

#save substitution matrix in a format readable by pastml
