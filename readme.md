# ACCURATE DETECTION OF CONVERGENT MUTATIONS IN PROTEIN ALIGNMENTS WITH CONDOR

### Marie MOREL, Frédéric LEMOINE, Anna ZHUKOVA and Olivier GASCUEL

Welcome to ConDor workflow repository ! 

In this repository you can find our ConDor workflow developed to detect convergent evolution in amino acid alignements. 
ConDor is available from a web service located at https://condor.pasteur.cloud/.

ConDor is adapted to the detection of evolutionary convergence at the resolution of a mutation especially in large datasets (several hundreds of sequences). To run ConDor you will need a multiple sequence alignment in fasta format inlcuding outgroup sequences, the corresponding tree in newick (with same sequence names as in the alignment), a file with outgroup sequence names and a file containing a list of sequences with the convergent phenotype.  

The output of ConDor is a list of mutations detected as convergent as they occur more often than expected under neutral evolution and correlate with the convergent phenotype (or a predictor of the phenotype). 


# Help
Please visit https://condor.pasteur.cloud/help if you have any questions regarding how to use ConDor and interpret the outputs. 



