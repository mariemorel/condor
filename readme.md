# ACCURATE DETECTION OF CONVERGENT MUTATIONS IN LARGE PROTEIN ALIGNMENTS WITH CONDOR

### Marie MOREL, Frédéric LEMOINE, Anna ZHUKOVA and Olivier GASCUEL

Welcome to ConDor workflow repository ! 

In this repository you can find our ConDor workflow developed to detect convergent evolution in amino acid alignments. 
ConDor is available from a web service located at https://condor.pasteur.cloud/.

ConDor is adapted to the detection of evolutionary convergence at the resolution of a mutation especially in large datasets (several hundreds of sequences). To run ConDor you will need a multiple sequence alignment in fasta format including outgroup sequences, the corresponding tree in newick (with same sequence names as in the alignment), a file with outgroup sequence names and a file containing a list of sequences with the convergent phenotype (can be optional).

ConDor is composed of two independent components Emergence and Correlation, which can be launched together or not.

The output of ConDor is a list of mutations detected as convergent as they occur more often than expected under neutral evolution (Emergence component) and/or correlate (Correlation component) with the convergent phenotype (or a predictor of the phenotype). 


# Usage
To launch a full analysis (condor = both Emergence and Correlation components) with the test data (sedge dataset) run : 

    nextflow run condor.nf --align test_data/cyp_coding.aa.coor_mays.fa --tree test_data/cyp_coding.phy_phyml_tree.txt --outgroup test_data/outgroup.txt --phenotype besnard2009_convergent_species.txt --resdir output --model best --nb_simu 100 --min_seq 2 --min_eem 2 --freqmode Fmodel --branches condor --correction holm --alpha 0.1 --bayes 2

You can also chose to run only the Emergence component. In this case you would not need to provide the phenotype file nor BayesTraits parameters : 

    nextflow run condor.nf --align test_data/cyp_coding.aa.coor_mays.fa --tree test_data/cyp_coding.phy_phyml_tree.txt --outgroup test_data/outgroup.txt --resdir output --model best --nb_simu 100 --min_seq 2 --min_eem 2 --freqmode Fmodel --branches emergence --correction holm --alpha 0.1

Finally, if you want to run only the Correlation component: 

    nextflow run condor.nf --align test_data/cyp_coding.aa.coor_mays.fa --tree test_data/cyp_coding.phy_phyml_tree.txt --outgroup test_data/outgroup.txt --phenotype besnard2009_convergent_species.txt --resdir output --model best --min_seq 2 --min_eem 2 --freqmode Fmodel --branches correlation --bayes 2

For larger datasets (>1000), we recommand to increase --min_seq 10 and --bayes 20

# prerequisite 
https://www.nextflow.io/docs/latest/getstarted.html#requirements


# Help
Please visit https://condor.pasteur.cloud/help for more details regarding how to use ConDor and interpret the outputs. 



