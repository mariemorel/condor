
path_file = "/pasteur/zeus/projets/p01/Evolbioinfo/users/mamorel/Projet_Convergence/Data/"

/*
//c3c4 sedges
params.align = path_file + "args/raw/C3C4/cyp_coding.aa.coor_mays.fa"
params.tree = path_file + "args/raw/C3C4/cyp_coding.phy_phyml_tree.txt"
params.outgroup = path_file+"args/raw/C3C4/outgroup.txt"
*/

/*
//HIV africa
params.align = path_file + "args/raw/africa/recomb_jphmm/align.noCRF.jphmm_outgroup.aa.fa"
params.tree = path_file + "args/raw/africa/recomb_jphmm/root.align.noCRF.jphmm_outgroup.fa.treefile"
params.outgroup = path_file+"args/raw/africa/recomb_jphmm/outgroup.txt"
*/


//synthetic africa
params.align = path_file+ "args/processed/synthetic_africa/synthetic_out_simulated_RTafrica_5.fasta"
params.tree = path_file+ "args/processed/africa/recomb_jphmm/root.align.noCRF.jphmm_outgroup.fa.treefile"
params.outgroup = path_file+"args/processed/africa/recomb_jphmm/outgroup.txt"
params.phenotype = path_file+"args/processed/africa/recomb_jphmm/id_phenotype.txt"


/*
//rhodopsin 
params.align = path_file + "args/raw/Rhodopsin/lessgappy.align_reroot.fa" 
params.tree = path_file + "args/raw/Rhodopsin/tree_rerooted.nhx"
params.outgroup = path_file+"args/raw/Rhodopsin/outgroup_bis.txt"
*/

params.resdir=path_file+"results/synthetic_africa/JTT_R4_nodrm5_bis/" //do not forget to change the resdir
params.model = "JTT+R4" // best: will choose the best model, other wise will take the given model
params.matrices = "$baseDir/assets/protein_model.txt"
params.correlation = "all"
params.nb_simu = 10000 //number of simulations to perform
params.min_seq = 10 //at least (11 for rhodopsine, 2 for c3c4, 10for synthetic, 10 for real HIV)
params.freqmode = "Fmodel" //if something else: FO. Need to be changed to allow other frequencies

params.correction = 'holm' // holm bonferroni correction , could be fdr_bh for benjamini hochberg
params.alpha = 0.1 //limit threshold (included)

align = file(params.align)
tree = file(params.tree)
outgroup = file(params.outgroup)
phenotype = file(params.phenotype)
matrices = file(params.matrices)
resdir=file(params.resdir)
resdir.with {mkdirs()}

nb_simu = params.nb_simu
min_seq = params.min_seq
model = params.model
alpha = params.alpha
freqmode = params.freqmode
correction = params.correction
correlation = params.correlation

//nb seqs and length align
//run iqtree model finder and take the best model or take the model given by user
process find_model {
    label 'iqtree2'

    input:
    file tree
    file align
    val model

    output:
    file "best_fit_model.txt" into ModelChannel, MatrixChannel
    shell:
    if ( model == 'best' )
    // modelfinder on user input tree and alignment
    '''
    iqtree -m MFP -s !{align} -te !{tree} -nt !{task.cpus}
    grep "Best-fit model" !{align}.iqtree | cut -d " " -f 6 > best_fit_model.txt
    '''
    else
    '''
    printf "!{model}" > best_fit_model.txt
    '''
}


process build_matrices {
    label 'python'

    input:
    file iqtree_modelrate from MatrixChannel
    file matrices

    output:
    file "*pastml_matrix" into PastmlMatrix
    file "*simulator_matrix.model" into SimulatorMatrix
    
    shell:
    //matrices of substitution rate readable by simulator or pastml
    '''
    matrix_name=`awk 'BEGIN{FS="+"} {print toupper($1) }' !{iqtree_modelrate}`
    matrix_pastml_format.py ${matrix_name} !{matrices}
    '''
}


process info_align {
    // retrieve length of alignment
    label 'goalign'

    input:
    file align
    output:
    stdout Length
    shell:
    '''
    len=`goalign stats -i !{align} | grep "length" | cut -f 2`
    printf $len
    '''
}

process info_align_nbseq {
    //retrieve Nb seqs in alignment
    label 'goalign'

    input:
    file align
    output:
    stdout Nb_seq
    shell:
    '''
    nb_seq=`goalign stats -i !{align} | grep "nseqs" | cut -f 2`
    printf $nb_seq
    '''
}

Stats_align = Length.merge(Nb_seq)
//lenght and nb seqs in alignment (I think we can remove nb seqs)


//reoptimize tree branch lengths and estimate rates and frequencies by ML 
process reoptimize_tree {
    label 'iqtree2'

    input:
    val freqmode
    file align
    file tree
    file iqtree_mode from ModelChannel
    tuple val(length), val(nb_seq) from Stats_align
    output:
    file "*.treefile" into TreeChannel 
    tuple val(length), file(align), file ("reestimated_rate"), file("frequencies.txt") into RatesparamsChannel
    file "reestimated_rate" into  RatesChannel
    file "frequencies.txt" into FreqChannel, FrequenciesChannel
    
    shell:  
    //te : fixed tree no tree search performed
    if ( freqmode == 'Fmodel' )
    //run iqtree with mode given by user. Reestimation rates. Frequencies of model retrieved from matrices file. 
    '''
    iqtreemode=`cat !{iqtree_mode}` 
    iqtree -m ${iqtreemode} -nt !{task.cpus} -s !{align} -te !{tree} -wsr -pre align
    tail -n+10 align.rate | cut -f 2 > reestimated_rate
    len=`wc -l reestimated_rate | cut -d " " -f 1`
    if [ "$len" -eq "0" ] ; then for i in {1..!{length}} ; do echo 1 >> reestimated_rate ;done ;  fi
    model=`awk 'BEGIN { FS="+" } {printf $1}' !{iqtree_mode}`
    freqs=`grep -i $model -A 20 protein_model.txt | tail -n 1 | sed 's/;//'` 
    AA="A R N D C Q E G H I L K M F P S T W Y V"
    paste <(tr ' ' '\n' <<< ${AA[*]}) <(tr ' ' '\n' <<< ${freqs[*]}) >  frequencies.txt
    '''
    else
    //optimized frequences. For now cannot work with empirical or given vector of frequencies
    //should work with no F or I 
    '''
    sed '/+F/!s/$/+FO/' !{iqtree_mode} | sed 's/+F[^+$]*/+FO/' > corrected_model
    iqtreemode=`cat corrected_model`
    iqtree -m ${iqtreemode} -nt !{task.cpus} -s !{align} -te !{tree} -wsr -pre align
    tail -n+10 align.rate | cut -f 2 > reestimated_rate
    len=`wc -l reestimated_rate | cut -d " " -f 1`
    if [ "$len" -eq "0" ] ; then for i in {1..!{length}} ; do echo 1 >> reestimated_rate ;done ;  fi
    for i in A R N D C Q E G H I L K M F P S T W Y V ; do grep "pi($i)" align.iqtree | awk -v var="$i" 'BEGIN{ORS="";print var"\\t"} {print $NF"\\n"}'; done > frequencies.txt  
    '''
}

process tree_rename{
    label 'gotree'

    publishDir "${resdir}", mode: 'copy'
    input:
    file tree from TreeChannel
    file outgroup
    output:
    file "named_tree" into SimulatorChannel, NamedtreeChannel
    file "rooted_*" into RootedtreeChannel, BayestreeChannel

    shell:
    //Remove branch length of root
    //remove outgroup ofter rerooting
    //give a name to internal nodes
    ''' 
    sed -i  '/);/!s/)[0-9]*.[0-9]*;/);/' !{tree}
    gotree reroot outgroup -r -i !{tree} -l !{outgroup} -o rooted_!{tree}
    gotree rename --internal --tips=false --auto -i rooted_!{tree} -o named_tree
    '''
}

process pars_align_file{
    label 'python'
 
    input : 
    tuple val(length), file(align), file(rates), file(freq) from RatesparamsChannel
    val min_seq

    output:
    tuple val(length), file(rates), file(freq), file('*pastml_input.tsv.gz') into Pastml_align mode flatten
    file "positions_to_test.txt" into PositionsChannel, ListPositionsChannel //positions numerotation from 1
  
    shell:
    //create a table of positions we test and for which we do acr
    '''
    pars_fasta_subset.py !{align} !{length} !{min_seq} refalign_
    '''
}

process acr_pastml{
    label 'pastml'

    publishDir "${resdir}", pattern: "work_pastml/named*.nw",  mode: 'copy'
    input:
    tuple val(length), file(rate), file(freq), file(input) from Pastml_align
    file tree from RootedtreeChannel
    file positions from PositionsChannel
    file matrix from PastmlMatrix

    output:
    tuple val(length),file(positions), file(rate), file("work_pastml/named.*nwk"), file("*pastml.ML.out.gz"), file("marginal_root.txt") into pastml_ML_out
    
     //rate sed numerotation from 1 ($line -1)
    //input numerotation from 0
    //create parameter file for pastml including rate per site (scaling factor) for each site and frequencies of amino acids
    //run pastml 
    //remove some temp files
    //retrive marginal proba for root (sed -n 2p)

    shell:
    '''
    align="!{input}"
    while read -r line; do 
    R=`sed -n "${line}"p !{rate}` 
    echo -e 'parameter\tvalue' > parameter_$((${line}-1)) 
    cat !{freq} >> parameter_$((${line}-1))
    echo -e "scaling_factor\t${R}" >> parameter_$((${line}-1)) ; done < !{positions}
      
    gunzip -c !{input} > ${align%.*.*}.input

    VAR=`while read -r line; do echo parameter_$((${line}-1)); done < !{positions}`
    ID=`while read -r line; do echo $((${line}-1)); done < !{positions}`
    MATRIX=`while read -r line; do echo !{matrix}; done < !{positions}`
    pastml --threads !{task.cpus} -t !{tree} -d ${align%.*.*}.input --prediction_method MAP -m CUSTOM_RATES --rate_matrix $MATRIX -c $ID -o ${align%.*.*}.pastml.ML.out --work_dir work_pastml --parameter $VAR
    gzip ${align%.*.*}.pastml.ML.out
    rm work_pastml/params*.tab
    rm parameter_*
    for i in `ls work_pastml/marginal_probabilities.character_*` ; do echo ${i//[^0-9]/} ; sed -n 2p $i ;  done > marginal_root.txt 
    '''
}


process pre_count{
    label 'python'

    publishDir "${resdir}", mode: 'copy'
    input : 
    tuple val(length), file(positions), file(rate), file(tree), file(pastml_acr), file(marginal_root) from pastml_ML_out
    val nb_simu

    output : 
    tuple file(positions), file(rate), file(tree), file("*pastml_acr.fasta") into python_count
    file "reconstructed_root" into Root_seq
    file "*marginal_posterior.txt" into Marginal_root

    shell:
    //transform pastml outpout in a fasta file and retrieve root with marginal proba
    '''
    pastml_fasta.py !{pastml_acr} !{positions} !{length} !{marginal_root} !{nb_simu} test_
    '''
}

//my simulator of sequence in python, using root, nb simulations and the ROOTED tree.
process simulator {
    //errorStrategy 'retry'
    //maxRetries 3
    label 'python'

publishDir "${resdir}", pattern: "count*.tsv.gz", mode: 'copy'
    input : 
    each x from ListPositionsChannel.readLines() //each tested positions (num from 1)
    file simulation_model from SimulatorMatrix //substitution matrix 
    file (rates) from RatesChannel //reestimated rates
    file (freq) from FreqChannel //frequencies (model or optimised)
    file (named_tree) from SimulatorChannel //rooted tree with named internal nodes
    file (root) from Marginal_root //root corresponding to marginal proba

    output : 
    file("count*npz") into MysimulationsChannel
    shell:
    //simulate and count EEMs from tips
    '''
    sed -n "/^$((!{x}-1))\t/p" !{root} | cut -f 2 > root.txt #root num from 0
    rate=`sed -n "!{x}p" !{rates}` # sed numerotation from 1
    
    output="!{x}_!{named_tree}_"
    simulator_counting_rates_from_root.py root.txt !{named_tree} ${output} ${rate} !{freq} !{simulation_model} 
    '''
}

Collect_simulations = MysimulationsChannel.collect()

//fasta file into a tab separated table understandable by pastml
//only input needed is align 
//56 different alignments 

process count_apparitions{
    label 'python'

    input : 
    tuple file(positions), file(rate), file(tree), file (align) from python_count
    output : 
    tuple file(positions), file(rate), file(align), file ("*substitutions_even_root.tsv"), file("*substitutions_aa_tips_per_base.tsv") into Subscribe_matrices, Ref_couting
    
    shell:
    //count EEMs from real data acr
    '''
    count_substitutions_from_tips.py !{align} !{tree} !{positions}
    '''
}

Subscribe_matrices.subscribe{positions, rate, align, freqs, substitutions ->  freqs.copyTo(file("${resdir}").resolve('ref_substitutions.txt'));}
 
// should 
process conclude_convergence{
    label 'python' //need to add statsmodels.api and statsmodels.stats in the python docker

    publishDir "${resdir}", mode: 'copy'
    
    input: 
    file simulation_model from SimulatorMatrix
    file align 
    file (freq) from FrequenciesChannel
    tuple file(positions), file(rate), file(acralign), file(ref_matrix), file(substitutions) from Ref_couting //rates for all positions
    file (counts) from Collect_simulations
    file(root) from Root_seq //only the interesting positions
    val nb_simu
    val min_seq
    val alpha //0.1
    val correction //holm or fdr_bh

    //named*.phy
     
    output:
    tuple file("detected_metrics.tsv"), file("all_results_metrics.tsv") into BayesChannel

    shell:
    '''
    R=`cat !{root}`
    convergent_substitutions_pvalue.py !{positions} ${R} !{rate} !{align} !{ref_matrix} !{substitutions} !{nb_simu} !{freq} !{simulation_model} !{min_seq} !{alpha} !{correction}
    '''
}

process prepare_BT {
    conda '/pasteur/appa/homes/mamorel/miniconda3/envs/jupyter-notebook'
    publishDir "${resdir}", mode: 'copy'
    //label 'python'
    input:
    file align
    tuple file(detected), file(tested) from BayesChannel
    file phenotype 
    val correlation

    output: 
    file "*binary_tested_sites.tsv" into BinaryTraits, BinaryBayes
    shell:
    if ( correlation == "all" )
    '''
    bayes_traits_preps.py !{align} !{tested} !{phenotype}
    '''
    else if ( correlation == "detected" )
    '''
    bayes_traits_preps.py !{align} !{detected} !{phenotype}
    '''
    else
    error "Invalid bayesTraits mode: ${correlation}"
}

process prepare_tree {
    label 'gotree'
    input: 
    file tree from BayestreeChannel 
    output:
    file "root_tree.nx" into NexusTree
    shell:
    '''
    gotree stats tips -i !{tree} | cut -f 4 | tail -n+2 > names.treefile.txt
    END=`wc -l names.treefile.txt | cut -d ' ' -f 1`
    for (( i=1; i<=${END}; i++ )); do echo $i >> id.treefile.txt ; done
    paste id.treefile.txt names.treefile.txt > map_file
    gotree rename -m map_file -r -i !{tree} -o !{tree}.rename
    gotree support clear -i !{tree}.rename | gotree reformat nexus -o !{tree}.rename.nx
    echo ';\n' >> map_file
    sed '/^BEGIN TREES;/a TRANSLATE\n' !{tree}.rename.nx > test_tree.nx
    sed -e '/^TRANSLATE/r map_file' test_tree.nx > root_tree.nx
    '''
}

process DataTraits {
    input: 
    file binary from BinaryTraits
    output:
    file ("data*") into BayesTraits mode flatten
    shell:
    ''' 
    END=`awk -F '\t' '{print NF}' !{binary} | sort -nu | tail -n 1`
    for (( i=2; i<=$((${END}-1)); i++ )); do cut -f 1,$i,$END !{binary} | tail -n+2 > data_$i.txt ; done
    '''
}

BayesTraits
    .map{it -> [it.getBaseName().split('_')[1] , it]}
    .map{it -> [Integer.parseInt(it[0])%50, it[1] ]}
    .groupTuple()
    .set{ SplittedBayesTraits }
    

process BayesTraits {
    errorStrategy 'retry'
    maxRetries 3
    memory "5G" 
    input: 
    tuple val(id), file(data_list) from SplittedBayesTraits
    file nx_tree from NexusTree
    output:
    file ("*Stones*") into Stones mode flatten
    shell:
    '''
    for data in `ls data*.txt` ; do 
    y=`echo $data` 
    x=${y//[!0-9]/} 
    echo -e  "3\n2\nPriorAll uniform 0 100\nStones 100 1000\nLogFile Dependent_MCMC_10_${x}\nRun" > cmd_MCMC$x.txt;
    BayesTraitsV3 !{nx_tree} $data < cmd_MCMC$x.txt; 
    echo -e  "2\n2\nPriorAll uniform 0 100\nStones 100 1000\nLogFile Independent_MCMC_10_${x}\nRun" > cmd_Independent_MCMC$x.txt; 
    BayesTraitsV3 !{nx_tree} $data < cmd_Independent_MCMC$x.txt; done
    '''
}

Collect_stones = Stones.collect()

process BayesFactor {
    input: 
    file stones from Collect_stones
    file nx_tree from NexusTree
    file binary from BinaryBayes
    output:
    tuple file(binary), file ("BayesFactor.txt") into CorrelationChannel
    shell:
    '''
    END=`awk -F '\t' '{print NF}' !{binary} | sort -nu | tail -n 1`
    for (( i=2; i<=$((${END}-1)); i++ )); do printf $i"\t" ; grep "Log marginal likelihood:" Dependent_MCMC_10_$i.Stones.txt | cut -f 2 ; done > Dependent_results.txt
    for (( i=2; i<=$((${END}-1)); i++ )); do printf $i"\t" ; grep "Log marginal likelihood:" Independent_MCMC_10_$i.Stones.txt | cut -f 2 ; done > Independent_results.txt
    paste Dependent_results.txt Independent_results.txt > Dep_Indep_results.txt
    while read p; do
        d=`echo $p | awk -F ' ' '{print $2}'`;
        i=`echo $p | awk -F ' ' '{print $4}'`;
        echo "2*($d- $i)" | bc >> tmp
    done <Dep_Indep_results.txt
    for i in `head -n 1 binary_tested_sites.tsv` ; do echo $i ; done | head -n -1 > names.txt
    paste Dep_Indep_results.txt tmp names.txt > BayesFactor.txt
    '''
}

process Correlation {
    conda '/pasteur/appa/homes/mamorel/miniconda3/envs/jupyter-notebook'
    publishDir "${resdir}", mode: 'copy'
    input: 
    tuple file (binary), file (bayesfactor) from CorrelationChannel
    output:
    file ("BayesFactor.txt")
    script:
    '''
    #!/usr/bin/env python
    
    import pandas as pd
    from collections import Counter
    
    binary_df = pd.read_csv("binary_tested_sites.tsv",sep='\t')
    bayes_df = pd.read_csv("BayesFactor.txt", header = None, sep='\t', names = ["x1", "log-dep", "x2", "log-indep", "BF", "posmut"])
    bayes_df.drop(["x1", "x2"], axis=1, inplace = True)
    corr = []
    for pos_mut in binary_df.columns[1:-1]:
        counts = Counter(binary_df[(binary_df[pos_mut]) == 1 ].phenotype)
        if counts[1] > counts[0]:
            corr.append("positive")
        else:
            corr.append("negative")
    
    bayes_df["correlation"] = corr
    bayes_df.to_csv("BayesFactor.txt", sep='\t', index=None)
    '''

}