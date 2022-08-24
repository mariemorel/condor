nextflow.enable.dsl=1
path_file = "/pasteur/zeus/projets/p01/Evolbioinfo/users/mamorel/Projet_Convergence/Data/"

/*
//c3c4 sedges
params.align = path_file + "args/raw/C3C4/cyp_coding.aa.coor_mays.fa"
params.tree = path_file + "args/raw/C3C4/cyp_coding.phy_phyml_tree.txt"
params.outgroup = path_file+"args/raw/C3C4/outgroup.txt"
params.phenotype = path_file+"args/raw/C3C4/besnard2009_convergent_species.txt"
*/

/*
//no c3c4 clade sedges adapted pheno
params.align = path_file + "args/processed/c3c4/no_c3c4_clade/cyp_coding_noc3c4.aa.coor_mays.fa"
params.tree = path_file + "args/processed/c3c4/no_c3c4_clade/cyp_coding_noc3c4.phy_phyml_tree.txt"
params.outgroup = path_file+"args/processed/c3c4/no_c3c4_clade/outgroup.txt"
params.phenotype = path_file+"args/processed/c3c4/no_c3c4_clade/besnard2009_c4pheno.txt" //phenotype could be optional ???
*/

/*
//HIV africa
params.align = path_file + "args/processed/africa/recomb_jphmm/align.noCRF.jphmm_outgroup.aa.fa"
params.tree = path_file + "args/processed/africa/recomb_jphmm/root.align.noCRF.jphmm_outgroup.fa.treefile"
params.outgroup = path_file+"args/processed/africa/recomb_jphmm/outgroup.txt"
params.phenotype = path_file+"args/processed/africa/recomb_jphmm/id_phenotype.txt"
*/

/*
//synthetic africa
params.align = path_file+ "args/processed/synthetic_africa/synthetic_out_simulated_RTafrica_5.fasta"
params.tree = path_file+ "args/processed/africa/recomb_jphmm/root.align.noCRF.jphmm_outgroup.fa.treefile"
params.outgroup = path_file+"args/processed/africa/recomb_jphmm/outgroup.txt"
params.phenotype = path_file+"args/processed/africa/recomb_jphmm/id_phenotype.txt"
*/

/*
//rhodopsin 
params.align = path_file + "args/processed/rhodopsin/lessgappy.align_reroot.fa" 
params.tree = path_file + "args/processed/rhodopsin/tree_rerooted.nhx"
params.outgroup = path_file+"args/processed/rhodopsin/outgroup_bis.txt"
params.phenotype = path_file+"args/processed/rhodopsin/id_phenotype_marine.txt" //id_phenotype_marine.txt id_phenotype.txt
*/

//declaration of parameters
params.help=false
params.resdir=path_file+"results/c3c4_phenotype/JTT+R3_noc3c4pheno/" //do not forget to change the resdir
params.model = "JTT+R3" // best: will choose the best model, other wise will take the given model
params.matrices = "$baseDir/assets/protein_model.txt"
params.nb_simu = 10000 //number of simulations to perform
params.min_seq = 2 //at least (11 for rhodopsine, 2 for c3c4, 10 for synthetic and HIV)
params.min_eem = 2 //strictly more than 2 EEMs
params.freqmode = "Fmodel" //Fmodel, if something else: FO. Need to be changed to allow other frequencies

/////OPTIONAL PARAMETERS
params.branches = "correlation" // condor, correlation, emergence. 
//params.branches_eem = "true"
//params.branches_corr = "true"
params.correlation = "all" //detected
params.correction = 'holm' // holm bonferroni correction , could be fdr_bh for benjamini hochberg
params.alpha = 0.1 //limit threshold (included)
params.bayes = 2

params.align="null"
params.tree="null"
params.outgroup="null"
params.phenotype="null"

def usage(){
    help="""CONDOR Usage:
    nextflow run condor.nf --align <input alignment (FASTA)> \\
                           --tree <input tree (NEWICK)> \\
                           --outgroup <outgroup file, one tip per line> \\
                           --phenotype <input tip phenotype data file, only with branches==(condor or correlation)> \\
                           --resdir <output directory> \\
                           --model <model or 'best'> \\
                           --matrices <directory where matrices are stored, default '$baseDir/assets/protein_model.txt'> \\
                           --nb_simu <number of simulations, default:10000> \\
                           --min_seq <min number of sequences having the mutation for convergence detection> \\
                           --min_eem <min number of EEMs> \\
                           --freqmode <amino acid frequencies: Fmodel or FO, default: Fmodel> \\
                           --branches <workflow run mode: condor, correlation or emergence, default: correlation> \\
                           --correction <multiple test correction, holm or fdr_bh, default: holm> \\
                           --alpha <alpha cutoff, default 0.1> \\
                           --bayes <log bayes factor threshold, default 2>
    """
    log.info(help.stripMargin())
}

if (params.help){
    usage()
    exit(0);
}

if (params.align==null || params.align == '' || !file(params.align).exists()){
    log.error("Error: Alignment file not defined or does not exist")
    usage()
    exit(1)
}
if (params.tree==null || params.tree == '' || !file(params.tree).exists()){
    log.error("Error: Tree file not defined or does not exist")
    usage()
    exit(1)
}

if (params.outgroup==null || params.outgroup == '' || !file(params.outgroup).exists()){
    log.error("Error: outgroup file not defined or does not exist")
    usage()
    exit(1)
}

if (params.matrices==null || params.matrices == '' || !file(params.matrices).exists()){
    log.error("Error: matrice file not defined or does not exist")
    usage()
    exit(1)
}

if (params.nb_simu <= 0){
    log.error("Error: number of simulations must be > 0")
    usage()
    exit(1)
}

if (!['Fmodel','FO'].contains(params.freqmode)){
    log.error("Error: wrong amino acid frequncy mode")
    usage()
    exit(1)
}

if (!["condor","correlation","emergence"].contains(params.branches)){
    log.error("Error: --branches must specify a valid run mode : condor, correlation, or emergence")
    usage()
    exit(1)
}

if (['condor','correlation'].contains(params.branches) && (params.phenotype==null || params.phenotype == '' || !file(params.phenotype).exists())){
    log.error("Error: phenotype file must be defined and must exist when --branches is condor or correlation")
    usage()
    exit(1)
}


//creation of parameters
align = file(params.align)
tree = file(params.tree)
outgroup = file(params.outgroup)
matrices = file(params.matrices)
model = params.model
nb_simu = params.nb_simu
min_seq = params.min_seq
min_eem = params.min_eem
freqmode = params.freqmode

//////////// Optional parameters
phenotype = file(params.phenotype)
correction = params.correction
alpha = params.alpha //risk 0.1 0.05
bayes = params.bayes //limit log Bayes Factor 2 20

//branches = params.branches
//branches_corr = params.branches_corr
//branches_eem = params.branches_eem

// create result directory
resdir=file(params.resdir)
resdir.with {mkdirs()}

//run iqtree model finder and take the best model or take the model given by user
process find_model {
    label 'iqtree'

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

//transform rate matrix in proper format
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
    label 'iqtree'

    input:
    val freqmode
    file align
    file tree
    file iqtree_mode from ModelChannel
    file matrices
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
    freqs=`grep -i $model -A 20 !{matrices} | tail -n 1 | sed 's/;//'`
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
    pastml_fasta.py !{pastml_acr} !{positions} !{length} !{marginal_root} !{nb_simu} ACR_
    '''
}

process count_apparitions{
    label 'python'

    input :
    val min_eem 
    tuple file(positions), file(rate), file(tree), file (align) from python_count
    output : 
    tuple file(rate), file(align), file ("*substitutions_even_root.tsv"), file("*substitutions_aa_tips_per_base.tsv") into Subscribe_matrices, Ref_couting
    file "positions_to_test_eem.txt" into Positions_simu, Positions_eem //num from 1
    file ("pos_mut_to_test_eem.txt") into Positions_correlation
    
    shell:
    //count EEMs from real data acr
    //min eem is strict >
    '''
    count_substitutions_from_tips.py !{align} !{tree} !{positions} !{min_eem} 
    '''
}

Subscribe_matrices.subscribe{rate, align, freqs, substitutions ->  freqs.copyTo(file("${resdir}").resolve('ref_substitutions.txt'));}

//////////END OF FIRST PART OF WORKFLOW : OBSERVED NB OF EEMS. 

//my simulator of sequence in python, using root, nb simulations and the ROOTED tree.
process simulator {
    //errorStrategy 'retry'
    //maxRetries 3
    label 'python'

    publishDir "${resdir}", pattern: "count*.tsv.gz", mode: 'copy'
    input : 
    //each x from ListPositionsChannel.readLines() //each tested positions (num from 1)
    each x from Positions_simu.flatMap{it.readLines()}
    file simulation_model from SimulatorMatrix //substitution matrix 
    file (rates) from RatesChannel //reestimated rates
    file (freq) from FreqChannel //frequencies (model or optimised)
    file (named_tree) from SimulatorChannel //rooted tree with named internal nodes
    file (root) from Marginal_root //root corresponding to marginal proba

    output : 
    file("count*npz") into MysimulationsChannel

    //when: branches_eem == "true" //("emergence" || "condor")
    when : branches == "emergence" || branches == 'condor'

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

process conclude_convergence{
    label 'python' //need to add statsmodels.api and statsmodels.stats in the python docker

    publishDir "${resdir}", mode: 'copy'
    
    input: 
    file simulation_model from SimulatorMatrix
    file align 
    file (freq) from FrequenciesChannel
    tuple file(rate), file(acralign), file(ref_matrix), file(substitutions) from Ref_couting //rates for all positions
    file (counts) from Collect_simulations
    file(root) from Root_seq //only the interesting positions starts from 1
    file(positions) from Positions_eem // starts from 1
    val nb_simu
    val min_seq // >= 
    val min_eem // > strict
    val alpha //0.1
    val correction //holm or fdr_bh

    //named*.phy
     
    output:
    tuple file("detected_metrics.tsv"), file("all_results_metrics.tsv") into EmergenceChannel, MergeChannel

    shell:
    '''
    convergent_substitutions_pvalue.py !{positions} !{root} !{rate} !{align} !{ref_matrix} !{substitutions} !{nb_simu} !{freq} !{simulation_model} !{min_seq} !{alpha} !{correction}
    '''
}


EmergenceChannel.subscribe{detected, tested ->  tested.copyTo(file("${resdir}").resolve('tested_results.tsv')); detected.copyTo(file("${resdir}").resolve('significant_results.tsv'));}

///////////END OF EMERGENCE PART OF WORKFLOW : EXPECTED NB OF EEMS. 

process prepare_BT {
    //conda '/pasteur/appa/homes/mamorel/miniconda3/envs/jupyter-notebook'
    publishDir "${resdir}", mode: 'copy'
    label 'python'
    input:
    file align
    file positions from Positions_correlation
    file phenotype 

    output: 
    file "*binary_tested_sites.tsv" into BinaryTraits, BinaryBayes

    when : branches == "correlation" || branches == 'condor'
    //when: branches_corr == 'true' //("correlation" || "condor")

    shell:
    '''
    bayes_traits_preps.py !{align} !{positions} !{phenotype}
    '''
}

/*
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
*/

process prepare_tree {
    label 'gotree'
    input: 
    file tree from BayestreeChannel 
    output:
    file "root_tree.nx" into NexusTree

    when : branches == "correlation" || branches == 'condor'
    
    //when: branches_corr == 'true' //("correlation" || "condor")

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

// BayesTraits
//     .map{it -> [it.getBaseName().split('_')[1] , it]}
//     .map{it -> [Integer.parseInt(it[0])%50, it[1] ]}
//     .groupTuple()
//     .set{ SplittedBayesTraits }

//could be replaced by maxForks 50   

process BayesTraits {
    errorStrategy 'retry'
    maxRetries 3
    memory "5G" 
    maxForks 50
    input: 
    file(data) from BayesTraits
    file nx_tree from NexusTree
    //tuple val(id), file(data_list) from SplittedBayesTraits
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
    tuple file(binary), file ("BayesFactor.txt") into CorrelationChannel, BayesChannel
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

process Conclude_BayesTraits {
    input: 
    val bayes
    tuple file(binary), file(bayesfactor) from BayesChannel
    output: 
    tuple file("detected_results.tsv"), file("tested_results.tsv") into Subscribe_BayesTraits
    
    when : branches == "correlation"     
    
    //when branches_corr == "true" && branches_eem == "false"  //"correlation"
    shell: 
    '''
    bayes_traits_filter.py !{bayesfactor} !{binary} !{bayes}
    '''

}

Subscribe_BayesTraits.subscribe{detected, tested ->  tested.copyTo(file("${resdir}").resolve('tested_results.tsv')); detected.copyTo(file("${resdir}").resolve('significant_results.tsv'));}

///////////END OF CORRELATION PART OF WORKFLOW. 

process Correlation {
    label 'python'
    publishDir "${resdir}", mode: 'copy'
    input: 
    val bayes
    tuple file (binary), file (bayesfactor) from CorrelationChannel
    tuple file (detected), file(tested) from MergeChannel
    output:
    tuple file ("BayesFactor.txt"), file("tested_results.tsv"), file("significant_results.tsv")
    shell:
    '''
    merge_results.py !{bayesfactor} !{binary} !{bayes} !{tested}
    '''
}