nextflow.enable.dsl=2

path_file = "$baseDir/test_data/"

/*
//c3c4 sedges
params.align = path_file + "C3C4/cyp_coding.aa.coor_mays.fa"
params.tree = path_file + "C3C4/cyp_coding.phy_phyml_tree.txt"
params.outgroup = path_file+"C3C4/outgroup.txt"
params.phenotype = path_file+"C3C4/besnard2009_convergent_species.txt"
*/

/*
//no c3c4 clade sedges adapted pheno
params.align = path_file + "no_c3c4_clade/cyp_coding_noc3c4.aa.coor_mays.fa"
params.tree = path_file + "no_c3c4_clade/cyp_coding_noc3c4.phy_phyml_tree.txt"
params.outgroup = path_file+"no_c3c4_clade/outgroup.txt"
params.phenotype = path_file+"no_c3c4_clade/besnard2009_c4pheno.txt" //phenotype could be optional ???
*/

/*
//HIV
params.align = path_file + "HIV/align.noCRF.jphmm_outgroup.aa.fa"
params.tree = path_file + "HIV/root.align.noCRF.jphmm_outgroup.fa.treefile"
params.outgroup = path_file+"HIV/outgroup.txt"
params.phenotype = path_file+"HIV/id_phenotype.txt"
*/

/*
//rhodopsin 
params.align = path_file + "rhodopsin/lessgappy.align_reroot.fa" 
params.tree = path_file + "rhodopsin/tree_rerooted.nhx"
params.outgroup = path_file+"rhodopsin/outgroup_bis.txt"
params.phenotype = path_file+"rhodopsin/id_phenotype_marine.txt" //id_phenotype_marine.txt id_phenotype.txt
*/

//declaration of parameters
params.help=false
params.resdir=path_file+"results/c3c4_phenotype/JTT+R3/" //do not forget to change the resdir
params.model = "JTT+R3" // best: will choose the best model, other wise will take the given model
params.matrices = "$baseDir/assets/protein_model.txt"
params.nb_simu = 10000 //number of simulations to perform
params.min_seq = 2 //at least (11 for rhodopsine, 2 for c3c4, 10 for synthetic and HIV) = 0.5% seq
params.min_eem = 3 // At least 3 EEMs
params.positions= "none" // if the file exists, then only the given positions (0-based coordinates) are analysed by condor and min_seq is useless, but the branches and rates are optimised on the full alignment.
params.freqmode = "Fmodel" //Fmodel, if something else: FO. Need to be changed to allow other frequencies

/////OPTIONAL PARAMETERS
params.branches = "condor" // condor, correlation, emergence. 
//params.branches_eem = "true"
//params.branches_corr = "true"
params.correction = 'holm' // holm bonferroni correction , could be fdr_bh for benjamini hochberg
params.alpha = 0.1 //limit threshold (included)
params.bayes = 2

params.align="null"
params.tree="null"
params.outgroup="null"
params.phenotype="null"

params.rates="rates.txt"

//creation of parameters
align = file(params.align)
tree = file(params.tree)
outgroup = file(params.outgroup)
matrices = file(params.matrices)
model = params.model
nb_simu = params.nb_simu
min_seq = params.min_seq
min_eem = params.min_eem
positions = file(params.positions)
freqmode = params.freqmode

//////////// Optional parameters
phenotype = file(params.phenotype)
correction = params.correction
alpha = params.alpha //risk 0.1 0.05
bayes = params.bayes //limit log Bayes Factor 2 20
branches = params.branches
//branches_corr = params.branches_corr
//branches_eem = params.branches_eem

// create result directory
resdir=file(params.resdir)
resdir.with {mkdirs()}



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
                           --branches <workflow run mode: condor, correlation or emergence, default: condor> \\
                           --correction <multiple test correction, holm (holm bonferroni) or fdr_bh, default: holm> \\
                           --alpha <alpha cutoff, default 0.1> \\
                           --bayes <log bayes factor threshold, default 2>
    """
    log.info(help.stripMargin())
}

//run iqtree model finder and take the best model or take the model given by user
process find_model {
    publishDir "${resdir}", mode: 'copy'
    
    label 'iqtree'

    input:
    path tree
    path align
    val model

    output:
    path "best_fit_model.txt"

    shell:
    if ( model == 'best' )
    // modelfinder on user input tree and alignment
    '''
    iqtree -m MFP -s !{align} -te !{tree} -nt !{task.cpus} -pre mfp_!{align}
    grep "Best-fit model" mfp_!{align}.iqtree | cut -d " " -f 6 > best_fit_model.txt
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
    path iqtree_modelrate
    path matrices

    output:
    path "*pastml_matrix"
    path "*simulator_matrix.model"
    
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
    path align

    output:
    stdout

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
    path align
    
    output:
    stdout

    shell:
    '''
    nb_seq=`goalign stats -i !{align} | grep "nseqs" | cut -f 2`
    printf $nb_seq
    '''
}

//reoptimize tree branch lengths and estimate rates and frequencies by ML 
process reoptimize_tree {
    label 'iqtree'

    publishDir "${resdir}", mode: 'copy'
    
    input:
    val freqmode
    path align
    path tree
    path iqtree_mode
    path matrices
    tuple val(length), val(nb_seq)

    output:
    path "align.treefile"
    tuple val(length), path(align), path("reestimated_rate"), path("frequencies.txt")
    path "reestimated_rate"
    path "frequencies.txt"
    
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
    path tree
    path outgroup

    output:
    path "named_tree"
    path "rooted_*"

    shell:
    //Remove branch length of root
    //remove outgroup ofter rerooting
    //give a name to internal nodes
    ''' 
    sed  '/);/!s/)[0-9]*.[0-9]*;/);/' !{tree} > !{tree}_tmp
    sed  's/[\\/\\|]/_/g' !{outgroup}  > outgroup_tmp
    gotree reroot outgroup -r -i !{tree}_tmp -l outgroup_tmp -o rooted_!{tree}
    gotree rename --internal --tips=false --auto -i rooted_!{tree} -o named_tree
    '''
}

process pars_align_file{
    label 'python'
 
    input : 
    tuple val(length), path(align), path(rates), path(freq)
    path positions
    val min_seq

    output:
    tuple val(length), path(rates), path(freq), path('*pastml_input.tsv.gz')
    path "positions_to_test.txt"
  
    shell:
    //create a table of positions we test and for which we do acr
    '''
    pars_fasta_subset.py !{align} !{length} !{min_seq} refalign_ !{positions}
    '''
}

process acr_pastml{
    label 'pastml'

    publishDir "${resdir}", pattern: "work_pastml/named*.nw",  mode: 'copy'

    input:
    tuple val(length), path(rate), path(freq), path(input)
    path tree
    path positions
    path matrix

    output:
    tuple val(length),path(positions), path(rate), path("work_pastml/named.*nwk"), path("*pastml.ML.out.gz"), path("marginal_root.txt")
    
    //rate sed numerotation from 1 ($line -1)
    //input numerotation from 0
    //create parameter file for pastml including rate per site (scaling factor) for each site and frequencies of amino acids
    //run pastml 
    //remove some temp files
    //retrieve marginal proba for root (sed -n 2p)

    shell:
    '''
    align="!{input}"
    while read -r line; do 
    R=`sed -n "${line}"p !{rate}` 
    echo -e 'parameter\tvalue' > parameter_$((${line}-1)) 
    cat !{freq} >> parameter_$((${line}-1))
    echo -e "scaling_factor\t${R}" >> parameter_$((${line}-1)) ; done < !{positions}
      
    gunzip -c !{input} | sed 's/[\\/\\|]/_/g' > ${align%.*.*}.input

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
    tuple val(length), path(positions), path(rate), path(tree), path(pastml_acr), path(marginal_root)
    val nb_simu

    output : 
    tuple path(positions), path(rate), path(tree), path("*pastml_acr.fasta"), path("reconstructed_root.txt")
    path "reconstructed_root.txt" //positions with min_seq //first column positions associated //root num from 1
    path "*marginal_posterior.txt" //still positions associated root num from 0

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
    val min_seq
    path align
    tuple path(positions), path(rate), path(tree), path (align_acr), path(root)

    output : 
    tuple path(rate), path(align_acr), path("*substitutions_even_root.tsv"), path("*substitutions_aa_tips_per_base.tsv")
    path "positions_to_test_eem.txt" //num from 1
    path "pos_mut_to_test_eem.txt"
    
    shell:
    //count EEMs from real data acr
    //min eem is >= so we subtract 1 to be strict >
    '''
    count_substitutions_from_tips.py !{align_acr} !{align} !{tree} !{positions} !{min_seq} !{min_eem-1} !{root} 
    '''
}

//////////END OF FIRST PART OF WORKFLOW : OBSERVED NB OF EEMS. 

//my simulator of sequence in python, using root, nb simulations and the ROOTED tree.
process simulator {
    //errorStrategy 'retry'
    //maxRetries 3
    label 'python'

    publishDir "${resdir}", pattern: "count*.tsv.gz", mode: 'copy'

    input : 
    each pos //each tested positions (num from 1) ListPositionsChannel.readLines()
    path simulation_model //substitution matrix 
    path rates //reestimated rates
    path freq  //frequencies (model or optimised)
    path named_tree //rooted tree with named internal nodes
    path root //root corresponding to marginal proba

    output : 
    path "count*npz"

    //when: branches_eem == "true" //("emergence" || "condor")
    when : branches == "emergence" || branches == 'condor'

    shell:
    //simulate and count EEMs from tips
    '''
    sed -n "/^$((!{pos}-1))\t/p" !{root} | cut -f 2 > root.txt #root num from 0
    rate=`sed -n "!{pos}p" !{rates}` # sed numerotation from 1
    
    output="!{pos}_!{named_tree}_"
    simulator_counting_rates_from_root.py root.txt !{named_tree} ${output} ${rate} !{freq} !{simulation_model} 
    '''
}

process conclude_convergence{
    label 'python' //need to add statsmodels.api and statsmodels.stats in the python docker

    publishDir "${resdir}", mode: 'copy'
    
    input: 
    path simulation_model
    path align 
    path freq
    tuple path(rate), path(acralign), path(ref_matrix), path(substitutions)//rates for all positions
    path counts
    path root //only the interesting positions starts from 1
    path positions // starts from 1
    val nb_simu
    val min_seq // >= 
    val min_eem // >= so we subtract 1 to be > strict
    val alpha //0.1
    val correction //holm or fdr_bh

    //named*.phy
    output:
    tuple path("detected_metrics.tsv"), path("all_results_metrics.tsv")

    shell:
    '''
    convergent_substitutions_pvalue.py !{positions} !{root} !{rate} !{align} !{ref_matrix} !{substitutions} !{nb_simu} !{freq} !{simulation_model} !{min_seq} !{min_eem-1} !{alpha} !{correction}
    '''
}

///////////END OF EMERGENCE PART OF WORKFLOW : EXPECTED NB OF EEMS. 

process prepare_BT {
    //conda '/pasteur/appa/homes/mamorel/miniconda3/envs/jupyter-notebook'
    publishDir "${resdir}", mode: 'copy'
    label 'python'

    input:
    path align
    path positions
    path phenotype 
    val branches

    output: 
    path "*binary_tested_sites.tsv"

    when : branches == "correlation" || branches == 'condor'
    //when: branches_corr == 'true' //("correlation" || "condor")

    shell:
    '''
    bayes_traits_preps.py !{align} !{positions} !{phenotype}
    '''
}


process prepare_tree {
    label 'gotree'

    input: 
    path tree
    val branches

    output:
    path "root_tree.nx"

    when : branches == "correlation" || branches == 'condor'
    
    //when: branches_corr == 'true' //("correlation" || "condor")

    shell:
    '''
    gotree support clear -i !{tree} | gotree rename -e ".*" -b "" --tips=false --internal | gotree reformat nexus --translate -o root_tree.nx
    '''
}

process DataTraits {

    input: 
    path binary

    output:
    path "data*"

    shell:
    ''' 
    END=`awk -F '\t' '{print NF}' !{binary} | sort -nu | tail -n 1`
    for (( i=2; i<=$((${END}-1)); i++ )); do cut -f 1,$i,$END !{binary} | tail -n+2 > data_$i.txt ; done
    '''
}

process BayesTraits {
    label 'bayestraits'

    input: 
    path data
    path nx_tree

    output:
    path "*Stones*"

    shell:
    '''
    y=`echo !{data}`
    x=${y//[!0-9]/}
    echo -e  "3\n2\nPriorAll uniform 0 100\nStones 100 1000\nLogFile Dependent_MCMC_10_${x}\nRun" > cmd_MCMC$x.txt;
    BayesTraitsV3 !{nx_tree} !{data} < cmd_MCMC$x.txt;
    echo -e  "2\n2\nPriorAll uniform 0 100\nStones 100 1000\nLogFile Independent_MCMC_10_${x}\nRun" > cmd_Independent_MCMC$x.txt;
    BayesTraitsV3 !{nx_tree} !{data} < cmd_Independent_MCMC$x.txt;
    '''
}


process BayesFactor {

    input: 
    path stones
    path nx_tree
    path binary

    output:
    tuple path(binary), path("BayesFactor_raw.txt")

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
    paste Dep_Indep_results.txt tmp names.txt > BayesFactor_raw.txt
    '''
}

process Conclude_BayesTraits {
    label 'python'
    
    input: 
    val bayes
    tuple path(binary), path(bayesfactor)
    path positions
    val branches

    output: 
    tuple path("bayes_detected_results.tsv"), path("bayes_tested_results.tsv")
    
    when : branches == "correlation"     
    
    //when branches_corr == "true" && branches_eem == "false"  //"correlation"
    shell: 
    '''
    bayes_traits_filter.py !{bayesfactor} !{binary} !{bayes} !{positions}
    '''

}

///////////END OF CORRELATION PART OF WORKFLOW. 

process Correlation {
    label 'python'
    publishDir "${resdir}", mode: 'copy'
    
    input: 
    val bayes
    tuple path(binary), path(bayesfactor)
    tuple path(detected), path(tested)
    
    output:
    tuple path("BayesFactor.txt"), path("tested_results.tsv"), path("significant_results.tsv")
    
    shell:
    '''
    merge_results.py !{bayesfactor} !{binary} !{bayes} !{tested}
    '''
}


workflow {

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

    // CONDOR Workflow 
    omodel = find_model(tree, align, model)
    omatrices = build_matrices(omodel, matrices)
    pastmlmatrix = omatrices[0]
    simulatormatrix = omatrices[1]
    length = info_align(align)
    nbseq = info_align_nbseq(align)
    //lenght and nb seqs in alignment (I think we can remove nb seqs)
    statsalign = length.combine(nbseq)

    otree = reoptimize_tree(freqmode, align, tree, omodel, matrices, statsalign)
    treechannel = otree[0]
    ratesparamchannel = otree[1]
    rateschannel = otree[2]
    freqchannel = otree[3]

    otreerename = tree_rename(treechannel, outgroup)
    namedtreechannel = otreerename[0]
    rootedtreechannel = otreerename[1]
    
    opars = pars_align_file(ratesparamchannel, positions, min_seq)
    pastmlalign = opars[0].transpose()
    positionschannel = opars[1]

    opastml = acr_pastml(pastmlalign, rootedtreechannel, positionschannel, pastmlmatrix)


    oprecount = pre_count(opastml, nb_simu)
    pythoncount = oprecount[0]
    rootseq = oprecount[1]
    marginalroot = oprecount[2]

    ocount = count_apparitions(min_eem, min_seq, align, pythoncount)
    refcounting = ocount[0]
    positionseem = ocount[1]
    positionseemfilter = ocount[2]
    refcounting.subscribe{rate, align, freqs, substitutions ->  freqs.copyTo(file("${resdir}").resolve('ref_substitutions.txt'));}

    simulchannel = simulator(positionseem.flatMap{it.readLines()}, simulatormatrix, rateschannel, freqchannel, namedtreechannel, marginalroot)
    collectsimu = simulchannel.collect()

    emergencechannel = conclude_convergence(simulatormatrix, align, freqchannel, refcounting, collectsimu, rootseq, positionseem, nb_simu, min_seq, min_eem, alpha, correction)
    emergencechannel.subscribe{detected, tested ->  tested.copyTo(file("${resdir}").resolve('tested_results.tsv')); detected.copyTo(file("${resdir}").resolve('significant_results.tsv'));}

    binarytraits = prepare_BT(align, positionseemfilter, phenotype, branches)

    nexustree = prepare_tree(rootedtreechannel, branches)

    bayestraitstmp = DataTraits(binarytraits)
    bayestraits = bayestraitstmp.flatten()

    stonestmp = BayesTraits(bayestraits, nexustree.first())
    stones = stonestmp.flatten().collect()

    bayeschannel = BayesFactor(stones, nexustree,binarytraits)

    obayestraits = Conclude_BayesTraits(bayes, bayeschannel, positionseemfilter, branches)
    obayestraits.subscribe{detected, tested ->  tested.copyTo(file("${resdir}").resolve('tested_results.tsv')); detected.copyTo(file("${resdir}").resolve('significant_results.tsv'));}

    occorelation = Correlation(bayes, bayeschannel, emergencechannel)
}
