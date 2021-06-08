
path_file = "/pasteur/zeus/projets/p01/Evolbioinfo/users/mamorel/Projet_Convergence/Data/"

params.align = path_file + "args/processed/naturepaper/B.pol.aa.woutgroup.norecomb.monophyletic.fasta"
params.tree = path_file + "trees/processed/naturepaper/root.B.pol.aa.woutgroup.norecomb.monophyletic.MLfreq.treefile"
params.outgroup = path_file+"args/processed/naturepaper/outgroup.losalamos.names"
params.resdir=path_file+"results/naturepaper/test_HIV/" //do not forget to change the resdir
params.model = "best" // best: will choose the best model, other wise will take the given model
params.matrices = "$baseDir/assets/protein_model.txt"

params.nb_simu = 10000
params.min_seq = 12

align = file(params.align)
tree = file(params.tree)
outgroup = file(params.outgroup)
matrices = file(params.matrices)
resdir=file(params.resdir)
resdir.with {mkdirs()}

nb_simu = params.nb_simu
min_seq = params.min_seq
model = params.model

//nb seqs and length align
//run iqtree model finder
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
    '''
    iqtree -m MFP -nt AUTO -s !{align} -te !{tree} -T !{task.cpus}
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
    '''
    matrix_name=`awk 'BEGIN{FS="+"} {print toupper($1) }' !{iqtree_modelrate}`
    matrix_pastml_format.py ${matrix_name} !{matrices}
    '''
}


process info_align {
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


//reoptimize tree branch lengths and estimate rates and frequencies by ML 
process reoptimize_tree {
    label 'iqtree'

    input:
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
//fixed tree no tree search performed
//should work with no F or I 
    '''
    sed -i '/+F/!s/+/+F+/' !{iqtree_mode}
    sed -i 's/+F/+FO/' !{iqtree_mode}
    iqtreemode=`cat !{iqtree_mode}`
    iqtree -m ${iqtreemode} -nt AUTO -s !{align} -te !{tree} -wsr -pre align
    
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
    file "rooted_*" into RootedtreeChannel
    //Remove the branch length for root

    shell:
    ''' 
    sed -i  '/);/!s/)[0-9]*.[0-9]*;/);/' !{tree}
    gotree reroot outgroup -i !{tree} -l !{outgroup} -o rooted_!{tree}
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
    
    //rate sed numerotation from 1
    //input numerotation from 0
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
    '''
    pastml_fasta.py !{pastml_acr} !{positions} !{length} !{marginal_root} !{nb_simu} test_
    '''
    //gunzip -c '1.out.gz' > pastml_output
    //for i in {2..!{length}};  do gunzip -c ${i}.out.gz | cut -f 2 > temp && paste pastml_output temp > test2 && mv test2 pastml_output ;  done
    //rm temp

}

//my simulator of sequence in python, using root, nb simulations and the ROOTED tree.
process simulator {
    //errorStrategy 'retry'
    //maxRetries 3
    label 'python'

    publishDir "${resdir}", pattern: "count*.tsv.gz", mode: 'copy'
    input : 
    each x from ListPositionsChannel.readLines()
    file simulation_model from SimulatorMatrix
    file (rates) from RatesChannel
    file (freq) from FreqChannel
    file (named_tree) from SimulatorChannel
    file (root) from Marginal_root

    output : 
    //tuple val(x), file("count*.tsv.gz") into MysimulationsChannel.collect() // We take all the counts into a single Channel
    file("count*npz") into MysimulationsChannel
    shell:
    '''
    sed -n "/^$((!{x}-1))\t/p" !{root} | cut -f 2 > root.txt
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
    '''
    count_substitutions_from_tips.py !{align} !{tree} !{positions}
    '''
}

Subscribe_matrices.subscribe{positions, rate, align, freqs, substitutions ->  freqs.copyTo(file("${resdir}").resolve('ref_substitutions.txt'));}
 
// should 
process conclude_convergence{
    label 'python'

    publishDir "${resdir}", mode: 'copy'
    
    input: 
    file simulation_model from SimulatorMatrix
    file align 
    file (freq) from FrequenciesChannel
    tuple file(positions),file(rate), file(acralign), file(ref_matrix), file(substitutions) from Ref_couting //rates for all positions
    file (counts) from Collect_simulations
    file(root) from Root_seq //only the interesting positions
    val nb_simu
    val min_seq

    //named*.phy
     
    output:
    tuple file("detected_metrics.tsv"), file("all_results_metrics.tsv") 

    shell:
    '''
    R=`cat !{root}`
    convergent_substitutions_pvalue.py !{positions} ${R} !{rate} !{align} !{ref_matrix} !{substitutions} !{nb_simu} !{freq} !{simulation_model} !{min_seq}
    '''
}
