Instructions
============

This is very much a work in progress and a personal exploration into
how pangenome alignment could work.  It is not to be considered a
robust and production ready suite of tools.

The main entry to running the simulations is the run_gfa_sim.sh
script.  The usage is

    Usage: run_gfa_sim.sh [options] [seed solver [out_prefix]]
    Options:
        -c,--config    FILE    Use FILE as configuration
        -s,--seed      INT     Specify random number seed [1]
        -p,--prefix    STR     Use STR as the output dir prefix [sim_]
           --solver    STR     Specify the solver [pathfinder]
        -a,--annotate  STR     GFA node weight algorithm [km]
           --shred_len INT     Shotgun read length
           --shred_err FLOAT   Shotgun read error rate (fraction)
        -t,--times     INT_LIST   Time limits provided to QUBO solvers
        -j,--jobs      INT     Number of runs of QUBO solvers
        -n,--training  INT     Number of strings to use as training set [10]
           --edge2node         Use edge2node version
           --trim-edges        Use trim_edges.pl
           --pathfinder        Use pathfinder to get subgraphs

Solver can be mqlib or pathfinder.

It uses a configuration file to control parameters such as the graph
complexity, but also which algorithms to use.  Use the CONFIG environment
variable or run_gfa_sim -c FILE to point to one of these configurations.
Premade configuration files are:

    config_base.sh                      Included by all other config files
    config_hifi.sh                      2kbp  long reads at 0.001 error
    config_illumina.sh                  200bp long reads at 0.001 error    
    config_hifi_{km,mg,ga,vg}.sh        Hifi using a specific --annotate option
    config_illumina_{km,mg,ga,vg}.sh    Ilumina using a specific --annotate opt

For example:

    run_gfa_sim.sh -c config_hifi_mg.sh -s 1 --solver pathfinder -p pf_hifi_mg_


Multiple runs can be launch via xargs:

    ( p=illumina; m=ga; seq 100 150 | xargs -I % -P 4 ./run_gfa_sim.sh -c config_${p}_$m.sh --pathfinder --solver mqlib --trim-edges -t 30 -j 1 -n 5 -p ~/lustre/tmp/mqlib-tmq10_${p}_${m}_ % )

Summaries from multiple runs with different solvers and annotations can be
produced with:

    ( for x in pf2-base sa-base ma-base mqlib-base mqlib-trim1.3;do for m in sa ma km mg ga vg;do p=illumina; if [ ! -e ~/lustre/tmp/${x}_${p}_${m}_00100 ];then continue;fi; printf "%-13s %3s " $x $m; awk '!/contig/ {for (i=2;i<11;i++) {a[i]+=$i}n++} END {for (i=2; i<11; i++) {printf("%7.1f ", a[i]/n)} print n}' ~/lustre/tmp/${x}_${p}_${m}_001*/*eval_cons* 2>/dev/null || echo;done;echo;done)


Other selected programs
=======================

run_syncasm_sim.sh
------------------

De novo assembly via SyncAsm.

run_miniasm_sim.sh
------------------

De novo assembly via Miniasm.


run_sim_create_gfa.sh
---------------------

Creates a synthetic population (via genome_create) and trains a GFA on
a subset of it (fofn.test).  Also produces a fofn.test list of sequences to
evaluate.

run_sim_add_gfa_weights_${annotate}.sh
--------------------------------------

A series of scripts to annotate the node weights in the population GFA by
aligning the fofn.test sequences.

Selected via the run_gfa_sim --annotate option.
This produces a series of seq*.gfa files with the annotated test sequence
graphs.

run_sim_solver_${solver}.sh
---------------------------

A series of scripts to solve the path.


gaf2nodeseq.pl
--------------

Merges a GFA graph, a GAF alignment file and the input FASTA used in
generating the GAF to identify the parts of the sequences that were aligned
against each GFA file (using the path and CIGAR fields).

These sequences are then written to a "nodeseq" file, which constitutes an
@node name and a series of sequences.  The first sequence is from the graph
while all subsequent sequences are from the input fasta (with KMER-1 prior
bases of context).

This is called by the run_sim_create_gfa.sh script if the solver is kmer2node.


kmer2node4
----------

A program which builds a kmer index from a nodeseq file and then compares
kmers in a set of query sequences to identify the depth of coverage for all
nodes in the graph.  This forms a rudimentary alignment.

The program produces a lot of debugging output, but piping the output to "grep
Node" will return the final answers.

See also merge_kmer2node.pl to allow running kmer2node4 with multiple kmers
and merge the results to form a single set of node weights.


pathfinder2seq.pl
-----------------

With an input GFA and pathfinder PATH information this creates a sequence by
concatenating (and complementing) nodes together to produce a candidate
assembly sequence.


run_sim_evaluate_path.sh
------------------------

Compares a true sequence to a candidate assembly sequence and reports how well
they match.  It compares A to B and B to A to get symmetric data on coverage
(how much of A is in B and vice versa).
