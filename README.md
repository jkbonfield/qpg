Musings on a PanGenome world
============================

This is very much a work in progress and a personal exploration into
how pangenome alignment could work.

The main software here is:

gaf2nodeseq.pl
--------------

Merges a GFA graph, a GAF alignment file and the input FASTA used in
generating the GAF to identify the parts of the sequences that were
aligned against each GFA file (using the path and CIGAR fields).

These sequences are then written to a "nodeseq" file, which
constitutes an @node name and a series of sequences.  The first
sequence is from the graph while all subsequent sequences are from the
input fasta (with KMER-1 prior bases of context).


kmer2node
---------

A program which builds a kmer index from a nodeseq file and then
compares kmers in a set of query sequences to identify the depth of
coverage for all nodes in the graph.

This forms an extremely rudimentary alignment.

The program produces a lot of debugging output, but piping the output
to "grep Node" will return the final answers.  A CSV file suitable for
loading into Bandage to plot depth can be produced with:

awk '{BEGIN {print "Node name,Freq"} /Node/ {printf("%s,%s\n",$2,$NF)}'

kmer2node2
----------

This has a -k INT option to specify the kmer size, which can be any
size unlike the kmer2node which is fixed hashing using memory 4^K (so
has a compile-time default of 14).

Note: remember to also specify the option last argument to
gaf2nodeseq.pl specifying the same kmer size so the indexing has
enough pre-node preamble to go on, particularly if you have a graph
with short nodes.


match_file.sh
-------------

A basic pipeline that aligns training data against the GFA, produces
the nodeseq file, runs kmer2node on query sequences against that, also
runs GraphAligner on the same query sequences, and finally summarises
node coverage from both kmer2node and GraphAligner >path counts.

As an example:

```
$ ./match_file.sh drb1+tangle1.gfa DRB-train.fa JAGYVH010000080.fa
=== JAGYVH010000080 1 ===
Node      aa-c1	len    603	exp  603.0	hit    587+1	ratio 1.00	1
Node    ae-a1c0	len     49	exp   49.0	hit     42+2	ratio 0.86	1
Node    af-a1i0	len     10	exp   10.0	hit      1+0	ratio 0.10	0
Node    ag-a1c1	len    211	exp  202.4	hit    190+1066	ratio 0.94	1
Node    ah-a1i1	len     47	exp   20.6	hit      5+42	ratio 0.24	1
Node    ai-a1c2	len    151	exp  136.9	hit    112+353	ratio 0.82	1
Node   aj-a1i2b	len    330	exp  143.3	hit      5+0	ratio 0.02	0
Node    ak-a1c3	len     65	exp   65.0	hit     62+0	ratio 0.95	1
Node      ba-c2	len   1967	exp 1938.4	hit   1911+211	ratio 0.99	1
Node      bb-i2	len     20	exp   20.0	hit     20+0	ratio 1.00	1
Node      ca-c3	len   1126	exp 1105.2	hit   1098+282	ratio 0.99	1
Node      da-c4	len    650	exp  648.7	hit    635+10	ratio 0.98	1
Node     db-rCT	len      2	exp    1.9	hit     10+0	ratio 4.98	8
Node     dc-rCA	len      2	exp    2.0	hit     25+0	ratio 12.44	20
Node      ea-c5	len    287	exp  286.7	hit    256+73	ratio 0.89	1
Node     ea-d5b	len     25	exp   25.0	hit     12+4	ratio 0.48	1
Node      ea-e5	len     59	exp   58.7	hit     34+19	ratio 0.58	1
Node      eb-i5	len     56	exp   56.0	hit     56+0	ratio 1.00	1
Node      fa-c6	len    258	exp  257.2	hit    258+2	ratio 1.00	1
Node      ga-c7	len    240	exp  239.2	hit    237+3	ratio 0.99	1
Node     la-c12	len    109	exp  108.5	hit     98+3	ratio 0.90	1
Node     ma-c13	len     81	exp   81.0	hit     54+0	ratio 0.67	1
Node     oa-c15	len    882	exp  881.4	hit    830+26	ratio 0.94	1
Node     ob-r16	len     16	exp   15.9	hit     47+1	ratio 2.96	3
Node     pa-c17	len   1008	exp 1005.7	hit    978+8	ratio 0.97	1
Node     pb-cnv	len    282	exp  212.1	hit    166+816	ratio 0.78	1
Node  pc-cnvstr	len      4	exp    2.8	hit     30+18	ratio 10.52	12
Node     qa-c19	len   1533	exp 1495.8	hit   1490+113	ratio 1.00	1
Node     ra-c20	len    886	exp  883.1	hit    882+3	ratio 1.00	1
Node     sa-c21	len    513	exp  508.9	hit      0+0	ratio 0.00	1
```

gaf2seq.pl
----------

A debugging tool which takes a GFA file plus GAF alignment path and
reports the sequence constructed from the GFA nodes.  This should be
the sequence from which the CIGAR string is applied.  It can be
compared to the fasta sequence (eg with dotter) to inspect how good
the alignment looks and whether there are nodes missing in the graph.


train8.fa
---------

A set of 8 sequences that jointly cover all nodes in the
drb1+tangle1.gfa graph.  These can be considered to be a minimum
training set for producing the nodeseq file.
