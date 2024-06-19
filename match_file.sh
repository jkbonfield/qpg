#!/bin/sh

# Compare a fasta file against a graph and training fasta.
# Usage: match_file.sh in.gfa train.fa query.fa
gfa=$1
train=$2
query=$3

ga=/nfs/sam_scratch/jkb/conda22/bin/GraphAligner

# Produce the training gaf
$ga -g $gfa -f $train -x vg -a train.gaf >/dev/null
./gaf2nodeseq.pl train.gaf $train $gfa > $gfa.nodeseq

# Compare each query against the nodeseq
samtools faidx $query
#for i in s1
for i in `cut -f 1 $query.fai`
do
    # Align
    samtools faidx $query $i > _.fa
    $ga -g $gfa -f _.fa -x vg -a query.gaf > /dev/null
    echo === $i `wc -l < query.gaf` ===

    # Reverse seq if reverse complemented as kmer2node is single-stranded atm
    if grep '<' query.gaf > /dev/null
    then
        samtools faidx -i $query $i|sed 's#/rc##' > _.fa
        $ga -g $gfa -f _.fa -x vg -a query.gaf > /dev/null
    fi

    # Produce kmer counts
    ./kmer2node $gfa.nodeseq _.fa | grep Node | sort -k2 -b | \

    # Report kmer2node counts along side GraphAligner counts
    perl -lane '
BEGIN {
    open(M,"query.gaf");
    while (<M>) {
        chomp($_);
        @F=split(/\s+/,$_);
        $F[5]=join(">",reverse split("<",$F[5]));
        if ($F[0] eq "'$i'") {
            foreach (split(">",$F[5])) {
                $h{$_}++;
            }
        }
    }
}

if ($F[-1] > 0 || $h{$F[1]}) {
    print "$_\t",$h{$F[1]}+0
}'
done
