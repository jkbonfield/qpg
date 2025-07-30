#!/bin/sh

gfa=$1
fasta=$2
shift 2
samtools faidx $fasta

qdir=${QDIR:-.}

# Index our query fastas
for q in $@
do
    bwa index $q 2>/dev/null
done


# Split the GFA on bubbles
rm $gfa.split*
$qdir/pangene.js call $gfa | ../bubble_split.pl $gfa $gfa.split

for i in $gfa.split.*
#for i in early.gfa.split.s91.s106
do
    # Foreach split, pull out the region in the original sequence,
    echo
    printf "%-35s" $i
    reg=`$qdir/gfa_region.pl $i`
    echo $reg
    samtools faidx $fasta $reg > _.fa
    len=`samtools fasta _.fa 2>/dev/null |awk '!/>/ {print length($0)}'`
    min_score=`expr $len "*" 3 / 10`

    # Align that region of fasta against each other fasta to create a series
    # of smaller bubble-specific problems.
    bwa index _.fa 2>/dev/null
    for q in $@
    do
	echo QUERY: $q
	bwa mem -T$min_score $q _.fa 2>/dev/null| 
	    awk '!/^@/ {print $6;last}'
    done
done

exit

_O.fa:83472-84729
_O.fa.245.fa.250.fa.263.fa.296.fa.392.fa

89954-90635
96754-97329
62731-62928


Or reverse bwa query & target
89955+680
96755+684
