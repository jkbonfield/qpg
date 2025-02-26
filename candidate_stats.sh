#!/bin/sh

# Wraps up candidate_stats.pl in a two way fashion
# (real vs eval and eval vs real).

echo "Per contig: (rname      len-t   len-q   covered %used   ncontig nbreaks nindel  ndiff   %identity)"
t=$1
q=$2
set -- `candidate_stats.pl $q $t | tail -1`
len=$2
cov=$3
set -- `candidate_stats.pl $t $q | tail -1`
printf "%-20s\t%d\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%s\n" $1 $2 $len $3 $cov $4 $5 $6 $7 $8
