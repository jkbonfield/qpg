#!/bin/sh

# Wraps up candidate_stats.pl in a two way fashion
# (real vs eval and eval vs real).

echo "Per contig: (rname      len-t   len-q   covered %used   ncontig nbreaks nindel  ndiff   %identity)"
t=$1
q=$2

set -- $(candidate_stats.pl $q $t | tail -n +2 | perl -lane '
$len += $F[1];
$cov += $F[1] * $F[2];
END { printf("%d %.2f%%", $len, $cov/$len); }')
len=$1
cov=$2

set -- `candidate_stats.pl $t $q | tail -1`
printf "%-20s\t%d\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%s\n" $1 $2 $len $3 $cov $4 $5 $6 $7 $8
