#!/usr/bin/bash

# I should have done this entirely in Perl

for i in $@
do
    # Filename is of the form seq_[0-9]+[0-9]+3-#[0-9]#[0-0].*gfa
    base=`echo $i | sed 's/\(.*-#[0-9]#[0-9]\).*/\1/'`

    # Perform a GraphAligner alignment to get the (hopefully) correct weights
    cp $base _tmp$$.fa
    GraphAligner --min-alignment-score 90 --seeds-minimizer-ignore-frequent 1e-2 -g $i -f _tmp$$.fa -x vg -a _tmp$$.gaf 2>/dev/null >/dev/null
    perl -lane 'foreach (split(/[<>]/, $F[5])) {$c{$_}++ if $_ ne ""} END {foreach (sort keys %c) {print "$_\t$c{$_}"}}' _tmp$$.gaf |sort > _tmp$$.weights1

    # Execute the pathfinder script to obtain its weights instead
    eval pathfinder ${PATHFINDER_OPTS} $i 2>/dev/null \
	| awk '/SUBGRAPH/ {x=1;next} /PATH/ {x=0} x {print $3,$6}' \
	| sort > _tmp$$.weights2

    # Report the difference
    # TODO: consider making this weighted by node length
    dist=$(perl -e '
open(FH, "<", shift(@ARGV));while (<FH>) {chomp($_);@F=split(/\s+/,$_);$w1{$F[0]}=$F[1]};close(FH);
open(FH, "<", shift(@ARGV));while (<FH>) {chomp($_);@F=split(/\s+/,$_);$w2{$F[0]}=$F[1]};close(FH);
foreach (keys(%w1)) {$w2{$_}=0 unless exists($w2{$_})};
foreach (keys(%w2)) {$w1{$_}=0 unless exists($w1{$_})};
my $dist = 0;
foreach (keys(%w1)) { $dist += abs($w1{$_} - $w2{$_}); }
print $dist
' _tmp$$.weights1 _tmp$$.weights2)

    echo $base $dist

    rm _tmp$$.*
    #tag_gfa_copy_numbers.pl 
done
