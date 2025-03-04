#!/usr/bin/perl -w
use strict;

# Turns our mqlib path output to the same format as a pathfinder path output
# eg [(0, 's36_-'), ... (283,'end')] to
# 
# PATH
# [ 1] s36+

while (<>) {
    next unless /^Path/;
    s/.*\[//;
    s/].*//;

    print "PATH\n";
    foreach (split(/\),/, $_)) {
	next if (/end/);
	m/(\d+),[ ']*(\S+)_([-+])/;
	printf("[%5d] %s%s\n", $1, $2, $3);
    }
}
