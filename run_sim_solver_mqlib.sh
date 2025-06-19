#!/bin/bash
mode=$2

. ${CONFIG:-$QDIR/config_illumina.sh} $mode

copy_numbers=$(perl -e '
use strict;
while (<>) {
    next unless /^S/;
    m/dc:f:([0-9.]*)/;
    print int($1/'$shred_depth' + .8), " ";
}
print "\n";' $1)

# Run mqlib
mqlib_oriented_qubo.py "foo" $1 "$copy_numbers" $timeout `dirname $1` | mqlib2path.pl

