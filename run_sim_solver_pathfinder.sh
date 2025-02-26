#!/bin/bash

. $QDIR/sim_path_config_hifi.sh

gfa=$1

$pathfinder_jkb -C40  $gfa 2>/dev/null
#  | tee $gfa.pf \
#  | sed -n '/PATH/,$p'|awk '/^\[/ {printf("%s ",$3)} END {print "\n"}'
