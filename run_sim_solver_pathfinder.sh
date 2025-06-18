#!/bin/bash
mode=$2

. ${CONFIG:-$QDIR/sim_path_config_hifi.sh} $mode

eval $pathfinder $pathfinder_opts $1 2>$1.pf.err
