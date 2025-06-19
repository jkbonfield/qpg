#!/bin/bash
mode=$2

. ${CONFIG:-$QDIR/config_illumina.sh} $mode

eval $pathfinder $pathfinder_opts $1 2>$1.pf.err
