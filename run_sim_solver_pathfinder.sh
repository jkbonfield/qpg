#!/bin/bash

. ${CONFIG:-$QDIR/config_illumina.sh}

eval $pathfinder $pathfinder_opts $1 2>$1.pf.err
