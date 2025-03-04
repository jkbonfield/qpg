#!/bin/bash

. ${CONFIG:-$QDIR/sim_path_config_hifi.sh}

eval $pathfinder $pathfinder_opts $1 2>$1.pf.err
