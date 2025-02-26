#!/bin/bash

. $QDIR/sim_path_config_hifi.sh
$pathfinder_jkb -C40 $1 2>$1.pf.err
