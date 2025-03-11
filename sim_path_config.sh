genome_len=1000000; # NB: unused by run_sim_gfa.sh

# Challenging: 400-800 nodes
#genome_opts="-P 50 -l 20000  -S 0.01 -C 0.001 -N 0.002 -n 0.002 -A 0.001 -L 0.0002 -T 0.0004"

# Complex: 100-200 nodes
#genome_opts="-P 50 -l 10000  -S 0.01 -C 0.0005 -N 0.005 -n 0.005 -A 0.0005 -L 0.0001 -T 0.00005 -I 0.0001"

# Medium: 40-80 nodes
genome_opts="-P 50 -l 5000 -S 0.001 -C 0.0003 -N 0.01 -n 0.01 -A 0.0005 -L 0.0001 -T 0.00005 -I 0.0001"

# Simple: 10-50 nodes
#genome_opts="-P 50 -l 5000 -S 0.001 -C 0.0002 -N 0.01 -n 0.01 -A 0.00005 -L 0.00001 -T 0.0002 -I 0.0001"

shred_len=${SHRED_LEN:-2000}
shred_err=${SHRED_ERR:-0.001}
shred_depth=${SHRED_DEPTH:-30}

mdbg_opts="-k 15 --density 0.01 -l 10 --minabund 4"
minimap2_opts="-k19 -Xw5 -m100"
miniasm_opts="-s1500 -F0.9 -c3 -m500 -e10"
syncasm_kmer=401
kmer2node_kmer=51

# Use minigraph
use_mg=1

pathfinder_cz3=/lustre/scratch127/qpg/cz3/QuantumTangle/pathfinder/pathfinder
pathfinder_jkb=/lustre/scratch125/ids/team117/sam/jkb/quantum/pathfinder/pathfinder
pathfinder=${PATHFINDER:-$pathfinder_jkb}
pathfinder_opts=${PATHFINDER_OPTS:--C40}

#mqlib timeout
timeout=${MQLIB_TIMEOUT:-10}

