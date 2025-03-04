genome_len=1000000

shred_len=2000
shred_err=0.001
shred_depth=30

mdbg_opts="-k 15 --density 0.01 -l 10 --minabund 4"
minimap2_opts="-k19 -Xw5 -m100"
miniasm_opts="-s1500 -F0.9 -c3 -m500 -e10"
syncasm_kmer=401
kmer2node_kmer=${KMER:-51}

pathfinder_cz3=/lustre/scratch127/qpg/cz3/QuantumTangle/pathfinder/pathfinder
pathfinder_jkb=/lustre/scratch125/ids/team117/sam/jkb/quantum/pathfinder/pathfinder
pathfinder=$pathfinder_jkb
pathfinder_opts="-a"

#mqlib timeout
timeout=30
