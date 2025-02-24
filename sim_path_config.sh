genome_len=1000000

shred_len=2000
shred_err=0.02
shred_depth=30

minimap2_opts="-k19 -Xw5 -m100"
miniasm_opts="-s1500 -F0.9 -c3 -m500 -e10"
syncasm_kmer=401
kmer2node_kmer=51

pathfinder=/lustre/scratch127/qpg/cz3/QuantumTangle/pathfinder/pathfinder
#pathfinder=/lustre/scratch125/ids/team117/sam/jkb/quantum/pathfinder/pathfinder
#pathfinder_opts="-a"
