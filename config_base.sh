genome_len=1000000; # NB: unused by run_sim_gfa.sh

#OLD Challenging: 400-800 nodes
#OLD genome_opts="-P 50 -l 20000  -S 0.01 -C 0.001 -N 0.002 -n 0.002 -A 0.001 -L 0.0002 -T 0.0004"

#OLD Complex: 100-200 nodes
#OLD genome_opts="-P 50 -l 10000  -S 0.01 -C 0.0005 -N 0.005 -n 0.005 -A 0.0005 -L 0.0001 -T 0.00005 -I 0.0001"

#OLD Medium: 40-80 nodes
#OLD genome_opts="-P 50 -l 5000 -S 0.001 -C 0.0003 -N 0.01 -n 0.01 -A 0.0005 -L 0.0001 -T 0.00005 -I 0.0001"

# Testing for high connectivity: 100 nodes, average 100 (???)
# genome_opts="-S 0.006 -C 0.0015 -N 0.01 -n 0.01 -A 0.0005 -L 0.0001 -T 0.005 -I 0.010 -l 10000 -P 100 -G 10"

# Challenging: 46-198 nodes, average 108
# genome_opts="-S 0.006 -C 0.0015 -N 0.01 -n 0.01 -A 0.0005 -L 0.0001 -T 0.0001 -I 0.0002 -l 10000 -P 100 -G 10"

# Medium: 25-150 nodes, average 79
# genome_opts="-S 0.005 -C 0.001 -N 0.01 -n 0.01 -A 0.0005 -L 0.0001 -T 0.0001 -I 0.0002 -l 10000 -P 100 -G 10"

#OLD Simple: 10-50 nodes
#OLD genome_opts="-P 50 -l 5000 -S 0.001 -C 0.0002 -N 0.01 -n 0.01 -A 0.00005 -L 0.00001 -T 0.0002 -I 0.0001"

# Simple: 12-76 nodes, average 44
genome_opts="-S 0.004 -C 0.0005 -N 0.01 -n 0.01 -A 0.0005 -L 0.0001 -T 0.0001 -I 0.0002 -l 10000 -P 100 -G 10"

# Very simple: 4-38 nodes, average 16
# genome_opts="-S 0.001 -C 0.0002 -N 0.01 -n 0.01 -A 0.0005 -L 0.0001 -T 0.0001 -I 0.0002 -l 10000 -P 100 -G 10"

shred_len=${SHRED_LEN:-2000}
shred_err=${SHRED_ERR:-0.001}
shred_depth=${SHRED_DEPTH:-30}

mdbg_opts="-k 15 --density 0.01 -l 10 --minabund 4"
minimap2_opts="-k19 -Xw5 -m100"
miniasm_opts="-s1500 -F0.9 -c3 -m500 -e10"
syncasm_kmer=401
kmer2node_kmer=51

# Annotation method
annotate=km; # kmer2node

#pathfinder_cz3=/lustre/scratch127/qpg/cz3/QuantumTangle/pathfinder/pathfinder
#pathfinder_jkb=/lustre/scratch125/ids/team117/sam/jkb/quantum/pathfinder/pathfinder
#pathfinder=${PATHFINDER:-$pathfinder_jkb}
pathfinder=/software/badger/opt/pangenomes/bin/pathfinder
#pathfinder_opts=${PATHFINDER_OPTS:--X40}

#Illumina
#ga6	9984.7	9870.1	94.6	95.3	1.4	1.9	1.6	0.3	96.7
#mg6	9984.7	9381.0	87.6	92.8	4.1	2.0	0.9	0.0	97.3
#km	9984.7	9180.9	88.5	95.2	3.2	1.3	1.1	0.1	97.1 old k2n
#km2	9984.7	9521.0	89.9	93.6	2.4	2.4	1.5	0.2	96.8 -m-s-cov 2
#km3	9984.7	10183.4	93.6	92.3	1.0	2.8	2.0	0.4	96.5 -m-s-cov 1
#km3	9984.7	10183.4	93.6	92.3	1.0	2.8	2.0	0.4	96.5 " -neighs 1
#km4    9984.7	9588.4	90.1	93.5	2.2	2.6	1.6	0.2	96.7 " " -c0
#km5    9984.7	9641.5	90.7	93.6	2.2	2.6	1.6	0.2	96.7 """ -bub-check
#km6    9984.7	9619.9	90.9	93.7	2.2	2.5	1.6	0.2	96.7 """" -X50
#km7    9984.7	9584.1	91.0	94.5	2.1	2.0	1.6	0.2	96.8 km6 r=r2 if 0
#km8    9984.7	8907.5	85.6	94.8	4.4	1.6	1.0	0.0	97.3 no -U
#km8    9984.7	9663.0	91.5	94.5	2.1	1.9	1.7	0.2	96.8 -U at 35,20
#km8    9984.7	9712.0	92.0	94.6	2.1	1.8	1.6	0.2	96.8 -U at 75,35,20
#km9    9984.7	9740.1	91.9	93.9	2.2	2.0	1.6	0.3	96.6 km8 poss+=.3
#km9    9984.7	9605.6	91.5	95.1	2.4	1.7	1.3	0.1	96.9 km8 poss+=.01
#km9    9984.7	9508.8	91.0	95.2	2.6	1.5	1.4	0.1	96.9 km8 count+=.01
#km9    9984.7	9508.8	91.0	95.2	2.6	1.5	1.4	0.1	96.9 also +=.1 +=.5
#kmA    9984.7	9297.5	89.7	95.6	2.9	1.3	1.2	0.1	97.0 " cnt+=.1 msc 2
#km2    9984.7	9335.3	90.0	95.6	2.9	1.3	1.2	0.1	97.0 " cnt+=.1 msc 2
#km2.50 9984.7	9335.3	90.0	95.6	2.9	1.3	1.2	0.1	97.0 ""
#km3.50 9984.7	10174.2	94.3	93.2	1.0	2.5	1.9	0.3	96.6 ""
#kmorig 9984.7	9486.7	90.6	95.3	2.6	1.5	1.3	0.1	96.9 with km9 prfle

#kmorig-with-km9 shows most gain is the options rather than kmer2node4 changes,
#but there's still some benefit to be had there.

#pathfinder_opts=${PATHFINDER_OPTS:--X40 --min-seq-cov 2}; # km2
#pathfinder_opts=${PATHFINDER_OPTS:--X50 --min-seq-cov 1}; # km3.50
#pathfinder_opts=${PATHFINDER_OPTS:--X40 --min-seq-cov 1}; # km3
#pathfinder_opts=${PATHFINDER_OPTS:--X40 --min-seq-cov 1 --neighbour-steps 1}; # km3
#pathfinder_opts=${PATHFINDER_OPTS:--X40 -c0 --min-seq-cov 1 --neighbour-steps 1}; # km4
#pathfinder_opts=${PATHFINDER_OPTS:--X40 -c0 --min-seq-cov 1 --neighbour-steps 1 --bub-check}; # km5
#pathfinder_opts=${PATHFINDER_OPTS:--X50 -c0 --min-seq-cov 1 --neighbour-steps 1 -v 3 --bub-check}; # km6, km7, km8, km9
pathfinder_opts=${PATHFINDER_OPTS:--X50 -c0 --min-seq-cov 1 --neighbour-steps 1 -v 3}
#--edge-to-seq appears to change nothing now.

#mqlib timeout
timeout=${MQLIB_TIMEOUT:-10}

