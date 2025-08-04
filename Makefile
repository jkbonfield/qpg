CC=gcc

#HTS=/nfs/users/nfs_j/jkb/work/samtools_master/htslib
#INCLUDES=-I$(HTS)
#LIBS=-L$(HTS) -Wl,-rpath,$(HTS) -lhts
INCLUDES = $(shell pkg-config htslib --cflags)
LIBS = $(shell pkg-config htslib --libs) -lhts -lm

CFLAGS=-g -O2
LDFLAGS=

PROGS=kmer2node4 genome_create genome_create
all: $(PROGS)

clean:
	-rm $(PROGS)

kmer2node4: kmer2node4.c
	$(CC) $(CFLAGS) $(LDFLAGS) $(INCLUDES) kmer2node4.c -o $@ -lz $(LIBS)

genome_create: genome_create.c
	$(CC) $(CFLAGS) $(LDFLAGS) $(INCLUDES) $< -o $@ $(LIBS) -lm

test check: test_pf test_mqlib_km

test_pf: test_pf_km test_pf_mg test_pf_ga test_pf_vg

test_pf_km:
	@-mkdir test.out
	./run_gfa_sim.sh --solver pathfinder -c config_illumina_km.sh -n 5 -p test.out/_$@_ -s 104
	sed -n '/=== Summary/,$$p' test.out/_$@_00104/sim.out > test.out/$@.txt
	diff -h test/$@.txt test.out/$@.txt

test_pf_mg:
	@-mkdir test.out
	./run_gfa_sim.sh --solver pathfinder -c config_illumina_mg.sh -n 5 -p test.out/_$@_ -s 104
	sed -n '/=== Summary/,$$p' test.out/_$@_00104/sim.out > test.out/$@.txt
	diff -h test/$@.txt test.out/$@.txt

test_pf_ga:
	@-mkdir test.out
	./run_gfa_sim.sh --solver pathfinder -c config_illumina_ga.sh -n 5 -p test.out/_$@_ -s 104
	sed -n '/=== Summary/,$$p' test.out/_$@_00104/sim.out > test.out/$@.txt
	diff -h test/$@.txt test.out/$@.txt

test_pf_vg:
	@-mkdir test.out
	./run_gfa_sim.sh --solver pathfinder -c config_illumina_vg.sh -n 5 -p test.out/_$@_ -s 104
	sed -n '/=== Summary/,$$p' test.out/_$@_00104/sim.out > test.out/$@.txt
	diff -h test/$@.txt test.out/$@.txt

# No fixed seed, so just test it worked.
test_mqlib_km:
	@-mkdir test.out
	./run_gfa_sim.sh --solver mqlib --pathfinder -t 5 -j 1 --trim-edges -c config_illumina_km.sh -n 5 -p test.out/_$@_ -s 104
	sed -n '/=== Summary/,$$p' test.out/_$@_00104/sim.out > test.out/$@.txt
	awk '/^seq_/ && ($$4 == 0 || $$5 == 0) {n++} END {if (n) {exit 1}}' test.out/$@.txt

# No fixed seed, so just test it worked.
# NB: test not automatically called at the moment to avoid licensing issues.
test_gurobi_km:
	@-mkdir test.out
	./run_gfa_sim.sh --solver gurobi --pathfinder -t 5 -j 1 --trim-edges -c config_illumina_km.sh -n 5 -p test.out/_$@_ -s 104
	sed -n '/=== Summary/,$$p' test.out/_$@_00104/sim.out > test.out/$@.txt
	awk '/^seq_/ && ($$4 == 0 || $$5 == 0) {n++} END {if (n) {exit 1}}' test.out/$@.txt

