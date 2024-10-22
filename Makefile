CC=gcc

#HTS=/nfs/users/nfs_j/jkb/work/samtools_master/htslib
#INCLUDES=-I$(HTS)
#LIBS=-L$(HTS) -Wl,-rpath,$(HTS) -lhts
INCLUDES = $(shell pkg-config htslib --cflags)
LIBS = $(shell pkg-config htslib --libs) -lhts

CFLAGS=-g
LDFLAGS=

PROGS=kmer2node kmer2node2 genome_create genome_create2
all: $(PROGS)

clean:
	-rm $(PROGS)

kmer2node: kmer2node.c buzhash.h
	$(CC) $(CFLAGS) $(LDFLAGS) $(INCLUDES) kmer2node.c -o $@ $(LIBS)

kmer2node2: kmer2node2.c buzhash.h
	$(CC) $(CFLAGS) $(LDFLAGS) $(INCLUDES) kmer2node2.c -o $@ $(LIBS)

genome_create: genome_create.c
	$(CC) $(CFLAGS) $(LDFLAGS) $(INCLUDES) $< -o $@ $(LIBS) -lm

genome_create2: genome_create2.c
	$(CC) $(CFLAGS) $(LDFLAGS) $(INCLUDES) $< -o $@ $(LIBS) -lm
