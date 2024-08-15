CC=gcc

#HTS=/nfs/users/nfs_j/jkb/work/samtools_master/htslib
#INCLUDES=-I$(HTS)
#LIBS=-L$(HTS) -Wl,-rpath,$(HTS) -lhts
INCLUDES = $(shell pkg-config htslib --cflags)
LIBS = $(shell pkg-config htslib --libs) -lhts

CFLAGS=-g
LDFLAGS=

PROGS=kmer2node genome_create
all: $(PROGS)

clean:
	-rm $(PROGS)

kmer2node: kmer2node.c
	$(CC) $(CFLAGS) $(LDFLAGS) $(INCLUDES) kmer2node.c -o kmer2node $(LIBS)

genome_create: genome_create.c
	$(CC) $(CFLAGS) $(LDFLAGS) $(INCLUDES) $< -o $@ $(LIBS) -lm
