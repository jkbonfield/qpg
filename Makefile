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
