CC=gcc

HTS=/nfs/users/nfs_j/jkb/work/samtools_master/htslib
INCLUDES=-I$(HTS)
LIBS=-L$(HTS) -Wl,-rpath,$(HTS) -lhts

CFLAGS=-g
LDFLAGS=

kmer2node: kmer2node.c
	$(CC) $(CFLAGS) $(LDFLAGS) $(INCLUDES) kmer2node.c -o kmer2node $(LIBS)
