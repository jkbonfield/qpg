#include <stdio.h>
#include <unistd.h>
#include <assert.h>

#include <htslib/sam.h>
#include <htslib/khash.h>
#include <htslib/kbitset.h>

#include "buzhash.h"

// kmer2node [-k kmersize] graph.nodeseq in.fasta

// nodeseq is a file consisting of "@name" and then one line per
// sequence that has been mapped to that node during a training run.
// The first sequence line is important as it is the graph node seq
// itself. We has this with the others, but the length matters as it
// determines the expected number of times the node occurs.

// Output is a list of node names and a weighting reflecting the
// number of relative kmers observed, normalised for the number we'd
// expect given the node length.

// kmer is the kmer size itself.
static int kmer = 50;
// kmer_idx is the size of kmer used in the nodeseq file
static int kmer_idx = 50;

// KMER_IDX is the hash table size for storing kmer entries.
#ifndef KMER_IDX
#define KMER_IDX 28
#endif

#define KSIZE (1<<KMER_IDX)
#define KMASK (KSIZE-1)

#define MAX_NODES (100*1024)

static int ndup = 0, nuniq = 0;

/* ----------------------------------------------------------------------
 * A node has a name, a length, and a bitvector of which kmers are present.
 */
typedef struct {
    int num; // Nth node number
    char *name;
    int length;

    // Maybe also an expected kmer count, as a fraction so small nodes will be
    // partial kmer.
    int kmer_unique, kmer_dup; // how many indexed uniquely/dups

    // Possible hits as we didn't switch node inbetween.
    // This is an upper-bound on the counting for cycles.
    int hit_possible, hit_count;
} node;

void node_free(node *n) {
    free(n->name);
}

node *node_create(char *name) {
    node *n = calloc(1, sizeof(*n));
    if (!n)
	return NULL;

    n->name = strdup(name);
    if (!n->name) {
	node_free(n);
	return NULL;
    }
    
    return n;
}

static int base4[256] = {0}, lookup_done=0;
static inline void base4_init(void) {
    if (!lookup_done) {
	lookup_done = 1;
	base4['a']=base4['A']=0;
	base4['c']=base4['C']=1;
	base4['g']=base4['G']=2;
	base4['t']=base4['T']=3;
    }
}

// From htslib/sam_internal.h
static inline void nibble2base(uint8_t *nib, char *seq, int len) {
    static const char code2base[512] =
        "===A=C=M=G=R=S=V=T=W=Y=H=K=D=B=N"
        "A=AAACAMAGARASAVATAWAYAHAKADABAN"
        "C=CACCCMCGCRCSCVCTCWCYCHCKCDCBCN"
        "M=MAMCMMMGMRMSMVMTMWMYMHMKMDMBMN"
        "G=GAGCGMGGGRGSGVGTGWGYGHGKGDGBGN"
        "R=RARCRMRGRRRSRVRTRWRYRHRKRDRBRN"
        "S=SASCSMSGSRSSSVSTSWSYSHSKSDSBSN"
        "V=VAVCVMVGVRVSVVVTVWVYVHVKVDVBVN"
        "T=TATCTMTGTRTSTVTTTWTYTHTKTDTBTN"
        "W=WAWCWMWGWRWSWVWTWWWYWHWKWDWBWN"
        "Y=YAYCYMYGYRYSYVYTYWYYYHYKYDYBYN"
        "H=HAHCHMHGHRHSHVHTHWHYHHHKHDHBHN"
        "K=KAKCKMKGKRKSKVKTKWKYKHKKKDKBKN"
        "D=DADCDMDGDRDSDVDTDWDYDHDKDDDBDN"
        "B=BABCBMBGBRBSBVBTBWBYBHBKBDBBBN"
        "N=NANCNMNGNRNSNVNTNWNYNHNKNDNBNN";

    int i, len2 = len/2;
    seq[0] = 0;

    for (i = 0; i < len2; i++)
        // Note size_t cast helps gcc optimiser.
        memcpy(&seq[i*2], &code2base[(size_t)nib[i]*2], 2);

    if ((i *= 2) < len)
        seq[i] = seq_nt16_str[bam_seqi(nib, i)];
}


/* ----------------------------------------------------------------------
 * A nodeset is a collection of nodes.  It is a node map indexed by name.
 */
KHASH_MAP_INIT_STR(node, node*)

typedef struct {
    khash_t(node) *nodes;
    int nnodes;
    int kmer[KSIZE];   // maps hash(kmer) to node number
    int unique[KSIZE]; // number of unique maps for this kmer
    node *num2node[MAX_NODES];
} nodeset;

nodeset *nodeset_create(void) {
    nodeset *ns = calloc(1, sizeof(*ns));
    if (!ns)
	return NULL;

    ns->nodes = kh_init(node);

    return ns;
}

void nodeset_free(nodeset *ns) {
    for (int i = 1; i <= ns->nnodes; i++)
	node_free(ns->num2node[i]);

    kh_destroy(node, ns->nodes);
}

node *nodeset_find(nodeset *ns, char *name, int create) {
    //printf("Processing %s", name);

    khiter_t k = kh_get(node, ns->nodes, name);
    if (k == kh_end(ns->nodes)) {
	if (!create) {
	    fprintf(stderr, "Could find node '%s'\n", name);
	    return NULL;
	}

	node *n = node_create(name);
	if (!n)
	    return NULL;

	n->num = ++ns->nnodes;
	ns->num2node[ns->nnodes] = n;

	int ret;
	k = kh_put(node, ns->nodes, strdup(name), &ret);
	kh_value(ns->nodes, k) = n;
    }

    return kh_value(ns->nodes, k);
}

#if 0
// A better hash with fewer collision than BuzHash, but slow as it's not
// rolling and we need to recompute every base kmer times.
uint32_t hash_seq(uint8_t *seq) {
    // FNV1a
    const uint32_t offset_basis = 2166136261;
    const uint32_t FNV_prime = 16777619;
    uint32_t h = offset_basis;
    for (int i = 0; i < kmer; i++)
	h = (h ^ (uint8_t) seq[i]) * FNV_prime;

    return h;
}
#endif

void nodeset_index_kmers(nodeset *ns, node *n, char *str, int bidir) {
    base4_init();

    int num = n->num;
    int i, j, len = strlen(str);

    // Assign "kmer" to node "num".  Duplicates get number -1
    // This is costly as we're recomputing hash every time rather than using
    // a rolling hash.
    uint32_t kh;
    uint8_t *str8 = (uint8_t *)str;
    for (i = 0; i < len - (kmer-1); i++) {
	kh = i
	    ? hash_shift(kh, str8[i-1], str8[i+kmer-1], kmer)
	    : hash_init(str8, kmer);

	//kh = hash_seq(str8+i); // FNV1a
	uint32_t k = kh & KMASK;

	int unique = 0;
	if (ns->kmer[k] && ns->kmer[k] != num) {
	    if (ns->kmer[k] != -1) {
		// Dups with an old node, we need to fix kmer_unique/dup.
		node *dup_n = ns->num2node[ns->kmer[k]];
		dup_n->kmer_unique -= ns->unique[k];
		dup_n->kmer_dup    += ns->unique[k];
	    }
	    ns->kmer[k] = -1;
	    n->kmer_dup++;
	} else {
	    ns->kmer[k] = num;
	    ns->unique[k]++;
	    n->kmer_unique++;
	    unique = 1;
	}
//	printf("%d Index %08x %.*s %s %s\n", bidir, k, kmer, str+i,
//	       n->name, unique?"":"dup");

	ndup  += !unique;
	nuniq +=  unique;
    }

    if (!bidir)
	return;

    // (Or we can do kmer^-1 to complement, but still needs reversal)

    // Reverse complement str in-situ and resubmit
    unsigned char c[256];
    c['A'] = c['a'] = 'T';
    c['C'] = c['c'] = 'G';
    c['G'] = c['g'] = 'C';
    c['T'] = c['t'] = 'A';
    for (i = 0, j = len-1; i < j; i++, j--) {
	char ci = c[str[i]];
	char cj = c[str[j]];
	str[i] = cj;
	str[j] = ci;
    }
    if (i == j)
	str[i] = c[str[i]];

    // Resubmit as other strand
    nodeset_index_kmers(ns, n, str, 0);
}

#define MAXLINE 1000000
static char line[MAXLINE];

nodeset *nodeset_load(char *fn, int bidir) {
    FILE *fp;
    nodeset *ns = NULL;

    fprintf(stderr, "Indexing nodeseq\n");
    if (!(fp = fopen(fn, "r"))) {
	perror(fn);
	goto err;
    }

    if (!(ns = nodeset_create()))
	goto err;

    node *n = NULL;
    int line_no = 0;
    while (fgets(line, MAXLINE, fp)) {
	size_t l = strlen(line);
	if (line[l-1] == '\n')
	    line[--l] = 0;

	if (*line == '@') {
	    // @name
	    n = nodeset_find(ns, line+1, 1);
	    if (!n)
		goto err;
	    line_no = 0;
	} else if (n) {
	    // seq
	    if (line_no == 0)
		n->length = l;
	    //node_add_kmers(n, line);
	    if (*line) {
		int offset = line_no && kmer_idx > kmer ? kmer_idx - kmer : 0;
		nodeset_index_kmers(ns, n, line+offset, bidir);
		line_no++;
	    }
	} else {
	    goto err;
	}
    }

    fclose(fp);
    return ns;

 err:
    if (fp)
	fclose(fp);
    if (ns)
	nodeset_free(ns);

    fprintf(stderr, "Indexing done\n");
}

void nodeset_report(nodeset *ns) {
    for (int i = 1; i <= ns->nnodes; i++) {
	node *n = ns->num2node[i];
	// First node, assuming sensible graph order, doesn't have a prefix
	// context sequence.
	//double expected = n->length - (i==1 ? KMER_IDX-1 : 0), expected2 = expected;
	double expected = n->length, expected2 = expected;
	// Account for expected unique vs dup hit rate.
	expected2 *= (n->kmer_unique+1.0) / (n->kmer_unique + n->kmer_dup+1.0);
	double ratio = (n->hit_count+0.01)/(expected2+0.01);
	// maximum possible based on length of node
	double ratio2 = (n->hit_count+n->hit_possible+0.01)/(n->length+0.01);
	//double ratio = (n->hit_count-expected+0.01)/(n->length+0.01);
	//double ratio = (n->hit_count+0.01)/(n->length-expected+0.01);
	if (ratio >= ratio2) // min
	    ratio = ratio2;
	// account for truncated nodes, eg at start and end of graph
	if (i==1 || i==ns->nnodes) {
	    ratio2 = (n->hit_count)/(n->hit_count+n->hit_possible+10.);
	    if (ratio < ratio2) // max
		ratio = ratio2;
	}
	printf("Node %10s\tlen %6d\texp %6.1f\thit %6d+%-6d\tratio %.2f\n",
	       n->name, n->length, expected2, n->hit_count,n->hit_possible,
	       ratio);
    }
}

/* ----------------------------------------------------------------------
 * The main application.
 */
void count_bam_kmers(nodeset *ns, bam1_t *b) {
#define _ 0
    //                             A C    G        T at 1, 2, 4, 8
    static int base16to4[16] = { _,0,1,_, 2,_,_,_, 3,_,_,_, _,_,_,_ };

    int i, len = b->core.l_qseq;
    uint8_t *seq = bam_get_seq(b);

//    printf("Seq %s\n", bam_get_qname(b));
    int last_node = -1, last_node_base = 0, last_node_poss = 0;
    int nposs_run = 0; // for node -1 and then changing node

    // Our hashing works on ASCII
    uint8_t *bases = malloc(b->core.l_qseq);
    nibble2base(seq, (char *)bases, b->core.l_qseq);
    
    uint32_t kh;
    for (i = 0; i < len-(kmer-1); i++) {
	kh = i
	    ? hash_shift(kh, bases[i-1], bases[i+kmer-1], kmer)
	    : hash_init(bases, kmer);
	//kh = hash_seq(bases+i); // FNV1a
	uint32_t k = kh & KMASK;
	
	int num = ns->kmer[k];
//	printf(" %2d %08x %.*s %s %d   %d %d %d\n", num, k, kmer, bases+i,
//	       num > 0 ? ns->num2node[num]->name : "?",
//	       i, nposs_run, last_node, last_node_poss);
	if (num > 0) {
	    ns->num2node[num]->hit_count++;
	    if (i > last_node_base+1 && last_node == num) {
		// Correct for missing kmers from SNPs
		int j;
		for (j = i-1; j > last_node_base && j > i-KMER_IDX; j--) {
//		    printf("Fix-up %d:%d  last_node_base %d\n",
//			   i, j, last_node_base);
		    ns->num2node[num]->hit_count++;
		}
//		if (j > last_node_base)
//		    printf("Possible hits from %d to %d\n", last_node_base, j);
		for (; j > last_node_base; j--)
		    ns->num2node[num]->hit_possible++;
	    }
	    if (nposs_run && last_node_poss > 0) {
//		printf("Possible %d hits in new node %s\n",
//		       nposs_run, ns->num2node[num]->name);
		ns->num2node[num]->hit_possible+=nposs_run;
	    }
//	} else if (last_node_poss > 0) {
//	    printf("Possible hit in node %d\n", last_node_poss);
//	    ns->num2node[last_node_poss]->hit_possible++;
	} else if (nposs_run && last_node_poss > 0) {
//		printf("Possible %d hits in NEW node %s\n",
//		       nposs_run, ns->num2node[last_node_poss]->name);
		ns->num2node[last_node_poss]->hit_possible+=nposs_run;
	}
	if (num) {
	    last_node_base = i; // records dup too so we only correct SNPs
	    last_node = num;
	    if (num > 0) {
		last_node_poss = num;
		nposs_run = 0;
	    } else {
		nposs_run++;
	    }
	}
    }
    free(bases);
    puts("");
}

void usage(int ret) {
    printf("kmer2node2 [options] graph.nodeseq in.fasta\n");
    printf("Options:\n");
    printf("   -h         Help\n");
    printf("   -k INT     Set kmer indexing / matching size\n");
    printf("   -K INT     Set kmer previously used in nodeseq creation\n");
    exit(ret);
}

int main(int argc, char **argv) {
    int ret = 1;
    samFile *sam = NULL;
    sam_hdr_t *hdr = NULL;
    bam1_t *b = NULL;
    int c;

    while ((c = getopt(argc, argv, "hk:K:")) != -1) {
	switch (c) {
	case 'k':
	    kmer = atoi(optarg);
	    break;

	case 'K':
	    kmer_idx = atoi(optarg);
	    break;

	case 'h':
	    usage(0);
	    break;

	default:
	    usage(1);
	    break;
	}
    }
    
    if (argc-optind != 2)
	usage(1);

    // Load the nodeseq file and mark kmer-to-node lookup table
    int bidir = 1;
    nodeset *ns = nodeset_load(argv[optind], bidir);
    if (!ns)
	return 1;

    if (!(sam = sam_open(argv[optind+1], "r"))) {
	perror(argv[optind+1]);
	goto err;
    }

    if (!(hdr = sam_hdr_read(sam)))
	goto err;

    // Read fasta as BAM as it's easy functionality we already have
    // (bar the annoying in-memory sequence encoding).
    if (!(b = bam_init1()))
	goto err;
    while (sam_read1(sam, hdr, b) >= 0) {
	count_bam_kmers(ns, b);
    }

    // Report node hit rates
    nodeset_report(ns);
    
    printf("%d / %d dups\n", ndup, ndup+nuniq);
    
    ret = 0;

 err:
    if (sam)
	ret = sam_close(sam) != 0;
    
    if (b)
	bam_destroy1(b);

    nodeset_free(ns);
    return ret;
}
