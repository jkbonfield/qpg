#include <stdio.h>

#include <htslib/sam.h>
#include <htslib/khash.h>
#include <htslib/kbitset.h>

/*
 * TODO:
 * Node changing from one node to another means kmers that fall between them.
 * eg: c1 c1 c1 ? ? ?...? ? a0c1 a0c1 a0c1.
 *
 * This is because we index c1 up to last base and a0c1 from 1st base,
 * but we have KMER-1 substrings in the query that overlap c1 and a0c1.
 *
 * We somehow need to add these in, base on graph edges.  If we can go from
 * node N to M, then we need a kmer constructed from the last few bases of N
 * and the first few of M to get the initial seeding correct.
 */

// kmer2node graph.nodeseq in.fasta

// nodeseq is a file consisting of "@name" and then one line per
// sequence that has been mapped to that node during a training run.
// The first sequence line is important as it is the graph node seq
// itself. We has this with the others, but the length matters as it
// determines the expected number of times the node occurs.

// Output is a list of node names and a weighting reflecting the
// number of relative kmers observed, normalised for the number we'd
// expect given the node length.

#ifndef KMER
#define KMER 14
#endif

#define KSIZE (1<<(2*KMER))
#define KMASK (KSIZE-1)

static int silent = 1;
static int ndup = 0, nuniq = 0;

// AND off the middle base from the KMER so we have BBBB-BBBB.
//#define KGAP  (KMASK & ~(3<<(2*(KMER/2))))

// AND off two bases within the KMER so we have BBB-BBB-BBB
//#define KGAP  (KMASK & ~((3<<(2*(KMER/3))) | (3<<(2*(2*KMER/3)))))

// AND off two bases within the KMER so we have BBBB(AC)BBBB
// If we're doing bi-directional these have to be central so we can reverse
//#define KGAP  ((KMASK & ~(15<<(2*((KMER-1)/2)))) | (1<<(2*((KMER-1)/2))))
#define KGAP  ((KMASK & ~(15<<(2*((KMER)/2)))) | (1<<(2*((KMER)/2))))

// AND off 3 bases
//#define KGAP  (KMASK & ~((3<<(2*(KMER/4)))|(3<<(2*(2*KMER/4)))|(3<<(2*(3*KMER/4)))))

// AND off 3 bases within the KMER so we have BBB--B-BBBB
//#define KGAP  (KMASK & ~(243<<(2*(KMER/2)-1)))

// Use the KMER as-is without gaps
//#define KGAP KMASK

#define MAX_NODES (100*1024)

// Debugging aid
static char *kmer2str(int kmer) {
    static char str[KMER+1];
    for (int i = 0; i < KMER; i++)
	str[KMER-1-i] = "ACGT"[(kmer >> (i*2))&3];
    str[KMER] = 0;
    return str;
}

/* ----------------------------------------------------------------------
 * A node has a name, a length, and a bitvector of which kmers are present.
 */
typedef struct {
    int num; // Nth node number
    char *name;
    int length;
    //kbitset_t *kmers; // REMOVE ME

    // Maybe also an expected kmer count, as a fraction so small nodes will be
    // partial kmer.
    //double kmer_dup; // expected fraction of kmers that are uniquely mapping
    int kmer_unique, kmer_dup; // how many indexed uniquely/dups

    int hit_count;
    // Possible hits as we didn't switch node inbetween.  This is an upper-bound on
    // the counting for cycles.
    int hit_possible;
} node;

void node_free(node *n) {
    free(n->name);
//    if (n->kmers)
//	kbs_destroy(n->kmers);
}

node *node_create(char *name) {
    node *n = calloc(1, sizeof(*n));
    if (!n)
	return NULL;

    n->name = strdup(name);
    //n->kmers = kbs_init(KSIZE);
    if (/*!n->kmers || */!n->name) {
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

/* ----------------------------------------------------------------------
 * A nodeset is a collection of nodes.  It is a node map indexed by name.
 */
KHASH_MAP_INIT_STR(node, node*)

typedef struct {
    khash_t(node) *nodes;
    int nnodes;
    int kmer[KSIZE];   // maps kmer to node number
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
//    for (khiter_t k = kh_begin(ns->nodes); k != kh_end(ns->nodes); k++) {
//	if (kh_exist(ns->nodes, k))
//	    node_free(kh_value(ns->nodes, k));
//    }
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

void nodeset_index_kmers(nodeset *ns, node *n, char *str, int bidir) {
    base4_init();

    int num = n->num;
    int kmer = 0, kmer2 = 0, i, j, len = strlen(str);
    for (i = 0; i < len && i < KMER-1; i++)
	kmer2 = (kmer2<<2) | base4[str[i]];

    // Assign "kmer" to node "num".  Duplicates get number -1
    for (; i < len; i++) {
	kmer2 = ((kmer2<<2) | base4[str[i]]) & KMASK;
	kmer = kmer2 & KGAP;
	//fprintf(stderr, "%s", kmer2str(kmer2));
	//fprintf(stderr, " %s\n", kmer2str(kmer));
	int unique = 0;
	if (ns->kmer[kmer] && ns->kmer[kmer] != num) {
	    if (ns->kmer[kmer] != -1) {
		// Dups with an old node, we need to fix kmer_unique/dup.
		node *dup_n = ns->num2node[ns->kmer[kmer]];
		dup_n->kmer_unique -= ns->unique[kmer];
		dup_n->kmer_dup    += ns->unique[kmer];
	    }
	    ns->kmer[kmer] = -1;
	    n->kmer_dup++;
	} else {
	    ns->kmer[kmer] = num;
	    ns->unique[kmer]++;
	    n->kmer_unique++;
	    unique = 1;
	}
	if (!silent)
	    printf("%d Index %08x %s %s %s\n", bidir, kmer, kmer2str(kmer),
		   n->name, unique?"":"dup");
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
		nodeset_index_kmers(ns, n, line, bidir);
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
}

void nodeset_report(nodeset *ns) {
    for (int i = 1; i <= ns->nnodes; i++) {
	node *n = ns->num2node[i];
	// First node, assuming sensible graph order, doesn't have a prefix
	// context sequence.
	//double expected = n->length - (i==1 ? KMER-1 : 0), expected2 = expected;
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

    int kmer = 0, kmer2 = 0, i, len = b->core.l_qseq;
    uint8_t *seq = bam_get_seq(b);
    for (i = 0; i < len && i < KMER-1; i++)
	kmer2 = (kmer2<<2) | base16to4[bam_seqi(seq, i)];

    if (!silent)
	printf("Seq %s\n", bam_get_qname(b));
    int last_node = -1, last_node_base = 0, last_node_poss = 0;
    int nposs_run = 0; // for node -1 and then changing node
    for (; i < len; i++) {
	kmer2 = ((kmer2<<2) | base16to4[bam_seqi(seq, i)]) & KMASK;
	kmer = kmer2 & KGAP;
	int num = ns->kmer[kmer];
	if (!silent)
	    printf(" %2d %08x %s %s %d   %d %d %d\n",
		   num, kmer, kmer2str(kmer),
		   num > 0 ? ns->num2node[num]->name : "?",
		   i, nposs_run, last_node, last_node_poss);
	if (num > 0) {
	    ns->num2node[num]->hit_count++;
	    if (i > last_node_base+1 && last_node == num) {
		// Correct for missing kmers from SNPs
		int j;
		for (j = i-1; j > last_node_base && j > i-KMER; j--) {
		    if (!silent)
			printf("Fix-up %d:%d  last_node_base %d\n",
			       i, j, last_node_base);
		    ns->num2node[num]->hit_count++;
		}
		if (!silent && j > last_node_base)
		    printf("Possible hits from %d to %d\n", last_node_base, j);
		for (; j > last_node_base; j--)
		    ns->num2node[num]->hit_possible++;
	    }
	    if (nposs_run && last_node_poss > 0) {
		if (!silent)
		    printf("Possible %d hits in new node %s\n",
			   nposs_run, ns->num2node[num]->name);
		ns->num2node[num]->hit_possible+=nposs_run;
	    }
//	} else if (last_node_poss > 0) {
//	    printf("Possible hit in node %d\n", last_node_poss);
//	    ns->num2node[last_node_poss]->hit_possible++;
	} else if (nposs_run && last_node_poss > 0) {
	    if (!silent)
		printf("Possible %d hits in NEW node %s\n",
		       nposs_run, ns->num2node[last_node_poss]->name);
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
    puts("");
}

int main(int argc, char **argv) {
    int ret = 1;
    samFile *sam = NULL;
    sam_hdr_t *hdr = NULL;
    bam1_t *b = NULL;
    
    if (argc != 3 && argc != 4) {
	fprintf(stderr, "Usage: kmer2node graph.nodeseq in.fasta [0/1]\n");
	return 1;
    }

    // Load the nodeseq file and mark kmer-to-node lookup table
    int bidir = argc == 4 ? atoi(argv[3]) : 0;
    nodeset *ns = nodeset_load(argv[1], bidir);
    if (!ns)
	return 1;

    if (!(sam = sam_open(argv[2], "r"))) {
	perror(argv[2]);
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
