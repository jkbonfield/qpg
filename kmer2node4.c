#include <stdio.h>
#include <unistd.h>
#include <assert.h>
#include <math.h>

#include <htslib/sam.h>
#include <htslib/khash.h>
#include <htslib/kbitset.h>

#include "buzhash.h"

// TODO: when indexing, our expected rate should be the maximum of any
// nodeseq rather than the average of all.  We want to match *a* node
// and not *all* nodes.  Or we record mean and s.d. somehow?

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
#define KMER_IDX_START 50
static int kmer_idx = KMER_IDX_START;
// Whether to use non-unique kmers in scoring
static int use_non_uniq = 0;
// Output edge metrics too;
static FILE *edge_fp = NULL;

static int small_node_fix = 1;

// KMER_IDX is the hash table size for storing kmer entries.
#ifndef KMER_IDX
// mg    = 1m56 13m47 0m55  85.1	91.2	4.7	2.3	1.0
// 28/20 = 3m19 13m22 10m26 88.3	93.2	2.8	1.9	1.2
// 24/20 = 2m48 12m56 7m5   88.2	93.1	2.8	1.9	1.3
// 24/10 = 2m31 11m52 5m53  88.3	93.1	2.8	2.0	1.3
// 22/10 = 2m00 11m2  3m1   88.5	93.0	2.7	2.0	1.3 <<<
// 19/10 = 1m45 10m54 1m19  88.6	92.4	2.8	2.0	1.3
//
// With rearrange nodeset struct
// 22/10 = 1m54 10m43 2m42  88.5	93.0	2.7	2.0	1.3 <<<
#define KMER_IDX 22
#endif

#define KSIZE (1<<KMER_IDX)
#define KMASK (KSIZE-1)

#define MAX_DUPS 10 // maximum number of nodes a kmer can be in.
#define MAX_NODES (15*1024)

static int ndup = 0, nuniq = 0;

//#define TRANSITION_COUNT

// Node->node edge counts
// FIXME: need to map ++, +-, -+ and -- transitions
KHASH_MAP_INIT_INT64(edge, int)
khash_t(edge) *edges = NULL;

KHASH_MAP_INIT_STR(gfa_edge, int)
khash_t(gfa_edge) *gfa_edges = NULL;

/* ----------------------------------------------------------------------
 * Loads edges from a GFA.  Minimal parsing.
 * Returns 0 on success, <0 on failure
 */
int gfa_load(char *fn) {
    FILE *fp = fopen(fn, "r");
    if (!fp) {
	perror(fn);
	return -1;
    }

    if (!(gfa_edges = kh_init(gfa_edge)))
	return -1;

    kstring_t ks = KS_INITIALIZE;
    while (ks.l = 0, kgetline(&ks, (kgets_func *)fgets, fp) >= 0) {
	if (ks.l && *ks.s == 'L') {
	    // Edge
	    char *n1, *n2, d1, d2, *s = ks.s+2;
	    n1 = s;
	    if (!(s = strchr(s, '\t')))
		continue;
	    *s++ = 0;

	    d1 = *s;
	    if (!(s = strchr(s, '\t')))
		continue;
	    n2 = ++s;
	    if (!(s = strchr(s, '\t')))
		continue;
	    *s++ = 0;
	    d2 = *s;

	    char *e = malloc(strlen(n1)+strlen(n2)+4);
	    if (!e)
		return -1;
	    sprintf(e, "%s%c,%s%c", n1,d1, n2,d2);

	    int ret;
	    khiter_t k = kh_put(gfa_edge, gfa_edges, e, &ret);
	    if (ret < -1)
		return -1;
	    if (ret == 0)
		free(e); // already present

	    kh_value(gfa_edges, k) = 0; // FIXME: should use a SET.
	}
    }
    fclose(fp);

    return 0;
}

/*
 * Checks whether node n1 and n2 are linked with directions d1 and d2, which
 * are '+' or '-'.
 *
 * TODO: permit distance N for N transitions away? Eg A and C are linked
 * if we have A->B->C?  This may also need to be adjusted for M basepair away
 * we can we track some expectation based on the offset between two kmers
 * found within a sequence.
 *
 * Returns 1 if edge exists,
 *         0 if not.
 */
int gfa_edge_exists(khash_t(gfa_edge) *g, char *n1, char d1, char *n2, char d2) {
    if (!g)
	return 0;

    char *e = malloc(strlen(n1)+strlen(n2)+4);
    if (!e)
	return 0;
    sprintf(e, "%s%c,%s%c", n1,d1, n2,d2);
    khiter_t k = kh_get(gfa_edge, g, e);
    //printf("Check %s %d,%d\n", e, k, k!=kh_end(g));

    if (k == kh_end(g)) {
	char dt = d1, *nt = n1;
	d1 = d2=='+'?'-':'+';
	n1 = n2;
	d2 = dt=='+'?'-':'+';
	n2 = nt;
	sprintf(e, "%s%c,%s%c", n1,d1, n2,d2);
	k = kh_get(gfa_edge, g, e);
	//printf("Check %s %d,%d\n", e, k, k!=kh_end(g));
    }
    free(e);

    return k != kh_end(g);
}

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
    int hit_count;
    double hit_possible;
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

//// A node list
//typedef struct node_list {
//    int node;
//    struct node_list *next;
//} node_list;

typedef struct {
    char     kdup;           // last used ele of kmer[].  >0 implies duplicated
    uint16_t kmer[MAX_DUPS]; // maps hash(kmer) to node number
    char     kdir[MAX_DUPS]; // 0=unknown, 1=+ 2=- 3=both (eg dup within node)
    int      unique;	     // number of unique maps for this kmer
} kc;

typedef struct {
    khash_t(node) *nodes;
    int nnodes;
    kc kc [KSIZE];
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

    // Remove homopolymers
    for (i=1;i<len;i++)
	if (str[i] != str[0])
	    break;
    if (i == len)
	return;

    //printf("Indexing %s\n", str);
    // Assign "kmer" to node "num".  Duplicates get number -1
    uint32_t kh;
    uint8_t *str8 = (uint8_t *)str;
    for (i = 0; i < len - (kmer-1); i++) {
	kh = i
	    ? hash_shift(kh, str8[i-1], str8[i+kmer-1], kmer)
	    : hash_init(str8, kmer);

	//kh = hash_seq(str8+i); // FNV1a
	uint32_t k = kh & KMASK;

	int unique = 0, kc;

	if (ns->kc[k].kmer[0] && ns->kc[k].kmer[0] != num) {
	    if (ns->kc[k].kdup == 0) {
		// Not previously observed as duplicated.
		// Fix up the first occurance kmer_unique/dup counts.
		node *dup_n = ns->num2node[ns->kc[k].kmer[0]];
		dup_n->kmer_unique -= ns->kc[k].unique;
		dup_n->kmer_dup    += ns->kc[k].unique;
	    }
	    int kc = ++ns->kc[k].kdup;

	    // While it's dup, it may not be a new dup to us.
	    for (kc = 1; kc < MAX_DUPS && ns->kc[k].kmer[kc]; kc++)
		if (ns->kc[k].kmer[kc] == num)
		    break;

	    if (kc < MAX_DUPS && ns->kc[k].kmer[kc] == 0) {
		// Not observed with this node before
		//printf("Set ns->kmer[%d][%d] = %d\n", k, kc, num);
		ns->kc[k].kdup = kc;
		ns->kc[k].kmer[kc] = num;
		int kdir=2-bidir; // 1(+)->1 0(-)->2
		ns->kc[k].kdir[kc] = (ns->kc[k].kdir[kc] && ns->kc[k].kdir[kc] != kdir)
		    ? 3 : kdir;
	    }
	    n->kmer_dup++;
	} else {	
	    // First time finding it or only ever found in this same node
	    ns->kc[k].kmer[0] = num;
	    ns->kc[k].unique++;
	    n->kmer_unique++;
	    int kdir=2-bidir; // 1(+)->1 0(-)->2
	    ns->kc[k].kdir[0] = (ns->kc[k].kdir[0] && ns->kc[k].kdir[0] != kdir)
		? 3 : kdir;
	    unique = 1;
	}
	//printf("%d Index %08x %.*s %s %s\n", bidir, k, kmer, str+i,
	//       n->name, unique?"":"dup");

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

    fprintf(stderr, "Indexing done\n");
    fclose(fp);
    return ns;

 err:
    fprintf(stderr, "Failed to index nodeseq\n");
    if (fp)
	fclose(fp);
    if (ns)
	nodeset_free(ns);

    return NULL;
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
	double ratio1 = (n->hit_count+0.01)/(expected2+0.01);
	// maximum possible based on length of node
	double ratio2 = (n->hit_count+n->hit_possible+0.01)/(n->length+0.01);
	//double ratio = (n->hit_count-expected+0.01)/(n->length+0.01);
	//double ratio = (n->hit_count+0.01)/(n->length-expected+0.01);
	double ratio = ratio1 < ratio2 ? ratio1 : ratio2;
	// account for truncated nodes, eg at start and end of graph
	if (i==1 || i==ns->nnodes) {
	    ratio2 = (n->hit_count)/(n->hit_count+n->hit_possible+10.);
	    if (ratio < ratio2) // max
		ratio = ratio2;
	}
	if (ratio == 0 && ratio2 >= 1)
	    ratio = ratio2;

//	if (use_non_uniq) {
//	    if (ratio < ratio2)
//		ratio = ratio2;
//
//	    printf("%f %f \n", expected2, n->length - (n->length - expected2)/1.5);
//	    expected2 = n->length - (n->length - expected2)/1.5;
//	    double ratio3 = (n->hit_count+n->hit_possible+0.01)/(expected2+.01);
//	    if (ratio < ratio3)
//		ratio = ratio3;
//	}

	// Also compensation for small nodes.
	if (small_node_fix) {
	    double r = sqrt(n->length)/10;
	    ratio *= r>1?1:r;
	}

	printf("Node %10s\tlen %6d\texp %6.1f %6d %6d\thit %6d+%-6d\tratio %.2f\n",
	       n->name, n->length, expected2, n->kmer_unique, n->kmer_dup, n->hit_count,(int)n->hit_possible,
	       ratio);
    }

    if (edge_fp) {
	fflush(stdout); // Work around mixing stdout with edge_fp
	khiter_t k;
	for (k = kh_begin(edges); k != kh_end(edges); k++) {
	    if (!kh_exist(edges, k))
		continue;

	    int64_t ee = kh_key(edges, k);
	    int n1 = ee>>32;
	    int k1 = n1&3; n1>>=2;
	    int n2 = ee & ((1ULL<<32)-1);
	    int k2 = n2&3; n2>>=2;

	    if (n1 && n2)
		fprintf(edge_fp, "Edge %8s%c %8s%c %8d\n",
			ns->num2node[n1]->name, "?+-b"[k1],
			ns->num2node[n2]->name, "?+-b"[k2],
			kh_value(edges, k));
	}
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

    int poss_nodes[MAX_NODES] = {0};
    int in_node_rep = 0;

#ifdef DEBUG
    printf("Seq %s\n", bam_get_qname(b));
#endif
    int last_node = -1, last_node_base = 0, last_node_poss = 0, last_dir = 0;
    int last_node_count = 0;
    int nposs_run = 0; // for node -1 and then changing node
    int nposs_dir = 0;

    // Our hashing works on ASCII
    uint8_t *bases = malloc(b->core.l_qseq);
    nibble2base(seq, (char *)bases, b->core.l_qseq);
    
    uint32_t kh, last_node_kmer;
    for (i = 0; i < len-(kmer-1); i++) {
	kh = i
	    ? hash_shift(kh, bases[i-1], bases[i+kmer-1], kmer)
	    : hash_init(bases, kmer);
	//kh = hash_seq(bases+i); // FNV1a
	uint32_t k = kh & KMASK;

	int dup = (ns->kc[k].kdup > 0), was_dup = dup;
	int num = ns->kc[k].kmer[0];
	int dir = ns->kc[k].kdir[0];

#ifdef DEBUG
	printf("Pos %d, dup=%d num=%d,%s dir=%d last=%d,%d x %d",
	       i, dup, num, num?ns->num2node[num]->name:"*", dir, last_node, last_node_poss, nposs_run);
	for (int kc=0; kc<MAX_DUPS && ns->kc[k].kmer[kc]; kc++) {
	    int n = ns->kc[k].kmer[kc];
	    printf(" %s", ns->num2node[n]->name);
	}
	printf("\n");
#endif

	// FIXME: we may start with a duplicate node and transition into
	// unique, but for now we only rescue the other way around.

	// Dup, but maybe graph unambiguates it for us.
	if (dup && last_node_poss > 0) {
	    int kc, kc1, kcn = 0;
	    for (kc = 0; kc < MAX_DUPS && ns->kc[k].kmer[kc]; kc++) {
		int n = ns->kc[k].kmer[kc];

		//if (n != last_node_poss)
		//    printf("k=%08x i=%d: last=%d,%d dup %d %d %s\n", k, i, last_node_poss, last_dir,  kc, n, ns->num2node[n]->name);
		if (n == last_node_poss ||
		    gfa_edge_exists(gfa_edges,
				    ns->num2node[last_node_poss]->name,
				    "?+-b"[last_dir],
				    ns->num2node[n]->name,
				    "?+-b"[ns->kc[k].kdir[kc]]) ||
		    gfa_edge_exists(gfa_edges,
				    ns->num2node[n]->name,
				    "?-+b"[ns->kc[k].kdir[kc]],
				    ns->num2node[last_node_poss]->name,
				    "?-+b"[last_dir])) {
		    kc1 = kc;
		    kcn++;
		}
	    }

	    //printf("kcn=%d last kc1=%d %s\n", kcn, kc1, ns->num2node[ns->kc[k].kmer[kc1]]->name);

// Attempt to track the number of copies of a node before/after a dup
// so we don't create transitions unless the number is significant.
//
// This doesn't appear to help, but overall it needs a big rewrite to.
// 1. Allocate all kmers to all nodes (unique or dups).
// 2. Stitch together uniquely mapping nodes to ensure consistency.
//    Maybe prune random hits so we have high scoring unique data only.
// 3. Stitch together non-uniquely mapping nodes between the unique ones
//    based on the GFA to see what routes are feasible.  This is a brute force
//    path enumeration between all pairs uniquely mapped nodes.
// 4. Allocate node weight to the paths found in 3, including completely
//    unobserved nodes.
// 5. Emit the full path for this sequence.
//
// The only difference between this and the whole thing is whether we do it
// per read or at the end.  At the end is basically rewriting pathfinder.
// Per read however turns it into a mapper and permits more complex analysis
// such as diploid path finding or creating newly observed routes.

#undef TRANSITION_COUNT

#ifdef TRANSITION_COUNT
	    if (kcn == 1 && last_node_count > 1)
#else
	    if (kcn == 1)
#endif
	    {
		kc = kc1;
		num = ns->kc[k].kmer[kc];
		dir = ns->kc[k].kdir[kc];
#ifdef DEBUG
		if (last_node_poss != num)
		    printf("Found transition %s %s %d (%d)\n",
			   ns->num2node[last_node_poss]->name,
			   ns->num2node[num]->name,
			   last_node_count, kc);
#endif
		dup = 0;
	    }
	}

//	printf(" %2d %08x %.*s %s %d   %d %d %d, %d %d\n", num, k, kmer, bases+i,
//	       num > 0 ? ns->num2node[num]->name : "?",
//	       i, nposs_run, last_node, last_node_poss, dir, last_dir);
	if (!dup && num > 0) {
	    if (!use_non_uniq && was_dup)
		ns->num2node[num]->hit_possible++;
	    else
		ns->num2node[num]->hit_count++;
	    in_node_rep++;
	    if (i > last_node_base+1 && last_node == num) {
		// Correct for missing kmers from SNPs
		in_node_rep+=last_node_base - (i-1);
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
		nposs_run = 0;
	    } else if (nposs_run && last_node) {
#ifdef TRANSITION_COUNT
		// Scan ahead to count how many of the new node we have
		int j;
		uint32_t kh2 = kh;
		printf("Scan ahead from %d of %d\n", i, len-(kmer-1));
		for (j=i+1; j < len-(kmer-1); j++) {
		    kh2 = hash_shift(kh2, bases[j-1], bases[j+kmer-1], kmer);
		    uint32_t k = kh2 & KMASK;
		    printf("%d %d, %c, %d ns->kmer[%u][0]=%d\n",
			   i, j, bases[j-1], num, k, ns->kc[k].kmer[0]);
		    if (ns->kc[k].kdup > 0 || ns->kc[k].kmer[0] != num)
			break;
		}
		printf("Run for %d more kmers\n", j-i-1);
#endif
		int nposs = 0, nposs_best;
		for (int n = 0; n < ns->nnodes; n++) {
		    if (!poss_nodes[n] || poss_nodes[n] < 0.5 * nposs_run)
			continue;
//		    printf("Poss_nodes[%d]=%d, nposs_run=%d, %s %s\t",
//			   n, poss_nodes[n], nposs_run,
//			   ns->num2node[n]->name, ns->num2node[num]->name);

		    // Possible candidate for the dup part.
		    int exists = 
			(gfa_edge_exists(gfa_edges,
					 ns->num2node[n]->name,
					 "?+-b"[nposs_dir],
					 ns->num2node[num]->name,
					 "?+-b"[dir]) ||
			 gfa_edge_exists(gfa_edges,
					 ns->num2node[num]->name,
					 "?-+b"[dir],
					 ns->num2node[n]->name,
					 "?-+b"[nposs_dir]));
//		    printf("    %d %d:  %d\n", n, poss_nodes[i], exists);
		    if (exists) {
			nposs_best = n;
			nposs++;
		    }
		}
		//if (nposs == 1) {
		//    printf("Possible %d hits from dup node %s\n", nposs_run,
		//	   ns->num2node[nposs_best]->name);
		//} else {
		//    printf("Possible %d hits from multiple dup nodes\n",
		//	   nposs_run);
		//}
#ifdef TRANSITION_COUNT
		if (nposs == 1 && j-i > 1)
#else
		if (nposs == 1 /*&& j-i > 1*/)
#endif
		{
//		    printf("    Add %d hits between %d and %d\n",
//			   nposs_run, nposs_best, num);
		    if (use_non_uniq) {
			ns->num2node[nposs_best]->hit_count += nposs_run;
		    } else {
			ns->num2node[nposs_best]->hit_possible += nposs_run;
		    }
		    nposs_run = 0;
		    if (edges) {
			int n1 = nposs_best, n2 = num, ret;
			int d1 = nposs_dir, d2 = dir;
			uint64_t ee = ((int64_t)(n1*4+d1)<<32) | (n2*4+d2);
			khiter_t k = kh_put(edge, edges, ee, &ret);
			if (ret > 0) {
			    // new
			    kh_value(edges, k) = 1;
			} else {
			    kh_value(edges, k)++;
			}
		    }
		}
	    }
//	} else if (last_node_poss > 0) {
//	    printf("Possible hit in node %d\n", last_node_poss);
//	    ns->num2node[last_node_poss]->hit_possible++;
	} else if (nposs_run && last_node_poss > 0) {
//		printf("Possible %d hits in NEW node %s\n",
//		       nposs_run, ns->num2node[last_node_poss]->name);
		ns->num2node[last_node_poss]->hit_possible+=nposs_run;
		in_node_rep += nposs_run;
		nposs_run = 0;
	} else if (dup && num) {
	    for (int kc=0; kc<MAX_DUPS && ns->kc[k].kmer[kc]; kc++) {
		int n = ns->kc[k].kmer[kc];
		if (use_non_uniq)
		    ns->num2node[n]->hit_count+=0.1;
		else
		    ns->num2node[n]->hit_possible+=0.1;
	    }
	}

	if (edge_fp && !dup && last_node_poss != num && last_node_poss > 0) {
	    khiter_t k;
	    int ret;
	    int n1 = last_node_poss, d1 = last_dir;
	    int n2 = num,       d2 = dir;
	    if (d2 == 2 && 0) {
		// reverse; n1- n2- => n2+ n1+
		d2 = 1;
		n2 = last_node_poss;
		d1 = dir==3 ? 3 : 3-dir;
		n1 = num;
	    }
	    int self_loop = gfa_edge_exists(gfa_edges,
					    ns->num2node[n1]->name, 1,
					    ns->num2node[n1]->name, 1);
	    if (self_loop && in_node_rep - ns->num2node[n1]->length > 0) {
		// Add self loop edge count
		int d1 = 1, d2 = 1;
		uint64_t ee = ((int64_t)(n1*4+d1)<<32) | (n1*4+d2);
		khiter_t k = kh_put(edge, edges, ee, &ret);
		int l = in_node_rep - ns->num2node[n1]->length;
		if (ret > 0) {
		    // new
		    kh_value(edges, k) = l;
		} else {
		    kh_value(edges, k) += l;
		}
	    }
	    //fprintf(stderr, "Switch node (len %d loop %d mat-run %d) %s%c -> %s%c -> %d\n",
	    //	    ns->num2node[n1]->length, self_loop, in_node_rep,
	    //	    ns->num2node[n1]->name, "?+-b"[d1],
	    //	    ns->num2node[n2]->name, "?+-b"[d2],
	    //	    gfa_edge_exists(gfa_edges,
	    //			    ns->num2node[n1]->name, "?+-"[d1],
	    //			    ns->num2node[n2]->name, "?+-"[d2]));
	    in_node_rep = 0;

	    uint64_t ee = ((int64_t)(n1*4+d1)<<32) | (n2*4+d2);
	    k = kh_put(edge, edges, ee, &ret);
	    if (ret > 0) {
		// new
		kh_value(edges, k) = 1;
	    } else {
		kh_value(edges, k)++;
	    }
	}

	if (num) {
	    last_node_base = i; // records dup too so we only correct SNPs
	    last_node = num ? num : (dup ? -1 : 0);
	    if (!dup) {
		if (num == last_node_poss)
		    last_node_count++;
		else
		    last_node_count = 1;
		last_node_poss = num;
		last_node_kmer = k;
		last_dir = dir;
		if (nposs_run)
		    memset(poss_nodes, 0, ns->nnodes * sizeof(*poss_nodes));
		nposs_run = 0;
	    } else {
		nposs_run++;
		for (int kc = 0; kc < MAX_DUPS && ns->kc[k].kmer[kc]; kc++)
		    poss_nodes[ns->kc[k].kmer[kc]]++;
		nposs_dir = dir;
	    }
	}
    }
    free(bases);
}

void usage(int ret) {
    printf("kmer2node2 [options] graph.nodeseq in.fasta\n");
    printf("Options:\n");
    printf("   -h         Help\n");
    printf("   -k INT     Set kmer indexing / matching size (also sets -K)\n");
    printf("   -K INT     Set kmer previously used in nodeseq creation\n");
    printf("   -U         Account for non-unique kmers in kmer ratio\n");
    printf("   -E FILE    Output edge weights too to FILE\n");
    printf("   -G FILE    Load GFA from FILE\n");
    printf("   -s         Do not adjust depth on small nodes\n");
    exit(ret);
}

int main(int argc, char **argv) {
    int ret = 1;
    samFile *sam = NULL;
    sam_hdr_t *hdr = NULL;
    bam1_t *b = NULL;
    int c;

    while ((c = getopt(argc, argv, "hk:K:UE:G:s")) != -1) {
	switch (c) {
	case 's':
	    small_node_fix = 0;
	    break;

	case 'k':
	    kmer = atoi(optarg);
	    // fall through

	case 'K':
	    kmer_idx = atoi(optarg);
	    break;

	case 'U':
	    use_non_uniq++;
	    break;

	case 'E':
	    edge_fp = fopen(optarg, "w");
	    if (!edge_fp) {
		perror(optarg);
		return 1;
	    }
	    break;

	case 'G':
	    if (gfa_load(optarg) < 0)
		return 1;
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
    if (edge_fp)
	edges = kh_init(edge);

    while (sam_read1(sam, hdr, b) >= 0) {
	count_bam_kmers(ns, b);
    }

    // Report node hit rates
    nodeset_report(ns);
    if (edge_fp) {
	kh_destroy(edge, edges);
	fclose(edge_fp);
    }

    if (gfa_edges)
	// FIXME: mem leak on keys.
	kh_destroy(gfa_edge, gfa_edges);
    
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
