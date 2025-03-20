#include <stdio.h>
#include <unistd.h>
#include <assert.h>
#include <math.h>

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
#define KMER_IDX_START 50
static int kmer_idx = KMER_IDX_START;
// Whether to use non-unique kmers in scoring
static int use_non_uniq = 0;
// Output edge metrics too;
static FILE *edge_fp = NULL;

// KMER_IDX is the hash table size for storing kmer entries.
#ifndef KMER_IDX
#define KMER_IDX 28
#endif

#define KSIZE (1<<KMER_IDX)
#define KMASK (KSIZE-1)

#define MAX_DUPS 10 // maximum number of nodes a kmer can be in.
#define MAX_NODES (100*1024)

static int ndup = 0, nuniq = 0;

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

//// A node list
//typedef struct node_list {
//    int node;
//    struct node_list *next;
//} node_list;

typedef struct {
    khash_t(node) *nodes;
    int nnodes;
    char kdup[KSIZE]; // last used ele of kmer[].  >0 implies duplicated
    int kmer[KSIZE][MAX_DUPS];   // maps hash(kmer) to node number
    char kdir[KSIZE][MAX_DUPS];  // 0=unknown, 1=+ 2=- 3=both (eg dup within node)
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

	if (ns->kmer[k][0] && ns->kmer[k][0] != num) {
	    if (ns->kdup[k] == 0) {
		// Not previously observed as duplicated.
		// Fix up the first occurance kmer_unique/dup counts.
		node *dup_n = ns->num2node[ns->kmer[k][0]];
		dup_n->kmer_unique -= ns->unique[k];
		dup_n->kmer_dup    += ns->unique[k];
	    }
	    int kc = ++ns->kdup[k];

	    // While it's dup, it may not be a new dup to us.
	    for (kc = 1; kc < MAX_DUPS && ns->kmer[k][kc]; kc++)
		if (ns->kmer[k][kc] == num)
		    break;

	    if (kc < MAX_DUPS && ns->kmer[k][kc] == 0) {
		// Not observed with this node before
		//printf("Set ns->kmer[%d][%d] = %d\n", k, kc, num);
		ns->kdup[k] = kc;
		ns->kmer[k][kc] = num;
		int kdir=2-bidir; // 1(+)->1 0(-)->2
		ns->kdir[k][kc] = (ns->kdir[k][kc] && ns->kdir[k][kc] != kdir)
		    ? 3 : kdir;
	    }
	    n->kmer_dup++;
	} else {	
	    // First time finding it or only ever found in this same node
	    ns->kmer[k][0] = num;
	    ns->unique[k]++;
	    n->kmer_unique++;
	    int kdir=2-bidir; // 1(+)->1 0(-)->2
	    ns->kdir[k][0] = (ns->kdir[k][0] && ns->kdir[k][0] != kdir)
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
	double r = sqrt(n->length)/10;
	ratio *= r>1?1:r;

	printf("Node %10s\tlen %6d\texp %6.1f\thit %6d+%-6d\tratio %.2f\n",
	       n->name, n->length, expected2, n->hit_count,n->hit_possible,
	       ratio);
    }

    if (edge_fp) {
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

    static int poss_nodes[MAX_NODES] = {0};

    printf("Seq %s\n", bam_get_qname(b));
    int last_node = -1, last_node_base = 0, last_node_poss = 0, last_dir = 0;
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

	int dup = (ns->kdup[k] > 0), was_dup = dup;
	int num = ns->kmer[k][0];
	int dir = ns->kdir[k][0];

	printf("Pos %d, dup=%d num=%d,%s dir=%d last=%d,%d x %d\n",
	       i, dup, num, num?ns->num2node[num]->name:"*", dir, last_node, last_node_poss, nposs_run);

	// FIXME: we may start with a duplicate node and transition into
	// unique, but for now we only rescue the other way around.

	// Dup, but maybe graph unambiguates it for us.
	if (dup && last_node_poss > 0) {
	    int kc, kc1, kcn = 0;
	    for (kc = 0; kc < MAX_DUPS && ns->kmer[k][kc]; kc++) {
		int n = ns->kmer[k][kc];

//		if (n != last_node_poss)
//		    printf("i=%d: last=%d, dup %d %d\n", i, last_node_poss, kc, n);
		if (n == last_node_poss ||
		    gfa_edge_exists(gfa_edges,
				    ns->num2node[last_node_poss]->name,
				    "?+-b"[last_dir],
				    ns->num2node[n]->name,
				    "?+-b"[ns->kdir[k][kc]]) ||
		    gfa_edge_exists(gfa_edges,
				    ns->num2node[n]->name,
				    "?-+b"[ns->kdir[k][kc]],
				    ns->num2node[last_node_poss]->name,
				    "?-+b"[last_dir])) {
		    kc1 = kc;
		    kcn++;
		}
	    }
	    if (kcn == 1) {
		kc = kc1;
		num = ns->kmer[k][kc];
		dir = ns->kdir[k][kc];
//		if (last_node_poss != num)
//		    printf("Found transition %d %d (%d)\n", last_node_poss, num, kc);
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
		nposs_run = 0;
	    } else if (nposs_run) {
		int nposs = 0, nposs_best;
		for (int n = 0; n < ns->nnodes; n++) {
		    if (!poss_nodes[n] || poss_nodes[n] < 0.5 * nposs_run)
			continue;

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
		if (nposs == 1) {
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
		nposs_run = 0;
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
//	    fprintf(stderr, "Switch node %s%c -> %s%c -> %d\n",
//	    	    ns->num2node[n1]->name, "?+-b"[d1],
//	    	    ns->num2node[n2]->name, "?+-b"[d2],
//		    gfa_edge_exists(gfa_edges,
//				    ns->num2node[n1]->name, "?+-"[d1],
//				    ns->num2node[n2]->name, "?+-"[d2]));

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
		last_node_poss = num;
		last_node_kmer = k;
		last_dir = dir;
		if (nposs_run)
		    memset(poss_nodes, 0, ns->nnodes * sizeof(*poss_nodes));
		nposs_run = 0;
	    } else {
		nposs_run++;
		for (int kc = 0; kc < MAX_DUPS && ns->kmer[k][kc]; kc++)
		    poss_nodes[ns->kmer[k][kc]]++;
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
    exit(ret);
}

int main(int argc, char **argv) {
    int ret = 1;
    samFile *sam = NULL;
    sam_hdr_t *hdr = NULL;
    bam1_t *b = NULL;
    int c;

    while ((c = getopt(argc, argv, "hk:K:UE:G:")) != -1) {
	switch (c) {
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
