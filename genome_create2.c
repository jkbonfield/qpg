#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <assert.h>

/*

TODO: output GFA.

We start with a random sequence and then we edit it.  These edits can be
translocations, CNVs, LINEs and SINEs, and STRs.  Pretty much it's all
repetition, either immediately adjacent or somewhere else.

We can model these in a grammar.  The grammar will initially simply be a
series of labelled sequence states, and later on we collapse repeated into
labels into multiple outputs for the same grammar rule, which translates
directly to GFA nodes and edges.

Essentially our sequence is a linked list of substrings (labelled states).
Only at the end does the linear linked list get deduplicated to turn into
a graph.

The initial sequence state is one long string with a single state S.

----------
Translocations:
We split state S into S1, T, S2, T, S3 where S1, S2 and S3 are the substrings
taken from the initial S and T are the copied component (or T' if reverse?).

If T spans more than one state then we may have to split partial components
first. Eg

[AAAA][BBBB][CCCC]
   <--trans--->
[A][a][BBBB][c][C]

So A became A and a halves and C becomes c and C.  This first step ensures
the boundaries of T always fall on a junction between states, making the
copy trivial.

----------
Repeats:
This generates a repeat and then mirrors it to multiple places.  It is a single
new sequence which is splattered all over the place.
Eg [SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS]
to [SSS][RR][SSSSS][RR][S][RR][RR][S]

This is iteratively done meaning we can partially overwrite existing repeats.
However each iteration is the same deal as before:
1. Turn [S] into [S1]^[S2] for the start and end coordinates (^) of R
2. Overwrite any internal states between start and end with the new R state.

----------
STRs / CNVs:
We're copying existing data into 2 or more (whole number) copies.

Eg
[SSSSSSSSSSSSSSSSSSS]
    <--STR----->
[SS][XX][XX][XX][SSS]

1. As we replicate rather than overwrite, we need to break the states at
the start and end of a repeat unit into [S1]^[S2].

2. Identify the list of states covered by the repeat unit.

3. Replicate that state list many times (removing anything previously there)


--------------------
Implementation:

A state should have {ID, position, seq, next_ptr}.
ID is just a number, as is position.
Seq is 1 or more bases obviously (so either C or pascal style strings).

A function to force a state split at a fixed position.

A function to insert a state at a fixed position (from a given start ptr)

A function to remove all states (starting at a specific state ptr) up to
a defined position.
(Hence we can insert followed by remove to perform overwrite.)

*/

// A sequence state.
typedef struct state {
    struct state *next;
    int id;
    int pos;
    int dir;
    int len;
    char *seq;
} state_t;

static double rep_snp_rate = 0.001;
static double STR_snp_rate = 0.02;

typedef struct {
    int length;
} opts;

static int id = 1;
state_t *state_new(int pos, char *seq, int len) {
    state_t *s = calloc(1, sizeof(*s));
    s->id = id++;
    s->pos = pos;
    s->dir = 0;
    s->len = len;
    if (len) {
	s->seq = malloc(len);
	if (seq) {
	    memcpy(s->seq, seq, len);
	} else {
	    memset(s->seq, 0, len);
	}
    }
    return s;
}

state_t *state_dup(state_t *s) {
    state_t *n = malloc(sizeof(*s));
    memcpy(n, s, sizeof(*s));
    n->seq = malloc(n->len);
    memcpy(n->seq, s->seq, n->len);
    return n;
}

void state_free(state_t *s, int recurse) {
    if (recurse && s->next)
	state_free(s->next, recurse);

    free(s->seq);
    free(s);
}

// Ensures position 'pos' is at the boundary of two states.
// The supplied state can be any state prior to pos.
//
// Returns the state prior to pos, or NULL if first position
state_t *state_split(state_t *s, int pos) {
    state_t *last = s;
    // Find closest state
    while (s->next && s->pos + s->len < pos)
	last = s, s = s->next;

    if (s->pos == pos)
	return last; // It's already at the boundary of two states
    if (pos >= s->pos + s->len)
	return s;    // At end of sequence

    // Split the current state
    state_t *sub = state_new(pos, s->seq + pos-s->pos, s->pos + s->len - pos);
    s->len = pos - s->pos;
    sub->next = s->next;
    s->next = sub;

    return s;
}

// As per state_split, but look for the same identifier and if found also
// split that in the same manner so the duplicated strings maintain their
// duplicated status.
state_t *state_split2(state_t *s, int pos) {
    state_t *s_orig = s;
    state_t *last = s;
    // Find closest state
    while (s->next && s->pos + s->len < pos)
	last = s, s = s->next;

    if (s->pos == pos)
	return last; // It's already at the boundary of two states
    if (pos >= s->pos + s->len)
	return s;    // At end of sequence

    // Split the current state; s is left(old), sub = right(new)
    state_t *sub = state_new(pos, s->seq + pos-s->pos, s->pos + s->len - pos);
    s->len = pos - s->pos;
    sub->next = s->next;
    s->next = sub;
    int sub_id = sub->id;

    // Look for dups.
    state_t *dup = s_orig;
    for (; dup; dup = dup->next) {
	if (dup->id == s->id && dup != s) {
	    int dpos = dup->pos + pos - s->pos;
	    fprintf(stderr, "Found dup of %d at pos %d (%d)\n", dup->id, dup->pos, dpos);
	    // Similarly duplicate that node too.
	    state_t *sub = state_new(dup->pos + pos - s->pos,
				     dup->seq + pos - s->pos,
				     dup->len - (pos - s->pos));
	    sub->id = sub_id;
	    dup->len = pos - s->pos;
	    sub->next = dup->next;
	    dup->next = sub;
	}
    }

    return s;
}

// Linearise the entire state
void state_to_seq(FILE *seq_fp, state_t *s) {
    fprintf(seq_fp, ">seq\n");
    while (s) {
	fprintf(seq_fp, "%.*s", s->len, s->seq);
	s = s->next;
    }
    fprintf(seq_fp, "\n");
}

// Generate GFA
#define MAX_NODES 10000
void state_to_gfa(state_t *s) {
    int used[MAX_NODES] = {0};
    state_t *orig_s = s, *last = NULL;
    while (s) {
	int edge_id = s->id | (last ? last->id<<16 : 0);
	// Wrong, need a better dedup
	if (last) // && (used[s->id] & edge_id) != edge_id)
	    printf("L\tN%d\t+\tN%d\t+\t0M\n", last->id, s->id);

	if (used[s->id] == 0) {
	    printf("S\tN%d\t%.*s\n", s->id, s->len, s->seq);
	    used[s->id] = edge_id;
	}

	last = s;
	s = s->next;
    }
}

#define MIN(a,b) ((a)<(b)?(a):(b))

// Allocate and fill out a buffer with a region of sequence.
// Caller to free() result after use.
char *state_substr(state_t *s, int pos, int len) {
    // Find first overlapping state
    while (s && s->pos + s->len < pos)
	s = s->next;

    // Fill out buffer
    char *seq = malloc(len+1), *seq_to = seq;
    seq[len] = 0;
    int end = pos+len;
    while (s && s->pos <= end) {
	int l = MIN(end, s->pos + s->len) - pos;
	memcpy(seq_to, s->seq + pos - s->pos, l);
	seq_to += l;
	pos += l;
	len -= l;
	s = s->next;
    }

    return seq;
}

// Debug
void dump_state(state_t *s) {
    fprintf(stderr, "===\n");
    while (s) {
	fprintf(stderr, "%p\t%d\t%d\t%.*s\n", s, s->id, s->pos, s->len, s->seq);
	s = s->next;
    }
}

// Removes states 'from' to 'to' inclusively
void state_remove(state_t *from, state_t *to) {
    state_t *next;
    while (from) {
	next = from->next;
	state_free(from, 0);
	from = (from == to ? NULL : next);
    }
}

// Replace position at 'pos' with 'seq' in a given state list
state_t *state_replace(state_t **s, int pos, char *seq, int len) {
    // Ensure state list is fragmented at pos and pos+len
    state_t *start = state_split(*s, pos);
    state_t *end   = state_split(start, pos+len);

    state_t *ins = state_new(pos, seq, len);
    ins->next = end ? end->next : NULL;

    // Remove between start->next and end inclusively
    if (start == end) {
	state_remove(start, end);
	start = NULL;
    } else {
	state_remove(start->next, end);
	start->next = ins;
    }

    // Special case of first pos being replaced
    if (pos == 0) {
	*s = ins;
	if (start)
	    state_free(start, 0);
    }

    return ins;
}

// Replicate nodes from->next to to inclusive, doubling them up and
// keeping the IDs.  Eg A|BCD|E to A|BCD|BCD|E where A=from and D=to.
// Rewrites the new replicated from and to so we can call it repeatedly
// to keep adding copies with their positions marching on.
void state_replicate(state_t **from, state_t **to) {
//    fprintf(stderr, "Rep %p/%d to %p/%d\n", (*from)->next, (*from)->next->pos,
//	    *to, (*to)->pos);

    // Replicate from..to
    state_t *last = NULL, *first = NULL;
    state_t *orig_from = *from;
    int pos = (*to)->pos + (*to)->len;
    do {
	*from = (*from)->next;
	state_t *n = state_dup(*from);
	if (last)
	    last->next = n;
	else
	    first = n;
	last = n;
	n->pos = pos;
	pos += n->len;
    } while (*from != *to);

    // Now link in the new first..last chain after the from..to chain
    last->next = (*to)->next;
    (*to)->next = first;

    *to = last;
}

// Add repeat elements, such as SINE and LINEs.
// For a given rep_len we use the same repeat sequence, but mutate it
// a little each time so the repeat elements are not identical.
void add_rep(state_t **s, int length, int count, int rep_len, int code) {
    static int rep_len_cache = 0;
    static char *rep = NULL, *repr = NULL;
    static int rep_id = MAX_NODES/2;
    rep_id = id++;

    if (count == 1)
	count = 2;

    // Create the repeat
    // One repeat in many places, which is different to the translocations
    // which is many items copies, but not many copies per item.
    if (rep_len_cache != rep_len) {
	char comp[256];
	comp['A']='T';
	comp['T']='A';
	comp['C']='G';
	comp['G']='C';
	rep = realloc(rep, rep_len);
	repr = realloc(repr, rep_len);
	fprintf(stderr, "New repeat element of length %d\n", rep_len);
	for (int i = 0; i < rep_len; i++)
	    rep[i] = "ACGT"[random()&3];
	for (int i = 0; i < rep_len; i++)
	    repr[rep_len-1-i] = comp[(unsigned)rep[i]];
    }

    if (length <= rep_len)
	return;

    for (int i = 0; i < count; i++) {
	int pos = random()%(length - rep_len);

	// Mutate the repeat
	//int dir = (random()/15551)&1;
	//int dir = i&1;
	int dir = 0; // Dir for now is fixed
	fprintf(stderr, "REP len %d dir %d at %d\n", rep_len, dir, pos);
	char *copy = dir ? rep : repr;

	for (int j = 0; j < rep_len; j++)
	    if ((random()&65535) < rep_snp_rate*65536)
		copy[j] = "ACGT"[random()&3];

	// Copy it
	//memcpy(seq+pos, copy, rep_len);

	state_t *ins = state_replace(s, pos, copy, rep_len);
	// NB: Not perfect as we can partially overwrite a previous copy,
	// giving both long and short components with the same rep_id.
	ins->id = rep_id;
    }

    rep_id++;
}

// Short tandem repeats, or larger copy-number variations.  Dist_l is the
// distribution of the repeat element length (eg 1 or 2 for
// homopolymer and dinucleotide STRs) and dist_n is the distribution
// of number of copies.
// "Clustered" is the likelihood that the next position will be
// correlated to the current one.
void add_STRs(state_t **s, int length, double count_f, int dist_l,
	      int dist_n, double clustered, double STR_trans_rate, char code) {
    char *seq;
    int pos = -1, sz=1;
    int count = count_f + (drand48() < (count_f - (int)count_f));
    fprintf(stderr, "STR count=%d\n", count);
    for (int i = 0; i < count; i++) {
	// repeat length and copy number distributions
	int rlen, copies, k;
	for (k = 1; k < 20000 && (random()&65535) < dist_l; k++);
	rlen = k;

	for (k = 1; k < 20000 && (random()&65535) < dist_n; k++);

	copies = k;
	copies = 1 + 3*k / log(rlen+1);
	    
	//char str[20000];
	//for (int k = 0; k < rlen; k++)
	//    str[k] = "ACGT"[random()&3];
	int local = 0;
	sz = rlen*copies;
	if ((random()&65535) < clustered*65535 && pos >= 0) {
	    int opos = pos;
	    do {
		pos = opos + (random()%(sz*4)) - (sz*2);
	    } while (pos + rlen >= length || pos < 0);
	    local = 1;
	} else {
	    pos = random()%(length-rlen);
	}

	if (pos == 0) pos = 1;

	while (pos + rlen*copies >= length)
	    copies--;
	if (copies <= 0)
	    continue;

	assert(pos + rlen*copies < length);

//	pos = 5;
//	rlen = 10;
//	copies = 3;

	char *str = state_substr(*s, pos, rlen);
	if (copies == 1) copies = 2;
	fprintf(stderr, "%d(%d)\t%d %d\t%.*s\n",
		pos, local, rlen, copies, rlen>30?30:rlen, str);
	free(str);
	
	// 1. Split at start_1, end_1 (beg/end of first rep)
	//    and end_N (end of last rep).
	// 2. Delete from end_1->next to end_N inclusive.
	// 3. Replicate nodes from start_1->next to end_1 'copies' times.
	// 4. Apply mutation rates too. (TODO)

	state_t *start_1 = state_split2(*s, pos);
	state_t *end_1   = state_split2(*s, pos+rlen);
	state_t *end_N   = state_split2(*s, pos+copies*rlen);

	state_t *del_from = end_1->next;
	state_t *del_to   = end_N;
	end_1->next = end_N->next;
	state_remove(del_from, del_to);

	for (int k = 0; k < copies-1; k++)
	    state_replicate(&start_1, &end_1);
    }
}

int genome_create(opts *o) {
    // Create initial state
    state_t *s = state_new(0, NULL, o->length);
    for (int i = 0; i < o->length; i++)
	s->seq[i] = "ACGT"[random()&3];

double SINE_rate = 0.0005;
double LINE_rate = 0.0002;
#define SINE_len 300
#define LINE_len 2000
char *name = "x";

    fprintf(stderr, "==== creating %s\n", name);
    fprintf(stderr, "LINEs: %d\n", (int)(o->length * LINE_rate));
    add_rep(&s, o->length, o->length * LINE_rate, LINE_len, 'l');

    fprintf(stderr, "SINEs: %d\n", (int)(o->length * SINE_rate));
    add_rep(&s, o->length, o->length * SINE_rate, SINE_len, 's');
    
     dump_state(s);

static double STR_rate = 0.001;
static double STR_snp_rate = 0.02;
static double CNV_rate = 0.0005;

    add_STRs(&s, o->length, 1, 30000, 60000, 0.95, 0.02, 'r'); // STR
    add_STRs(&s, o->length, 1, 65400, 40000, 0.60, 0.10, 'r'); // CNV

    int id = 50;
    state_t *ins;

    // Copy states from start..end to pos, copying state numbers too.
    //state_copy(...);

    dump_state(s);

    FILE *seq_fp = fopen("_.fa", "w");
    state_to_seq(seq_fp, s);
    fclose(seq_fp);

    state_to_gfa(s);
    state_free(s, 1);

    return 0;
}


int main(int argc, char **argv) {
    int c;
    opts o;

    srandom(0);

    while ((c = getopt(argc, argv, "l:s:")) != -1) {
	switch (c) {
	case 'l':
	    o.length = atoi(optarg);
	    break;
	case 's':
	    srandom(atoi(optarg));
	    break;

	default:
	    fprintf(stderr, "Usage: genome_create [options]\n\n");
	    return 1;
	}
    }

    return genome_create(&o);
}
