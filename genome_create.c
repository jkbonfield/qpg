#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

static int length = 100000;
static double STR_rate = 0.001;
static double STR_snp_rate = 0.02;
static double CNV_rate = 0.0005;
#define SINE_len 300
#define LINE_len 2000
static double SINE_rate = 0.00005;
static double LINE_rate = 0.00002;
static double rep_snp_rate = 0.001;
static double trans_rate = 0.0002;
static FILE *seq_out = NULL;
static FILE *meta_out = NULL;

// Rate of new items vs edit existing items
static double STR_new_rate = 1;
static double STR_edit_rate = 0;
static double CNV_new_rate = 1;
static double CNV_edit_rate = 0;
static double SNP_edit_rate = 0.005;
//static double SNP_edit_rate = 0.05; // extreme!
static double INS_edit_rate = 0.0005;
static double DEL_edit_rate = 0.0005;

// Add repeat elements, such as SINE and LINEs.
// For a given rep_len we use the same repeat sequence, but mutate it
// a little each time so the repeat elements are not identical.
void add_rep(char *seq, char *meta, int length, int count, int rep_len,
	     int code) {
    static int rep_len_cache = 0;
    static char *rep = NULL, *repr = NULL;

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
	int dir = i&1;
	fprintf(stderr, "REP len %d dir %d at %d\n", rep_len, dir, pos);
	char *copy = dir ? rep : repr;

	for (int j = 0; j < rep_len; j++)
	    if ((random()&65535) < rep_snp_rate*65536)
		copy[j] = "ACGT"[random()&3];

	// Copy it
	memcpy(seq+pos, copy, rep_len);
	memset(meta+pos, code, rep_len);
	meta[pos] = toupper(code);
    }
}

// Short tandem repeats, or larger copy-number variations.  Dist_l is the
// distribution of the repeat element length (eg 1 or 2 for
// homopolymer and dinucleotide STRs) and dist_n is the distribution
// of number of copies.
// "Clustered" is the likelihood that the next position will be
// correlated to the current one.
void add_STRs(char *seq, char *meta, int length, double count_f, int dist_l,
	      int dist_n, double clustered, double STR_trans_rate, char code) {
    //printf("Count=%d\n", count);
    int pos = -1, sz=1;
    int count = count_f + (drand48() < (count_f - (int)count_f));
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
	char *str = seq + pos;
	//memcpy(str,  seq + pos, rlen);

	fprintf(stderr, "%d(%d)\t%d %d\t%.*s\n",
		pos, local, rlen, copies, rlen>30?30:rlen, str);

	pos += rlen;
	for (int k = 0; k < copies; k++) {
	    if (pos + rlen >= length)
		break;
	    memcpy(seq+pos, str, rlen);
	    memset(meta+pos, code, rlen);
	    meta[pos] = toupper(code);
	    // Mutate at random
	    for (int l = 0; l < rlen; l++)
		if ((random()&65535) < STR_snp_rate*65536)
		    seq[pos+l] = "ACGT"[random()&3];

	    // random translocations induced at STRs/CNVs
	    pos = (random()&65535) < STR_trans_rate*65536
		? random()%(length-rlen)
		: pos + rlen;
	    //pos += rlen;
	}
    }
}

// Like add_STRs, but we're editing existing STRs, extending or shrinking
void edit_STRs(char *seq, char *meta, int length, double count_f, int dist_l,
	       int dist_n, double clustered, double STR_trans_rate, char
	       code) {
    int pos = -1, sz=1;
    int count = count_f + (drand48() < (count_f - (int)count_f));
    if (!count)
	return;

    // Identify known repeat starting points
    int *rep_start = malloc(length * sizeof(*rep_start));
    int *rep_len = malloc(length * sizeof(*rep_len));
    int nreps = 0;
    for (int i = 0; i < length; i++) {
	if (meta[i] == toupper(code)) {
	    rep_start[nreps] = i;
	    while (++i < length && meta[i] == code)
		;
	    rep_len[nreps] = i-- - rep_start[nreps];
	    //fprintf(stderr, "rep %d len %d\n", rep_start[nreps], rep_len[nreps]);
	    nreps++;
	}
    }

    if (nreps == 0)
	return;

    for (int i = 0; i < count; i++) {
	// Find a repeat element in meta.
	int rnum = random()%nreps;

	int rstart = rep_start[rnum];
	int rlen = rep_len[rnum];

	// FIXME: can't change length.
	// Maybe ignore and trim all to fixed later?
	if (random()%2) {
	    // Incr; lose bases at end of seq
	    fprintf(stderr, "INCR STR at %d+%d\n", rstart, rlen);
	    memmove(seq+rstart+rlen, seq+rstart, length - (rstart+rlen));
	    memmove(meta+rstart+rlen, meta+rstart, length - (rstart+rlen));
	    if (rstart+rlen+rlen < length) {
		memcpy(seq+rstart+rlen, seq+rstart, rlen);
		memcpy(meta+rstart+rlen, seq+rstart, rlen);
	    }
	} else {
	    // Decr; stutter at end of seq
	    fprintf(stderr, "DECR STR at %d+%d\n", rstart, rlen);
	    memmove(seq+rstart, seq+rstart+rlen, length - (rstart+rlen));
	    memmove(meta+rstart, meta+rstart+rlen, length - (rstart+rlen));
	}

	// Mutate at random
	for (int l = 0; l < rlen; l++)
	    if ((random()&65535) < STR_snp_rate*65536)
		seq[rstart+l] = "ACGT"[random()&3];
    }

    free(rep_start);
    free(rep_len);
}

// Translocations
void add_trans(char *seq, char *meta, int length, double count_f, char code) {
    char comp[256];
    comp['A']='T';
    comp['T']='A';
    comp['C']='G';
    comp['G']='C';

    int count = count_f + (drand48() < (count_f - (int)count_f));
    for (int i = 0; i < count; i++) {
	int tv_len = random()%256;
	while (random()%3)
	    tv_len*=2;
	while (tv_len > length/8)
	    tv_len/=3;

	int pos1 = random()%(length-tv_len);
	int pos2 = random()%(length-tv_len);
	fprintf(stderr, "Translocation %d from %d to %d\n",
		tv_len, pos1, pos2);

	if (i%2) {
	    // rev/comp
	    char *copy = malloc(tv_len);
	    for (int j = 0; j < tv_len; j++)
		copy[tv_len-1-j] = comp[(unsigned)seq[pos1+j]];
	    memcpy(seq+pos2, copy, tv_len);
	    free(copy);
	} else {
	    memmove(seq+pos2, seq+pos1, tv_len);
	}
	memset(meta+pos2, code, tv_len);
	meta[pos2] = toupper(code);
    }
}

void genome_create(char *bases, char *meta, char *name) {
    fprintf(stderr, "==== creating %s\n", name);
    fprintf(stderr, "LINEs: %d\n", (int)(length * LINE_rate));
    add_rep(bases, meta, length, length * LINE_rate, LINE_len, 'l');
    fprintf(stderr, "SINEs: %d\n", (int)(length * SINE_rate));
    add_rep(bases, meta, length, length * SINE_rate, SINE_len, 's');

    fprintf(stderr, "STRs\n");
    //add_STRs(bases, length, length * STR_rate, 30000, 55000);
    add_STRs(bases, meta, length, length * STR_rate * STR_new_rate,
	     30000, 60000, 0.95, 0.02, 'r');
    edit_STRs(bases, meta, length, length * STR_rate * STR_edit_rate,
	      30000, 60000, 0.95, 0.02, 'r');

    fprintf(stderr, "CNVs\n");
    add_STRs(bases, meta, length, length * CNV_rate * CNV_new_rate,
	     65400, 40000, 0.6, 0.1, 'c');
    edit_STRs(bases, meta, length, length * CNV_rate * CNV_edit_rate,
	      65400, 40000, 0.6, 0.1, 'c');

    fprintf(stderr, "Translocations\n");
    add_trans(bases, meta, length, length * trans_rate, 't');

    // Random SNP, INS and DEL mutations, for when we do derived sequences.
    for (int i = 0; i < length; i++)
	if (drand48() < SNP_edit_rate)
	    bases[i] = "ACGT"[random()&3];

    // TODO: indels, use a tmp copy so we can insert/del without memmoves.

    fprintf(seq_out, ">%s\n%s\n", name, bases);
    if (meta_out)
	fprintf(meta_out, ">%s\n%s\n", name, meta);

}

void population_create(int count) {
    // Initial pass;
    char **bases = malloc(count * sizeof(*bases));
    char **meta = malloc(count * sizeof(*meta));
    if (!bases || !meta) abort();

    bases[0] = malloc(length+1);
    meta[0] = malloc(length+1);
    if (!bases[0] || !meta[0]) abort();

    for (int i = 0; i < length; i++)
	bases[0][i] = "ACGT"[random()&3];
    memset(meta[0], '.', length);
    bases[0][length] = 0;
    meta[0][length] = 0;

    genome_create(bases[0], meta[0], "seq_0000#1#1");

    // proportion of the standard counts for new vs editing
    STR_new_rate = 0.01;
    STR_edit_rate = 0.5;
    CNV_new_rate = 0.01;
    CNV_edit_rate = 0.3;
    LINE_rate /= 50;
    SINE_rate /= 50;
    trans_rate /= 10;

    char name[100];
    for (int i = 1; i < count; i++) {
	// FIXME: make diploid
	int prev = random()%i;
	bases[i] = malloc(length+1);
	meta[i] = malloc(length+1);
	if (!bases[i] || !meta[i]) abort();

	memcpy(bases[i], bases[prev], length+1);
	memcpy(meta[i], meta[prev], length+1);

	sprintf(name, "seq_%04d#1#1", i);
	genome_create(bases[i], meta[i], name);
    }

    for (int i = 0; i < count; i++) {
	free(bases[i]);
	free(meta[i]);
    }
    free(bases);
    free(meta);
}

int main(int argc, char **argv) {
    int seed = 0, count =1;
    int opt;
    seq_out = stdout;

    while ((opt = getopt(argc, argv, "l:s:S:C:N:n:A:L:T:o:O:P:E:")) != -1) {
	switch (opt) {
	case 'P':
	    count = atoi(optarg);
	    break;
	case 'l':
	    length = atoi(optarg);
	    break;
	case 's':
	    seed = atoi(optarg);
	    break;

        // STRs and CNVs
	case 'S':
	    STR_rate = atof(optarg);
	    break;
	case 'C':
	    CNV_rate = atof(optarg);
	    break;
	case 'N':
	    STR_snp_rate = atof(optarg);
	    break;
	case 'n':
	    rep_snp_rate = atof(optarg);
	    break;
	case 'E':
	    SNP_edit_rate = atof(optarg);
	    break;

	// translocations
	case 'T':
	    trans_rate = atof(optarg);
	    break;

	// repeats
	case 'A':
	    SINE_rate = atof(optarg);
	    break;
	case 'L':
	    LINE_rate = atof(optarg);
	    break;

	// output locations
	case 'o':
	    seq_out = fopen(optarg, "w");
	    if (!seq_out) {
		perror(optarg);
		exit(1);
	    }
	    break;
	case 'O':
	    meta_out = fopen(optarg, "w");
	    if (!meta_out) {
		perror(optarg);
		exit(1);
	    }
	    break;

	default:
	    fprintf(stderr, "Usage: genome_create [options]\n\n"
		    "Options:\n"
		    "    -l length      Genome length [%d]\n"
		    "    -s seed        Random number seed [%d]\n"
		    "    -S fraction    Rate of STR occurance [%f]\n"
		    "    -C fraction    Rate of CNV occurance [%f]\n"
		    "    -N fraction    SNP rate inside STRs/CNVs [%f]\n"
		    "    -n fraction    SNP rate inside repeats [%f]\n"
		    "    -A fraction    Rate of 500bp repeat element [%f]\n"
		    "    -L fraction    Rate of 2000bp repeat element [%f]\n"
		    "    -T fraction    Rate of translocations [%f]\n"
		    "    -o FILE        Filename for fasta sequence [stdout]\n"
		    "    -O FILE        Filename for fasta meta-data [/dev/null]\n",
		    length, seed, STR_rate, CNV_rate, STR_snp_rate,
		    rep_snp_rate, SINE_rate, LINE_rate, trans_rate
		    );
	    return 1;
	}
    }

    srandom(seed);

    population_create(count);
    //genome_create();

    int err = 0;
    if (seq_out)  err |= fclose(seq_out) < 0;
    if (meta_out) err |= fclose(meta_out) < 0;

    return err;
}
