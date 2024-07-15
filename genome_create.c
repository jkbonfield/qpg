#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

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

// Add repeat elements, such as SINE and LINEs.
// For a given rep_len we use the same repeat sequence, but mutate it
// a little each time so the repeat elements are not identical.
void add_rep(char *seq, int length, int count, int rep_len) {
    static int rep_len_cache = 0;
    static char *rep = NULL, *repr = NULL;

    // Create the repeat
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
    }
}

// Short tandem repeats, or larger copy-number variations.  Dist_l is the
// distribution of the repeat element length (eg 1 or 2 for
// homopolymer and dinucleotide STRs) and dist_n is the distribution
// of number of copies.
// "Clustered" is the likelihood that the next position will be
// correlated to the current one.
void add_STRs(char *seq, int length, int count, int dist_l, int dist_n,
	      double clustered, double STR_trans_rate) {
    //printf("Count=%d\n", count);
    int pos = -1, sz=1;
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

// Translocations
void add_trans(char *seq, int length, int count) {
    char comp[256];
    comp['A']='T';
    comp['T']='A';
    comp['C']='G';
    comp['G']='C';

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
    }
}

void genome_create(void) {
    char *bases = malloc(length+1);
    if (!bases) abort();

    // Initial pass
    int i;
    for (i = 0; i < length; i++)
	bases[i] = "ACGT"[random()&3];
    bases[length] = 0;

    fprintf(stderr, "LINEs: %d\n", (int)(length * LINE_rate));
    add_rep(bases, length, length * LINE_rate, LINE_len);
    fprintf(stderr, "SINEs: %d\n", (int)(length * SINE_rate));
    add_rep(bases, length, length * SINE_rate, SINE_len);

    fprintf(stderr, "STRs\n");
    //add_STRs(bases, length, length * STR_rate, 30000, 55000);
    add_STRs(bases, length, length * STR_rate, 30000, 60000, 0.95, 0.02);

    fprintf(stderr, "CNVs\n");
    add_STRs(bases, length, length * CNV_rate, 65400, 40000, 0.6, 0.1);

    fprintf(stderr, "Translocations\n");
    add_trans(bases, length, length * trans_rate);

    printf(">seq\n%s\n", bases);
}

int main(int argc, char **argv) {
    int seed = 0;
    int opt;
    while ((opt = getopt(argc, argv, "l:s:S:C:N:n:A:L:")) != -1) {
	switch (opt) {
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
		    "    -T fraction    Rate of translocations [%f]\n",
		    length, seed, STR_rate, CNV_rate, STR_snp_rate,
		    rep_snp_rate, SINE_rate, LINE_rate, trans_rate
		    );
	    return 1;
	}
    }

    srandom(seed);

    genome_create();
    return 0;
}
