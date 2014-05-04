/*
 * bic-seq.c
 *
 * Concept by Ruibin Xi, algorithm by Joe Luquette.
 * Modified by Ruibin Xi to get it work in R (March 22 2010)
 */

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ll.h"
#include "rbtree.h"
#include "bin.h"
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

double N;
double lambda = 1.0;

BIN_LIST bins = ll_new();

BIN *
bin_new(int tumor, int total, double freq, int start, int end)
{
	BIN *nb = malloc(sizeof(BIN));

	if (nb == NULL) {
		printf("bin_new: malloc: %s\n", strerror(errno));
		exit(1);
	}

	nb->tumor = tumor;
	nb->total = total;
	nb->freq = freq;
	nb->start = start;
	nb->end = end;
	nb->prev = NULL;
	nb->next = NULL;

	return nb;
}

void
bin_print(BIN *b)
{
	printf("%d\t%d\t%lf\t%d\t%d\n",
		b->tumor, b->total, b->freq, b->start, b->end);
}

/*
 * Windows are a doubly linked list for fast next() and prev()
 * and are stored in a red-black tree keyed by BIC difference.
 *
 * Some operations on windows:
 *     - bic_diff(window): compute the difference between the
 *       current BIC of the window and the BIC after running
 *       merge() on it.
 *     - merge(window): merge window.  Merging a window requires
 *       extensive modification of the (b - 1) windows to the left
 *       and right of the merged window.  See merge() for details.
 *     - rbinsert(index, window): insert window into a fast index
 *       keyed by the BIC difference.
 *     - rbdelete(index, window): removes window from the fast index
 *     - ll_next(window): next window in the global window list.  The
 *       global list is ordered by increasing genomic coordinates.
 *     - ll_prev(window): previous window in the global list.
 *     - window_print(window)
 *     - window_new(b, start): a new window with b bins, starting at
 *       bin 'start'.
 */
typedef struct window {
	int id;               /* XXX: hack to allow duplicate keys */
	double bic_diff;
	int b;                /* Number of bins in window */
	BIN *bin_start;       /* The first bin in this window */
	BLOCK *idx_entry;     /* The entry stored in the BIC idx */
	struct window *next;
	struct window *prev;
} WINDOW;

typedef LL(WINDOW) WINDOW_LIST;
WINDOW_LIST windows = ll_new();

void
window_print(WINDOW *w)
{
	int i, start, end;
	BIN *b;

	/* Find start/end */
	b = w->bin_start;
	start = b->start;
	end = b->end;
	for ( ; ll_next(b); b = ll_next(b))
		if (b->end > 0)
			end = b->end;

	printf("WINDOW (b=%d): (%d,%d), merge gain: %lf\n",
		w->b, start, end, w->bic_diff);
	b = w->bin_start;
	for (i = 0; i < w->b; ++i) {
		bin_print(b);
		b = ll_next(b);
	}
}

/*
 * Compute w's BIC.  Let:
 *    k be the number of parameters in w's model,
 *    N be the number of Bernoulli trials observed in w,
 *    L be the value of the likelihood function evaluated at MLE
 * then BIC = k*ln(N) - 2*ln(L).
 *
 * bins is a pointer to the first bin in a window, b is the number
 * of bins to be included in this window.
 */
double
compute_bic(BIN *bins, int b)
{
	int i;
	double logL, logN;
	BIN *x;

	logL = 0;  /* ln(P(mle)) */
	x = bins;
	for (i = 0; i < b; ++i) {
		if (x->tumor != 0 && x->tumor != x->total) {
			/* Compute log(L(mle)), L=likelihood function */
			logL += x->tumor * log(x->freq) +
				(x->total - x->tumor) * log(1.0 - x->freq);
		} /* else logL is 0 */

		x = x->next;
	}

	/* Some windows are empty.  In that case, log(N) is set to 0 */
	logN = N == 0 ? 0 : log(N);

	return b*logN*lambda - 2*logL;
}

/* Compute the difference in BIC between w as is and merge(w) */
double
bic_diff(WINDOW *w)
{
	int i;
	BIN *x, merged_bin;
	WINDOW merged_w;

	/* Create a fake bin */
	merged_bin.tumor = 0;
	merged_bin.total = 0;
	merged_bin.next = NULL;
	merged_bin.prev = NULL;
	x = w->bin_start;
	for (i = 0; i < w->b; ++i) {
		merged_bin.tumor += x->tumor;
		merged_bin.total += x->total;
		x = ll_next(x);
	}
	merged_bin.freq = merged_bin.total > 0 ?
		(double)merged_bin.tumor / (double)merged_bin.total : 0;

	merged_w.b = 1;  /* Fake window representing merge(w) */

	return compute_bic(&merged_bin, 1) - compute_bic(w->bin_start, w->b);
}

WINDOW *
window_new(int b, BIN *bin_start)
{
	static int last_id;
	WINDOW *nw = malloc(sizeof(WINDOW));

	if (nw == NULL) {
		printf("window_new: malloc: %s\n", strerror(errno));
		exit(1);
	}

	nw->id = last_id++;
	nw->b = b;
	nw->bin_start = bin_start;
	nw->idx_entry = mkblock(nw);
	nw->prev = NULL;
	nw->next = NULL;

	nw->bic_diff = bic_diff(nw);

	return nw;
}

WINDOW *
window_free(WINDOW *w)
{
	rmblock(w->idx_entry);
	free(w);
	return NULL;
}

/* Must return < 0 if a < b; > 0 if a > b; and =0 if a == b */
/* XXX: hack for equal keys: break tie by a unique ID assigned to */
/* each window. */
int
bcomp(BLOCK *a, BLOCK *b)
{
	double d;
	WINDOW *x = (WINDOW *)DATA(a);
	WINDOW *y = (WINDOW *)DATA(b);

	/* break ties on bic_diff using an ID */
	d = x->bic_diff - y->bic_diff;
	d = d == 0 ? x->id - y->id : d;
	return d < 0 ? -1 : (d > 0 ? 1 : 0);
}

/*
 * After deleting bins in merge(), a window on the right-most edge
 * of the chromosome may no longer be valid.  This function checks
 * the bin list to ensure that from w's start position, there are
 * at least w->b bins for w to represent.
 */
int
window_can_exist(WINDOW *w)
{
	int i;
	BIN *bin;

	bin = w->bin_start;
	for (i = 0; i < w->b && bin != NULL; ++i)
		bin = ll_next(bin);

	return i == w->b;
}

void
merge(RBTREE *idx, WINDOW *w)
{
	int i;
	BIN *b, *new_bin, *tmpb;
	WINDOW *z, *tmpw;

//if (w->bin_start->start > 191550000 && w->bin_start->end < 191750000) window_print(w);

	/*
	 * Step 1: merge the (b-1) bins right of the first bin in window
	 * w into the first bin.
	 * NOTE: Ruibin's R script outputs start=0 and end=0 for bins
	 * with 0 total reads.  So we only extend the bin's end point
	 * when the end is not 0.
	 */
	new_bin = w->bin_start;
	b = ll_next(w->bin_start);
	for (i = 1; i < w->b && b != NULL; ++i) {
		if (new_bin->start == 0)
			new_bin->start = b->start;
		new_bin->tumor += b->tumor;
		new_bin->total += b->total;
		if (b->end > 0)
			new_bin->end = b->end;
		tmpb = ll_next(b);
		ll_delete(bins, b);
		free(b);
		b = tmpb;
	}
	new_bin->freq = new_bin->total > 0 ?
		(double)new_bin->tumor / (double)new_bin->total : 0;

	/*
	 * Step 2: delete the (b-1) windows right of w.  They no longer
	 * exist since the (b-1) starting bins they correspond to no
	 * longer exist.  These windows are not merged; they're simply
	 * deleted.
	 * NOTE: there may be fewer than (b-1) windows right of w.
	 */
	z = ll_next(w);
	for (i = 1; i < w->b && z != NULL; ++i) {
		rbdelete(idx, z->idx_entry);
		tmpw = ll_next(z);
		ll_delete(windows, z);
		z = window_free(z);
		z = tmpw;
	}

	/*
	 * Step 3: recompute the BIC for w and the (b-1) windows left
	 * of w.  Recomputing the BIC implies adjusting the BIC index:
	 * delete the entry in the index for the old BIC and insert a
	 * new entry for the newly computed BIC.
	 * NOTE: there may be fewer than (b-1) windows left of w.
	 * MOREOVER: w and some of the (b-1) windows left of it may not
	 * have reason to exist after the bin merge.
	 * XXX: for right now, we'll use a non-optimal O(b) solution to
	 * determine which windows can stick around.  we'll fix this later
	 * if it's a performance problem.
	 */
	z = w;
	for (i = 0; i < w->b && z != NULL; ++i) {
		tmpw = z;
		rbdelete(idx, z->idx_entry);

		if (window_can_exist(z)) {
			z->bic_diff = bic_diff(z);
			if (rbinsert(idx, z->idx_entry) < 0) {
				printf("merge: key already exists in index\n");
				exit(1);
			}
		} else {
			ll_delete(windows, z);
			z = window_free(z);
		}

		z = ll_prev(tmpw);
	}
}

BLOCK *
rbmin_exhaustive(RBTREE *t, BLOCK *x)
{
	BLOCK *m, *z;
	extern BLOCK *NIL;

	m = x;
	if (LEFT(x) != NIL) {
		z = rbmin_exhaustive(t, LEFT(x));
		m = BCOMP(t, m, z) < 0 ? m : z;
	}

	if (RIGHT(x) != NIL) {
		z = rbmin_exhaustive(t, RIGHT(x));
		m = BCOMP(t, m, z) < 0 ? m : z;
	}

	return m;
}

WINDOW *
window_min()
{
	WINDOW *w, *m;

	m = ll_head(windows);
	for (w = ll_head(windows); w != NULL; w = ll_next(w))
		if (w->bic_diff < m->bic_diff)
			m = w;

	return m;
}
 
/*
 * Run bic_seq on the current list of bins, attempting to merge
 * b bins at a time.  Continues merging until no bin merge will
 * yield a lower BIC.
 *
 * Returns the number of successful merges.  If 0 merges occur,
 * no further levels should be attempted.
 */
int
bic_seq_level(int b)
{
	int i, p, q;
	int nmerged;
	BIN *x;
	BLOCK *m;
	RBTREE *idx;
	WINDOW *w;

	/* Compute the current BIC (this is just for display) */
	printf("\t%d-bins: total bins: %d\n", b, ll_length(bins));
	printf("\t%d-bins: initial BIC: %lf\n",
		b, compute_bic(ll_head(bins), ll_length(bins)));

	/*
	 * Compute a list of windows (based on b) and index them.  To
	 * index, we compute the BIC of the window before and after a
	 * a simulated merge() and record the difference.  The diff-
	 * erence keys the index.  The most negative value is the best.
	 */
	idx = rbcreate("BIC tree", bcomp, NULL);  /* No block-key comparator */
	x = ll_head(bins);
	p = q = 0;
	for (i = 0; i < ll_length(bins) - b+1; ++i) {
		w = window_new(b, x);
		w->bic_diff < 0 ? ++p : ++q;      /* Just for display */
		ll_append(windows, w);

		if (rbinsert(idx, w->idx_entry) < 0) {
			printf("rbinsert: key already exists in tree\n");
			exit(1);
		}
		x = ll_next(x);
	}
	//printf("\t%d-bins: %d windows, %d with BIC improvement, %d without\n"
	//	"\t%d-bins: index on BIC improvement: %d nodes, bl-height=%d\n",
	//	b, ll_length(windows), p, q,
	//	b, NODENO(idx), BH(idx));

	/* Merge until no window can lower the BIC of the model */
	nmerged = 0;
	for ( ; ; ) {
		/* m == NULL means idx is empty */
		m = rbmin(idx);
		if (m == NULL)
			break;

		w = (WINDOW *)DATA(m);
		if (w != NULL && w->bic_diff < 0) {
			merge(idx, w);
			++nmerged;
		} else
			break;
	}

	ll_dealloc(windows);
	idx = rbclose(idx);
	return nmerged;
}

void
bic_seq()
{
	int nmerged = 0;
	int level = 2;
        int flag = 0;
	
	while(flag==0){
		do {
			//printf("Merging bins by %ds..\n", level);
			nmerged = bic_seq_level(level);
			printf("Merged %d bins of size %d; %d remaining.\n",
				nmerged, level, ll_length(bins));
			++level;
		} while (nmerged>0);

		level = level -1;
		if(level<=2&&nmerged==0) flag = 1; else level = level-1;
	}
}

/*
 * Input file is the result of preprocessing by Ruibin's R
 * script.  One line encodes one bin.  Format is as follows:
 *	tumor reads | total reads | freq (tumor/total) | start | end
 */

SEXP bic_seq_Rinterface(SEXP bin_file_name, SEXP lmbda_SEXP)
{
	int tumor, total, start, end, i,nrow,ncol;
	double freq,*seg, *lambda_r=REAL(lmbda_SEXP);
	const char *bin_file;
	FILE *f;
	BIN *b=NULL;
	SEXP segment;

       if(!isString(bin_file_name)||length(bin_file_name)!=1){
                error("Multiple bin files specified.\n");
        }
        else bin_file = CHAR(STRING_ELT(bin_file_name,0));

	if(length(lmbda_SEXP)>1){
		warning("mutiple lambda specified, only the first one is used.\n");
	}


	lambda = lambda_r[0];
	if(lambda<=0){
		error("Parameter misspecification: the penalty lambda must be a positive number.\n");
	}

        f = fopen(bin_file, "r");
        if (f == NULL) {
                printf("fopen: %s: %s\n", bin_file, strerror(errno));
                exit(1);
        }

	N = 0;
        while (fscanf(f, "%d %d %lf %d %d",
                &tumor, &total, &freq, &start, &end) != EOF) {
                ll_append(bins, bin_new(tumor, total, freq, start, end));
                N += total;
        }
        printf("Initial bins detected: %d\n", ll_length(bins));
        printf("lambda is %g\n",lambda);

        bic_seq(bins);

	/*Now write the segments into a matrix*/
	nrow = ll_length(bins);
	ncol = 5;
	PROTECT(segment = allocMatrix(REALSXP, nrow, ncol));
	seg = REAL(segment);
	i = 0;
	for (b = ll_head(bins); b != NULL; b = ll_next(b)) {
		seg[i] = b->tumor;
		seg[i+nrow] = b->total;
		seg[i+2*nrow] =  b->freq;
		seg[i+3*nrow] = b->start;
		seg[i+4*nrow] = b->end;
		i++;	
	}
	ll_dealloc(bins);
	UNPROTECT(1);

	return(segment);
}
