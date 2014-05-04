#ifndef BIN_H

/* The smallest data type the computation handles */
typedef struct bin {
	int tumor;
	int total;
	double freq;
	int start;
	int end;
	struct bin *prev;
	struct bin *next;
} BIN;

typedef LL(BIN) BIN_LIST;

#endif
