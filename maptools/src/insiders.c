# include <R.h>
# include <Rdefines.h>

#define ROFFSET 1

int pipbb(double pt1, double pt2, double *bbs);

int between(double x, double low, double up); 

SEXP mtInsiders(SEXP n1, SEXP bbs) {

	int n, pc=0;
	int i, j, k, k1;
	double bbi[4], bbj[4];
	int *yes, jhit[4], hsum;
	SEXP ans;

	n = INTEGER_POINTER(n1)[0];
	PROTECT(ans = NEW_LIST(n)); pc++;
	yes = (int *) R_alloc((size_t) n, sizeof(int));
	for (i=0; i < n; i++) {
		for (j=0; j < n; j++) yes[j] = 0;
		bbi[0] = NUMERIC_POINTER(bbs)[i];
		bbi[1] = NUMERIC_POINTER(bbs)[i + n];
		bbi[2] = NUMERIC_POINTER(bbs)[i + 2*n];
		bbi[3] = NUMERIC_POINTER(bbs)[i + 3*n];
		k = 0;
		for (j=0; j < n; j++) {
			if (i != j) {
				hsum = 0;
				bbj[0] = NUMERIC_POINTER(bbs)[j];
				bbj[1] = NUMERIC_POINTER(bbs)[j + n];
				bbj[2] = NUMERIC_POINTER(bbs)[j + 2*n];
				bbj[3] = NUMERIC_POINTER(bbs)[j + 3*n];
				for (k1=0; k1 < 4; k1++) jhit[k1] = 0;
				jhit[0] = pipbb(bbi[2], bbi[3], bbj);
    				jhit[1] = pipbb(bbi[0], bbi[1], bbj);
				jhit[2] = pipbb(bbi[0], bbi[3], bbj);
				jhit[3] = pipbb(bbi[2], bbi[1], bbj);

				for (k1=0; k1 < 4; k1++) 
					hsum = hsum + jhit[k1];
				if (hsum == 4) { yes[j] = 1;
					k = k + yes[j];
				}
			}
		}
		
		if (k != 0) {
			SET_VECTOR_ELT(ans, i, NEW_INTEGER(k));
			for (j=0, k1=0; j < n; j++) {
				if (yes[j] > 0)
					INTEGER_POINTER(VECTOR_ELT(ans, 
						i))[k1++] = j + ROFFSET;
			}
		}
	}
	UNPROTECT(pc); /* ans */
	return(ans);
}


int between(double x, double low, double up) {
	if (x >= low && x <= up) return(1);
	else return(0);
}

int pipbb(double pt1, double pt2, double *bbs) {
	if ((between(pt1, bbs[0], bbs[2]) == 1) && 
		(between(pt2, bbs[1], bbs[3]) == 1)) return(1);
	else return(0);
} 

