#include <stdlib.h>
#include <string.h>
#include "shapefil.h"
#include <R.h>
#include <Rinternals.h>

int 	SHPRingDir_2d ( SHPObject *psCShape, int Ring );
void	RFindCG( int *n, double *x, double *y, double *xc, double *yc, 
		double *area );
SEXP mtInsiders(SEXP n1, SEXP bbs);
SEXP R_point_in_polygon_mt(SEXP px, SEXP py, SEXP polx, SEXP poly);
SEXP Rgshhs(SEXP fn, SEXP mode, SEXP dolim, SEXP lim, SEXP level, SEXP minarea);
SEXP Rshapeget(SEXP shpnm, SEXP repair);
SEXP Rshapeinfo1(SEXP shpname);
SEXP shpwritepoint(SEXP fname, SEXP shapes, SEXP ncol);
SEXP shpwritepolys(SEXP fname, SEXP shapes);
SEXP shpwritelines(SEXP fname, SEXP shapes);

