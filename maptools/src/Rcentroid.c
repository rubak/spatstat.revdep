
/* 
 * Modified for R by Nicholas Lewin-Koh also made modifications for
 * multipart polygons. Sept 29, 2000
 */


#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include "maptools.h"

SEXP RshpCentrd_2d (SEXP);
SEXP R_RingCentrd_2d (int , SEXP, double *);


/* **************************************************************************
 * RshpCentrd_2d
 *
 * Return the single mathematical / geometric centroid of a potentially 
 * complex/compound RShapeObject
 *
 * reject non area SHP Types
 * 
 * **************************************************************************/
SEXP RshpCentrd_2d (SEXP call) {
    int		ring, ringPrev, ring_nVertices, rStart, nprts;
    int         i,j,totvert;
    double	Area, ringArea;
    SEXP ringCentrd, Cent, shape, flag, ringVerts;
    shape = CADR(call);
    flag = CADDR(call);     
   
/*     if ( !(SHPDimension(INTEGER(getAttrib(shape,install("shp.type")))[0])  */
/*            & SHPD_AREA) )   */
/*         	 error("Not a class of shape with defined 2d area"); */

   nprts = INTEGER(getAttrib(shape, install("nParts")))[0];
   Area = 0;
   if(INTEGER(flag)[0]==0 ||nprts==1){
     PROTECT(Cent=allocVector(REALSXP, 2));
     REAL(Cent)[0] = 0.0;
     REAL(Cent)[1] = 0.0;
   }
   else{
     PROTECT(Cent=allocMatrix(REALSXP, nprts, 2));
   }
   /* for each ring in compound / complex object calc the ring cntrd	*/
   
   ringPrev = INTEGER(getAttrib(shape, install("nVerts")))[0];
   totvert = INTEGER(getAttrib(shape, install("nVerts")))[0];

   if(nprts==0) nprts=1;
   for ( ring = nprts-1; ring >= 0; ring-- ) {
     rStart = INTEGER(VECTOR_ELT(shape,0))[ring];
     ring_nVertices = ringPrev - rStart;
/*  Rprintf("ringPrev= %d, rStart=%d, ring_nVertices=%d \n", */
/*          ringPrev, rStart, ring_nVertices); */

   PROTECT(ringVerts=allocMatrix(REALSXP, ring_nVertices, 2));
   for(i=rStart,j=0;i<ringPrev ;i++,j++){
     REAL(ringVerts)[j]=REAL(VECTOR_ELT(shape,1))[i];
     REAL(ringVerts)[j+ring_nVertices]=REAL(VECTOR_ELT(shape,1))[i+totvert];
   }
/*  Rprintf(" matrix begin %f, matrix end: %f \n", */
/*  	     REAL(ringVerts)[0],REAL(ringVerts)[(2*ring_nVertices)-1]); */   
     PROTECT(ringCentrd = 
             R_RingCentrd_2d (ring_nVertices, ringVerts, &ringArea));  
/*  Rprintf("xcent: %f, ycent: %f, area: %f\n ",  */
/*              REAL(ringCentrd)[0],REAL(ringCentrd)[1],ringArea ); */

     /* use Superposition of these rings to build a composite Centroid	*/
     /* sum the ring centrds * ringAreas,  at the end divide by total area */
     if(INTEGER(flag)[0]==0 ||nprts==1){
       REAL(Cent)[0] +=  REAL(ringCentrd)[0] * ringArea;
       REAL(Cent)[1] +=  REAL(ringCentrd)[1] * ringArea; 
     }
     else{
       REAL(Cent)[ring]= REAL(ringCentrd)[0];
       REAL(Cent)[ring+nprts]= REAL(ringCentrd)[1]; 
     }
     Area += ringArea; 
     ringPrev = rStart;
     UNPROTECT(2);
    }    

     /* hold on the division by AREA until were at the end */
   if(INTEGER(flag)[0]==0 ||nprts==1){
     REAL(Cent)[0] = REAL(Cent)[0] / Area;
     REAL(Cent)[1] = REAL(Cent)[1] / Area;
     UNPROTECT(1);   
     return ( Cent );
   }
   else{
     UNPROTECT(1);   
     return ( Cent );
   }
}


/* **************************************************************************
 * RingCentroid_2d
 * Copyright (c) 1999, Carl Anderson
 *
 * This code is based in part on the earlier work of Frank Warmerdam
 * 
 *
 * Return the mathematical / geometric centroid of a single closed ring
 *
 * **************************************************************************/
SEXP R_RingCentrd_2d (int nVert, SEXP xy, double *Area ) {
  int		iv /*, jv */;
/*  int		sign_x, sign_y; */
  double	/* dy_Area, */ dx_Area, Cx_accum, Cy_accum, ppx, ppy;
  double 	x_base, y_base, x, y;
  SEXP          RingCent;
/* the centroid of a closed Ring is defined as
 *
 *      Cx = sum (cx * dArea ) / Total Area
 *  and
 *      Cy = sum (cy * dArea ) / Total Area
 */      
   
  x_base = REAL(xy)[0];
  y_base = REAL(xy)[nVert];
  
  Cy_accum = 0.0;
  Cx_accum = 0.0;

  ppx = REAL(xy)[1] - x_base;
  ppy = REAL(xy)[nVert + 1] - y_base;
  *Area = 0;

/* Skip the closing vector */
  for ( iv = 2; iv <= nVert - 2; iv++ ) {
    x = REAL(xy)[iv] - x_base;
    y = REAL(xy)[nVert + iv] - y_base;

    /* calc the area and centroid of triangle built out of an arbitrary  */
    /* base_point on the ring and each successive pair on the ring  */
    
    /* Area of a triangle is the cross product of its defining vectors	 */
    /* Centroid of a triangle is the average of its vertices		 */

    dx_Area =  ((x * ppy) - (y * ppx)) * 0.5;
    *Area += dx_Area;
    
    Cx_accum += ( ppx + x ) * dx_Area;       
    Cy_accum += ( ppy + y ) * dx_Area;
/*  #ifdef DEBUG2 */
/*      printf("(ringcentrd_2d)  Pp( %f, %f), P(%f, %f)\n", ppx, ppy, x, y); */
/*      printf("(ringcentrd_2d)    dA: %f, sA: %f, Cx: %f, Cy: %f \n",  */
/*  		dx_Area, *Area, Cx_accum, Cy_accum); */
/*  #endif   */  
    ppx = x;
    ppy = y;
  }

/*  #ifdef DEBUG2 */
/*    printf("(ringcentrd_2d)  Cx: %f, Cy: %f \n",  */
/*    	( Cx_accum / ( *Area * 3) ), ( Cy_accum / (*Area * 3) )); */
/*  #endif */

  /* adjust back to world coords 
  */
  PROTECT(RingCent=allocVector(REALSXP,2));
  REAL(RingCent)[0] = ( Cx_accum / ( *Area * 3)) + x_base;
  REAL(RingCent)[1] = ( Cy_accum / ( *Area * 3)) + y_base;
  UNPROTECT(1);   
  return (RingCent);
}


/*	Modified 24 May 2005 Roger S. Bivand for maptools
	Written by Joseph O'Rourke
	orourke@cs.smith.edu
	October 27, 1995

	Computes the centroid (center of gravity) of an arbitrary
	simple polygon via a weighted sum of signed triangle areas,
	weighted by the centroid of each triangle.
	Reads x,y coordinates from stdin.  
	NB: Assumes points are entered in ccw order!  
	E.g., input for square:
		0	0
		10	0
		10	10
		0	10
	This solves Exercise 12, p.47, of my text,
	Computational Geometry in C.  See the book for an explanation
	of why this works. Follow links from
		http://cs.smith.edu/~orourke/

*/

#define DIM     2               /* Dimension of points */
typedef double  tPointd[DIM];   /* type double point */

/*#define PMAX    1000     	 Max # of pts in polygon */
/* typedef tPointd *tPolygond; */ /* type double polygon */

double  Area2( tPointd a, tPointd b, tPointd c );
void    FindCG( int n, tPointd *P, tPointd CG, double *Areasum2 );
void    Centroid3( tPointd p1, tPointd p2, tPointd p3, tPointd c );

void	RFindCG( int *n, double *x, double *y, double *xc, double *yc, 
		double *area ) {

	int i, nn;
	tPointd *P;
	tPointd CG;
	double Areasum2;
	nn = n[0];
	P = (tPointd *) R_alloc((size_t) nn, sizeof(tPointd));
	for (i=0; i<nn; i++) {
		P[i][0] = x[i];
		P[i][1] = y[i];
	}
	FindCG(nn, P, CG, &Areasum2);
	xc[0] = CG[0];
	yc[0] = CG[1];
	area[0] = Areasum2/2;
	return;
}

/*      
        Returns the cg in CG.  Computes the weighted sum of
	each triangle's area times its centroid.  Twice area
	and three times centroid is used to avoid division
	until the last moment.
*/
void     FindCG( int n, tPointd *P, tPointd CG, double *Areasum2)
{
        int     i;
        double  A2;        /* Partial area sum */    
	tPointd Cent3;

	CG[0] = 0;
	CG[1] = 0;
        Areasum2[0] = 0;
	for (i = 1; i < n-1; i++) {
	        Centroid3( P[0], P[i], P[i+1], Cent3 );
	        A2 =  Area2( P[0], P[i], P[i+1]);
		CG[0] += A2 * Cent3[0];
		CG[1] += A2 * Cent3[1];
		Areasum2[0] += A2;
	      }
        CG[0] /= 3 * Areasum2[0];
        CG[1] /= 3 * Areasum2[0];
	return;
}
/*
	Returns three times the centroid.  The factor of 3 is
	left in to permit division to be avoided until later.
*/
void    Centroid3( tPointd p1, tPointd p2, tPointd p3, tPointd c )
{
        c[0] = p1[0] + p2[0] + p3[0];
        c[1] = p1[1] + p2[1] + p3[1];
	return;
}
/* 
        Returns twice the signed area of the triangle determined by a,b,c,
        positive if a,b,c are oriented ccw, and negative if cw.
*/
double     Area2( tPointd a, tPointd b, tPointd c )
{
	double area;
	area = (b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1]);
	return(area);
}

