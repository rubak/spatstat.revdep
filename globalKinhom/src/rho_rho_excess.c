/* This file contains modified code from the spatstat package,
which is distributed as free software under the GNU Public Licence >=2,
and Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018

The modified code is distributed under the same licence:
Copyright (C) Thomas Shaw
Licence: GNU Public Licence >= 2 */

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Utils.h>
#include <R_ext/Rdynload.h>

#define TWOPI 6.2831853071795
void rho_rho_excess(nquery, xq, yq, ndata, xd, yd, nsep, xh, yh, rmaxi, sig, result) 
  /* inputs */
  int *nquery;            /* number of locations to be interrogated */
  double *xq, *yq;    /* (x,y) coordinates to be interrogated */
  int *ndata;            /* number of data points */
  double *xd, *yd;    /* (x,y) coordinates of data */
  int *nsep;
  double *xh, *yh;  /* (x,y) coordinates of h (one for now) */
  double *rmaxi;    /* maximum distance at which points contribute */
  double *sig;      /* Gaussian sd */
  /* output */
  double *result;   /* vector of computed density values */
{
  double coef, resulti; 
  double sigma, twosig2; 
  int i, j, jleft, h, hleft;
  double rmax, r2max;
  int nq,nd,nh;
  double xpi, ypi, xqi, yqi, d2, d2p, dx, dy, dxp, dyp;
  double x1left, hxleft;

  sigma = *sig;
  twosig2 = 2.0 * sigma * sigma;
  coef = 1.0/(TWOPI * sigma * sigma);
  coef = coef*coef;

  rmax = *rmaxi;
  r2max = rmax*rmax;

  nq = *nquery;
  nd = *ndata;
  nh = *nsep;

  if(nq == 0 || nd == 0 || nh == 0)
    return;

  /* jleft[i] <= jleft[i+1], likewise kleft, so only initialize once */
  /* i indexes the query points xq,yq */
  jleft=0;
  hleft = nh-1;

  for (i = 0; i < nq; i++) {
    xqi = xq[i];
    yqi = yq[i];

    /* given xq, range of h and range of x1 are determined. */
    /* index for h is h, index for x1 is j */
    /* get relevant range of j for this xq */
    x1left = xqi - rmax;
    while((xd[jleft] < x1left) && (jleft + 1 < nd)) ++jleft;

    /* get leftmost h s.t. xpi = xqi + xh[h] is in the window (> 0) */
    hxleft = -xqi;
    while ((xh[hleft] > hxleft) && (hleft > 0)) --hleft;

    /* TODO: only check ||h|| < rmax? */
    for (h = hleft; h < nh; h++) { /* h loop */
      xpi = xqi + xh[h];
      if (xpi > 1) break;
      ypi = yqi + yh[h];
      /* NOTE: this assumes window is unit square */
      if (ypi > 1 || ypi < 0) continue;

      resulti = 0;

      for (j = jleft; j < nd; j++) { /* first data loop */
        dx = xd[j] - xqi;
        if (dx > rmax) break;
        dy = yd[j] - yqi;

        d2 = dx*dx + dy*dy;
        if (d2 <= r2max) {
          dxp = dx - xh[h];
          if (dxp > rmax) break;
          dyp = dy - yh[h];

          d2p = dxp * dxp + dyp * dyp;
          if (d2p < r2max) {
            /* now contribute */
            resulti += exp(-(d2 + d2p)/twosig2);
          }
        } /* first data loop */
      }
      result[nh*i + h] = resulti*coef;
    } /* end h loop */
  }
}

static R_NativePrimitiveArgType rre_fn_t[] = {
    INTSXP, REALSXP, REALSXP,
    INTSXP, REALSXP, REALSXP,
    INTSXP, REALSXP, REALSXP,
    REALSXP, REALSXP,
    REALSXP
};
    
static const R_CMethodDef cMethods[] = {
    {"rho_rho_excess", (DL_FUNC) &rho_rho_excess, 12, rre_fn_t},
    {NULL, NULL, 0, NULL}
};

void R_init_globalKinhom(DllInfo *info) {
    R_registerRoutines(info, cMethods, NULL, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
    R_forceSymbols(info, TRUE);
}

