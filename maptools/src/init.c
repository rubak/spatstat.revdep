/* Copyright 2017 by Roger S. Bivand. */

#include "maptools.h"
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

static const R_CMethodDef CEntries[]  = {
    {"RFindCG", (DL_FUNC) &RFindCG, 6},
    {NULL, NULL, 0}

};

static R_CallMethodDef CallEntries[] = {
    {"mtInsiders", (DL_FUNC) &mtInsiders, 2},
    {"R_point_in_polygon_mt", (DL_FUNC) &R_point_in_polygon_mt, 4},
    {"Rgshhs", (DL_FUNC) &Rgshhs, 6},
    {"Rshapeget", (DL_FUNC) &Rshapeget, 2},
    {"Rshapeinfo1", (DL_FUNC) &Rshapeinfo1, 1},
    {"shpwritepoint", (DL_FUNC) &shpwritepoint, 3},
    {"shpwritepolys", (DL_FUNC) &shpwritepolys, 2},
    {"shpwritelines", (DL_FUNC) &shpwritelines, 2},
    {NULL, NULL, 0}

};


void 
#ifdef HAVE_VISIBILITY_ATTRIBUTE
__attribute__ ((visibility ("default")))
#endif
R_init_maptools(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);

}



