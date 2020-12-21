/* Copyright (c) 2004, Nicholas J. Lewin-Koh and Roger Bivand */

#include "maptools.h"

#include <R.h>
#include <Rdefines.h>

SEXP shpwritepolys(SEXP fname, SEXP shapes)
{
    SHPHandle   hSHP;
    SHPObject   *psShape;
    int         nShapeType = SHPT_POLYGON, i, j, k, kk;
    int		nShapes, maxnParts=0, maxnVerts=0, pc=0, nDims;
    int		*nParts, *nVerts, *panPartStart, *from, *to;
    double      *padfX, *padfY, *padfZ=NULL;
    SEXP	SnParts, Spstart, SnDims;

    PROTECT(SnDims = NEW_CHARACTER(1)); pc++;
    SET_STRING_ELT(SnDims, 0, COPY_TO_USER_STRING("nDims"));

    nDims = INTEGER_POINTER(GET_ATTR(shapes, SnDims))[0];
 
    if (nDims == 2) nShapeType = SHPT_POLYGON;
    else if (nDims == 3) nShapeType = SHPT_POLYGONZ;
    else error("Invalid dimension");

/* -------------------------------------------------------------------- */
/*      Create the requested layer.                                     */
/* -------------------------------------------------------------------- */

    hSHP = SHPCreate(R_ExpandFileName(CHAR(STRING_ELT(fname,0))), nShapeType );

    if( hSHP == NULL )
    {
         error("Unable to create:%s\n", CHAR(STRING_ELT(fname,0)) );
    }

    nShapes = LENGTH(shapes);
    nParts = (int *) R_alloc((size_t) nShapes, sizeof(int));
    nVerts = (int *) R_alloc((size_t) nShapes, sizeof(int));
    PROTECT(SnParts = NEW_CHARACTER(1)); pc++;
    SET_STRING_ELT(SnParts, 0, COPY_TO_USER_STRING("nParts"));
    PROTECT(Spstart = NEW_CHARACTER(1)); pc++;
    SET_STRING_ELT(Spstart, 0, COPY_TO_USER_STRING("pstart"));

    for (i = 0; i < nShapes; i++) {
      nParts[i] = INTEGER_POINTER(GET_ATTR(VECTOR_ELT(shapes, i), SnParts))[0];
      if (nParts[i] > maxnParts) maxnParts = nParts[i];
      nVerts[i] = INTEGER_POINTER(VECTOR_ELT(GET_ATTR(VECTOR_ELT(shapes, i), 
		    Spstart), 1))[(nParts[i]-1)] - (nParts[i]-1);
      if (nVerts[i] > maxnVerts) maxnVerts = nVerts[i];
    } 
    panPartStart = (int *) R_alloc((size_t) maxnParts, sizeof(int));
    from = (int *) R_alloc((size_t) maxnParts, sizeof(int));
    to = (int *) R_alloc((size_t) maxnParts, sizeof(int));
    if (maxnVerts > 1000000 || maxnVerts < 1)
      error("Old polylist object cannot be exported");
    padfX = (double *) R_alloc((size_t) maxnVerts, sizeof(double));
    padfY = (double *) R_alloc((size_t) maxnVerts, sizeof(double)); 
    if (nShapeType == SHPT_POLYGONZ)
        padfZ = (double *) R_alloc((size_t) maxnVerts, sizeof(double));

    for (i = 0; i < nShapes; i++) {
      kk = 0;
      for (j = 0; j < nParts[i]; j++) {
        from[j] = INTEGER_POINTER(VECTOR_ELT(GET_ATTR(VECTOR_ELT(shapes, i), 
		    Spstart), 0))[j] - 1;
        panPartStart[j] = from[j] - j;
        to[j] = INTEGER_POINTER(VECTOR_ELT(GET_ATTR(VECTOR_ELT(shapes, i), 
		    Spstart), 1))[j] - 1;
        for (k=from[j]; k<=to[j]; k++) {
          padfX[kk] = NUMERIC_POINTER(VECTOR_ELT(shapes, i))[k];
          padfY[kk] = NUMERIC_POINTER(VECTOR_ELT(shapes,
                        i))[k+nVerts[i]+(nParts[i]-1)];
          if (nShapeType == SHPT_POLYGONZ)
              padfZ[kk] = NUMERIC_POINTER(VECTOR_ELT(shapes,
                        i))[k+2*(nVerts[i]+(nParts[i]-1))];

          kk++;
        }
      }
      if (kk != nVerts[i]) error("wrong number of vertices in polylist");

      if (nShapeType == SHPT_POLYGONZ)
          psShape = SHPCreateObject(nShapeType, -1, nParts[i], panPartStart, 
                    NULL, nVerts[i], padfX, padfY, padfZ, NULL);
      else psShape = SHPCreateObject(nShapeType, -1, nParts[i], panPartStart, 
                    NULL, nVerts[i], padfX, padfY, NULL, NULL);

      SHPWriteObject( hSHP, -1, psShape );
      SHPDestroyObject( psShape );
    } 



    SHPClose( hSHP );
    UNPROTECT(pc);

    return R_NilValue;
}


SEXP shpwritelines(SEXP fname, SEXP shapes)
{
    SHPHandle   hSHP;
    SHPObject   *psShape;
    int         nShapeType, i, j, k, kk;
    int		nShapes, maxnParts=0, maxnVerts=0, pc=0;
    int		*nParts, *nVerts, *panPartStart, *from, *to;
    double      *padfX, *padfY;
    SEXP	SnParts, Spstart;
 
    nShapeType = SHPT_ARC;

/* -------------------------------------------------------------------- */
/*      Create the requested layer.                                     */
/* -------------------------------------------------------------------- */

    hSHP = SHPCreate(R_ExpandFileName(CHAR(STRING_ELT(fname,0))), nShapeType );

    if( hSHP == NULL )
    {
         error("Unable to create:%s\n", CHAR(STRING_ELT(fname,0)) );
    }

    nShapes = GET_LENGTH(shapes);
    nParts = (int *) R_alloc((size_t) nShapes, sizeof(int));
    nVerts = (int *) R_alloc((size_t) nShapes, sizeof(int));
    PROTECT(SnParts = NEW_CHARACTER(1)); pc++;
    SET_STRING_ELT(SnParts, 0, COPY_TO_USER_STRING("nParts"));
    PROTECT(Spstart = NEW_CHARACTER(1)); pc++;
    SET_STRING_ELT(Spstart, 0, COPY_TO_USER_STRING("pstart"));

    for (i = 0; i < nShapes; i++) {
      nParts[i] = INTEGER_POINTER(GET_ATTR(VECTOR_ELT(shapes, i), SnParts))[0];
      if (nParts[i] > maxnParts) maxnParts = nParts[i];
      nVerts[i] = INTEGER_POINTER(VECTOR_ELT(GET_ATTR(VECTOR_ELT(shapes, i), 
		    Spstart), 1))[(nParts[i]-1)] - (nParts[i]-1);
      if (nVerts[i] > maxnVerts) maxnVerts = nVerts[i];
    } 
    panPartStart = (int *) R_alloc((size_t) maxnParts, sizeof(int));
    from = (int *) R_alloc((size_t) maxnParts, sizeof(int));
    to = (int *) R_alloc((size_t) maxnParts, sizeof(int));
/*    for (i = 0; i < nShapes; i++) {
      nVerts[i] = INTEGER_POINTER(GET_DIM(VECTOR_ELT(shapes, i)))[0];
      if (nVerts[i] > maxnVerts) maxnVerts = nVerts[i];
    } */
    if (maxnVerts < 1)
      error("list object cannot be exported");
    padfX = (double *) R_alloc((size_t) maxnVerts, sizeof(double));
    padfY = (double *) R_alloc((size_t) maxnVerts, sizeof(double)); 

    for (i = 0; i < nShapes; i++) {
      kk = 0;
      for (j = 0; j < nParts[i]; j++) {
        from[j] = INTEGER_POINTER(VECTOR_ELT(GET_ATTR(VECTOR_ELT(shapes, i), 
		    Spstart), 0))[j] - 1;
        panPartStart[j] = from[j] - j;
        to[j] = INTEGER_POINTER(VECTOR_ELT(GET_ATTR(VECTOR_ELT(shapes, i), 
		    Spstart), 1))[j] - 1;
        for (k=from[j]; k<=to[j]; k++) {
          padfX[kk] = NUMERIC_POINTER(VECTOR_ELT(shapes, i))[k];
          padfY[kk] = NUMERIC_POINTER(VECTOR_ELT(shapes,
                        i))[k+nVerts[i]+(nParts[i]-1)];
          kk++;
        }
      }
      if (kk != nVerts[i]) error("wrong number of vertices in polylist");

      psShape = SHPCreateObject(nShapeType, -1, nParts[i], panPartStart, 
                    NULL, nVerts[i], padfX, padfY, NULL, NULL);
/*      for (j = 0; j < nVerts[i]; j++) {
          padfX[j] = NUMERIC_POINTER(VECTOR_ELT(shapes, i))[j];
          padfY[j] = NUMERIC_POINTER(VECTOR_ELT(shapes, i))[j+nVerts[i]];
      }

      psShape = SHPCreateObject(nShapeType, -1, 0, NULL, NULL, nVerts[i], 
        padfX, padfY, NULL, NULL);*/

      SHPWriteObject( hSHP, -1, psShape );
      SHPDestroyObject( psShape );
    } 

    SHPClose( hSHP );
    UNPROTECT(pc);

    return R_NilValue;
}


