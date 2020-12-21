/* Copyright (c) 2004-6, Nicholas J. Lewin-Koh and Roger Bivand */

#include "maptools.h"

#include <R.h>
#include <Rdefines.h>

SEXP shpwritepoint(SEXP fname, SEXP shapes, SEXP ncol)

{
    SHPHandle   hSHP;
    SHPObject   *psShape;
    int         nShapeType, i, nShapes;
 
    if (INTEGER_POINTER(ncol)[0] == 2) nShapeType = SHPT_POINT;
    else nShapeType = SHPT_POINTZ;

/* -------------------------------------------------------------------- */
/*      Create the requested layer.                                     */
/* -------------------------------------------------------------------- */

    hSHP = SHPCreate(R_ExpandFileName(CHAR(STRING_ELT(fname,0))), nShapeType);

    if( hSHP == NULL )
    {
         error("Unable to create:%s\n", CHAR(STRING_ELT(fname,0)) );
    }

    nShapes = LENGTH(shapes)/INTEGER_POINTER(ncol)[0];
    if (nShapeType == SHPT_POINT) {
      for (i = 0; i < nShapes; i++) {
        psShape = SHPCreateObject(nShapeType, -1, 0, NULL, NULL, 1, 
          &NUMERIC_POINTER(shapes)[i], &NUMERIC_POINTER(shapes)[i + nShapes], 
	  NULL, NULL);

        SHPWriteObject(hSHP, -1, psShape);
        SHPDestroyObject(psShape);
      }
    } else {
      for (i = 0; i < nShapes; i++) {
        psShape = SHPCreateObject(nShapeType, -1, 0, NULL, NULL, 1, 
          &NUMERIC_POINTER(shapes)[i], &NUMERIC_POINTER(shapes)[i + nShapes], 
	  &NUMERIC_POINTER(shapes)[i + (2*nShapes)], NULL);

        SHPWriteObject(hSHP, -1, psShape);
        SHPDestroyObject(psShape);
      }
    }

    SHPClose(hSHP);

    return R_NilValue;
}


