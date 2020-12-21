#include "shapefil.h"
#include <R.h>
#include <Rdefines.h>

/* INTERFACE */ void Rshapeinfo(char **, int *, int *, double *, double *);

void Rshapeinfo(char **shpnm, int *Shapetype, int *Entities, double
		*MinBound, double *MaxBound)

{
    SHPHandle	hSHP;
    int    nShapeType, nEntities, i;
    double  adfMinBound[4], adfMaxBound[4];
/*      const char 	*pszPlus; */
/* -------------------------------------------------------------------- */
/*      Open the passed shapefile.                                      */
/* -------------------------------------------------------------------- */
    hSHP = SHPOpen( shpnm[0] , "rb" );

    if( hSHP == NULL )
    {
/*	REprintf( "Unable to open:%s\n", shpnm[0] );
	exit( 1 ); */
	error("No such file");
    }

/* -------------------------------------------------------------------- */
/*      Print out the file bounds.                                      */
/* -------------------------------------------------------------------- */
      SHPGetInfo( hSHP, &nEntities, &nShapeType, adfMinBound, adfMaxBound ); 
    
      *Entities = nEntities; 
      *Shapetype = nShapeType;
      for (i=0;i<4;i++){ 
        MinBound[i]=adfMinBound[i];
        MaxBound[i]=adfMaxBound[i];
    }

/*    Rprintf ("Info for %s\n", shpnm[0]);
    Rprintf("Shapefile Type: %s(%d)   # of Shapes: %ld\n\n",
            SHPTypeName( nShapeType ), nShapeType, nEntities );
    Rprintf("File Bounds: (%15.10lg,%15.10lg)\n\t(%15.10lg,%15.10lg)\n",
	    MinBound[0], MinBound[1], MaxBound[0], MaxBound[1] );*/
    

    SHPClose( hSHP );   
    return;
}





SEXP Rshapeinfo1(SEXP shpname)

{
    SEXP res, nms;
    SHPHandle	hSHP;
    int    nShapeType, nEntities, i, pc=0;
    double  adfMinBound[4], adfMaxBound[4];
    PROTECT(res = NEW_LIST(5)); pc++;
    PROTECT(nms = NEW_CHARACTER(5)); pc++;
    SET_STRING_ELT(nms, 0, COPY_TO_USER_STRING("fname"));
    SET_STRING_ELT(nms, 1, COPY_TO_USER_STRING("type"));
    SET_STRING_ELT(nms, 2, COPY_TO_USER_STRING("entities"));
    SET_STRING_ELT(nms, 3, COPY_TO_USER_STRING("minbounds"));
    SET_STRING_ELT(nms, 4, COPY_TO_USER_STRING("maxbounds"));
    setAttrib(res, R_NamesSymbol, nms);
    SET_VECTOR_ELT(res, 0, NEW_CHARACTER(1));
    SET_VECTOR_ELT(res, 1, NEW_INTEGER(1));
    SET_VECTOR_ELT(res, 2, NEW_INTEGER(1));
    SET_VECTOR_ELT(res, 3, NEW_NUMERIC(4));
    SET_VECTOR_ELT(res, 4, NEW_NUMERIC(4));
    SET_STRING_ELT(VECTOR_ELT(res, 0), 0, STRING_ELT(shpname, 0));
    
/*      const char 	*pszPlus; */
/* -------------------------------------------------------------------- */
/*      Open the passed shapefile.                                      */
/* -------------------------------------------------------------------- */
    hSHP = SHPOpen(CHAR(STRING_ELT(shpname, 0)), "rb" );

    if( hSHP == NULL ) error("Error opening SHP file");

/* -------------------------------------------------------------------- */
/*      Print out the file bounds.                                      */
/* -------------------------------------------------------------------- */
    SHPGetInfo( hSHP, &nEntities, &nShapeType, adfMinBound, adfMaxBound ); 
    INTEGER_POINTER(VECTOR_ELT(res, 1))[0] = nShapeType;
    INTEGER_POINTER(VECTOR_ELT(res, 2))[0] = nEntities;
    for (i=0; i<4; i++) { 
        NUMERIC_POINTER(VECTOR_ELT(res, 3))[i] = adfMinBound[i];
        NUMERIC_POINTER(VECTOR_ELT(res, 4))[i] = adfMaxBound[i];
    }

    SHPClose( hSHP );   
    UNPROTECT(pc);
    return(res);
}
