/*	$Id: Rgshhs.c 321 2017-02-23 09:29:12Z rsbivand $
 *
 *	Copyright (c) 1996-2011 by P. Wessel and W. H. F. Smith
 *	See COPYING file for copying and redistribution conditions.
 *
 * PROGRAM:	gshhs.c
 * AUTHOR:	Paul Wessel (pwessel@hawaii.edu)
 * CREATED:	JAN. 28, 1996
 * PURPOSE:	To extract ASCII data from binary shoreline data
 *		as described in the 1996 Wessel & Smith JGR Data Analysis Note.
 * VERSION:	1.1 (Byte flipping added)
 *		1.2 18-MAY-1999:
 *		   Explicit binary open for DOS systems
 *		   POSIX.1 compliant
 *		1.3 08-NOV-1999: Released under GNU GPL
 *		1.4 05-SEPT-2000: Made a GMT supplement; FLIP no longer needed
 *		1.5 14-SEPT-2004: Updated to deal with latest GSHHS database (1.3)
 *		1.6 02-MAY-2006: Updated to deal with latest GSHHS database (1.4)
 *		1.7 11-NOV-2006: Fixed bug in computing level (&& vs &)
 *		1.8 31-MAR-2007: Updated to deal with latest GSHHS database (1.5)
 *		1.9 27-AUG-2007: Handle line data as well as polygon data
 *		1.10 15-FEB-2008: Updated to deal with latest GSHHS database (1.6)
 *		1.12 15-JUN-2009: Now contains information on container polygon,
 *				the polygons ancestor in the full resolution, and
 *				a flag to tell if a lake is a riverlake.
 *				Updated to deal with latest GSHHS database (2.0)
 *		1.13 15-JUL-2011: Now contains improved area information (2.2.0),
 *				 and revised greenwhich flags (now 2-bit; see gshhs.h).
 *				 Also added -A and -G as suggested by José Luis García Pallero,
 *				 as well as -Qe|i to control river-lake output.
 *
 *	Copyright (c) 1996-2004 by P. Wessel and W. H. F. Smith
 *	See COPYING file for copying and redistribution conditions.
 *
 *	This program is free software; you can redistribute it and/or modify
 *	it under the terms of the GNU General Public License as published by
 *	the Free Software Foundation; version 2 of the License.
 *
 *	This program is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *	GNU General Public License for more details.
 *
 *	Contact info: www.soest.hawaii.edu/pwessel */

/*
 * This modification of gshhs.c is Copyright (c) 2005-2011 Roger Bivand
 * Modification to swap function taken from Rsystat.c in foreign 071117
*/

#include "Rgshhs.h"
#include <R.h>
#include <Rdefines.h>
#include <R_ext/Applic.h>
#include <Rconfig.h>

static void swapb(void *result, int size)
{
#ifndef WORDS_BIGENDIAN
    int i;
    char *p = result, tmp;

    if (size == 1) return;
    for (i = 0; i < size/2; i++) {
	tmp = p[i];
	p[i] = p[size - i - 1];
	p[size - i - 1] = tmp;
    }
#endif
}



int getNpols(FILE *);

int gshhs_pipbb(double pt1, double pt2, double *bbs);

int gshhs_between(double x, double low, double up); 


SEXP Rgshhs(SEXP fn, SEXP mode, SEXP dolim, SEXP lim, SEXP level, SEXP minarea) 
{
	FILE *fp;
	double w, e, s, n, area, lon, lat, scale = 10.0;
/*	char source;*/
	char msg[255];
	const char *name[2] = {"polygon", "line"};
	int k, line, max_east = 270000000, n_read, /*flip,*/ Level, version, greenwich, src, m/*, river*/;
	struct POINT p;
	struct GSHHS h;
	int npols, pc=0;
	SEXP res, resnames, plist, choice, chosen, clip, subset;
	int i, ipols;
	signed int fpos;
	double bb[4], bbi[4];
	int k1[4], k2[4], j, j1, j2;

	fp = fopen (CHAR(STRING_ELT(fn, 0)), "rb");
	if (fp == NULL ) {
		sprintf(msg, "Could not find file %s", CHAR(STRING_ELT(fn, 0)));
		error(msg);

	}

	npols = getNpols(fp);
        line = 0;
        res = R_NilValue;
	if (INTEGER_POINTER(mode)[0] == 0) {
		fclose (fp);

		PROTECT(res = NEW_INTEGER(1)); pc++;
		INTEGER_POINTER(res)[0] = npols;
		UNPROTECT(pc); /* res */
		return(res);

	} else if (INTEGER_POINTER(mode)[0] > 0) {

		rewind(fp);

		PROTECT(res = NEW_LIST(14)); pc++;

		PROTECT(resnames = NEW_CHARACTER(14)); pc++;
		SET_STRING_ELT(resnames, 0, COPY_TO_USER_STRING("id"));
		SET_STRING_ELT(resnames, 1, COPY_TO_USER_STRING("n"));
		SET_STRING_ELT(resnames, 2, COPY_TO_USER_STRING("level"));
		SET_STRING_ELT(resnames, 3, COPY_TO_USER_STRING("source"));
		SET_STRING_ELT(resnames, 4, COPY_TO_USER_STRING("greenwich"));
		SET_STRING_ELT(resnames, 5, COPY_TO_USER_STRING("fpos"));
		SET_STRING_ELT(resnames, 6, COPY_TO_USER_STRING("area"));
		SET_STRING_ELT(resnames, 7, COPY_TO_USER_STRING("west"));
		SET_STRING_ELT(resnames, 8, COPY_TO_USER_STRING("east"));
		SET_STRING_ELT(resnames, 9, COPY_TO_USER_STRING("south"));
		SET_STRING_ELT(resnames, 10, COPY_TO_USER_STRING("north"));
		SET_STRING_ELT(resnames, 11, COPY_TO_USER_STRING("line"));
		SET_STRING_ELT(resnames, 12, COPY_TO_USER_STRING("container"));
		SET_STRING_ELT(resnames, 13, COPY_TO_USER_STRING("ancestor"));
		setAttrib(res, R_NamesSymbol, resnames);

		SET_VECTOR_ELT(res, 0, NEW_INTEGER(npols));
		SET_VECTOR_ELT(res, 1, NEW_INTEGER(npols));
		SET_VECTOR_ELT(res, 2, NEW_INTEGER(npols));
		SET_VECTOR_ELT(res, 3, NEW_INTEGER(npols));
		SET_VECTOR_ELT(res, 4, NEW_INTEGER(npols));
		SET_VECTOR_ELT(res, 5, NEW_INTEGER(npols));
		SET_VECTOR_ELT(res, 6, NEW_NUMERIC(npols));
		SET_VECTOR_ELT(res, 7, NEW_NUMERIC(npols));
		SET_VECTOR_ELT(res, 8, NEW_NUMERIC(npols));
		SET_VECTOR_ELT(res, 9, NEW_NUMERIC(npols));
		SET_VECTOR_ELT(res, 10, NEW_NUMERIC(npols));
		SET_VECTOR_ELT(res, 11, NEW_INTEGER(npols));
		SET_VECTOR_ELT(res, 12, NEW_INTEGER(npols));
		SET_VECTOR_ELT(res, 13, NEW_INTEGER(npols));

		fpos =  (signed int) ftell(fp);
		n_read = (int) fread ((void *)&h, (size_t)sizeof (struct GSHHS), 		    (size_t)1, fp);
/*		version = (h.flag >> 8) & 255;
		flip = (version != GSHHS_DATA_VERSION);	 Take as sign that byte-swabbing is needed */
/*		flip = (! (h.level > 0 && h.level < 5));	
 Take as sign that byte-swabbing is needed */
		i = 0;
		while (n_read == 1) {
/*		    if (flip) {*/
			swapb (&h.id, sizeof(int));
			swapb (&h.n, sizeof(int));
		/*	h.level = swapb ((unsigned int)h.level); */
			swapb (&h.west, sizeof(int));
			swapb (&h.east, sizeof(int));
			swapb (&h.south, sizeof(int));
			swapb (&h.north, sizeof(int));
			swapb (&h.area, sizeof(int));
		/*	h.version  = swapb ((unsigned int)h.version);
			h.greenwich = swabi2 ((unsigned int)h.greenwich);
			h.source = swabi2 ((unsigned int)h.source);*/
			swapb (&h.flag, sizeof(int));
                        swapb (&h.ancestor, sizeof(int));
                        swapb (&h.container, sizeof(int));
/*		    }*/
		    Level = h.flag & 255;
		    version = (h.flag >> 8) & 255;
		    if (!(version >= 9)) 
			warning("Data not same version as software %d", version);
		    greenwich = (h.flag >> 16) & 3;			/* Greenwich is 0-3 */
		    src = (h.flag >> 24) & 1;			/* Greenwich is 0 (WDBII) or 1 (WVS) */
/*		    river = (h.flag >> 25) & 1;			 River is 0 (not river) or 1 (is river) */
		    w = h.west  * GSHHS_SCL;	
/* Convert from microdegrees to degrees */
		    e = h.east  * GSHHS_SCL;
		    s = h.south * GSHHS_SCL;
		    n = h.north * GSHHS_SCL;
/*		    source = (src == 1) ? 'W' : 'C';	 Either WVS or CIA (WDBII) pedigree */
		    line = (h.area) ? 0 : 1;		/* Either Polygon (0) or Line (1) (if no area) */
		    m = h.flag >> 26;
		    scale = pow (10.0, (double)m);		/* Area scale */
		    area = h.area / scale;				/* Now im km^2 */
		    INTEGER_POINTER(VECTOR_ELT(res, 0))[i] = (signed int) h.id;
		    INTEGER_POINTER(VECTOR_ELT(res, 1))[i] = (signed int) h.n;
		    INTEGER_POINTER(VECTOR_ELT(res, 2))[i] = 
			(signed int) Level;
		    INTEGER_POINTER(VECTOR_ELT(res, 3))[i] = 
			(signed int) src;
		    INTEGER_POINTER(VECTOR_ELT(res, 4))[i] = 
			(signed int) greenwich;
		    INTEGER_POINTER(VECTOR_ELT(res, 5))[i] = (signed int) fpos;
		    NUMERIC_POINTER(VECTOR_ELT(res, 6))[i] = area;
		    NUMERIC_POINTER(VECTOR_ELT(res, 7))[i] = w;
		    NUMERIC_POINTER(VECTOR_ELT(res, 8))[i] = e;
		    NUMERIC_POINTER(VECTOR_ELT(res, 9))[i] = s;
		    NUMERIC_POINTER(VECTOR_ELT(res, 10))[i] = n;
		    INTEGER_POINTER(VECTOR_ELT(res, 11))[i] = (signed int) line;
		    INTEGER_POINTER(VECTOR_ELT(res, 12))[i] = (signed int) h.container;
		    INTEGER_POINTER(VECTOR_ELT(res, 13))[i] = (signed int) h.ancestor;

		    fseek (fp, (long) (h.n * (int) sizeof(struct POINT)), SEEK_CUR);

		    fpos =  (signed int) ftell(fp);
		    n_read = (int) fread((void *)&h, (size_t) sizeof (struct GSHHS), (size_t)1, fp);
		    i++;
		}

	}
	if (INTEGER_POINTER(mode)[0] == 1) {
		fclose (fp);
		UNPROTECT(pc);
		return(res);
	} else {
		if (INTEGER_POINTER(mode)[0] > 1) {

		    PROTECT(subset = NEW_INTEGER(npols)); pc++;
		    for (i=0; i<npols; i++) {
			INTEGER_POINTER(subset)[i] = 1;
			if (INTEGER_POINTER(VECTOR_ELT(res, 2))[i] > 
			    INTEGER_POINTER(level)[0]) 
			    INTEGER_POINTER(subset)[i] = 0;
			if (NUMERIC_POINTER(VECTOR_ELT(res, 6))[i] < 
			    NUMERIC_POINTER(minarea)[0]) 
			    INTEGER_POINTER(subset)[i] = 0;
		    }

		    if (LOGICAL_POINTER(dolim)[0] == TRUE) {
			PROTECT(choice = NEW_LIST(2)); pc++;
			SET_VECTOR_ELT(choice, 0, NEW_INTEGER(npols));
			SET_VECTOR_ELT(choice, 1, NEW_INTEGER(npols));
			for (i=0; i<npols; i++) {
			    INTEGER_POINTER(VECTOR_ELT(choice, 0))[i] = 0;
			    INTEGER_POINTER(VECTOR_ELT(choice, 1))[i] = 0;
			}
			ipols = 0;
			bb[0] = NUMERIC_POINTER(lim)[0];
			bb[1] = NUMERIC_POINTER(lim)[1];
			bb[2] = NUMERIC_POINTER(lim)[2];
			bb[3] = NUMERIC_POINTER(lim)[3];
			for (i=0; i<npols; i++) {
			  if (INTEGER_POINTER(subset)[i] == 1) {
			    j = 0;
			    bbi[0] = NUMERIC_POINTER(VECTOR_ELT(res, 7))[i];
			    bbi[1] = NUMERIC_POINTER(VECTOR_ELT(res, 8))[i];
			    bbi[2] = NUMERIC_POINTER(VECTOR_ELT(res, 9))[i];
			    bbi[3] = NUMERIC_POINTER(VECTOR_ELT(res, 10))[i];
			    for (k=0; k<4; k++) k1[k] = 0;
			    for (k=0; k<4; k++) k2[k] = 0;
			    k1[0] = gshhs_pipbb(bb[0], bb[2], bbi);
			    k1[1] = gshhs_pipbb(bb[0], bb[3], bbi);
			    k1[2] = gshhs_pipbb(bb[1], bb[2], bbi);
			    k1[3] = gshhs_pipbb(bb[1], bb[3], bbi);
			    k2[0] = gshhs_pipbb(bbi[0], bbi[2], bb);
			    k2[1] = gshhs_pipbb(bbi[0], bbi[3], bb);
			    k2[2] = gshhs_pipbb(bbi[1], bbi[2], bb);
			    k2[3] = gshhs_pipbb(bbi[1], bbi[3], bb);
			    for (k=0, j1=0; k<4; k++) j1+= k1[k];
			    for (k=0, j2=0; k<4; k++) j2+= k2[k];
			    INTEGER_POINTER(VECTOR_ELT(choice, 0))[i] = j1;
			    INTEGER_POINTER(VECTOR_ELT(choice, 1))[i] = j2;
			    if (j1 != 0 || j2 != 0) ipols++;
			  }
			} /* npols */
			PROTECT(chosen = NEW_INTEGER(ipols)); pc++;
			PROTECT(clip = NEW_INTEGER(ipols)); pc++;
			for (i=0, j=0; i<npols; i++) {
			    if (INTEGER_POINTER(VECTOR_ELT(choice, 0))[i] != 0
				|| INTEGER_POINTER(VECTOR_ELT(choice, 1))[i] 
				!= 0) {
				INTEGER_POINTER(chosen)[j] = i;
				INTEGER_POINTER(clip)[j] = 0;
				if (INTEGER_POINTER(VECTOR_ELT(choice, 1))[i] 
				    != 4) INTEGER_POINTER(clip)[j] = 1;
				j++;
			    }
			}


		        if (INTEGER_POINTER(mode)[0] == 4) {
		            fclose (fp);

		            UNPROTECT(pc);
		            return(choice);
		        }

		        if (INTEGER_POINTER(mode)[0] == 3) {
		            fclose (fp);

		            UNPROTECT(pc);
		            return(clip);
		        }

		    } /* dolim */ else {
			for (i=0, ipols=0; i<npols; i++)
			    ipols += INTEGER_POINTER(subset)[i];
			PROTECT(chosen = NEW_INTEGER(ipols)); pc++;
			for (i=0, j=0; i<npols; i++) {
			    if (INTEGER_POINTER(subset)[i] == 1) {
			        INTEGER_POINTER(chosen)[j] = i;
				j++;
			    }
			}
		    }

		    if (INTEGER_POINTER(mode)[0] == 2) {
		        fclose (fp);

		        UNPROTECT(pc);
		        return(chosen);
		    }

		    rewind(fp);
		    PROTECT(plist = NEW_LIST(ipols)); pc++;
		    for (i=0; i < ipols; i++) {
			j = INTEGER_POINTER(chosen)[i];
			if (j > 1) max_east = 180000000;
			j1 = INTEGER_POINTER(VECTOR_ELT(res, 1))[j];
			j2 = INTEGER_POINTER(VECTOR_ELT(res, 5))[j];
			fseek (fp, (long)((int) sizeof(struct GSHHS) + j2), SEEK_SET);
			SET_VECTOR_ELT(plist, i, allocMatrix(REALSXP, j1, 2));
			for (k = 0; k < j1; k++) {
			    if (fread ((void *)&p, 
				(size_t) sizeof(struct POINT), 
				(size_t) 1, fp) != 1) {
					sprintf (msg, 
			"Error reading file %s for %s %d, point %d.\n", 
			CHAR(STRING_ELT(fn, 0)), name[line], 
			INTEGER_POINTER(VECTOR_ELT(res, 0))[j], k);
					error(msg);
			    }
/*			    if (flip) {*/
				swapb (&p.x, sizeof(int));
				swapb (&p.y, sizeof(int));
/*			    }*/
			    lon = (INTEGER_POINTER(VECTOR_ELT(res, 4))[j] 
			    	&& p.x > max_east) ? 
				p.x * GSHHS_SCL - 360.0 : p.x * GSHHS_SCL;
			    lat = p.y * GSHHS_SCL;
			    NUMERIC_POINTER(VECTOR_ELT(plist, i))[k] =  lon;
			    NUMERIC_POINTER(VECTOR_ELT(plist, i))[k+j1] =  lat;
			}
		    }
		    fclose (fp);

		    UNPROTECT(pc);
		    return(plist);
		}
	}
        error("shouldn't be here");
        return(R_NilValue);
}

int getNpols(FILE *fp) {
	struct GSHHS h;
	int n_read/*, flip, version*/;
	int n;

	n_read = (int) fread ((void *)&h, (size_t)sizeof (struct GSHHS), 
		(size_t)1, fp);
/*	version = (h.flag >> 8) & 255;
	flip = (version != GSHHS_DATA_VERSION);	 Take as sign that byte-swabbing is needed */
/*	flip = (! (h.level > 0 && h.level < 5));	
 Take as sign that byte-swabbing is needed */
	
	n=0;
	while (n_read == 1) {
/*		if (flip) {*/
			swapb (&h.n, sizeof(int));
/*		}*/
		fseek (fp, (long)(h.n * (int) sizeof(struct POINT)), SEEK_CUR);
		n_read = (int) fread((void *)&h, (size_t)sizeof (struct GSHHS), 
			(size_t)1, fp);
		n++;
	}

	return(n);
}

int gshhs_between(double x, double low, double up) {
	if (x >= low && x <= up) return(1);
	else return(0);
}

int gshhs_pipbb(double pt1, double pt2, double *bbs) {
	if ((gshhs_between(pt1, bbs[0], bbs[1]) == 1) && 
		(gshhs_between(pt2, bbs[2], bbs[3]) == 1)) return(1);
	else return(0);
} 

