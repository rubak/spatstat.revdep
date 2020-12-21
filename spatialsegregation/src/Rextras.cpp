/***********************************
 * stuff needed by other functions:
 *
 *
 * Tuomas Rajala <tuomas@sokkelo.net>
 *
 *   */
#include "Rextras.h"
/*******************************************************/
SEXP getListElement(SEXP list, const char *str)
// same as R  list$'str'. from the online manual.
{
       SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
       int i;

       for (i = 0; i < length(list); i++)
         if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
           elmt = VECTOR_ELT(list, i);
           break;
         }
       return elmt;
}

/*******************************************************/

// turn an R-vector-of-int-vectors to C-vector-of-int-vectors
void VectsxpToVector(SEXP nodelistR, std::vector<std::vector<int> > &result)
{
//	std::vector<std::vector<int> > result(length(nodelistR));
//	std::vector<int> *resnode;
	SEXP node;
	PROTECT(nodelistR = coerceVector(nodelistR,VECSXP) );
	int i,j;
	result.resize(length(nodelistR));
	for(i=0;i<length(nodelistR);i++)
	{
		node = VECTOR_ELT(nodelistR, i);
		PROTECT(node = coerceVector(node, INTSXP));
		for(j=0; j<length(node);j++)
			result[i].push_back(INTEGER(node)[j] );
		UNPROTECT(1);
	}
	UNPROTECT(1);
}

/**********************************************************************************/

SEXP vectorToSEXP(std::vector<std::vector<int> > nodelist)
//transform a std::vector<std::vector<int> > to SEXP, desctructive
{
	SEXP graph, *node;

	PROTECT(graph = allocVector(VECSXP, nodelist.size()));
	int i,j, *p, n;
	for(i=0;i< (int)nodelist.size();i++)
	{
		node = new SEXP;
		PROTECT(*node = allocVector(INTSXP, nodelist[i].size() ) );
		p = INTEGER(*node);
		n = nodelist[i].size();
		if(n<1) ;// p[0]=NULL;
		else
			for(j=0;j<n;j++)
			{
				p[j] = (int) nodelist[i].at(j);
			};
		nodelist[i].clear();
		SET_VECTOR_ELT(graph, i, *node);
		UNPROTECT(1);
	}
	UNPROTECT(1);
	return graph;
}
