#include "Fun.h"

Fun::Fun()
{
}

Fun::~Fun()
{
}


void Fun::Init(Graph *g0, double *par0, int *parn, int *gt, int *ft, double *fpar0, int *trans0, int *included0, int *dbg0)
{
	int i;
	graph = g0;
	std::vector<double> * pvec;
	for(i=0;i<*parn;i++)
	{
		parvec.push_back(par0[i]);
		pvec = new std::vector<double>;
		value.push_back(*pvec);
		value.at(i).push_back(0.0);
		delete pvec;
	}
	included = included0;
	autoborder = 0;
	if(included[0]<0) {
		g0->pp->calcEdgeDists(); // precompute the edge distances for adaptive minus-sampling.
		autoborder = 1;
	}
	gtype = gt;
	ftype = ft;
	fpar = fpar0;
	dbg = dbg0;

	trans = trans0;
	if(*trans>0){
		graph->pp->calcTransWeights();
	}
}


SEXP Fun::toSEXP(SEXP pp)
//transform a std::vector<double > to SEXP
{
	SEXP res, *onevalue;
	double *p, *m;

	PROTECT(res = allocVector(VECSXP, this->value.size()));
	int i,j;
	for(i=0;i < (int)this->value.size();i++)
	{
		onevalue = new SEXP;
		PROTECT(*onevalue = allocVector(REALSXP, this->value.at(i).size()));
		p = REAL(*onevalue);

		for(j=0;j<(int)this->value.at(i).size();j++)
			p[j] = this->value.at(i).at(j);
		SET_VECTOR_ELT(res, i, *onevalue);
		UNPROTECT(1);
		delete onevalue;
	}
//	set the mass2 element of pp
	m = REAL(getListElement(pp,"mass2"));
	for(int i=0; i< graph->pp->size();i++)
		m[i] = graph->pp->getMass2(&i);
	UNPROTECT(1);
	return res;
}


void Fun::calculate()
{
	int i;
	std::vector<double> resvec;
	for(i=parvec.size()-1 ; i >= 0 ; i--)
	{
		if(*this->dbg)Rprintf("Fun %i/%i: graph[",(int)parvec.size()-i,(int)parvec.size());

		// update graph
		graph->par = &parvec[i];
		graph->sg_calc();
		*graph->oldpar = *graph->par;
		if(autoborder){
			if(*this->dbg)Rprintf("][minus]");
			updateInclude();
		}
		if(*this->dbg)Rprintf("Value[ ");
		// calc index
		if(*ftype == 1)
			resvec = mingling(graph, fpar, dbg, included);
		if(*ftype == 2)
			resvec = shannon(graph, fpar, dbg, included);
		if(*ftype == 3)
			resvec = simpson(graph, fpar, dbg, included);
		if(*ftype == 4)
			resvec = isar(graph, fpar, dbg, included);
		if(*ftype == 5)
			resvec = mci(graph, fpar, dbg, included);
		if(*ftype == 6)
			resvec = biomass(graph, fpar, dbg, included);
				value.at(i) = resvec;
		if(*this->dbg)Rprintf(" ]\n");
	}
}

void Fun::re_calculate()
{
		std::vector<double> resvec;
		if(autoborder) updateInclude();
		if(*this->dbg)Rprintf(" Value[ ");
		// calc index
		if(*ftype == 1)
			resvec = mingling(graph, fpar, dbg, included);
		if(*ftype == 2)
			resvec = shannon(graph, fpar, dbg ,included);
		if(*ftype == 3)
			resvec = simpson(graph, fpar, dbg, included);
		if(*ftype == 4)
			resvec = isar(graph, fpar, dbg, included);
		if(*ftype == 5)
			resvec = mci(graph, fpar, dbg, included);
		value.at(0) = resvec;
		if(*this->dbg)Rprintf(" ]\n");
}

void Fun::updateInclude() {
	for(int i=0; i<graph->pp->size();i++) {
		included[i]=1;
		if(graph->pp->getEdgeDist(&i) < *graph->par) included[i]=0;
	}
}
