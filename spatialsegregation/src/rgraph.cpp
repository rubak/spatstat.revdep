#include <R.h>
#include "Graph.h"
#include "Pp.h"


extern "C" {

SEXP graph_c(SEXP Args)
{
	Pp pp;
	double *prepR, *par, *d0, *dm1;
	int *gtype, *toroidal, *dbg, *incl, *i0;
	Graph graph;
	d0 = new double;
	dm1 = new double;
	d0[0] =0.0;
	dm1[0] = -1.0;
	i0 = new int;
	i0[0]=0;
//start parsing the args
	Args = CDR(Args);
	pp.Init(CAR(Args)); // init pp

	Args = CDR(Args);
	gtype = INTEGER(CAR(Args)); //what type of graph

	Args = CDR(Args);
	par = REAL(CAR(Args)); // graph par

	Args = CDR(Args);
	prepR = REAL(CAR(Args)); // if preprocessing

	Args = CDR(Args);
	toroidal = INTEGER(CAR(Args)); // if toroidal correction

	Args = CDR(Args);
	incl = INTEGER(CAR(Args)); // inclusion vector

	Args = CDR(Args);
	dbg = INTEGER(CAR(Args)); // if debug messages

//	Pp *pp0, int *gtype0, double *par0, double *prepR0, int *doDists0, double *preDists, int *toroidal0, int *inc0, double *wMatrix, int *dbg0
//	graph.Init(&pp, gtype, par, prepR , doDists, d0, toroidal, incl, weightMatrix, dbg);

	graph.Init(&pp, gtype, par, prepR,       i0, dm1, toroidal, incl,          dm1, dbg);
	graph.sg_calc();

	if(*dbg)Rprintf("\n");
	return graph.toSEXP();

}


}
