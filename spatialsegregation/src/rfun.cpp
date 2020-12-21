#include <R.h>
#include "Fun.h"



extern "C" {

//	res<-.External("fun_c", as.integer(dbg), X, as.numeric(funpars),
//					as.integer(ntypei), as.numeric(parvec),
//					as.integer(funi), as.integer(toroidal),
//					as.numeric(prepRange), as.integer(doDists),
//					as.integer(included), prepGraph,
//					as.integer(prepGraphIsTarget),
//					as.numeric(weightMatrix),
//					PACKAGE="spatialsegregation")

SEXP fun_c(SEXP Args)
{

	Pp pp;
	Graph graph;
	Fun fun;
	double *prepR, *fpar, *par, *parvec, *d0, *weightMatrix;
	int *gtype, *doDists, *toroidal, *ftype, *dbg, parn, *incl, prepG=0, *prepGraphIsTarget, *translate;
	SEXP prepGraph, origpp;

	d0 = new double;
	d0[0] =-1.0;
//start parsing the args
	Args = CDR(Args);
	dbg = INTEGER(CAR(Args)); // if debug messages
	if(*dbg)Rprintf("Parsing parameter:");

	Args = CDR(Args);
	origpp = CAR(Args);
	pp.Init(CAR(Args)); // init pp

	if(*dbg)Rprintf(".");
	Args = CDR(Args);
	fpar = REAL(CAR(Args)); // additional function parameters.

	if(*dbg)Rprintf(".");
	Args = CDR(Args);
	gtype = INTEGER(CAR(Args)); //what type of graph

	if(*dbg)Rprintf(".");
	Args = CDR(Args);
	parvec = REAL(CAR(Args)); // graph parvec
	parn = length(CAR(Args));

	if(*dbg)Rprintf(".");
	Args = CDR(Args);
	ftype = INTEGER(CAR(Args)); //what type of function

	if(*dbg)Rprintf(".");
	Args = CDR(Args);
	toroidal = INTEGER(CAR(Args)); // if toroidal correction

	if(*dbg)Rprintf(".");
	Args = CDR(Args);
	prepR = REAL(CAR(Args)); // if preprocessing R given

	if(*dbg)Rprintf(".");
	Args = CDR(Args);
	doDists = INTEGER(CAR(Args)); // if the point-to-point distances should be precalculated and stored

	if(*dbg)Rprintf(".");
	Args = CDR(Args);
	translate = INTEGER(CAR(Args)); // if translation weights should be used in mingling&simpson
	if(*dbg)Rprintf(".");

	Args = CDR(Args);
	incl = INTEGER(CAR(Args)); // the point inclusion vector

	if(*dbg)Rprintf(".");
	Args = CDR(Args);
	prepGraph = CAR(Args); // possibly precalculated graph
	prepG = 1- INTEGER(getListElement(prepGraph,"isnull"))[0];

	if(*dbg)Rprintf(".");
	Args = CDR(Args);
	prepGraphIsTarget = INTEGER(CAR(Args)); // use only the precalculated graph to get a value

	if(*dbg)Rprintf(".");
	Args = CDR(Args);
	weightMatrix = REAL(CAR(Args)); // weightMatrix


	if(*dbg)Rprintf("done.\n");
	par = &parvec[parn-1];



	if(*dbg)Rprintf("Init graph...");
	graph.Init(&pp, gtype, par, prepR , doDists, d0, toroidal, incl, weightMatrix, dbg);
	if(prepG)// if precalculated graph, set the edges, overrule the prepR-parameter
	{
		graph.setNodelist(prepGraph);
		*graph.oldpar = *par+1;
	}

	if(*dbg)Rprintf("Init fun...");
	fun.Init(&graph, parvec, &parn, gtype, ftype, fpar,translate, incl, dbg);
	if(*dbg)Rprintf("done.\n");


	if(*dbg)Rprintf("Calculating:\n");
	if(*prepGraphIsTarget){ fun.re_calculate(); }
	else{  fun.calculate();  }

	if(*dbg)Rprintf("done.\n");
	return fun.toSEXP(origpp);
}



} //extern
