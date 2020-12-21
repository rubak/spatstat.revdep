/* Functional class for the segregation measures
 * ftype=
 *        1   mingling
 *        2   shannon
 *        3   simpson
 *        4   ISAR
 *        5   MCI
 *        6   biomass sum
 *
 * Supports Geometric and k-nn graphs, and toroidal correction
 * TODO: border correction
 * by: Tuomas Rajala
 *
 * 	280410
 *
 */
#include <R.h>
#include <vector>
#include "Graph.h"
#include "mingling.h"
#include "shannon.h"
#include "simpson.h"
#include "isar.h"
#include "mci.h"
#include "mean_sd.h"
#include "biomass.h"
#ifndef FUN_H_
#define FUN_H_


class Fun
{
	Graph *graph;
	std::vector<std::vector <double> > value;
	std::vector<double> parvec;

	int *gtype; // 0 = geometric, 1 = knn
	int *ftype;
	int *included;
	int autoborder; // should we adapt the minus-border correction
	int *trans; // translation weights should be computed
	double *fpar;

	int *dbg;
	void updateInclude();
public:
	Fun();
	virtual ~Fun();
	void Init(Graph *g0, double *par0, int *parn, int *gt, int *ft, double *fpar, int *trans0, int *included0, int *dbg0);
	void calculate();
	void re_calculate();
	SEXP toSEXP(SEXP);
};

#endif /*FUN_H_*/
