#ifndef ISAR_H_
#define ISAR_H_

#include "Graph.h"
#include <vector>
std::vector<double> isar(Graph *graph, double *fpar, int *dbg, int *included);
std::vector<double> isar_normal(Graph *graph, double *fpar, int *dbg, int *included);
std::vector<double> isar_normal_empty(Graph *graph, double *fpar, int *dbg, int *included);
std::vector<double> isar_wdeg(Graph *graph, double *fpar, int *dbg, int *included);
std::vector<double> isar_markweighted(Graph *graph, double *fpar, int *dbg, int *included);


#endif /*ISAR_H_*/
