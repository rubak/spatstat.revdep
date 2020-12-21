#ifndef SHANNON_H_
#define SHANNON_H_

#include "Graph.h"
#include <vector>

std::vector<double> shannon0(Graph *graph, double *fpar, int *dbg, int *included);
std::vector<double> shannon(Graph *graph, double *fpar, int *dbg, int *included);
std::vector<double> piitauf(Graph *graph, double *fpar, int *dbg, int *included);
std::vector<double> shannon_v2(Graph *graph, double *fpar, int *dbg, int *included);
#endif /*SHANNON_H_*/
