#ifndef SIMPSON_H_
#define SIMPSON_H_

#include "Graph.h"
#include <vector>
std::vector<double> simpson_typewise(Graph *graph, double *fpar, int *dbg, int *included);
std::vector<double> simpson0(Graph *graph, double *fpar, int *dbg, int *included);
std::vector<double> simpson(Graph *graph, double *fpar, int *dbg, int *included);

#endif /*SIMPSON_H_*/
