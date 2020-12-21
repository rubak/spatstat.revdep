 
#ifndef MORPHO_H_
#define MORPHO_H_
  
#include <R.h>
#include <vector>
#include "Graph.h"
  
std::vector<std::vector<double> >  morphoArcs(Graph *graph);
void morphoArcsMinus(std::vector<std::vector<double> > *arcs, std::vector<double> arc);

#endif /* MORPHO_H_ */
