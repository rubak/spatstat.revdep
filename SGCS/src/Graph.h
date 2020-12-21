#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <vector>
#include "Pp.h"

#ifndef GRAPH_H_
#define GRAPH_H_

class Graph
{
public:
	double par;
	double oldpar;
	double prepR;
	int    gtype;
	int    given;
  int    dbg;
  int    saveMemory;
	
  Pp *pp; //the point pattern
  std::vector<std::vector<int> > nodelist;
  Graph(Pp *pp0, int typ, double par, double prepR, int saveMemory, int dbg);
	virtual ~Graph();
	
  double Dist(int *, int *);
  int size();
	
  void setNodelist(SEXP);
  
  void sg_calc();
	void sg_geometric();
	void sg_geometric(double );
	void sg_shrink_geometric(double );
	void sg_knn();
	void sg_shrink_knn();
};

#endif /*GRAPH_H_*/


int compare_doubles(const void *a, const void *b);

