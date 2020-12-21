/*
 * Components.h
 *
 *  Created on: 19.9.2008
 *      Author: tarajala
 */
#include <vector>
#include <R.h>
#include "Graph.h"

#ifndef COMPONENTS_H_
#define COMPONENTS_H_

class Components {
  int n;
  std::vector <int> connectionTriangle;
  int (Components::*pconn)(int*, int*);
  int connectionsPrecalculated(int *, int *);
  int computeConnected(int *, int *);
public:
	Components();
	virtual ~Components();
	std::vector<std::vector<int> > componentlist;
	void calculate(Graph *graph);
	int connected(int *, int *);
  void preComputeConnected();
  SEXP toSEXP();

};

#endif /* COMPONENTS_H_ */
