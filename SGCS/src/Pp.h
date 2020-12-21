// point pattern class

#include <Rdefines.h>
#include <Rinternals.h>
#include <R.h>
#include <vector>
#include "Rextras.h"
#ifndef PP_H_
#define PP_H_

class Pp {
  double *coordinates;
  double *masses;
  int *types;
  int n;
  int dim;
  double windowArea;
  double *bbox, *xlim, *ylim, *zlim;
  
  bool toroidal;
  
  double lambda;
  std::vector<double > lambdas; // for multitype patterns
  std::vector<int > typevec; // for multitype patterns
  int ntypes;
  bool type_given, mass_given;
  
  // for precomputed distances and weights
  double *edgeDistanceVector;
  double *pairwiseDistanceTriangle;
  double *pairwiseWeightTriangle;
  
  double (Pp::*pdist)(int*, int*);
  double distEuclidian(int*, int*);
	double distEuclidian3(int*, int*);
  double distToroidal(int *i, int *j);
  double distToroidal3(int *i, int *j);
  double distPrecalculated(int*, int*);
  
  double (Pp::*pedgedist)(int *);
  double computeEdgeDistance(int *);
  double edgeDistancePrecalculated(int *);
	
  double (Pp::*pweight)(int*, int*);
	double weightAll1(int *, int *);
	double weightPrecalculated(int *, int *);

public:
  Pp(SEXP);
  virtual ~Pp();
  
  double getCoord(int *, int *);
  double getX(int *);
	double getY(int *);
	double getZ(int *);
  
	int    getType(int *);
  double getMass(int *);

	int	   getTypevec(int *);
	void   setToroidal(bool *);
	int    size();
	int    getDim();
  int    getNtypes();
  
  bool    included(int *, double *r); // check border 
  
  double getBoundingBoxExtent(int, int);
  
	double getDistance(int *, int *);
	void   setDistances(double *);
  void   calculateDistances();

  double getEdgeDistance(int *i);
  void   setEdgeDistances(double *);
	
  double getWeight(int *, int *);
	void   setWeights(double *);
  
	SEXP toSEXP();
  
};


#endif /*PP_H_*/
