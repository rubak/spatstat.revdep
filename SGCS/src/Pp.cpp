//spatgraphs
#include <Rmath.h>
#include "Pp.h"
/********************************************************************************************/
Pp::~Pp()
{
}
/********************************************************************************************/
Pp::Pp(SEXP Argspp) {
	int i,j;
	  
  SEXP temp; 
  
  /// constants ///
  n = INTEGER(getListElement(Argspp, "n"))[0];
  dim = INTEGER(getListElement(Argspp, "dim"))[0];
  
  ///	set-up location data ///
  coordinates = REAL(getListElement(Argspp, "coord"));
  
  /// parse marks ////
  // type mark
  temp = getListElement(Argspp, "type");
  if(!Rf_isNull(temp)){
    types = INTEGER(temp);
    type_given = true;
  }else type_given = false;
  //  If point types given, collect them to a vector
  if(type_given){
    int old;
    typevec.clear();
    for(i=0;i < n ;i++) { // also attach the type to each point
    old = 0;
    for(j=0;j<(int)typevec.size();j++)
    if(typevec.at(j)==types[i]){ old = 1;break;}
    if(!old)
    typevec.push_back(types[i]);
    }
    ntypes = typevec.size();
  }else ntypes = 1;
  // real valued mark
  temp = getListElement(Argspp,"mass");
  if(!Rf_isNull(temp)){
    masses = REAL(temp);
    mass_given = true;
  } else mass_given = false;
  /// window ///
  // bounding box should be always present.
  bbox = REAL(getListElement(Argspp, "bbox"));
  // area
  temp = getListElement(Argspp, "area");
  if(Rf_isNull(temp)){
    windowArea = 1.0;
    for(i=0; i < dim ; i++) windowArea *= bbox[1+i*2]-bbox[i*2];
  }else windowArea = REAL(temp)[0]; // Precomputed window area, in case of more complicated than rectangle.
  
  // split bbox to variables for easy access
  xlim = new double(2); xlim[0] = bbox[0]; xlim[1] = bbox[1];
  ylim = new double(2); ylim[0] = bbox[2]; ylim[1] = bbox[3];
  if(dim==3){  zlim = new double(2); zlim[0] = bbox[4]; zlim[1] = bbox[5];  }  
  
  ///	intensities ///
	lambda = 0;
	if(ntypes > 1){
	  for(i=0; i < ntypes; i++) {
	    lambdas.push_back(0.0);
	    for(j=0; j < n; j++)
	    if(types[j]==i+1)
	    lambdas[i]=lambdas[i]+1.0;
	    lambdas[i]=lambdas[i]/windowArea;
	    lambda += lambdas[i];
	  }
	} else lambda = (double)n / windowArea;

  /// distance metric ///
  toroidal = 0; // for past.
	setToroidal(&toroidal);
  /// if distances precomputed
  temp = getListElement(Argspp, "pairwise_distances");
  if(!Rf_isNull(temp)){
    setDistances(REAL(temp));
  }
  /// translation weights ///
  temp = getListElement(Argspp, "weights");
  if(!Rf_isNull(temp)){  
    setWeights(REAL(temp));
  }else pweight = &Pp::weightAll1;
  // edge distance
  temp = getListElement(Argspp, "edgeDistances");
  if(!Rf_isNull(temp)){  
    setEdgeDistances(REAL(temp));
    pedgedist = &Pp::edgeDistancePrecalculated;
  }else pedgedist = &Pp::computeEdgeDistance;
}
/********************************************************************************************/
double Pp::distEuclidian(int *i, int *j)
{
  if(*i==*j) return 0.0;
		if(*i>*j) return distEuclidian(j, i);
			return 	sqrt(
					pow( getX(i)- getX(j)  , 2) +
					pow( getY(i)- getY(j)  , 2));
}
/********************************************************************************************/
double Pp::distToroidal(int *i, int *j)
{
  if(*i==*j) return 0.0;
		if(*i>*j) return distToroidal(j, i);
		return	sqrt(
	  pow( fmin2( xlim[1]-xlim[0]-fabs(getX(i)-getX(j)) , fabs(getX(i)-getX(j)) ) , 2) +
		  pow( fmin2( ylim[1]-ylim[0]-fabs(getY(i)-getY(j)) , fabs(getY(i)-getY(j)) ) ,2)   );
}
/********************************************************************************************/
double Pp::distEuclidian3(int *i, int *j)
{
	if(*i==*j) return 0.0;
		if(*i>*j) return distEuclidian3(j, i);
			return 	sqrt(
					pow( getX(i)- getX(j)  , 2) +
					pow( getY(i)- getY(j)  , 2) +
					pow( getZ(i)- getZ(j)  , 2)   );
}
/********************************************************************************************/
double Pp::distToroidal3(int *i, int *j)
{
  if(*i==*j) return 0.0;
		if(*i>*j) return distToroidal3(j, i);
		return	sqrt(
		  pow( fmin2( xlim[1]-xlim[0]-fabs(getX(i)-getX(j)) , fabs(getX(i)-getX(j)) ) , 2) +
		  pow( fmin2( ylim[1]-ylim[0]-fabs(getY(i)-getY(j)) , fabs(getY(i)-getY(j)) ) , 2) +
		  pow( fmin2( zlim[1]-zlim[0]-fabs(getZ(i)-getZ(j)) , fabs(getZ(i)-getZ(j)) ) , 2)   );
}
/********************************************************************************************/
double Pp::distPrecalculated(int *i, int *j)
{
	if(*i==*j) return 0.0;
	if(*i>*j) return distPrecalculated(j, i);
	return pairwiseDistanceTriangle[ *j-*i -1 + (int)((*i)*n-(*i)*(*i+1)/2) ];
}
/********************************************************************************************/
void Pp::calculateDistances()
{
	int i,j, k=0;
  pairwiseDistanceTriangle = new double [n*(n-1)/2];
  
	for(i=0; i < n-1;i++)
		for(j=i+1; j<n;j++)
		{
			pairwiseDistanceTriangle[k] = (this->*pdist)(&i, &j);
      k++;
		}
	pdist = &Pp::distPrecalculated;
}
/********************************************************************************************/
double Pp::getDistance(int *i, int *j)
{
	return (this->*pdist)(i,j);
}
/********************************************************************************************/
void Pp::setDistances(double *dvec)
{
  pairwiseDistanceTriangle = dvec;
  pdist = &Pp::distPrecalculated;
}

/********************************************************************************************/
bool Pp::included(int *i, double *r){
  return (getEdgeDistance(i) >= *r);
} // check border 

double Pp::getEdgeDistance(int *i){
  return (this->*pedgedist)(i);
}

double Pp::edgeDistancePrecalculated(int *i){
	return edgeDistanceVector[*i];
}

double Pp::computeEdgeDistance(int *i){
  double bx,by,bz;
  bx = fmin2(getX(i)-xlim[0], xlim[1]-getX(i));
  by = fmin2(getY(i)-ylim[0], ylim[1]-getY(i)); 
  bx = fmin2(bx,by);
  if(dim==3) {
    bz = fmin2(getZ(i)-zlim[0], zlim[1]-getZ(i));
    bx = fmin2(bx, bz);
  }
  return bx;
}

void Pp::setEdgeDistances(double *evec){
  edgeDistanceVector = evec;
  pedgedist = &Pp::edgeDistancePrecalculated;
}


/********************************************************************************************/
double Pp::weightAll1(int *i, int *j)
{
	return 1.0;
}
/********************************************************************************************/
double Pp::weightPrecalculated(int *i, int *j)
{
	if(*i==*j) return windowArea;
	if(*i>*j){ return weightPrecalculated(j, i); }
	else return pairwiseWeightTriangle[ *j-*i -1 + (int)((*i)*n-(*i)*(*i+1)/2) ];
}
/********************************************************************************************/
void Pp::setWeights(double *wvec)
{
  pairwiseWeightTriangle = wvec;
  pweight = &Pp::weightPrecalculated;
}
/********************************************************************************************/
double Pp::getWeight(int *i, int *j)
{
	return (this->*pweight)(i,j);
}
/********************************************************************************************/
double  Pp::getCoord(int *i, int *j){return coordinates[*i + (*j)*n];}
double  Pp::getX(int *i) {return coordinates[*i];}
double  Pp::getY(int *i) {return coordinates[*i+n];}
double  Pp::getZ(int *i) {return coordinates[*i+2*n];}
int     Pp::getType(int *i) {return types[*i];}
double  Pp::getMass(int *i){return masses[*i];}
int	    Pp::getTypevec(int *i){return this->typevec.at(*i);}
int 	  Pp::size()      {return this->n;   }
int     Pp::getDim()       {return this->dim;   }
int     Pp::getNtypes(){return this->ntypes;}
double  Pp::getBoundingBoxExtent(int d, int k){return bbox[k + 2*d];}
/********************************************************************************************/
void  Pp::setToroidal(bool *i){
  this->toroidal = *i;
  if(*i){
    pdist = &Pp::distToroidal;
    if(dim == 3) pdist = &Pp::distToroidal3;
  }
  else{
    pdist = &Pp::distEuclidian;
    if(dim == 3) pdist = &Pp::distEuclidian3;
  }
}

