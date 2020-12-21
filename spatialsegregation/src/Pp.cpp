//spatgraphs
#include <Rmath.h>
#include "Pp.h"
/********************************************************************************************/
Pp::Pp()
{
}
/********************************************************************************************/
Pp::~Pp()
{
}
/********************************************************************************************/
void Pp::Init(SEXP Argspp)
{
	int i,j, *type, old;
	double *x, *y, *z, *mass, *la;
	Point *p;

	tor = 0;
	toroidal = &tor;

	m = length(getListElement(Argspp, "x"));
	x = REAL(getListElement(Argspp, "x"));
	y = REAL(getListElement(Argspp, "y"));
	z = REAL(getListElement(Argspp, "z"));
	type = INTEGER(getListElement(Argspp, "types"));
	la = REAL(getListElement(Argspp, "area"));
	mass = REAL(getListElement(Argspp,"mass"));
	bdist = REAL(getListElement(Argspp,"bdist"));
	windowArea = *la;

//	set points
	this->points.clear();
	for(i=0; i < this->size(); i++)
	{
		p = new Point(x[i], y[i], z[i]);
		p->setT(&type[i]);
		p->setMass(&mass[i]);
		points.push_back(*p);
		delete p;
	}
	m = points.size();

//	collect types into a vector
	typevec.clear();
	for(i=0;i < m ;i++)
	{
		old = 0;
		for(j=0;j<(int)typevec.size();j++)
			if(typevec.at(j)==type[i]){ old = 1;break;}
		if(!old)
			typevec.push_back(type[i]);
	}
	ntypes = typevec.size();

// window
	xlim = REAL(getListElement(getListElement(Argspp, "window") ,"x"));
	ylim = REAL(getListElement(getListElement(Argspp, "window") ,"y"));
	zlim = REAL(getListElement(getListElement(Argspp, "window") ,"z"));
//	intensities
	lambda = 0;
	for(i=0;i<ntypes; i++)
	{
		lambdas.push_back(0.0);
		for(j=0;j<m;j++)
			if(type[j]==i+1)
				lambdas[i]=lambdas[i]+1.0;
		lambdas[i]=lambdas[i]/windowArea;
		lambda += lambdas[i];
	}

// distance
	dist = &Pp::distEuclidian;
//	edge distances for border correction: For rectangular window only.
	edgeDistp = &Pp::edgeDist;
//	weights
	weight = &Pp::weightAll1;
}
/********************************************************************************************/
double Pp::distEuclidian(int *i, int *j)
{
	if(*i==*j) return 0.0;
		if(*i>*j) return distEuclidian(j, i);
		if(*toroidal)
			return	sqrt(
						pow( fmin2( xlim[1]-xlim[0]-fabs(points.at(*i).getX()-points.at(*j).getX()) , fabs(points.at(*i).getX()-points.at(*j).getX()) ) ,2) +
						pow( fmin2( ylim[1]-ylim[0]-fabs(points.at(*i).getY()-points.at(*j).getY()) , fabs(points.at(*i).getY()-points.at(*j).getY()) ) ,2) +
						pow( fmin2( zlim[1]-zlim[0]-fabs(points.at(*i).getZ()-points.at(*j).getZ()) , fabs(points.at(*i).getZ()-points.at(*j).getZ()) ) ,2)   );
		else
			return 	sqrt(
					pow( points.at(*i).getX()- points.at(*j).getX()  ,2) +
					pow( points.at(*i).getY()- points.at(*j).getY()  ,2) +
					pow( points.at(*i).getZ()- points.at(*j).getZ()  ,2)   );
}
/********************************************************************************************/
double Pp::distPrecalculated(int *i, int *j)
{
	if(*i==*j) return 0.0;
	if(*i>*j) return distPrecalculated(j, i);
	return distTriangle.at( *j-*i -1 + (int)((*i)*m-(*i)*(*i+1)/2) );
}
/********************************************************************************************/
void Pp::calcDists()
{
	int i,j;
	for(i=0; i < m-1;i++)
		for(j=i+1; j<m;j++)
		{
			distTriangle.push_back(distEuclidian(&i, &j));
		}
	dist = &Pp::distPrecalculated;
}
/********************************************************************************************/
double Pp::getDist(int *i, int *j)
{
	return (this->*dist)(i,j);
}
/********************************************************************************************/
double Pp::getEdgeDist(int *i)
{
	return (this->*edgeDistp)(i);
}
/********************************************************************************************/
void Pp::setDist(int *i, int *j, double d)
{
	if(*i>*j){ setDist(j, i, d);}
	else if(*i!=*j) distTriangle.at( *j-*i -1 + (int)((*i)*m-(*i)*(*i+1)/2) ) = d;
}
/********************************************************************************************/
void Pp::setDists(double *dvec)
{
	int i;
	distTriangle.resize(m*(m-1)/2);
	for(i=0; i < (int)distTriangle.size(); i++)
			distTriangle.at(i) = dvec[i];

	dist = &Pp::distPrecalculated;
}
/********************************************************************************************/
double Pp::weightAll1(int *i, int *j)
{
	return 1.0;
}
/********************************************************************************************/
double Pp::weightTrans(int *i, int *j)
{
	if(*i==*j) return windowArea;
	if(*i>*j){ return weightTrans(j, i); }
	else return weightTriangle.at( *j-*i -1 + (int)((*i)*m-(*i)*(*i+1)/2) );
}
/********************************************************************************************/

void Pp::setAllTransWeights(double d)
{
	int i;
	weightTriangle.resize(m*(m-1)/2);
	for(i=0; i < (int)weightTriangle.size(); i++)
			weightTriangle.at(i) = d;

	weight = &Pp::weightTrans;
}
/********************************************************************************************/
void Pp::calcTransWeights()
{
	int i,j;
	double d;
	weightTriangle.resize(m*(m-1)/2);
	for(i=0; i<(int)this->m-1; i++)
		for(j=i+1; j<(int)this->m; j++){
			d = (xlim[1]-fabs(getX(&i)-getX(&j)))*(ylim[1]-fabs(getY(&i)-getY(&j)))*(zlim[1]-fabs(getZ(&i)-getZ(&j)));
			setWeight(&i, &j, d);
		}
	weight = &Pp::weightTrans;
}
/********************************************************************************************/
//void Pp::calcEdgeDists() {
//	distEdge.clear();
//	for(int i=0; i<this->m;i++)
//		distEdge.push_back(edgeDist(&i));
//	edgeDistp = &Pp::edgeDistPrecalculated;
//}
void Pp::calcEdgeDists() {
	distEdge.clear();
	for(int i=0; i<this->m;i++)
		distEdge.push_back(bdist[i]);
	edgeDistp = &Pp::edgeDistPrecalculated;
}
/********************************************************************************************/
double Pp::edgeDist(int *i) {
	return (double)fmin2(fmin2(xlim[1]-getX(i), getX(i)-xlim[0]), fmin2(ylim[1]-getY(i),getY(i)-ylim[0]));
}
double Pp::edgeDistPrecalculated(int *i) {
	return distEdge.at(*i);
}
/********************************************************************************************/
void Pp::setWeight(int *i, int *j, double d)
{
	if(*i>*j){ setWeight(j, i, d);}
	else weightTriangle.at( *j-*i -1 + (int)((*i)*m-(*i)*(*i+1)/2) ) = d;
}
/********************************************************************************************/
double Pp::getWeight(int *i, int *j)
{
	return (this->*weight)(i,j);
}
/********************************************************************************************/
double  Pp::getX(int *i) {return this->points[*i].getX();}
double  Pp::getY(int *i) {return this->points[*i].getY();}
double  Pp::getZ(int *i) {return this->points[*i].getZ();}
int     Pp::getT(int *i) {return this->points[*i].getT();}
int	    Pp::getTypevec(int *i){return this->typevec.at(*i);}
void 	Pp::setToroidal(int *i){this->tor = *i;}
double  Pp::getMass(int *i){return this->points[*i].getMass();}
void    Pp::setMass(int *i, double *x){this->points[*i].setMass(x);}
int 	Pp::size()      {return this->m;   }
int     Pp::nsize(int *i){return this->points[*i].nsize();}
int     Pp::getCluster(int *i){return this->points[*i].getCluster();}
double  Pp::getMass2(int *i){return this->points[*i].getMass2();}
void    Pp::setMass2(int *i, double *x){this->points[*i].setMass2(x);}
int     Pp::getNtypes(){return this->ntypes;}
/********************************************************************************************/

int Pp::Empty(int *i, int *j, int *k)
// check if the circumcircle of three point triangle is empty of other points.
// See: http://mathworld.wolfram.com/Circumcircle.html
{
	int l;
	double x0,y0,R2,bx,by,a,c,d2,xxyy1,xxyy2,xxyy3,x13,x23,x21,y13,y21,y23;
	xxyy1 = points.at(*i).getX()*points.at(*i).getX()+points.at(*i).getY()*points.at(*i).getY();
	xxyy2 = points.at(*j).getX()*points.at(*j).getX()+points.at(*j).getY()*points.at(*j).getY();
	xxyy3 =	points.at(*k).getX()*points.at(*k).getX()+points.at(*k).getY()*points.at(*k).getY();
	y23 = points.at(*j).getY()-points.at(*k).getY();
	y13 = points.at(*i).getY()-points.at(*k).getY();
	y21 = points.at(*j).getY()-points.at(*i).getY();
	x23 = points.at(*j).getX()-points.at(*k).getX();
	x13 = points.at(*i).getX()-points.at(*k).getX();
	x21 = points.at(*j).getX()-points.at(*i).getX();
	bx = -( xxyy1*y23-xxyy2*y13-xxyy3*y21 );
	by =  ( xxyy1*x23-xxyy2*x13-xxyy3*x21 );
	a  = points.at(*i).getX()*y23-points.at(*j).getX()*y13-points.at(*k).getX()*y21;
	c  = -(xxyy1*(points.at(*j).getX()*points.at(*k).getY()-points.at(*j).getY()*points.at(*k).getX())
			-xxyy2*(points.at(*i).getX()*points.at(*k).getY()-points.at(*i).getY()*points.at(*k).getX())
			-xxyy3*(points.at(*j).getX()*points.at(*i).getY()-points.at(*i).getX()*points.at(*j).getY()));
	R2 = (bx*bx+by*by-4.0*a*c)/(4.0*a*a);
	x0 = -bx/(2.0*a);
	y0 = -by/(2.0*a);
	for(l=0;l<m;l++)
	{
		if( (l!=*i) & (l!=*j) & (l!=*k))
		{
			d2 = pow(x0-points.at(l).getX(),2)+pow(y0-points.at(l).getY(),2);
			if(d2<R2) return 0;
		}
	}
	return 1;
}
/********************************************************************************************/
int Pp::EmptyConstrained(int *i, int *j, int *k, std::vector<int> * node)
// check if the circumcircle of three point triangle is empty of other points in the neighbourhood
// See: http://mathworld.wolfram.com/Circumcircle.html
{
	int l, no;
	double x0,y0,R2,bx,by,a,c,d2,xxyy1,xxyy2,xxyy3,x13,x23,x21,y13,y21,y23;
	xxyy1 = points.at(*i).getX()*points.at(*i).getX()+points.at(*i).getY()*points.at(*i).getY();
	xxyy2 = points.at(*j).getX()*points.at(*j).getX()+points.at(*j).getY()*points.at(*j).getY();
	xxyy3 =	points.at(*k).getX()*points.at(*k).getX()+points.at(*k).getY()*points.at(*k).getY();
	y23 = points.at(*j).getY()-points.at(*k).getY();
	y13 = points.at(*i).getY()-points.at(*k).getY();
	y21 = points.at(*j).getY()-points.at(*i).getY();
	x23 = points.at(*j).getX()-points.at(*k).getX();
	x13 = points.at(*i).getX()-points.at(*k).getX();
	x21 = points.at(*j).getX()-points.at(*i).getX();
	bx = -( xxyy1*y23-xxyy2*y13-xxyy3*y21 );
	by =  ( xxyy1*x23-xxyy2*x13-xxyy3*x21 );
	a  = points.at(*i).getX()*y23-points.at(*j).getX()*y13-points.at(*k).getX()*y21;
	c  = -(xxyy1*(points.at(*j).getX()*points.at(*k).getY()-points.at(*j).getY()*points.at(*k).getX())
			-xxyy2*(points.at(*i).getX()*points.at(*k).getY()-points.at(*i).getY()*points.at(*k).getX())
			-xxyy3*(points.at(*j).getX()*points.at(*i).getY()-points.at(*i).getX()*points.at(*j).getY()));
	R2 = (bx*bx+by*by-4.0*a*c)/(4.0*a*a);
	x0 = -bx/(2.0*a);
	y0 = -by/(2.0*a);
	for(no=0; no<(int)node->size(); no++)
	{
		l = node->at(no)-1;
		if( (l!=*i) & (l!=*j) & (l!=*k))
		{
			d2 = pow(x0-points.at(l).getX(),2)+pow(y0-points.at(l).getY(),2);
			if(d2<R2) return 0;
		}
	}
	return 1;
}



// old init
/********************************************************************************************/
void Pp::Init(double *x, double *y, double *z, int *type, double *mass, int *n, double *xlim0, double *ylim0, double *zlim0)
{
	int i,j,old;
	Point *p;
	m = *n;
	std::vector<int> temp;
	this->points.clear();
	for(i=0;i<m;i++)
	{
		old = 0;
		for(j=0;j<(int)temp.size();j++)
			if(temp.at(j)==type[i]){ old = 1;break;}
		if(!old)
			temp.push_back(type[i]);
		p = new Point(x[i], y[i], z[i]);
		p->setT(&type[i]);
		p->setMass(&mass[i]);
		points.push_back(*p);
		delete p;

	}
	m = points.size();
	ntypes = temp.size();

	xlim = xlim0;
	ylim = ylim0;
	zlim = zlim0;
}
/********************************************************************************************/
