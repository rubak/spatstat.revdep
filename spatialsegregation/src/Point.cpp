#include "Point.h"

/********************************************************************************************/
Point::Point()
{
}
/********************************************************************************************/
Point::Point(double x0, double y0)
{
	this->x = x0;
	this->y = y0;
	this->z = 0;
	this->type = 0;
	this->neighbours.clear();
}
/********************************************************************************************/
Point::Point(double x0, double y0, double z0)
{
	this->x = x0;
	this->y = y0;
	this->z = z0;
	this->type = 0;
	this->neighbours.clear();
}
/********************************************************************************************/
Point::~Point()
{
}
/********************************************************************************************/

double  Point::getX() {return this->x;}
double  Point::getY() {return this->y;}
double  Point::getZ() {return this->z;}
int     Point::getT() {return this->type;}
void    Point::setT(int *t) {this->type = *t;}
void    Point::setMass(double *x){this->mass = *x;}
double  Point::getMass(){return this->mass;}
/********************************************************************************************/
int 	Point::nsize() {return this->neighbours.size();}
int     Point::getCluster(){return this->cluster;}
void    Point::clearNeighbourhood(){this->neighbours.clear();}
int 	Point::getId(){return this->id;}
void 	Point::setId(int *i){this->id = *i;}
void    Point::setCluster(int *i){this->cluster = *i;}
void    Point::setMass2(double *x){this->mass2 = *x;}
double  Point::getMass2(){return this->mass2;}
/********************************************************************************************/
void    Point::move(double *x, double *y){this->x =*x;this->y =*y;}
void    Point::move(double *x, double *y, double *z)
{
	this->x =*x;this->y =*y;
	this->z =*z;
}
/********************************************************************************************/
int 	Point::getNeighbour(int *i)
{
	if(*i<this->nsize())
		return this->neighbours.at(*i);
	else
		return -1;
}

/********************************************************************************************/
void    Point::addNeighbour(int *i)
{
	int j,isnew=1;
	for(j=0;j<(int)this->neighbours.size();j++)
		if(this->neighbours.at(j)==*i)
		{
			isnew = 0;
			break;
		}
	if(isnew & (*i != this->id))this->neighbours.push_back(*i);
}
/********************************************************************************************/
void    Point::removeNeighbour(int *i)
{
	int j, loc = -1;
	std::vector<int> *newneighs;

	for(j=0; j < (int)this->neighbours.size(); j++)
		if(this->neighbours.at(j)==*i)
		{
			loc = j;
			break;
		}
	if(loc > -1)
	{
		newneighs = new std::vector<int>;
		newneighs->clear();
		for(j=0;j<loc;j++)
			newneighs->push_back(this->neighbours.at(j));

		for(j=loc+1;j<(int)this->neighbours.size();j++)
			newneighs->push_back(this->neighbours.at(j));
		this->neighbours.swap(*newneighs);
		delete newneighs;
	}
}
/********************************************************************************************/

