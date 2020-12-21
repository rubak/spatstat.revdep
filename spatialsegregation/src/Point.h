#include <vector>
#ifndef POINT_H_
#define POINT_H_


class Point
{
	double x;
	double y;
	double z;
	double mass;
	double mass2;
	int id;
	int type;

	int cluster;
	std::vector<int> neighbours;

public:
	Point();
	Point(double , double);
	Point(double , double, double);
	virtual ~Point();

	double  getX();
	double  getY();
	double  getZ();
	int     getT();
	void	setT(int *);
	double  getMass();
	void    setMass(double *);

	int 	getId();
	void 	setId(int *);
	int		nsize();
	int 	getCluster();
	double  getMass2();
	void    setMass2(double *);
	int     getNeighbour(int *);
	void    addNeighbour(int *i);
	void    removeNeighbour(int *i);
	void    clearNeighbourhood();
	void    setCluster(int *i);

	void    move(double*, double*);
	void	move(double*, double*, double*);
};

#endif /*POINT_H_*/
