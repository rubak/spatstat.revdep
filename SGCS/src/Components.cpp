/*
 * Components.cpp
 *
 *  Created on: 19.9.2008
 *      Author: tarajala
 */

#include "Components.h"

Components::Components() {
  pconn = &Components::computeConnected;
}

Components::~Components() {

}

void Components::calculate(Graph *graph)
{
	int h, i, j, k, l, g, isnew, loop, sizei, sizeg, d, x;
	n = graph->size(); // internal n

	std::vector<std::vector<int> > complist;
	complist.resize(n);

	if(graph->dbg) Rprintf("grouping, ");
	for(i=0; i < n; i++)
	{
		complist.at(i).clear();
		for(j=0; j < (int)graph->nodelist.at(i).size(); j++)
		{

			complist.at(i).push_back( graph->nodelist.at(i).at(j)-1 );
		}
	}

	loop=1;
	d=0;
	if(graph->dbg) Rprintf("sorting, ");
	i=0;
	while(loop)
	{
		sizei = complist.at(i).size();

		if(sizei - d > 0 ) //unvisited neigh's left
		{
			for(k=d; k < sizei; k++ )
			{
				g = complist.at(i).at(k);
				sizeg = complist.at(g).size();
				for(j = sizeg-1; j>=0; j--) //start union of c(i) U c(g) => c(i)
				{
					isnew=1;
					h = complist.at(g).at(j);//from rear, for .pop_back()
					if( h != i)
					{
						for(l = 0; l< (int)complist.at(i).size(); l++) //is it new?
						{

							if(complist.at(i).at(l) == h) isnew=0;
						}
						if(isnew > 0)
							complist.at(i).push_back(h);
					}
					complist.at(g).pop_back();
				}// end of union
				complist.at(g).clear();
				d++;
			}
		}
		else // all neigh's of neighs i collected to c(i)
		{
			d=0;
			complist.at(i).push_back(i);
			i++;
		}
		if(i>=n)loop=0;
	}; // eo clustering loop

	std::vector<int>  *p;

	componentlist.clear();

	if(graph->dbg) Rprintf("cleaning");
	x=0;
	for(i=0; i< n ; i++)
	{
		p = new std::vector<int>;
		p->resize(0);
		if(complist.at(i).size()>0)
		{
			x++;
			for(j=0; j < (int)complist.at(i).size(); j++)
			{
				p->push_back(complist.at(i).at(j));
				complist.at(complist.at(i).at(j)).clear();

			}
			componentlist.push_back(*p);
		}
//		delete p;
	}
	if(graph->dbg) Rprintf(" (x=%i)ok ",x);
  
  
}


int Components::computeConnected(int *i, int *j)
{
	int k,l,m;
	for(k = 0; k < (int) componentlist.size(); k++)
	{
		for(l=0; l < (int) componentlist.at(k).size();l++)
			if( componentlist.at(k).at(l)== (*i) )
				for (m = 0; m < (int)componentlist.at(k).size(); ++m)
				{
					if(componentlist.at(k).at(m) == (*j))	{return 1;}
				}
	}
	return 0;
}
/********************************************************************************************/
int Components::connectionsPrecalculated(int *i, int *j){
  if(*i==*j) return 0.0;
	if(*i>*j) return connectionsPrecalculated(j, i);
	return connectionTriangle[ *j-*i -1 + (int)((*i)*n-(*i)*(*i+1)/2) ];  
}

/********************************************************************************************/
int Components::connected(int *i, int *j){
  return (this->*pconn)(i,j);
}

/********************************************************************************************/
void Components::preComputeConnected() {
  connectionTriangle.resize( n*(n-1)/2);
  int i,j;
  for(i=0; i < n-1;i++)
    for(j=i+1; j < n; j++) connectionTriangle[ j-i -1 + (int)((i)*n-(i)*(i+1)/2) ] = connected(&i,&j);
  pconn = &Components::connectionsPrecalculated;
}


/********************************************************************************************/
SEXP Components::toSEXP()
//transform a std::vector<std::vector<int> > to SEXP, desctructive
{
	SEXP comps, *node;

	PROTECT(comps = allocVector(VECSXP, this->componentlist.size()));
	int i,j, *p, n;
	for(i=0;i< (int)this->componentlist.size();i++)
	{
		node = new SEXP;
		PROTECT(*node = allocVector(INTSXP, this->componentlist.at(i).size()) );
		p = INTEGER(*node);
		n = this->componentlist.at(i).size();
		if(n<1) ;//p[0]=NULL;
		else
			for(j=0;j<n;j++)
			{
				p[j] = (int) this->componentlist.at(i).at(j);
			};
		this->componentlist.at(i).clear();
		SET_VECTOR_ELT(comps, i, *node);
		UNPROTECT(1);
	}
	UNPROTECT(1);
	return comps;
}
