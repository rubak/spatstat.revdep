#include "Graph.h"

/********************************************************************************************/
Graph::~Graph()
{
}
/********************************************************************************************/
Graph::Graph(Pp *pp0, 
                 int typ,
                 double par0,
                 double prepR0,
                 int saveMemory0,
                 int dbg0)
{
	pp = pp0;
	par = par0;
  oldpar = par;
  saveMemory = saveMemory0;
  
  prepR=prepR0;
	
	nodelist.resize(pp->size());

  given = 0;
  
  dbg = dbg0;

	gtype = typ;
}
/********************************************************************************************/
void Graph::setNodelist(SEXP preGraph)
{
  std::vector<std::vector<int> > prepNodelist;
  VectsxpToVector(getListElement(preGraph,"edges"), prepNodelist);
  if(dbg) Rprintf("Restoring given edges...");
	nodelist.clear();
	for(int i=0;i<(int) prepNodelist.size();i++)
		nodelist.push_back(prepNodelist.at(i));
	par = REAL(getListElement(preGraph,"parameters"))[0];
  given = 1;
	if(dbg) Rprintf("ok. ");
}
/********************************************************************************************/
double Graph::Dist(int *i, int *j) {return pp->getDistance(i,j);}
int    Graph::size()               {return nodelist.size();}
/********************************************************************************************/
//The graph methods
/********************************************************************************************/

void Graph::sg_calc()
{
	// preprocess if requested
	if(prepR>0 && oldpar<= par )
	{
		if(dbg) Rprintf("Preprocessing[");
		this->sg_geometric(prepR);
		if(dbg) Rprintf("] ok. ");
	}
	if(gtype==0)
	{
		if(oldpar > par)
			this->sg_shrink_geometric(par);
		else
			this->sg_geometric();
	}
	else if(gtype==1)
	{
		if(oldpar > par)
			this->sg_shrink_knn();
		else
			this->sg_knn();
	}
}

/********************************************************************************************/
void Graph::sg_geometric()
{
 Graph::sg_geometric(par);
}

void Graph::sg_geometric(double R)
{
	if(dbg) Rprintf("Geometric (R=%f):",R);
	int i,j;
	double dist;
	for(i=0;i<(pp->size()-1);i++)
		for(j=i+1;j<pp->size();j++)
		{
			dist = Dist(&i,&j);
			if(dist<=R){
				nodelist[i].push_back(j+1);
				nodelist[j].push_back(i+1);
			}
		}
	if(dbg) Rprintf(" Ok.");
}

void Graph::sg_shrink_geometric(double R)
{
	if(dbg) Rprintf("Geometric (R=%f) (shrinking):",R);
	int i,j,j0;
	double dist;
	std::vector<int> *node;
	for(i=0; i < pp->size() ; i++)
	{
		node = new std::vector<int>;
		for(j=0;j < (int)this->nodelist[i].size() ; j++)
		{
			j0 = nodelist[i][j]-1;
			dist = Dist(&i,&j0);
			if(dist<=R)
				node->push_back(j0+1);
		}
		nodelist[i].clear();
		for (j = 0; j < (int)node->size(); ++j) this->nodelist[i].push_back(node->at(j));
		delete node;
	}
	if(dbg) Rprintf(" Ok.");
}
/********************************************************************************************/
void Graph::sg_knn()
{

	int i,j,l, k;
	k = (int) par;
	std::vector<int> *node;

	if(prepR==0)// if not preprocessed
	{
		if(dbg) Rprintf("%i-nn: ", k);
		int n = pp->size();
		double *dists2_i, *dists2_i2;
		dists2_i = new double[n];
		dists2_i2 = new double[n];
		for(i=0;i<pp->size();i++) //for each point
		{
			for(j=0;j<pp->size();j++) dists2_i2[j]=dists2_i[j]= Dist(&i,&j); //gather the distances to others
			qsort( dists2_i, pp->size(), sizeof(double), compare_doubles); // sort distances, rising
			for(j=1;j <= k;j++) // find the k nearest
				for(l=0;l<pp->size();l++)
					if( dists2_i[j] == dists2_i2[l] ) //with distance comparison
					{
						nodelist[i].push_back(l+1);
						break;
					}

		}
	}
	else{ //preprocessed
		if(dbg) Rprintf("%i-nn (shrinking):", k);
		double *dists2_i, *dists2_i2;
		for(i=0;i<pp->size();i++) //for each point
		{
			node = new std::vector<int>;
			dists2_i = new double [nodelist[i].size()];
			dists2_i2 = new double [nodelist[i].size()];
			if((int)nodelist[i].size()<k){ Rprintf("\n preprocessing R too small, not enough neighbours (point #%i)!!\n",i+1); return;}
			for(l=	0;l< (int)nodelist[i].size();l++)
			{
				j = nodelist[i][l]-1;
				dists2_i2[l]=Dist(&i,&j); //gather the distances to others, given preprocessing
				dists2_i[l]=dists2_i2[l];
			}
			qsort( dists2_i, nodelist[i].size() , sizeof(double),compare_doubles); // sort distances, rising
			for(j=0;j<k;j++) // find the k nearest
				for(l=0;l<(int)nodelist[i].size();l++)
					if( dists2_i[j] == dists2_i2[l] ) //with distance comparison
					{
						node->push_back(nodelist[i][l]);
						break;
					}
			nodelist[i].clear();nodelist[i].resize(0);
			for(j=0;j < (int)node->size();j++) nodelist[i].push_back( (*node)[j] );
			delete node;
			delete[] dists2_i;
			delete[] dists2_i2;
		}
	}
	 if(dbg) Rprintf(" Ok.");
}

void Graph::sg_shrink_knn()
{

	double R0=prepR, R=1;
	prepR = R;
	this->sg_knn();
	prepR = R0;
}

/********************************************************************************************/


int compare_doubles(const void *a, const void *b)
{
  const double *da = (const double *) a;
  const double *db = (const double *) b;

  return (*da > *db) - (*da < *db);
}

// EOF
