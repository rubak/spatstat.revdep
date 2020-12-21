#include <Rmath.h>
#include "Graph.h"

#ifndef MAX_DOUBLE
const double MAX_DOUBLE = 9999999;
#endif

Graph::Graph()
{
}
/********************************************************************************************/
Graph::~Graph()
{
}
/********************************************************************************************/
void Graph::Init(Pp *pp0, int *gtype0, double *par0, double *prepR0, int *doDists0,
		double *preDists, int *toroidal0, int *inc0, double *wMatrix, int *dbg0)
{
	int typein,i,j;
	if(*dbg0)Rprintf("intializing graph-object... ");

	pp = pp0;
	par=par0;
	prepR=prepR0;
	preEdges = 0;
	opar = *par;
	oldpar = &opar;
	doDists=doDists0;
	dbg=dbg0;
	inc = inc0;
	nodelist.resize(pp->size());
	gtype = gtype0;
	mdeg = 0.0;
	pp->setToroidal(toroidal0);
	typeIncluded.clear();
	weightMatrix = wMatrix;

	if(weightMatrix[0]<0)getTypeToTypeWeightp = &Graph::getTypeToTypeWeight_all1; // if no weights given
	else{
		if(*dbg)Rprintf("Type-to-Type weight matrix received.\n");
		getTypeToTypeWeightp = &Graph::getTypeToTypeWeight_weighted; // weights given
	}

	for(i=0; i< pp->getNtypes(); i++)//need to ignore some types based on the inclusion vector
	{
		typein=0;
		for(j=0;j< pp->size(); j++)
			if( (pp->getT(&j) == pp->getTypevec(&i)) & inc[j])
			{
				typein=1;
				break;
			}
		typeIncluded.push_back(typein);
	}

	if(preDists[0]>=0)
	{
		if(*dbg)Rprintf("Setting precalculated distances...");
		pp->setDists(preDists);
		if(*dbg)Rprintf("ok. ");
	}
	else if(*doDists)  // precalculate distance triangle
	{
		if(*dbg)Rprintf("Precalculating distances...");
		pp->calcDists();
		if(*dbg)Rprintf("ok. ");
	}
	if(*dbg)Rprintf(" done.\n");

}

/********************************************************************************************/
void Graph::setNodelist(std::vector<std::vector<int> > *nodelist_new)
{
	int i;
	nodelist.clear();nodelist.resize(0);
	for(i=0;i<(int)nodelist_new->size();i++)
		nodelist.push_back(nodelist_new->at(i));
	preEdges = 1;
}
/********************************************************************************************/
void Graph::setNodelist(SEXP prepGraph)
{
	if(*dbg)Rprintf("setting precalculated edges...");
	nodelist.clear();nodelist.resize(0);
	VectsxpToVector(getListElement(prepGraph,"edges"), nodelist);
	preEdges = 1;
	if(*dbg)Rprintf("ok.");
}
/********************************************************************************************/
SEXP Graph::toSEXP()
//transform a std::vector<std::vector<int> > to SEXP, desctructive
{
	SEXP graph, *node;
	PROTECT(graph = allocVector(VECSXP, this->nodelist.size()));
	int i,j, *p, n;
	for(i=0;i< (int) this->nodelist.size();i++)
	{
		node = new SEXP;
		PROTECT(*node = allocVector(INTSXP, this->nodelist[i].size()) );
		p = INTEGER(*node);
		n = (int) this->nodelist[i].size();
		if(n<1) ;//p[0]=NULL;
		else
			for(j=0;j<n;j++)
			{
				p[j] = (int) this->nodelist[i][j];
			};
		this->nodelist[i].clear();
		SET_VECTOR_ELT(graph, i, *node);
		UNPROTECT(1);
		delete node;
	}
	UNPROTECT(1);
	return graph;
}
/********************************************************************************************/
void Graph::remove_duplicates()
// remove duplicates in edge lists, especially the simple delaunay method creates them
{
	int i,j,k,isnew;
	std::vector<int> *node;
	for(i=0; i < (int) this->nodelist.size() ;i++)
	{
		node = new std::vector<int>;
		node->resize(0);
		for(j=0; j < (int) this->nodelist.at(i).size() ; j++)
		{
			isnew=1;
			for(k=0;k<(int)node->size();k++)
			{
				if(node->at(k) == this->nodelist.at(i).at(j))
				{
					isnew=0;
					break;
				}
			}
			if(isnew)
				node->push_back(nodelist.at(i).at(j));
		}
		nodelist.at(i).swap(*node);
		delete node;
	}
}
/********************************************************************************************/
void Graph::addNew(int i, int j)
// add j to i
{
	int k,isnew=1;
	for(k=0; k< (int) nodelist.at(i).size();k++)
	{
		if(nodelist.at(i).at(k) == j)
		{
			isnew=0;
			break;
		}
	}
	if(isnew)
		nodelist.at(i).push_back(j);
}
/********************************************************************************************/
/********************************************************************************************/
// Weight matrix functions
double Graph::getTypeToTypeWeight(int *t1, int *t2){
	return (this->*getTypeToTypeWeightp)(t1,t2);
}
// if truly weighted
double Graph::getTypeToTypeWeight_weighted(int *t1, int *t2){
	return weightMatrix[(*t2-1)*pp->getNtypes()+(*t1-1)];
}
// otherwise
double Graph::getTypeToTypeWeight_all1(int *t1, int *t2){
	return 1;
}
/********************************************************************************************/
//The graph methods
/********************************************************************************************/
void Graph::sg_calc()
{
	// preprocess if requested
	if(prepR[0]>0 && *oldpar<= *par )
	{
		if(*dbg)Rprintf("Preprocessing[");
		this->sg_geometric(prepR);
		preEdges = 1;
		if(*dbg)Rprintf("] ok.\n ");
	}
	//start the calculation
	if(*gtype==0) //geometric
	{
		if(preEdges)
			this->sg_shrink_geometric(par);
		else
			this->sg_geometric();
	}
	else if(*gtype==1) //knn
	{
		if(preEdges)
			this->sg_shrink_knn();
		else
			this->sg_knn();
	}
	else if(*gtype==2) this->sg_mass_geometric();
	else if(*gtype==3) this->sg_gabriel();
	else if(*gtype==4) this->sg_delaunay();
	else if(*gtype==5) this->sg_MST();
	else if(*gtype==6) this->sg_markcross();
	else if(*gtype==7) this->sg_SIG();
	else if(*gtype==8) this->sg_RST();
	else if(*gtype==9) this->sg_RNG();
	else if(*gtype==10) this->sg_CCC();
	else if(*gtype==11) this->sg_STIR();
	else if(*gtype==12) this->sg_big_geometric();
	preEdges = 1;
}

/********************************************************************************************/
void Graph::sg_geometric()
{
 Graph::sg_geometric(par);
}

void Graph::sg_geometric(double *R)
{
	if(*dbg)Rprintf("Geometric (R=%f):",*R);
	int i,j;
	double dist;
	for(i=0;i<(pp->size()-1);i++)
		for(j=i+1;j<pp->size();j++)
		{
			dist = pp->getDist(&i, &j);
			if(dist<*R){
				nodelist[i].push_back(j+1);
				nodelist[j].push_back(i+1);
			}
		}
	mdeg = this->pp->lambda*PI*(*R)*(*R);
	if(*dbg)Rprintf(" Ok.");
}

void Graph::sg_big_geometric()
{
	if(*dbg)Rprintf("Big geometric (R=%f):",*par);
	int i,j;
	double dist;
	for(i=0;i<pp->size();i++)
		if(inc[i])
			for(j=0;j<pp->size();j++)
				if(i!=j)
				{
					dist = pp->getDist(&i, &j);
					if(dist<*par){
						nodelist[i].push_back(j+1);
					}
				}
	mdeg = this->pp->lambda*PI*(*par)*(*par);
	if(*dbg)Rprintf(" Ok.");
	*gtype = 0;
}

void Graph::sg_shrink_geometric(double *R)
{
	if(*dbg)Rprintf("Geometric (R=%f) (shrinking):",*R);
	int i,j,j0;
	double dist;
	std::vector<int> *node;
	for(i=0; i < pp->size() ; i++)
		if(inc[i])
		{
			node = new std::vector<int>;
			for(j=0;j < (int) this->nodelist[i].size() ; j++)
			{
				j0 = nodelist[i][j]-1;
				dist = pp->getDist(&i,&j0);
				if(dist<*R)
					node->push_back(j0+1);
			}
			nodelist[i].clear();nodelist[i].resize(0);
			for (j = 0; j < (int)node->size(); ++j) this->nodelist[i].push_back(node->at(j));
			delete node;
	}
	mdeg = this->pp->lambda*PI*(*R)*(*R);
	if(*dbg)Rprintf(" ok.");
}
/********************************************************************************************/
void Graph::sg_mass_geometric()
{
	if(*dbg)Rprintf("Mass-geometric:");
	int i,j;
	double dist;
	for(i=0;i<pp->size();i++)
		for(j=0;j<pp->size();j++)
		{
			if(i!=j)
			{
				dist = pp->getDist(&i, &j);
				if(dist< pp->getMass(&i)){
					nodelist[i].push_back(j+1);
	//				nodelist[j].push_back(i+1);
				}
			}
		}
	if(*dbg)Rprintf(" Ok.");
}
/********************************************************************************************/
void Graph::sg_knn()
{

	int i,j,l,*k,kk, mink;
	kk = (int) par[0];
	k = &kk;
	std::vector<int> *node;
	if(preEdges==0)// if not preprocessed
	{
		if(*dbg)Rprintf("%i-nn:",*k);
		double *dists2_i = new double[pp->size()], *dists2_i2 = new double[pp->size()];
		for(i=0;i<pp->size();i++) //for each point
		if(inc[i])
		{
			for(j=0;j<pp->size();j++) dists2_i2[j]=dists2_i[j]= pp->getDist(&i, &j); //gather the distances to others
			qsort( dists2_i, pp->size(), sizeof(double),compare_doubles); // sort distances, rising
			for(j=1;j<=*k;j++) // find the k nearest
				for(l=0;l<pp->size();l++)
					if( dists2_i[j] == dists2_i2[l] ) //with distance comparison
					{
						nodelist[i].push_back(l+1);
						break;
					}
		}
	}
	else{ //preprocessed
		if(*dbg)Rprintf("%i-nn (shrinking):",*k);
		double *dists2_i, *dists2_i2;
		for(i=0;i<pp->size();i++) //for each point
		if(inc[i])
		{
			node = new std::vector<int>;
			dists2_i = new double [nodelist[i].size()];
			dists2_i2 = new double [nodelist[i].size()];
			mink = *k;
			if( (int) nodelist[i].size()<*k )
			{
				mink = nodelist[i].size();
				Rprintf("\n preprocessing R too small, not enough neighbours (point #%i)!!\n",i+1);
			}

			for(l=0;l< (int) nodelist[i].size();l++)
			{
				j = nodelist[i][l]-1;
				dists2_i2[l]=pp->getDist(&i, &j); //gather the distances to others, given preprocessing
				dists2_i[l]=dists2_i2[l];
			}
			qsort( dists2_i, nodelist[i].size() , sizeof(double),compare_doubles); // sort distances, rising
			for(j=0; j<mink ;j++) // find the k nearest
				for(l=0;l< (int) nodelist[i].size();l++)
					if( dists2_i[j] == dists2_i2[l] ) //with distance comparison
					{
						node->push_back(nodelist[i][l]);
						break;
					}
			nodelist[i].clear();nodelist[i].resize(0);
			for(j=0;j < (int) node->size();j++) nodelist[i].push_back( node->at(j) );
			delete node;
			delete[] dists2_i;
			delete[] dists2_i2;
		}
	}
	mdeg = kk;
	if(*dbg)Rprintf(" Ok.");
}

void Graph::sg_shrink_knn()
{
	double *R0=prepR, R=1;
	prepR = &R;
	this->sg_knn();
	prepR = R0;
}

/********************************************************************************************/
void Graph::sg_gabriel()
{
	int kk = (int) par[0];
	if(*dbg & (kk>0) )Rprintf("%i-",kk);
	if(*dbg)Rprintf("Gabriel:");
	int i,j,k, empty,m,l,h;
	double x0,y0,R2, d;
	std::vector<int> *node;

	if(preEdges == 0) // no preprocessing done,a heavy looping
	  for(i=0;i<(pp->size()-1);i++)
	  {
		  for(j=i+1;j<pp->size();j++)
		  {
			  x0 = fabs(pp->getX(&i)-pp->getX(&j))/2.0+fmin2(pp->getX(&i),pp->getX(&j));
			  y0 = fabs(pp->getY(&i)-pp->getY(&j))/2.0+fmin2(pp->getY(&i),pp->getY(&j));
			  R2 = ( pow(pp->getX(&i)-pp->getX(&j),2) + pow(pp->getY(&i)-pp->getY(&j),2) )/4.0;
			  //		brute force
			  empty = 1+kk;
			  for(k=0;k<pp->size();k++)
			  {
				  if(k != i)
					  if( k != j)
					  {
						  d = pow(x0-pp->getX(&k),2) + pow(y0-pp->getY(&k),2);
						  if( d<R2 )
						  {
							  empty = empty - 1;
							  if(empty == 0) break;
						  }
					  }
			  }
			  if(empty)
			  {
				  this->nodelist[i].push_back(j+1);this->nodelist[j].push_back(i+1);
			  }
		  }
	  }
	else{ // preprocessed: nodelist has the restricted neighbourhoods to look trough
		if(*dbg)Rprintf("(prepd): ");
		for(i = 0 ; i< pp->size() ;  i++)
		{
			if(inc[i])
			{
				node = new std::vector<int>;
				for( l=0 ; l < (int) this->nodelist[i].size(); l++ )
				{
					j = this->nodelist[i][l]-1;
					x0 = fabs(this->pp->getX(&i)-this->pp->getX(&j))/2.0+fmin2(this->pp->getX(&i),this->pp->getX(&j));
					y0 = fabs(this->pp->getY(&i)-this->pp->getY(&j))/2.0+fmin2(this->pp->getY(&i),this->pp->getY(&j));
					R2 = (pow(this->pp->getX(&i)-this->pp->getX(&j),2) + pow(this->pp->getY(&i)-this->pp->getY(&j),2) )/4.0;
					empty = 1+kk;
					for(m=0; m <(int) this->nodelist[i].size();m++) // the small ball is included in the preprocessing ball
					{
						k =  (int) this->nodelist[i][m]-1;
						if(k != i)
							if( k != j)
							{
								d = pow(x0-pp->getX(&k),2) + pow(y0-pp->getY(&k),2);
								if( d<R2 )
								{
									empty = empty - 1;
									if(empty == 0) break;
								}
							}
					}
					if(empty)
					{
						node->push_back(j+1);
					}
				}
				nodelist[i].clear();nodelist[i].resize(0);
				for(h=0;h<(int)node->size();h++) nodelist[i].push_back(node->at(h));
				delete node;
			}
		}
	}
	mdeg = 4;
	if(*dbg)Rprintf(" Ok.");
}

/********************************************************************************************/
void Graph::sg_delaunay()
{
//Naive algorithm, checks the interiors of triangle circumcircles.
//For 2D patterns

	if(*dbg)Rprintf("Delaunay: ");
	int i,j,k,l,h;
	std::vector<int> *node;
	if(preEdges==0) // no preprocessing done, heavy looping
	{
		if(*dbg)Rprintf("(raw):");
		for(i = 0 ; i< pp->size()-2 ; i++ )
			for(j = i+1 ; j < pp->size()-1 ; j++ )
				for(k = j+1 ; k < pp->size() ; k++ )
					if( pp->Empty(&i,&j,&k) )
					{
						addNew(i,j+1); addNew(i,k+1);
						addNew(j,i+1); addNew(j,k+1);
						addNew(k,i+1); addNew(k,j+1);
					}
	}
	else{ // preprocessed: nodelist has the restricted neighbourhoods to look trough for triangles
		if(*dbg)Rprintf("(prepd): ");
		for(i = 0 ; i< pp->size() ;  i++)
		{
			if(inc[i])
			{
				node = new std::vector<int>;
				for( l=0 ; l < (int)nodelist[i].size()-1; l++ )
				{
					j = nodelist[i][l]-1;
					for(h=l+1; h < (int)nodelist[i].size();h++ )
					{
						k = nodelist[i][h]-1;
						if( pp->EmptyConstrained(&i, &j, &k, &nodelist.at(i)) )
						{
							node->push_back(j+1);
							node->push_back(k+1);
						}
					}
				}
				nodelist[i].clear();nodelist[i].resize(0);
				for(h=0;h<(int)node->size();h++) nodelist[i].push_back((*node).at(h));
				delete node;
			}
		}
		this->remove_duplicates();
	}
	mdeg = 6;
	if(*dbg)Rprintf(" Ok.");

}
/********************************************************************************************/
void Graph::sg_MST()
{
  if(*this->dbg) Rprintf("MST:");
  int i,j,k=0,l=0,zz,k0=0,l0=0;
  int *done = new int[pp->size()],dn;
  double apu0,apu1,apu2;
  done[0] = 0;
  dn = 1;
  int left=pp->size()-dn;
  while( left > 0 )
  {
    apu2 = MAX_DOUBLE;
    for(i=1; i<pp->size();i++){
      zz = 1;
      apu0=apu2;
      for(j=0; j<dn;j++){
        if(i == done[j] ) {
          zz=0;
          break;
        }
        apu1 = pp->getDist(&i,&done[j]);
        if( apu1<apu0 ){
          apu0=apu1;
          k0=i;
          l0=done[j];
        }
      }
      if(zz){
        if(apu0<apu2){
          apu2=apu0;
          k=k0;
          l=l0;
        }
      }
    }
    done[dn] = k;
    dn++;
    left--;
    this->nodelist[l].push_back(k+1);
  }
  if(*this->dbg)Rprintf(" Ok.");
}
/********************************************************************************************/
void Graph::sg_markcross()
{
	if(*dbg)Rprintf("Markcross: ");
	int i,j;
	double dist;
	for(i=0;i<(pp->size()-1);i++)
		for(j=i+1;j<pp->size();j++)
		{
			dist = pp->getDist(&i, &j);
			if(dist< (pp->getMass(&i)+pp->getMass(&j))){
				nodelist[i].push_back(j+1);
				nodelist[j].push_back(i+1);
			}
		}
	if(*dbg)Rprintf(" Ok.");
}
/********************************************************************************************/
void Graph::sg_SIG()
{
	if(*dbg)Rprintf("Spheres-of-Influence:");
	int i,j,dbg0=*dbg;
	double dist;
	for(i=0;i<pp->size();i++)
	{
		dist = MAX_DOUBLE;
		for(j=0;j<pp->size();j++)
			if(i!=j) dist = fmin2(dist, pp->getDist(&i, &j));
		pp->setMass(&i,&dist);
	}
	*dbg=0;
	sg_markcross();
	*dbg=dbg0;
	if(*dbg)Rprintf(" Ok.");
}
/********************************************************************************************/
void Graph::sg_RST()
{
  if(*dbg) Rprintf("Radial Spanning Tree (o=(%f,%f,%f)): ",par[0],par[1],par[2]);
  nodelist.resize(pp->size()-1);
  int i,j,k,foc_i=pp->size()-1;
  double apu0,apu1,apu2,apu3;

  for(i=0;i<pp->size()-1;i++)
  {
    apu0 = pp->getDist(&i,&foc_i);//dists[i*(*n+1)+*n];
    apu3=MAX_DOUBLE;
    k=-1;
    for(j=0;j<pp->size()-1;j++)
    {
      if(j!=i)
      {
        apu1 = pp->getDist(&j,&foc_i);//dists[j*(*n+1)+*n];
        if(apu1 < apu0 )
        {
          apu2 = pp->getDist(&i, &j);//dists[i*(*n+1)+j];
          if( apu2 < apu3 )
          {
            apu3 = apu2;
            k = j;
          }
        }
      }
    }
    if(k>-1) addNew(k,i+1);//e[k*(*n)+i] = 1;
  }
  if(*dbg) Rprintf(" Ok.");
}
/********************************************************************************************/
void Graph::sg_RNG()
{
	if(*dbg) Rprintf("Relative neighbourhood: ");
	int i,j,k,isempty;
    for(i=0;i<(pp->size()-1);i++)
    {
        for(j=i+1;j<pp->size();j++)
        {
        	isempty = 1;
        	for(k=0;k<pp->size();k++)
        		if( (k!=i) & (k!=j) )
        			if(pp->getDist(&i,&k) < pp->getDist(&i, &j))
        				if(pp->getDist(&j,&k) < pp->getDist(&j,&i))
        				{isempty=0;break;}
        	if(isempty)
        	{
        		addNew(i,j+1);
        		addNew(j,i+1);
        	}
        }
    }
    if(*dbg) Rprintf(" Ok.");
}
/********************************************************************************************/
void Graph::sg_CCC()
{
	double m=MAX_DOUBLE, mm = -MAX_DOUBLE, apu;
	if(*dbg) Rprintf("Class Cover Catch for type=%i: ",(int)par[0]);
	int i,j, type0=(int)par[0];
	for(i=0; i<pp->size();i++)
	{
		pp->setMass(&i,&mm);
		if(pp->getT(&i)==type0)
		{
			pp->setMass(&i, &m);
			for(j=0;j<pp->size();j++)
				if( (j!=i) & (pp->getT(&j)!=type0) )
				{
					apu = fmin2(pp->getMass(&i),pp->getDist(&i, &j));
					pp->setMass(&i, &apu);
				}
		}
	}
	for(i=0;i<pp->size();i++) //TODO: optimize this
		if(pp->getT(&i)==type0)
			for(j=0;j<pp->size();j++)
				if(i!=j)
					if(pp->getT(&j)==type0)
						if(pp->getDist(&i, &j)< pp->getMass(&i))
							addNew(i,j+1);
	if(*dbg) Rprintf(" Ok.");
}
/********************************************************************************************/
double Attenuate(double r, double alpha)
{
	return pow(1.0+r,-alpha);
}

void Graph::sg_STIR()
{
	if(*dbg) Rprintf("Signal-To-Noise-Ratio graph, noise=%f,alpha=%f,beta=%f,gamma=%f: ",par[0],par[1],par[2],par[3]);
	int i,j;
	double noise0 = par[0], alpha=par[1],beta=par[2],gamma=par[3], *noise = new double[pp->size()], sij,sji,s;
	//first we compute the interference for each location including the sending transmitter
	for(i=0;i<pp->size();i++)
	{
		noise[i]=0.0;
		for(j=0;j<pp->size();j++)
			if(i!=j) noise[i] = noise[i]+ pp->getMass(&j)*Attenuate(pp->getDist(&i, &j),alpha);
	}

	for(i=0;i<(pp->size()-1);i++)
		for(j=i+1;j<pp->size();j++)
		{
			sij = pp->getMass(&i)*Attenuate(pp->getDist(&i, &j),alpha);
			sji = pp->getMass(&j)*Attenuate(pp->getDist(&i, &j),alpha);
			sij = sij/(noise0 + gamma* (noise[j]-sij));
			sij = sji/(noise0 + gamma* (noise[i]-sji));
			s = fmin2(sij,sji);
			if(s >= beta )
			{
				addNew(i,j+1);
				addNew(j,i+1);
			}
		}
	if(*dbg) Rprintf(" Ok.");
}


/********************************************************************************************/
// cut (remove) the graph edges longer than R
void Graph::sg_cut(double *R)
{
	int i,j,k, count=0;
	if(*dbg)Rprintf("Cutting the graph (R=%f):",*R);
	std::vector<int > *pnode;
	for(i=0;i < pp->size();i++)
	{
		pnode = new std::vector<int>;
		pnode->resize(0);
		for(j=0; j < (int)nodelist.at(i).size();j++)
		{
			k = nodelist.at(i).at(j)-1;
			if( pp->getDist(&i, &k) < *R )
				pnode->push_back(k+1);
			else
				count++;
		}
		nodelist.at(i).swap(*pnode);
		delete pnode;
	}
	if(*dbg)Rprintf(" ok (%i edges cut). ",count);
}
/********************************************************************************************/
// prune branches less than lev hops long

void Graph::sg_prune(double *lev)
{
	int level = (int) *lev , i, leaf, count=0, prev, next;
	std::vector<int> left;
	std::vector<int> branch, *pnode;
	left.resize(0);
	branch.resize(0);
	if(*dbg)Rprintf("Pruning the graph (level=%i):",level);

	for(i=0; i < (int)nodelist.size(); i++ ) // get the leaves
	{
		if( (int)nodelist.at(i).size() == 1 )
			left.push_back(i+1);
	}
	if(*dbg)Rprintf("found %i leaves, pruning...",(int)left.size());

	while(!left.empty()) // go each branch trough starting from the leaf
	{
		leaf = left.back();
		branch.push_back(leaf);
		prev = leaf;
		next = nodelist.at(leaf-1).at(0);
		while((int) nodelist.at(next-1).size()==2)
		{
			branch.push_back(next);
			if(nodelist.at(next-1).at(0) != prev)
			{
				prev = next;
				next = nodelist.at(next-1).at(0);
			}
			else
			{
				prev = next;
				next = nodelist.at(next-1).at(1);
			}
//			for(i=0;i < branch.size();i++)
//				if(j == branch.at(i)){notnew=1; break;}
//			if(notnew)break;
		}
//		Rprintf("leaf:%i, branch length: %i\n",left.back(),branch.size());
		if((int)branch.size() <= level) // if short enough branch, cut it.
		{
			pnode = new std::vector<int>;
			pnode->resize(0);
			for(i=0; i < (int)branch.size();i++)
				nodelist.at(branch.at(i)-1).clear();

			for(i=0; i < (int)nodelist.at(next-1).size();i++)
				if(nodelist.at(next-1).at(i) != prev)
					pnode->push_back(nodelist.at(next-1).at(i));
			nodelist.at(next-1).swap(*pnode);
			delete pnode;
			count++;
		}
		left.pop_back();
		branch.clear();
		branch.resize(0);
	}
	if(*dbg)Rprintf(" Ok (%i branches pruned).",count);
}
/********************************************************************************************/


/**********************************************************************************/
int compare_doubles(const void *a, const void *b)
{
  const double *da = (const double *) a;
  const double *db = (const double *) b;

  return (*da > *db) - (*da < *db);
}

// EOF
