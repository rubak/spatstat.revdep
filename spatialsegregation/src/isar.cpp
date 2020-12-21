#include "isar.h"
/**********************************************************************************/
std::vector<double> isar(Graph *graph, double *fpar, int *dbg, int *included)
{
	if(fpar[1]==1) return isar_wdeg(graph, fpar, dbg, included);
	else if(fpar[1]==2) return isar_markweighted(graph, fpar, dbg, included);
	else if(fpar[1]==3) return isar_normal_empty(graph, fpar, dbg, included);
	return isar_normal(graph, fpar, dbg, included);
}
/*********************************************************************************/
std::vector<double> isar_normal(Graph *graph, double *fpar, int *dbg, int *included)
{
	if(*dbg)Rprintf("isar[type=%i",(int) fpar[0]);
	int i, j, k,l, dbg0,nt=0;
	double a1[1], m;
	std::vector<double> value;
	value.clear();

	if((int)fpar[0]<=0)
	{
		dbg0 = *dbg;
		*dbg = 0;
		for(i=0;i< graph->pp->getNtypes();i++)
		{
			if(graph->typeIncluded.at(i))
			{
				nt++;
				a1[0] = (double) graph->pp->getTypevec(&i);
				value.push_back(isar_normal(graph,a1,dbg,included).at(0));
			}
		}
//		Rprintf("nt=%i,",nt);
		*dbg = dbg0;
	}
	else
	{
		int target_type = (int) fpar[0];
		value.push_back(0.0);
		int n=0,typej;

		for(i=0;i< (int)graph->nodelist.size();i++) //sum over target species...
		{
			if(included[i] && graph->pp->getT(&i)==target_type) //...
			{
				n++; // one more of target_type
				m=0; // how many types around i?
				for(j=0; j< graph->pp->getNtypes();j++) // sum over all species
				{
					typej = graph->pp->getTypevec(&j);
					for(k=0;k < (int)graph->nodelist[i].size();k++)  // check if type j present
					{
						l = graph->nodelist[i][k]-1;
						if(graph->pp->getT(&l)==typej) // if i's neighbour k is of type j
						{
							m = m + graph->getTypeToTypeWeight(&target_type, &typej);
							break;
						}
					}
				}
				value.at(0)=value.at(0)+ m;
				graph->pp->setMass2(&i,&m);
			}
		}
		if(n>0) value.at(0) = value.at(0)/(double)n; // mean value over the species
	}

	if(*dbg)Rprintf("]");
	return value;
}

/*********************************************************************************/
std::vector<double> isar_normal_empty(Graph *graph, double *fpar, int *dbg, int *included)
{
	if(*dbg)Rprintf("isar(empty)[type=%i",(int) fpar[0]);
	int i, j, k,l, dbg0;
	double a1[1], m;
	std::vector<double> value;
	value.clear();

	if((int)fpar[0]<=0)
	{
		dbg0 = *dbg;
		*dbg = 0;
		for(i=0;i< graph->pp->getNtypes();i++)
		{
			if(graph->typeIncluded.at(i))
			{
				a1[0] = (double) graph->pp->getTypevec(&i);
				value.push_back(isar_normal_empty(graph,a1,dbg,included).at(0));
			}
		}
		*dbg = dbg0;
	}
	else
	{
		int target_type = (int) fpar[0];
		value.push_back(0.0);
		int n, typej, empty;
		for(j=0; j< graph->pp->getNtypes();j++) // sum P_tj over j
		{
			n=0;
			m=0.0;
			typej = graph->pp->getTypevec(&j);
			for(i=0;i< (int)graph->nodelist.size();i++) //sum over target species to estimate P_tj...
			{
				if(included[i] && graph->pp->getT(&i)==target_type) //...
				{
					n++; // one more of target_type
					empty = 1;
					for(k=0;k < (int)graph->nodelist[i].size();k++)  // check if type j present
					{
						l = graph->nodelist[i][k]-1;
						if(graph->pp->getT(&l)==typej) // if i's neighbour k is of type j
						{
							empty = 0;
							break;
						}
					}
					m = m + empty;
				}
			}
			if(n>0) m = m/(double)n;
			value.at(0)=value.at(0)+ (1.0-m);
			graph->pp->setMass2(&i,&m);
		}
	}
	if(*dbg)Rprintf("]");
	return value;
}

/**********************************************************************************/
std::vector<double> isar_wdeg(Graph *graph, double *fpar, int *dbg, int *included)
{
	if(*dbg)Rprintf("isar (degree weighted)[type=%i",(int) fpar[0]);
	int i, j, n, k,l, dbg0;
	double m;
	double a1[1],vi;
	std::vector<double> value;
	value.clear();

	if((int)fpar[0]<=0)
	{
		dbg0 = *dbg;
		*dbg = 0;
		for(i=0;i< graph->pp->getNtypes();i++)
		{
			if(graph->pp->lambdas[i]>0)
			{
				a1[0] = (double) graph->pp->getTypevec(&i);
				value.push_back(isar_wdeg(graph,a1,dbg,included).at(0));
			}
		}
		*dbg = dbg0;
	}
	else
	{
		int target_type = (int) fpar[0];
		int typej=0;
		value.push_back(0.0);
		n=0;
		for(i=0;i< (int)graph->nodelist.size();i++) //sum over target species...
		{
			if(included[i] && graph->pp->getT(&i)==target_type) //...
			{
				m=0; // how many types around i?
				n++;
				for(j=0; j< graph->pp->getNtypes();j++) // sum over all species
				{
					typej = graph->pp->getTypevec(&j);
					for(k=0;k < (int)graph->nodelist[i].size();k++)  // check if type j present
					{
						l = graph->nodelist[i][k]-1;
						if(graph->pp->getT(&l)==typej) // if i's neighbour k is of type j
						{
							m=m+graph->getTypeToTypeWeight(&target_type, &typej);
							break;
						}
					}
				}
				if(graph->nodelist[i].size()>0){
					vi = (double)m/graph->nodelist[i].size();
					value.at(0)=value.at(0)+ vi;
					graph->pp->setMass2(&i, &vi);
				}else n--;
			}
		}
		if(n>0) value.at(0) = value.at(0)/(double)n;
	}

	if(*dbg)Rprintf("]");
	return value;
}


/**********************************************************************************/
// not complete.
std::vector<double> isar_markweighted(Graph *graph, double *fpar, int *dbg, int *included)
{
	if(*dbg)Rprintf("isar (mass weighted)[type=%i",(int) fpar[0]);
	int i, j, n, k, l, dbg0;
	double a1[1], m;
	std::vector<double> value;
	value.clear();

	if((int)fpar[0]<=0)
	{
		dbg0 = *dbg;
		*dbg = 0;
		for(i=0;i< graph->pp->getNtypes();i++)
		{
			if(graph->pp->lambdas[i]>0)
			{
				a1[0] = (double) graph->pp->getTypevec(&i);
				value.push_back(isar_markweighted(graph,a1,dbg,included).at(0));
			}
		}
		*dbg = dbg0;
	}
	else
	{
		int target_type = (int) fpar[0];
		int isnew, typej;
		value.push_back(0.0);
		n=0;
		std::vector<int> p;
		std::vector<int> counts;
		std::vector<double> M;

		for(i=0;i< (int)graph->nodelist.size();i++) //sum over target species...
		{
			m=0;
			if(included[i] && graph->pp->getT(&i)==target_type) //...
			{
				p.clear();
				counts.clear();
				M.clear();
				n++;
				for(j=0; j < (int)graph->nodelist.at(i).size();j++)  // check if type j present
				{
					isnew=1;
					l = graph->nodelist[i][j]-1;
					typej = graph->pp->getT(&l);
					for(k=0; k < (int) p.size(); k++)
					{
						if(typej == p.at(k))
						{
							counts.at(k)++;
							M.at(k) = M.at(k) + graph->pp->getMass(&l)*graph->getTypeToTypeWeight(&target_type, &typej);
							isnew = 0;
							break;
						}
					}
					if(isnew)
					{
						counts.push_back(1);
						l = graph->nodelist.at(i).at(j)-1;
						M.push_back(graph->pp->getMass(&l));
						p.push_back(typej);
					}
				}
				for(k=0; k < (int)p.size();k++)
					m = m + (double) M.at(k)/(double)counts.at(k);
				value.at(0)=value.at(0)+(double)m;
			}
			graph->pp->setMass2(&i, &m);

		}
		value.at(0) = value.at(0)/(double)n;

	}

	if(*dbg)Rprintf("]");
	return value;
}
