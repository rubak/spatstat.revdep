#include "mingling.h"
/**********************************************************************************/

std::vector<double> mingling(Graph *graph, double *fpar, int *dbg, int *included)
{
	if(*dbg)Rprintf("mingling[");
	int target_type;
	int i,j,k,n=0,m, dbg0;
	double a1[2], neq, mdegi, mdeg, vtemp;
	std::vector<double> value;
	value.clear();

	if(*dbg)Rprintf("(type=%i,ratio=%i)",(int)fpar[0],(int)fpar[1]);
	if((int)fpar[0]==0)
	{
		dbg0 = *dbg;
		*dbg = 0;
		a1[1] = fpar[1];
		for(i=0;i< graph->pp->getNtypes();i++)
		{
			if(graph->typeIncluded.at(i))
			{
				a1[0] = (double) graph->pp->getTypevec(&i);
				value.push_back(mingling(graph, a1, dbg ,included).at(0));
			}

		}
		*dbg = dbg0;
	}
	else // target type given
	{
		target_type = (int) fpar[0];
		value.push_back(0.0);
		mdeg=0;
		for(i=0;i< (int)graph->nodelist.size() ;i++)
			if(included[i] &&  graph->pp->getT(&i) == target_type)
			{
				m = graph->nodelist[i].size();
				if(m>0)
				{
					n++;
					neq = 0;
					mdegi = 0;
					for(j=0;j<m;j++)
					{
						k = graph->nodelist[i][j]-1;
						if(target_type != graph->pp->getT(&k)) neq= neq+1.0/graph->pp->getWeight(&i,&k);
						mdegi = mdegi +  1.0/graph->pp->getWeight(&i,&k);
					}

					value.at(0) = value.at(0)+ (double) neq;
					mdeg = mdeg +  mdegi;
					vtemp = neq/(double)m;
					graph->pp->setMass2(&i, &vtemp);
				}
			}
		if(n>0){
			value.at(0) = value.at(0)/mdeg;
		}
		else value.at(0)=NA_REAL;

		if(fpar[1]>0) // ratio version (1-M)/ (lambda_t/lambda)
		{
			if(*dbg)Rprintf("M=%1.3f -> ",value.at(0));
			double ala;
			for(i=0;i < graph->pp->getNtypes(); i++) if(graph->pp->getTypevec(&i) == target_type) break;
			ala = graph->pp->lambdas[i] / graph->pp->lambda;
			value.at(0) = (1.0-value.at(0))/(double)ala;
			if(*dbg)Rprintf("%1.3f",value.at(0));
		}
	}
	if(*dbg)Rprintf("]");
	return value;
}

/**********************************************************************************/
