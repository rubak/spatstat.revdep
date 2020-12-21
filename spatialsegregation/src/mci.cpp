#include "mci.h"
/**********************************************************************************/

std::vector<double> mci(Graph *graph, double *fpar, int *dbg, int *included)
{
	if(*dbg)Rprintf("MCI[%i",(int)fpar[0]);
	int target_type;
	int i,j,k,l,n=0,in, dbg0;
	double a1, *pa, pt, mi;
	pa = &a1;
	std::vector<double> value;
	value.clear();

	if((int)fpar[0]==0) // target type not specified
	{
		dbg0 = *dbg;
		*dbg = 0;
		for(i=0;i< graph->pp->getNtypes();i++)
		{
			if(graph->typeIncluded.at(i))
			{
				a1 = (double) graph->pp->getTypevec(&i);
				value.push_back(mci(graph, pa, dbg ,included).at(0));
			}

		}
		*dbg = dbg0;
	}
	else // target type given
	{
		target_type = (int) fpar[0];
		value.push_back(0.0);
		for(i=0;i< (int)graph->nodelist.size(); i++)//go through points
		{
			if(included[i] &&  graph->pp->getT(&i) == target_type) // for each target type point...
			{
				n++;
				mi = 0.0;
				for(j=0; j< graph->pp->getNtypes(); j++) // sum over all species
				{
					in=0;
					for(k=0;k < (int)graph->nodelist.at(i).size(); k++)  // check if type j present
					{
						l = graph->nodelist.at(i).at(k)-1;
						if(graph->pp->getT(&l)==graph->pp->getTypevec(&j))
						{
							in = 1;
							break;
						}
					}
					pt = exp(-1.0*graph->pp->lambdas[j]*graph->mdeg/graph->pp->lambda);
					if(in) // species j present in neighborhood of i
						mi += log(1.0-pt);
					else
						mi += log(pt);
				}
//				Rprintf("\n");
				value.at(0) = value.at(0) - mi;
			}
		}
		if(n>0) value.at(0) = value.at(0)/(double)n; //I-value of target type
		for(j=0; j<graph->pp->getNtypes(); j++) // CSR substraction
		{
			pt = exp(-1.0*graph->mdeg*graph->pp->lambdas[j]/graph->pp->lambda);
			value.at(0) += pt*log(pt)+(1-pt)*log(1-pt);
		}
	}
	if(*dbg)Rprintf("]");
	return value;
}

/**********************************************************************************/
