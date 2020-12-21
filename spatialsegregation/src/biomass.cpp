#include "biomass.h"
/**********************************************************************************/

std::vector<double> biomass(Graph *graph, double *fpar, int *dbg, int *included)
{
	if(*dbg)Rprintf("biomass[");
	int target_type;
	int i,j,k,n=0,m,dbg0;
	double a1[2], vtemp;
	std::vector<double> value;
	value.clear();
	if(*dbg)Rprintf("(type=%i, mean=%i)",(int)fpar[0], (int) fpar[1]);
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
				value.push_back(biomass(graph, a1, dbg ,included).at(0));
			}

		}
		*dbg = dbg0;
	}
	else // target type given
	{
		target_type = (int) fpar[0];
		value.push_back(0.0);
		// double ko = 99;
		double vi;
		for(i=0;i< (int)graph->nodelist.size() ;i++)
			if(included[i] &&  graph->pp->getT(&i) == target_type)
			{
			  //Rprintf("%f", graph->pp->getMass(&i));
				vi = 0.0;
				m = graph->nodelist[i].size();
				if(m>0)
				{
					n++;
					for(j=0;j<m;j++)
					{
						k = graph->nodelist[i][j]-1;
						vi = vi + graph->pp->getMass(&k);
					}

					if(fpar[1]>0) // if mean version
					{
						vi = vi/(double)m;
						//graph->pp->setMass2(&i, &vtemp);
					}
					//else graph->pp->setMass2(&i, &vi);
					value.at(0) = value.at(0) + vi;//graph->pp->getMass2(&i);
				}
				// graph->pp->setMass2(&i, &ko);
				//Rprintf(" - %f\n", graph->pp->getMass(&i));
			}
		if(n>0) value.at(0) = value.at(0)/(double)n;
	}

	if(*dbg)Rprintf("]");
	return value;
}

/**********************************************************************************/
