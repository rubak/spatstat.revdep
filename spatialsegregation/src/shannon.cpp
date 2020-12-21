#include "shannon.h"



std::vector<double> shannon0(Graph *graph, double *fpar, int *dbg, int *included)
// this is obsolete: we compute only piitau:s in C and do the scaling in R
{
	if(*dbg)Rprintf("shannon[");
	int i,S=graph->pp->getNtypes();
	double Eglobal=0.0, Elocal=0.0, ratio, lambda0=0, X;
	std::vector<double> value;


	for(i=0;i<S;i++)// compute overall intensity&format piitau
		lambda0= lambda0 + graph->pp->lambdas[i];

	for(i=0;i<S;i++)// compute the global entropy
	{
		ratio = graph->pp->lambdas[i]/lambda0;
		X = 0.0;
		if(ratio>0) X =  graph->pp->lambdas[i] * (log((double)ratio) / log((double)S) );
		Eglobal= Eglobal + X;
	}
	Eglobal = - Eglobal/lambda0;//=~E(N)

	value = piitauf(graph, fpar, dbg, included); // compute the pi_tau's, most demanding part

	for(i=0;i<S;i++)//~pi_tau
	{
		ratio = value.at(i);
		if(ratio>0)
			Elocal = Elocal + ratio* (log((double)ratio)/log((double)S));
	}
	//and finally the estimate
	value.clear();
	value.push_back(0.0);

	value.at(0) =  -Elocal /(double) Eglobal;

	if(*dbg)Rprintf("]");
	return value;
}

std::vector<double> shannon(Graph *graph, double *fpar, int *dbg, int *included)
{
	if(fpar[0]==0) return piitauf(graph, fpar, dbg, included);
	else return shannon_v2(graph, fpar, dbg, included);
}


std::vector<double> piitauf(Graph *graph, double *fpar, int *dbg, int *included)
{
	if(*dbg)Rprintf("piitau[");
	int i,j,k,l,m, S=graph->pp->getNtypes();
	double piitauj, *piitau = new double[S];
	std::vector<double> value;
	value.clear();

	for(i=0;i<S;i++)//compute the pii_tau's
	{
//		if(graph->pp->lambdas[i]>0) // not right, should not exclude any targets
		{
			m = 0;
			piitau[i]=0.0;

			for(j=0;j< graph->pp->size();j++) // first the piitau for each node...
			 if(included[j]) // exclude the minus-sampling sources
			 {
				m=m+1;
				piitauj = 0.0;
				if(graph->nodelist[j].size()>0)
				{
					for(k=0;k<(int)graph->nodelist[j].size(); k++)
					{
						l = graph->nodelist[j][k]-1;
						if(graph->pp->getT(&l)==graph->pp->getTypevec(&i))
							piitauj = piitauj + 1;
					}
					piitauj = piitauj / (double)graph->nodelist[j].size();
				}
				piitau[i] = piitau[i] + piitauj; // sum for the arith. mean...
			 }
			if(m>0) piitau[i] = piitau[i]/(double) m;//then the arith. mean over all sampled nodes
			value.push_back(piitau[i]);
		}
	}
	if(*dbg)Rprintf("]");
	return value;
}


// Second version of local entropy: Use the real number of local species as the log-base,
// then calculate the mean of entropy per type. Returns per type values
std::vector<double> shannon_v2(Graph *graph, double *fpar, int *dbg, int *included)
{
	if(*dbg)Rprintf("local entropies[");
	int i,j,k, S=graph->pp->getNtypes(), *Htau = new int[S];
	double v, pii;
	std::vector<int> counts;
	std::vector<double> locpitau;
	std::vector<double> value;


	for(k=0;k<S;k++){counts.push_back(0);value.push_back(0.0);} // zero the counts and values

//	main loop over points
	for(i=0;i<graph->pp->size();i++)//compute the local entropy for each included node
	{
		if(included[i]) // focus point included in calculation
		{
			for(k=0; k< S; k++) Htau[k]=0; // set local type counts 0
			locpitau.clear();

			counts[graph->pp->getT(&i)-1]++; // add one to global count vector

			Htau[graph->pp->getT(&i)-1]++;   //TODO: Should the point be aware of itself's type?

			// add the counts of neighbour types to local count vector
			for(j=0; j< (int)graph->nodelist[i].size(); j++)
			{
				k = graph->nodelist[i][j]-1; // type index of neighbour j
				Htau[ graph->pp->getT(&k) ] = 1 + Htau[ graph->pp->getT(&k) ]; // add one to local count vector
			}

			// count which types are present
			for(k=0; k< S; k++)
			{
				if(Htau[k]>0) locpitau.push_back(Htau[k]);
			}

			v = 0;
			// do the local entropy calculation
			if(locpitau.size()>1) // if we have multitype neighbourhood (others than focal point type), otherwise entropy =0
				for(k=0;k < (int)locpitau.size();k++)
				{
					pii = locpitau[k]/(double)graph->nodelist[i].size(); // the local relative frequency of type k
					v = v + pii * (log((double)pii)/log((double)locpitau.size())); // sum over type pi*log(pi)
				}

			value.at(graph->pp->getT(&i)-1) = value.at(graph->pp->getT(&i)-1) + v; // contributes to its own type's mean entropy
		}
	}

	// now calculate the mean entropy per type
	for(k=0;k<S;k++)
	{
		if(counts.at(k)>0)
			value.at(k) = value.at(k)/(double)counts.at(k);
	}

	if(*dbg)Rprintf("]");
	return value;

}


//EOF
