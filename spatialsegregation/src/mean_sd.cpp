#include "mean_sd.h"


std::vector<double> mean(std::vector<double> v)
{
	std::vector<double> res;
	if(v.size()<1)
	{
		res.push_back(0.0);
		return res;
	}

	for(int i=1; i < (int)v.size();i++)
		v.at(0) = v.at(0)+v.at(i);
	res.push_back( v.at(0)/ (double)v.size() );
	return res;
}

std::vector<double>  sd(std::vector<double> v)
{
	std::vector<double> res;
	if(v.size()<2)
	{
		res.push_back(0.0);
		return res;
	}

	double m=mean(v).at(0),S=0.0;

	for(int i=0;i < (int)v.size();i++)
		S = S + pow(m-v.at(i),2);
	res.push_back(sqrt(S));
	return res;
}

