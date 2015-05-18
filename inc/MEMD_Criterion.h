#ifndef MEMD_CRITERION_H
#define MEMD_CRITERION_H

#include <vector>

bool stopEMD( 
	double *r
	, double *seq
	);

bool stopSift( 
	double *m
	, double *seq
	, std::vector<double> &env_mean 
	);

#endif