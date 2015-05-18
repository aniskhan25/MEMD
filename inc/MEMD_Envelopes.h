#ifndef MEMD_ENVELOPES_H
#define MEMD_ENVELOPES_H

#include <vector>

void envelopeMean(
	double* m
	, double* seq	
	, std::vector<double>& env_mean
	, std::vector<unsigned int>& ner
	, std::vector<unsigned int>& nzm
	, std::vector<double>& ampp
	);

#endif