#ifndef MEMD_PEAKS_H
#define MEMD_PEAKS_H

#include <bitset>
#include <vector>

#include "UTY_Types.h"

template<unsigned int T>
void peaks ( 
	std::vector<double> X
	, std::bitset<T>& locs 
	);

template<unsigned int T>
void peaks ( 
	std::vector<double> X
	, const std::bitset<T>& a
	, std::bitset<T>& locs
	, int type = MAXIMA 
	);

template<unsigned int T>
void localPeaks ( 
	std::vector<double> X
	, std::bitset<T>& loc_min
	, std::bitset<T>& loc_max 
	);

#include "MEMD_Peaks.tpp"

#endif