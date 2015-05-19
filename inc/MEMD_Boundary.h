#ifndef MEMD_BOUNDARY_H
#define MEMD_BOUNDARY_H

#include <vector>
#include <bitset>

template<unsigned int T>
bool boundaryConditions( 
std::vector<double> X
, const std::bitset<T> indmin
, const std::bitset<T> indmax	
, std::bitset<3*T>& tmin
, std::bitset<3*T>& tmax
);

#include "MEMD_Boundary.tpp"

#endif