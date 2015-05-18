#ifndef UTY_COMMON_H
#define UTY_COMMON_H

#include <vector>
#include <bitset>

#include "UTY_Types.h"

inline int prime ( int n );

inline double myRound( double x );

template<unsigned int T>
void zeroCrossing( std::vector<double> X, std::bitset<T>& zerCross );

inline bool getSign( double data );

template<unsigned int T>
int getLength( std::bitset<T> bitvec );

template<unsigned int T>
int getIndex( std::bitset<T> bitvec, int direction );

template<unsigned int T>
std::vector<int> getIndex( std::bitset<T> bitvec, unsigned int startidx, unsigned int endidx, int direction = FORWARD );

template<unsigned int T>
std::bitset<T> flipLR( std::bitset<T> bitvec, unsigned int startidx, unsigned int endidx );

#include "UTY_Common.tpp"

#endif