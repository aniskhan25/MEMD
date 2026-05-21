#ifndef UTY_COMMON_H
#define UTY_COMMON_H

#include <vector>
#include <bitset>
#include <cstddef>

#include "UTY_Types.h"

inline int prime ( int n );

inline double myRound( double x );

template<size_t T>
void zeroCrossing( std::vector<double> X, std::bitset<T>& zerCross );

inline bool getSign( double data );

template<size_t T>
int getLength( std::bitset<T> bitvec );

template<size_t T>
int getIndex( std::bitset<T> bitvec, int direction );

template<size_t T>
std::vector<int> getIndex( std::bitset<T> bitvec, unsigned int startidx, unsigned int endidx, int direction = FORWARD );

template<size_t T>
std::bitset<T> flipLR( std::bitset<T> bitvec, unsigned int startidx, unsigned int endidx );

#include "UTY_Common.tpp"

#endif