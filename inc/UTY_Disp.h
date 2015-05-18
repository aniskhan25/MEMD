#ifndef UTY_DISP_H
#define UTY_DISP_H

#include <vector>
#include <bitset>

inline void printBitVector( std::vector<bool> bitvec );

template<unsigned int T>
inline void printBitVector( std::bitset<T> bitvec );

template<typename T>
inline void printVector( std::vector<T> vec );

template<typename T>
inline void printVector( std::vector<T> vec, unsigned int n );

template<typename T>
inline void printVector( T* vec, unsigned int n );

template<typename T>
inline void printVector( std::vector<T> vec, unsigned int rows, unsigned int cols );

template<typename T>
inline void printVector( T* vec, unsigned int rows, unsigned int cols );

template<typename T>
inline void writeFile( std::string filename, std::vector<T> vec, unsigned int rows, unsigned int cols );

template<typename T>
inline void writeFile( std::string filename, std::vector<T> vec, unsigned int rows, unsigned int cols );

template<typename T>
inline void writeFileTabbed( std::string filename, std::vector<T> vec );

template<typename T>
inline void writeFileTabbed( std::string filename, std::vector<T> vec, unsigned int n );

template<typename T>
inline void writeFileTabbed( std::string filename, std::vector<T> vec, unsigned int rows, unsigned int cols );

template<typename T>
inline void writeFileTabbed( std::string filename, T* vec, unsigned int rows, unsigned int cols );

template<typename T>
inline void readFile( std::string filename, std::vector<T>& vec );

template<typename T>
inline bool compareVectors( std::vector<T>& vec1, std::vector<T>& vec2 );

#include "UTY_Disp.tpp"

#endif