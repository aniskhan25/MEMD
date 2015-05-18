
#include <vector>
#include <bitset>
#include <cmath>

#include "MEMD_Config.h"

template<unsigned int T>
void peaks ( std::vector<double> X, std::bitset<T>& locs )
{
	double diff1, diff2;

	for ( unsigned int j = 0; j < N - 2; j++ )
	{
		diff1 = X[j+1] - X[j];
		diff2 = X[j+2] - X[j+1];
		if ( diff1 > 0 &  diff2 < 0 )
		locs[j] = true;
		else
		locs[j] = false;
	}
}

template<unsigned int T>
void peaks ( std::vector<double> X, const std::bitset<T>& a, std::bitset<T>& locs, int type )
{
	double diff1, diff2;

	double temp[3];
	
	unsigned int i = 0;	
	for ( unsigned int j = 0; j < N; j++ )
	{
		if ( a[j] )
		{
			temp[i] = (type) * X[j]; 
			++i;
		}
		if ( i > 2 )
		{	
			diff1 = temp[1] - temp[0];
			diff2 = temp[2] - temp[1];
			if ( (diff1 > 0) &  (diff2 < 0) )
			locs[j-1] = true;
			else
			locs[j-1] = false;
			
			--i;
			
			temp[i-2] = temp[i-1];
			temp[i-1] = temp[i];
		}
	}
}

template<unsigned int T>
void localPeaks ( std::vector<double> X, std::bitset<T>& loc_min, std::bitset<T>& loc_max )
{	
	std::bitset<N> a;

	double diff;
	for ( unsigned int j = 0; j < N - 1; j++ )
	{
		diff = X[j+1] - X[j];
		a[j] = ( diff != 0 ) ? true : false;
	}

	bool flag = false;
	unsigned int d = 1;				
	for ( unsigned int i = 0; i < N ; i++ )
	{	
		if ( a[i] ) flag = true;
		
		if ( !a[i] & flag) ++d;
		
		if ( (d > 1) & a[i] )
		{
			a[i] = false;
			a[i - (int) floor( (double) d/2 )] = true;
			
			d = 1;
			flag = false;
		}
	}
	a[N-1] = true;
	
	flag = false;
	for ( unsigned int i = 0; i<a.size(); ++i)
	{
		flag = ( a[i] ? true : false );
	}

	if ( flag )
	{
		// Maxima		
		peaks( X, a, loc_max, MAXIMA );
		// Minima
		peaks( X, a, loc_min, MINIMA );
	}
}