
#include <cmath>

#include "MEMD_Criterion.h"

#include "MEMD_Config.h"
#include "MEMD_Envelopes.h"
#include "MEMD_Peaks.h"

#include "UTY_Types.h"
#include "UTY_Common.h"
#include "UTY_Disp.h"

bool stopSift( 
	double *m
	, double *seq
	, std::vector<double> &env_mean
	)
{
	std::vector<unsigned int> nem(NDIR);
	std::vector<unsigned int> nzm(NDIR);
	std::vector<double> ampp(N);
	
	envelopeMean( m, seq, env_mean, nem, nzm, ampp);
	
	//printVector<double>(env_mean,N,NDIM);
	//printVector( nem );
		
	double sx;	
	unsigned int counter = 0;
	bool flag1 = false;
	for ( unsigned int i = 0; i < N; i++ )
	{
		sx = 0.0;
		for ( unsigned int j = 0; j < NDIM; j++ )
		{
			sx += pow( env_mean[j + i*NDIM], 2.0 );
		}
		sx = sqrt( sx );

		if( ampp[i] != 0.0 ) sx /= ampp[i];

		counter += ( ( sx > SD ) ? 1 : 0 );

		flag1 |= ( ( sx > SD2 ) ? true : false );
	}

	bool flag2 = false;
	for ( unsigned int i = 0; i < NDIR; i++ )
	{
		flag2 |= ( ( nem[i] > 2) ? true : false );
	}

	return !( ( ( (counter / N) > TOL ) | flag1 ) & flag2 ); // TODO. verify parenthesis
}

bool stopEMD( 
	double *r
	, double *seq 
	)
{
	std::vector<double> y(N);
	
	std::bitset<N> loc_min;
	std::bitset<N> loc_max;
	
	std::vector<double> tht(NDIM - 1);
	std::vector<double> dir_vec(NDIM);
		
	double b1, b2, sum, prod;
	
	unsigned int ner = 0;
	
	bool stp = true;
	
	for ( unsigned int it = 0; it < NDIR; ++it )
	{
		if ( NDIM != 3 ) // Multivariate signal (for NDIM ~=3) with hammersley sequence
		{
			// Linear normalisation of hammersley sequence in the range of -1.00 - 1.00
			
			sum = 0.0;
			b2 = 2*seq[(NDIM - 1) + it*NDIM] - 1;
			for ( int j = (NDIM - 2); j >=0; --j )
			{
				b1 = 2*seq[j + it*NDIM] - 1;			
				sum += pow( b2, 2.0 );
				tht[j] = atan2( sqrt( sum ), b1 );
				b2 = b1;
			}
			
			// Find angles corresponding to the normalised sequence
			
			prod = 1.0;			
			for ( unsigned int j = 0; j < (NDIM - 1); ++j )
			{
				dir_vec[j] = cos(tht[j]) * prod;
				prod *= sin( tht[j] );				
			}
			dir_vec[NDIM-1] = prod;
						
			//printVector(dir_vec);
		}
		else // Trivariate signal with Hammersley sequence
		{ // TODO. not tested
			// Linear normalisation of Hammersley sequence in the range of -1.0 - 1.0
			
			double tt = 2*seq[0 + it*NDIM] - 1;
			if ( tt < -1 ) tt = -1;
			if ( tt > +1 ) tt = +1;
									
			double phirad = seq[1 + it*NDIM] * (2 * PI);
			
			// Normalize angle from 0 - 2*pi
			
			double st = sqrt( 1.0 - pow( tt, 2.0 ) );
						
			dir_vec[0] = st * cos( phirad );
			dir_vec[1] = st * sin( phirad );
			dir_vec.push_back( tt );
		}
		
		// Projection of input signal on nth (out of total ndir) direction vectors
		
		for ( unsigned int i = 0; i < N; ++i )
		{
			sum = 0.0;
			for ( unsigned int j = 0; j < NDIM; ++j )
			{
				sum += ( r[j + i*NDIM] * dir_vec[j] );
			}
			y[i] = sum;
		}
		
		//printVector(y);
		
		// Calculates the extrema of the projected signal

		localPeaks ( y, loc_min, loc_max );
		
		//printBitVector(loc_min);
		//printBitVector(loc_max);
	
		// Stops if the all projected signals have less than 3 extrema	
		
		ner = getLength(loc_min) + getLength(loc_max);
		stp &= ( ner < 3 ) ? true : false;
	}
	
	//std::cout << "Stop = " << (stp ? '1' : '0');
	
	return stp;
}


