
#include <bitset>
#include <cmath>

#include "MEMD_Config.h"
#include "MEMD_Boundary.h"
#include "MEMD_Peaks.h"

#include "UTY_Common.h"
#include "UTY_Spline.h"
#include "UTY_Disp.h"

// computes the mean of the envelopes and the mode amplitude estimate

void envelopeMean(
	double* m
	, double* seq
	, std::vector<double>& env_mean
	, std::vector<unsigned int>& ner
	, std::vector<unsigned int>& nzm
	, std::vector<double>& ampp
	)
{			
	std::bitset<N> loc_min;
	std::bitset<N> loc_max;
	std::bitset<N> zerCross;
	
	std::bitset<3*N> tmin;
	std::bitset<3*N> tmax;
	
	std::vector<double> tht(NDIM - 1);
	std::vector<double> dir_vec(NDIM);
	
	std::vector<double> env_min (NDIM*N);
	std::vector<double> env_max (NDIM*N);
	
	std::vector<double> x(N);
	std::vector<double> y(N);
	std::vector<double> u(N);
	std::vector<double> v(N);

	std::vector<double> ypp(N);
	
	std::vector<double> a1(N);
	std::vector<double> a2(N);
	std::vector<double> a3(N);
	std::vector<double> a4(N);
	std::vector<double> a5(N);
	std::vector<double> b(N);
	
	double b1, b2, sum, prod;
	
	unsigned int count = 0;
	
	for ( unsigned int j = 0; j < N*NDIM; ++j )
	{
		env_mean[j] = 0.0;
	}
	
	for ( unsigned int j = 0; j < N; ++j )
	{
		ampp[j] = 0.0;
	}
	
	for ( unsigned int it = 0; it < NDIR; ++it )
	//for ( unsigned int it = 0; it < 1; ++it )
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
				sum += ( m[j + i*NDIM] * dir_vec[j] );
			}
			y[i] = sum;
		}
		
		//printVector( y );
				
		// Calculates the extrema of the projected signal

		localPeaks ( y, loc_min, loc_max );
				
		//printBitVector(loc_min);
		//printBitVector(loc_max);
		
		ner[it] = getLength(loc_min) + getLength(loc_max);

		zeroCrossing( y, zerCross );
		
		//printBitVector(zerCross);
		
		nzm[it] = getLength(zerCross);
		
		bool mode = boundaryConditions( y, loc_min, loc_max, tmin, tmax );
		
		//printBitVector(tmin);
		
		// Calculate multidimensional envelopes using spline interpolation
		// Only done if number of extrema of the projected signal exceed 3		
		if( mode )
		{			
			unsigned int N1;
			for ( unsigned int i = 0; i < NDIM; ++i )
			{
				// interpolate minima
				
				for ( unsigned int j = 0; j < N; ++j )
				{
					x[j] = 0.0;
					y[j] = 0.0;
					u[j] = j;
					v[j] = 0.0;
				}
				
				N1 = 0;
				for ( unsigned int j = 0; j < 3*N; ++j )
				{					
					if ( tmin[j] )
					{	
						int idx1;
						
						if ( j < N ) 
							idx1 = j - N;
						else if ( j >= (2*N) ) 
							idx1 = (j - 2*N) + N;
						else idx1 = j % N;
						
						x[N1] = idx1;

						unsigned int idx2;

						if ( j < N ) 
							idx2 = N - j;
						else if ( j >= (2*N) )							
							idx2 = N - (j%N) - 4;
						else idx2 = j % N;
						
						y[N1] = m[i + ( idx2 )*NDIM];
						++N1;
					}
				}

				//writeFileTabbed<double>( "F:/TEST/env_min.txt", m, N, NDIM );

				//cout << endl;
				
				//cout << N1 << " : ";
				
				//printBitVector( tmin );
				
				//printVector<double>( x, N1 );
				//printVector<double>( y, N1 );
				//printVector<double>( u );
				
				//chermiteCatmullRom( x, y, N1, u, v, N );				
				//chermiteMonotone( x, y, N1, u, v, N );
				splineCubicSpecial ( x, y, ypp, a1, a2, a3, a4, a5, b, N1, u, v, N );				
				//splineCubicTaha( x, y, ypp, a1, a2, a3, a4, a5, b, N1, u, v, N );
				
				//printVector<double>( b, N1 );
				//printVector<double>( c, N1 );
				//printVector<double>( b, N1 );
				//printVector<double>( v );
				
				for ( unsigned int j = 0; j < N; ++j )
				{
					env_min[i + j*NDIM] = v[j];
				}
				
				// interpolate maxima
				
				
				for ( unsigned int j = 0; j < N; ++j )
				{
					x[j] = 0.0;
					y[j] = 0.0;
					u[j] = j;
					v[j] = 0.0;
				}
				
				N1 = 0;
				for ( unsigned int j = 0; j < 3*N; ++j )
				{					
					if ( tmax[j] )
					{						
						int idx1;
						
						if ( j < N ) 
							idx1 = j - N;
						else if ( j >= (2*N) ) 
							idx1 = (j - 2*N) + N;
						else idx1 = j % N;
						
						x[N1] = idx1;

						unsigned int idx2;

						if ( j < N ) 
							idx2 = N - j;
						else if ( j >= (2*N) )							
							idx2 = N - (j%N) - 4;
						else idx2 = j % N;
						
						y[N1] = m[i + ( idx2 )*NDIM];
						++N1;
					}
				}

				//printVector<double>( x );
				//printVector<double>( y );
				//printVector<double>( u );
				
				//chermiteCatmullRom( x, y, N1, u, v, N );
				//chermiteMonotone( x, y, N1, u, v, N );
				splineCubicSpecial ( x, y, ypp, a1, a2, a3, a4, a5, b, N1, u, v, N );
				//splineCubicTaha( x, y, ypp, a1, a2, a3, a4, a5, b, N1, u, v, N );
								
				//printVector<double>( v );
				
				for ( unsigned int j = 0; j < N; ++j )
				{
					env_max[i + j*NDIM] = v[j];
				}
			}
			
			//readFile<double>( "F:/TEST/env_min.csv", env_min );
			//readFile<double>( "F:/TEST/env_max.csv", env_max );

			//writeFileTabbed<double>( "F:/TEST/env_min.txt", env_min, N, NDIM );
			//writeFileTabbed<double>( "F:/TEST/env_max.txt", env_max, N, NDIM );
			
			for ( unsigned int j = 0; j < N; ++j )
			{
				sum = 0.0;
				for ( unsigned int i = 0; i < NDIM; ++i )
				{
					sum += pow( ( env_max[i + j*NDIM] - env_min[i + j*NDIM] ), 2.0 );
					
					env_mean[i + j*NDIM] += ( ( env_max[i + j*NDIM] + env_min[i + j*NDIM]) / 2.0 );
				}
				ampp[j] += ( sqrt( sum ) / 2.0 );
			}
		
			//writeFileTabbed<double>( "F:/TEST/env_mean.txt", env_mean, N, NDIM );

			//printVector<double>(ampp);
			//printVector<double>(env_min,N,NDIM);
			//printVector<double>(env_max,N,NDIM);
		}
		else // if the projected signal has inadequate extrema
		{
			++count;
		}
	}
	
	//printVector<unsigned int>(ner);

	//printVector<double>(env_mean,N,NDIM);
	
	if ( count < NDIR )
	{
		for ( unsigned int j = 0; j < N; ++j )
		{			
			for ( unsigned int i = 0; i < NDIM; ++i )
			{
				env_mean[i + j*NDIM] /= ( NDIR - count );
			}
			ampp[j] /= ( NDIR - count );
		}
	}

	//writeFileTabbed<double>( "F:/TEST/env_mean.txt", env_mean, N, NDIM );

	//printVector<double>(env_mean,N,NDIM);
	//printVector(ampp);

	/*
	std::vector<double> env_mean_ref (NDIM*N);
	readFile<double>( "F:/TEST/env_mean.csv", env_mean_ref );			
	std::cout << ( ( compareVectors( env_mean, env_mean_ref ) ) ? "Passed" : "Failed" ) << std::endl;

	std::vector<double> ampp_ref (NDIM*N);
	readFile<double>( "F:/TEST/ampp.csv", ampp_ref );
	std::cout << ( ( compareVectors( ampp, ampp_ref ) ) ? "Passed" : "Failed" ) << std::endl;
	*/
}
