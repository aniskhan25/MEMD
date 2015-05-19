
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <cstdlib>

//#define LOGGING_LEVEL_1

#include "UTY_Time.h"
#include "UTY_Disp.h"
#include "UTY_Hammersley.h"

#include "MEMD_Config.h"
#include "MEMD_Criterion.h"
#include "MEMD_Envelopes.h"

int main()
{	
	unsigned int nbit, n_imf;
	double *x = new double[N*NDIM];
	double *r = new double[N*NDIM];
	double *m = new double[N*NDIM];
	double *seq = new double[NDIM*NDIR];
	std::vector<double> env_mean(NDIM*N);
   
	//std::ifstream data1("F:/TEST/xdata.csv");
	std::ifstream data1("F:/TEST/ma.csv");
	std::string line;
	int i = 0;
	while (getline(data1, line))
	{
		std::stringstream  lineStream(line);
		std::string  cell;
		while (getline(lineStream, cell, ','))
		{
			double ds = atof(cell.c_str());
			x[i] = ds;
			++i;
		}
	}
	
	int64 starttime = GetTimeMs64();

	hammersley( seq );
	
	//printVector( seq,  NDIR, NDIM );

	for ( unsigned int j = 0; j < N*NDIM; ++j )
	{
			r[j] = x[j];
	}
	
	bool stop_sift;
	nbit = 0;
	n_imf = 0;
	while ( ! stopEMD( r, seq ) )
	//if ( ! stopEMD( r, seq ) )
	{
		for ( unsigned int j = 0; j < N*NDIM; ++j )
		{
			m[j] = r[j];
		}
		
		//writeFileTabbed<double>( "F:/TEST/mdata.txt", m, N, NDIM );
		
		// computation of mean and stopping criterion		
		stop_sift = stopSift( m, seq, env_mean );
						
		//printVector<double>( m, N, NDIM );
		
		// In case the current mode is so small that machine precision can cause
		// spurious extrema to appear
		
		double mx1 = m[0], mx2 = x[0];
		for ( unsigned int j = 1; j < N*NDIM; ++j )
		{
			if ( fabs( m[j] ) > mx1 ) mx1 = m[j];
			if ( fabs( x[j] ) > mx2 ) mx2 = x[j];
		}
		
		//std::cout << "mx1 = " << mx1 << " : mx2 = " << mx2 << std::endl;
		
		if ( mx1 < ( 1e-10 ) * mx2 )
		{
			if ( !stop_sift )
			std::cout << "emd:warning, forced stop of EMD : too small amplitude" << std::endl;
			else
			std::cout << "forced stop of EMD : too small amplitude" << std::endl;
			//break;
		}
		
		nbit = 0;
		// sifting loop
		while ( ( !stop_sift ) & ( nbit < MAXITERATIONS ) )
		//if ( ( !stop_sift ) & ( nbit < MAXITERATIONS ) )
		{	

			//writeFileTabbed<double>( "F:/TEST/mdata.txt", m, N, NDIM );

			//sifting
			for ( unsigned int j = 0; j < N*NDIM; ++j )
			{
				m[j] -= env_mean[j];
			}
		
			//printVector<double>( env_mean, N, NDIM );
			//printVector<double>( m, N, NDIM );

			//writeFileTabbed<double>( "F:/TEST/env_mean.txt", env_mean, N, NDIM );
			//writeFileTabbed<double>( "F:/TEST/mdata.txt", m, N, NDIM );			
						
			// computation of mean and stopping criterion
			stop_sift = stopSift( m, seq, env_mean );
			
			++nbit;
			std::cout << nbit << std::endl;
			
			if ( ( nbit == (MAXITERATIONS - 1) ) &  ( nbit > 100 ) )
			std::cout << "emd:warning','forced stop of sifting : too many iterations" << std::endl;						
		}
		
		//q(:,n_imf,:) = m';
		//printVector( m, N, NDIM );		
	
		std::cout << nbit << std::endl;
		//printVector( m, N, NDIM );
	
		writeFileTabbed<double>( "F:/TEST/mdata.txt", m, N, NDIM );

		++n_imf;
						
		for ( unsigned int j = 0; j < N*NDIM; ++j )
		{
			r[j] -= m[j];
		}
		
		nbit = 0;
	}
	
	// Stores the residue
	//q(:, n_imf, :)=r';

	int64 endtime= GetTimeMs64();
	std::cout << "Elapsed time: " << (endtime - starttime)/1000 << " seconds.\n";

	//printVector<double>( r, N, NDIM );

	x = NULL; free(x);
	r = NULL; free(r);
	m = NULL; free(m);
	seq = NULL; free(seq);

	getchar();
	return 0;
}