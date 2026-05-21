
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <iomanip>

//#define LOGGING_LEVEL_1

#include "UTY_Time.h"
#include "UTY_Disp.h"
#include "UTY_Hammersley.h"

#include "MEMD_Config.h"
#include "MEMD_Criterion.h"
#include "MEMD_Envelopes.h"

int main(int argc, char* argv[])
{	
	unsigned int nbit, n_imf;
	const std::string input_file = (argc > 1) ? argv[1] : "data/demo_signal.csv";
	const std::string output_prefix = (argc > 2) ? argv[2] : "output/memd";

	double *x = new double[N*NDIM]();
	double *r = new double[N*NDIM]();
	double *m = new double[N*NDIM]();
	double *seq = new double[NDIM*NDIR]();
	std::vector<double> env_mean(NDIM*N);

	std::ifstream data1(input_file.c_str());
	if (!data1) {
		std::cerr << "Unable to open input file: " << input_file << std::endl;
		return 1;
	}

	std::string line;
	unsigned int i = 0;
	while (getline(data1, line) && i < N*NDIM)
	{
		std::stringstream  lineStream(line);
		std::string  cell;
		while (getline(lineStream, cell, ',') && i < N*NDIM)
		{
			double ds = atof(cell.c_str());
			x[i] = ds;
			++i;
		}
	}
	if (i != N*NDIM) {
		std::cerr << "Expected " << (N*NDIM) << " values (" << N << " rows x " << NDIM
		          << " columns), but read " << i << " from " << input_file << std::endl;
		return 1;
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
			stop_sift = true;
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
	
		std::ostringstream imf_file;
		imf_file << output_prefix << "_imf" << std::setw(2) << std::setfill('0') << (n_imf + 1) << ".txt";
		writeFileTabbed<double>( imf_file.str(), m, N, NDIM );

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

	writeFileTabbed<double>( output_prefix + "_residue.txt", r, N, NDIM );
	std::cout << "Wrote " << n_imf << " IMF file(s) and residue with prefix " << output_prefix << std::endl;

	delete[] x;
	delete[] r;
	delete[] m;
	delete[] seq;

	return 0;
}