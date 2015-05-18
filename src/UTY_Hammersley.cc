
#include <cmath>

#include "UTY_Hammersley.h"

#include "MEMD_Config.h"

#include "UTY_Common.h"

//
// Computes N elements of a Hammersley subsequence.
//
void hammersley( double *r )
{
	int j, seed2, digit;
	double base_inv, temp;

#define FIDDLE 0.5
#define STEP 1

	int base = -NDIR;
	unsigned int it=0;
	do {
		//
		//  Calculate the data.
		//
		if ( 1 < base )
		{			
			for ( j = 0; j < NDIR; j++ )
			{
				base_inv = 1.0 / ( ( double ) base );

				temp = 0.0;
				seed2 = ( j + STEP );
				while ( seed2 != 0 )
				{
					digit = seed2 % base;
					temp += ( ( double ) digit ) * base_inv;
					base_inv = base_inv / ( ( double ) base );
					seed2 /= base;					
				}				
				r[it + j*NDIM] = temp;
			}
		}
		//
		//  In the following computation, the value of FIDDLE can be:
		//
		//    0,   for the sequence 0/N, 1/N, ..., N-1/N
		//    1,   for the sequence 1/N, 2/N, ..., N/N
		//    1/2, for the sequence 1/(2N), 3/(2N), ..., (2*N-1)/(2N)
		//
		else
		{
			for ( j = 0; j < NDIR; j++ )
			{				
				temp = ( j + STEP ) % (-base + 1);
				temp = ( ( double ) ( temp ) + FIDDLE ) / ( double ) ( -base );
				r[it + j*NDIM] = temp;
			}
		}

		it++;
		base = prime( it );
	} while( it < NDIM );
}
