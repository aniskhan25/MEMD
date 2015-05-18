
#include <iostream>
#include <vector>
#include <bitset>
#include <algorithm>

#include "MEMD_Config.h"

#include "UTY_Types.h"
#include "UTY_Common.h"
#include "UTY_Disp.h"

template<unsigned int T>
bool boundaryConditions( 
std::vector<double> X
, const std::bitset<T> indmin
, const std::bitset<T> indmax	
, std::bitset<3*T>& tmin
, std::bitset<3*T>& tmax
)
{
	//printf("In\n");

	unsigned int 
	lenmin, lenmax
	, lsym, rsym;
	
	int 	
	firstmin, firstmax	
	, lastmin, lastmax
	, first_tlmin, first_tlmax
	, last_trmin, last_trmax;
	
	std::vector<int> lmax;
	std::vector<int> lmin;
	std::vector<int> rmax;
	std::vector<int> rmin;
	std::vector<int> tlmax;
	std::vector<int> tlmin;
	std::vector<int> trmax;
	std::vector<int> trmin;
	
	std::vector<int> temp;

	lenmin = getLength( indmin );
	lenmax = getLength( indmax );

	if ( (lenmin + lenmax) < 3 ) return false;

	firstmin = getIndex( indmin, FORWARD );	
	firstmax = getIndex( indmax, FORWARD );
	
	lastmin = getIndex( indmin, BACKWARD );
	lastmax = getIndex( indmax, BACKWARD );

	//std::cout << "min : (" << firstmin << "," << lastmin << ")" << std::endl;
	//std::cout << "max : (" << firstmax << "," << lastmax << ")" << std::endl;
	
	// boundary conditions for interpolations
	if ( firstmax < firstmin )
	{
		//std::cout << "in 1 ..." << std::endl;
		if ( X[0] > X[firstmin] )
		{
			//std::cout << "in 1.1 ..." << std::endl;
			lmax = getIndex( indmax, 1, NBSYM, FORWARD );
			lmin = getIndex( indmin, 0, (NBSYM - 1), FORWARD );
			lsym = firstmax;
		}
		else
		{
			//std::cout << "in 1.2 ..." << std::endl;
			
			lmax = getIndex( indmax, 0, (NBSYM - 1), FORWARD );
			lmin = getIndex( indmin, 0, (NBSYM - 2), FORWARD ); lmin.push_back(0);
			lsym = 0;
			
			//std::cout << "lmax : "; printVector( lmax );
			//std::cout << "lmin : "; printVector( lmin );
		}
	}
	else
	{   
		//std::cout << "in 2 ..." << std::endl;	
		if ( X[0] < X[firstmax] )
		{
			//std::cout << "in 2.1 ..." << std::endl;	
			lmax = getIndex( indmax, 0, (NBSYM - 1), FORWARD );
			lmin = getIndex( indmin, 1, NBSYM, FORWARD );
			lsym = firstmin;

			//std::cout << "lmax : "; printVector( lmax );
			//std::cout << "lmin : "; printVector( lmin );
		}
		else
		{
			//std::cout << "in 2.2 ..." << std::endl;	
			lmax = getIndex( indmax, 0, (NBSYM - 2), FORWARD ); lmax.push_back(0);
			lmin = getIndex( indmin, 0, (NBSYM - 1), FORWARD );
			lsym = 0;
		}
	}

	if ( lastmax < lastmin )
	{
		//std::cout << "in 3 ..." << endl;	
		if ( X[N-1] < X[lastmax] )
		{
			//std::cout << "in 3.1 ..." << std::endl;
			rmax = getIndex( indmax, 0, (NBSYM - 1), BACKWARD );
			rmin = getIndex( indmin, 1, NBSYM, BACKWARD );
			rsym = lastmin;
			
			//std::cout << "rmax : "; printVector( rmax );
			//std::cout << "rmin : "; printVector( rmin );
		}
		else
		{
			//std::cout << "in 3.2 ..." << std::endl;	
			//temp = getIndex( indmax, 0, (NBSYM - 2), BACKWARD ); rmax.push_back( N - 1 ); for(unsigned int i = 0; i < temp.size(); ++i) rmax.push_back(temp[i]);
			rmax = getIndex( indmax, 0, (NBSYM - 2), BACKWARD ); rmax.push_back( N - 1 );
			rmin = getIndex( indmin, 0, (NBSYM - 1), BACKWARD );			
			rsym = N - 1;
		}
	}
	else
	{
		//std::cout << "in 4 ..." << std::endl;	
		if ( X[N-1] > X[lastmin] )
		{
			//std::cout << "in 4.1 ..." << std::endl;	
			rmax = getIndex( indmax, 1, NBSYM, BACKWARD );
			rmin = getIndex( indmin, 0, (NBSYM - 1), BACKWARD );
			rsym = lastmax;
		}
		else
		{
			//std::cout << "in 4.2 ..." << std::endl;	
			rmax = getIndex( indmax, 0, (NBSYM - 1), BACKWARD );
			//temp = getIndex( indmin, 0, (NBSYM - 2), BACKWARD ); rmin.push_back( N - 1 ); for(unsigned int i = 0; i < temp.size(); ++i) rmin.push_back(temp[i]);
			rmin = getIndex( indmin, 0, (NBSYM - 2), BACKWARD ); rmin.push_back( N - 1 );
			rsym = N - 1;
		}
	}

	//std::cout << "lsym : " << lsym << std::endl;	
	//std::cout << "rsym : " << rsym << std::endl;	
	
	for ( unsigned int j = 0; j < lmin.size(); j++ ) {
		tlmin.push_back( 2*lsym - lmin[j] );		
	}
	for ( unsigned int j = 0; j < lmax.size(); j++ ) {
		tlmax.push_back( 2*lsym - lmax[j] );		
	}
	for ( unsigned int j = 0; j < rmin.size(); j++ ) {
		trmin.push_back( 2*rsym - rmin[j] );		
	}
	for ( unsigned int j = 0; j < rmax.size(); j++ ) {
		trmax.push_back( 2*rsym - rmax[j] );
	}

	/*
	for ( unsigned int j = 0; j < NBSYM; j++ )
	{
		//std::cout << "in 5.1 ..." << std::endl;
		tlmin.push_back( 2*lsym - lmin[j] );
		tlmax.push_back( 2*lsym - lmax[j] );
		trmin.push_back( 2*rsym - rmin[j] );
		trmax.push_back( 2*rsym - rmax[j] );
	}
	*/
	
	//std::cout << "tlmax : "; printVector( tlmax );
	//std::cout << "tlmin : "; printVector( tlmin );
	//std::cout << "trmax : "; printVector( trmax );
	//std::cout << "trmin : "; printVector( trmin );
	
	// in case symmetrized parts do not extend enough
		
	first_tlmin = *std::min_element( tlmin.begin(), tlmin.end() );	
	first_tlmax = *std::min_element( tlmax.begin(), tlmax.end() );
	
	//std::cout << "min = " << first_tlmin << std::endl;
	//std::cout << "max = " << first_tlmax << std::endl;
	
	if ( (first_tlmin > 0) | (first_tlmax > 0) ) // TODO. verify
	{
		//std::cout << "in 6 ..." << std::endl;	
		
		if ( lsym == firstmax )
		lmax = getIndex( indmax, 0, (NBSYM - 1), FORWARD );
		else
		lmin = getIndex( indmin, 0, (NBSYM - 1), FORWARD );
		
		if ( lsym == 0 ) std::cout << "error" << std::endl;
		
		lsym = 0;
		
		//std::cout << "in 6.1 ..." << std::endl;

		for ( unsigned int j = 0; j < lmax.size(); j++ ) {
			if ( j > tlmax.size()-1 )
				tlmax.push_back( 2*rsym - lmax[j] );
			else
				tlmax[j] = 2*lsym - lmax[j];			
		}

		for ( unsigned int j = 0; j < lmin.size(); j++ ) {
			if ( j > tlmin.size()-1 )
				tlmin.push_back( 2*rsym - lmin[j] );
			else
				tlmin[j] = 2*lsym - lmin[j];			
		}
		
		/*
		for ( unsigned int j = 0; j < NBSYM; j++ )
		{			
			tlmax[j] = 2*lsym - lmax[j];
			tlmin[j] = 2*lsym - lmin[j];
		}
		*/

		//std::cout << "tlmax : "; printVector( tlmax );
		//std::cout << "tlmin : "; printVector( tlmin );
	}
	
	last_trmin = *std::max_element(trmin.begin(), trmin.end() );
	last_trmax = *std::max_element(trmax.begin(), trmax.end() );
	
	//std::cout << "min = " << last_trmin << std::endl;
	//std::cout << "max = " << last_trmax << std::endl;
		
	if ( (last_trmin < (N - 1)) | (last_trmax < (N - 1)) ) // TODO. verify
	{
		//std::cout << "in 7 ..." << std::endl; // TODO. enters this condition for -O3 ???
		
		if ( rsym == lastmax )		
			rmax = getIndex( indmax, 0, (NBSYM -1), BACKWARD );			
		else		
			rmin = getIndex( indmin, 0, (NBSYM - 1), BACKWARD );
		
		if ( rsym == (N - 1) ) std::cout << "error" << std::endl;
		
		//std::cout << "rsym = " << rsym << std::endl;	
		
		rsym = N - 1;
		
		//std::cout << "rmax : "; printVector<int>( rmax );
		//std::cout << "rmin : "; printVector<int>( rmin );
		
		//std::cout << "in 7.1 ..." << std::endl;

		for ( unsigned int j = 0; j < rmax.size(); j++ ) {
			if ( j > trmax.size()-1 )
				trmax.push_back( 2*rsym - rmax[j] );
			else
				trmax[j] = 2*rsym - rmax[j];			
		}

		for ( unsigned int j = 0; j < rmin.size(); j++ ) {
			if ( j > trmin.size()-1 )
				trmin.push_back( 2*rsym - rmin[j] );
			else
				trmin[j] = 2*rsym - rmin[j];
		}

		/*
		for ( unsigned int j = 0; j < NBSYM; j++ )
		{	
			trmax[j] = 2*rsym - rmax[j];
			trmin[j] = 2*rsym - rmin[j];
		}
		*/

		//std::cout << "trmax : "; printVector<int>( trmax );
		//std::cout << "trmin : "; printVector<int>( trmin );
	}

	//std::cout << "trmax : "; printVector( trmax );
	//std::cout << "trmin : "; printVector( trmin );

	//std::cout << "tlmax : "; printVector( tlmax );
	//std::cout << "tlmin : "; printVector( tlmin );

	for ( unsigned int j = 0; j < 3*N; j++ )
	{
		tmin[j] = false;
		tmax[j] = false;
	}

	for ( unsigned int j = 0; j < N; j++ )
	{
		tmin[N + j] = indmin[j];
		tmax[N + j] = indmax[j];		
	}
	
	for ( unsigned int j = 0; j < tlmin.size(); j++ ) {
		tmin[N + tlmin[j]] = true;
	}

	for ( unsigned int j = 0; j < tlmax.size(); j++ ) {
		tmax[N + tlmax[j]] = true;
	}
	
	for ( unsigned int j = 0; j < trmin.size(); j++ ) {	
		tmin[N + trmin[j]] = true;
	}

	for ( unsigned int j = 0; j < trmax.size(); j++ ) {
		tmax[N + trmax[j]] = true;
	}

	/*
	for ( unsigned int j = 0; j < NBSYM; j++ )
	{
		//std::cout << j << " : " << ( N - 1 ) << " : " << tlmax[j] << std::endl;
		 
		tmin[( N  ) + tlmin[j]] = true;
		tmax[( N  ) + tlmax[j]] = true;
		
		tmin[( N  ) + trmin[j]] = true;
		tmax[( N  ) + trmax[j]] = true;
	}
	*/

	//std::cout << "tmin : "; printBitVector( tmin ); cout << "size = " << tmin.size() << std::endl;
	//std::cout << "tmax : "; printBitVector( tmax ); cout << "size = " << tmax.size() << std::endl;
	
	//printf("Out\n");

	return true; // the projected signal has inadequate extrema
}
