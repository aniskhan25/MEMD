
#include <iostream>
#include <cmath>

#include "UTY_Spline.h"

#include "UTY_Common.h"

//void chermiteCatmullRom( 
//	std::vector<double> x
//	, std::vector<double> y
//	, unsigned int N1
//	, std::vector<double> xi
//	, std::vector<double>& yi
//	, unsigned int N2 
//	)
//{
//	double t, t2, t3, h, h2, h00, h10, h01, h11, a, b;
//
//	int k, klo, khi;
//
//	for( unsigned int i = 0 ; i < N2 ; i++ )
//	{
//		// Find the right place in the table by means of a bisection.
//		klo = 0;
//		khi = (N1-1);
//		while ( (khi - klo) > 1 )
//		{
//			k = (int) myRound( ( khi + klo) / 2.0 );
//			if ( x[k] > xi[i] )
//				khi = k;
//			else
//				klo = k;
//		}
//
//		h = x[khi] - x[klo];
//		if ( h == 0.0 )	std::cout << "Bad x input to chermite ==> x values must be distinct" << std::endl;      
//
//		// Catmull-Rom spline
//		if ( klo == 0 )
//		{
//			a = ( y[khi] - y[klo] ) / h;
//			b = ( y[khi+1] - y[klo] ) / ( x[khi+1] - x[klo] );
//		}
//		else if ( khi == (N1-1) )
//		{
//			a = ( y[khi] - y[klo-1] ) / ( x[khi] - x[klo-1] );
//			b = ( y[khi] - y[klo] ) / h;
//		}
//		else
//		{
//			a = ( y[khi] - y[klo-1]) / (x[khi] - x[klo-1] );
//			b = ( y[khi+1] - y[klo]) / (x[khi+1] - x[klo] );
//		}
//
//		// Evaluate cubic Hermite polynomial
//		t = ( xi[i] - x[klo] ) / h;
//		t2 = t*t;
//		t3 = t2*t;
//		h2 = h*h;
//		h00 = 2*t3 - 3*t2 + 1;
//		h10 = t3 - 2*t2 + t;
//		h01 = -2*t3 + 3*t2;
//		h11 = t3 - t2;
//
//		yi[i] = h00*y[klo] + h10*h*a + h01*y[khi] + h11*h*b;   
//	}
//}
//
//// Monotone interpolation
//
//void chermiteMonotone( std::vector<double> x, std::vector<double> y, unsigned int N1, std::vector<double> xi, std::vector<double>& yi, unsigned int N2 )
//{
//	double 
//		t, t2, t3, t4
//		, h, ih, h2
//		, h00, h10, h01, h11
//		, a, b, alpha, beta, tau
//		, xlo, ylo, xhi, yhi;
//
//	int k, klo, khi;
//
//	for( unsigned int i = 0 ; i < N2 ; i++ )
//	{
//		// Find the right place in the table by means of a bisection.
//		klo = 0;
//		khi = (N1 - 1);
//		while ( (khi - klo) > 1 )
//		{
//			//k = (int) myRound( ( khi + klo) * 0.5 );
//			k = ( khi + klo) / 2;
//			//k = static_cast<int>( ( ( khi + klo) * 0.5 ) + 0.5 );
//
//			if ( x[k] > i )
//				khi = k;
//			else
//				klo = k;
//		}
//
//		xlo = x[klo];
//		ylo = y[klo];
//		xhi = x[khi];
//		yhi = y[khi];
//
//		h = xhi- xlo;
//		if ( h == 0.0 )	std::cout << "Bad x input to chermite ==> x values must be distinct" << std::endl;
//
//		ih = 1 / h;
//
//		t4 = ( yhi - ylo ) * ih;
//
//		// Monotone interpolation
//		if ( klo == 0 )
//		{
//			a = t4;
//			b = ( t4 + ( y[khi+1] - yhi ) / ( x[khi+1] - xhi ) ) * 0.5;
//		}
//		else if ( khi == (N1 - 1) )
//		{
//			a = ( ( ylo - y[klo-1]) / ( xlo - x[klo-1] ) + t4 ) * 0.5;
//			b = t4;
//		}
//		else
//		{
//			a = ( ( ylo - y[klo-1] ) / ( xlo - x[klo-1] ) + t4 ) * 0.5;
//			b = ( t4 + ( y[khi+1] - yhi ) / ( x[khi+1] - xhi ) ) * 0.5;
//		}
//
//		if ( yhi == ylo )
//		{
//			a = 0;
//			b = 0;
//		}
//		else
//		{
//			alpha = a / t4;
//			beta  = b / t4;
//			if ( (alpha < 0) | (beta < 0) )
//			{
//				a = 0;
//				b = 0;
//			}
//			else
//			{
//				if ( (alpha*alpha + beta*beta) > 9 )
//				{
//					tau = 3/sqrt(alpha*alpha + beta*beta);
//					a = tau*alpha * t4;
//					b = tau*beta  * t4;
//				}
//
//			}
//		}
//
//		// Evaluate cubic Hermite polynomial
//		t = ( i - xlo ) * ih;
//		t2 = t*t;
//		t3 = t2*t;
//		h2 = h*h;
//		h00 = 2*t3 - 3*t2 + 1;
//		h10 = t3 - 2*t2 + t;
//		h01 = -2*t3 + 3*t2;
//		h11 = t3 - t2;
//
//		yi[i] = h00*ylo + h10*h*a + h01*yhi + h11*h*b;
//	}
//}
//
///*
//*	Natural Splines
//*	---------------
//*	Here the end-conditions are determined by setting the second
//*	derivative of the spline at the end-points to equal to zero.
//*
//*	There are n-2 unknowns (y[i]'' at x[2], ..., x[n-1]) and n-2
//*	equations to determine them.  Either Choleski or Gaussian
//*	elimination could be used.
//*/
//
//void naturalSpline( 
//	std::vector<double> x
//	, std::vector<double> y
//	, unsigned int n
//	, std::vector<double>& b
//	, std::vector<double>& c
//	, std::vector<double>& d
//	)
//{
//	int nm1;
//	double t;
//
//	//x--; y--; b--; c--; d--;
//
//	if(n < 2) return;
//
//	if(n < 3) 
//	{
//		t = ( y[2] - y[1] );
//		b[1] = t / ( x[2] - x[1] );
//		b[2] = b[1];
//		c[1] = c[2] = d[1] = d[2] = 0.0;
//		return;
//	}
//
//	nm1 = n - 1;
//
//	/* Set up the tridiagonal system */
//	/* b = diagonal, d = offdiagonal, c = right hand side */
//
//	d[1] = x[2] - x[1];
//	c[2] = ( y[2] - y[1] ) / d[1];
//	for( unsigned int i = 2 ; i < n ; i++ )
//	{
//		d[i] = x[i+1] - x[i];
//		b[i] = 2.0 * ( d[i-1] + d[i] );
//		c[i+1] = ( y[i+1] - y[i] ) / d[i];
//		c[i] = c[i+1] - c[i];
//	}
//
//	/* Gaussian elimination */
//
//	for( unsigned int i = 3 ; i < n ; i++ )
//	{
//		t = d[i-1] / b[i-1];
//		b[i] = b[i] - t*d[i-1];
//		c[i] = c[i] - t*c[i-1];
//	}
//
//	/* Backward substitution */
//
//	c[nm1] = c[nm1] / b[nm1];
//	for( unsigned int i = n-2 ; i > 1 ; i-- )
//		c[i] = ( c[i] - d[i]*c[i+1] ) / b[i];
//
//	/* End conditions */
//
//	c[1] = c[n] = 0.0;
//
//	/* Get cubic coefficients */
//
//	b[1] = (y[2] - y[1]) / d[1] - d[1] * c[2];
//	c[1] = 0.0;
//	d[1] = c[2] / d[1];
//	b[n] = ( y[n] - y[nm1] ) / d[nm1] + d[nm1] * c[nm1];
//	for( unsigned int i = 2 ; i < n ; i++ )
//	{
//		b[i] = ( y[i+1] - y[i] ) / d[i] - d[i] * ( c[i+1] + 2.0*c[i] );
//		d[i] = ( c[i+1] - c[i] ) / d[i];
//		c[i] = 3.0 * c[i];
//	}
//	c[n] = 0.0;
//	d[n] = 0.0;
//
//	return;
//}
//
//void splineEval( 
//	std::vector<double> x
//	, std::vector<double> y
//	, std::vector<double> b
//	, std::vector<double> c
//	, std::vector<double> d
//	, unsigned int n1
//	, std::vector<double> u
//	, std::vector<double>& v
//	, unsigned int n2
//	)
//{
//	/* 
//	* Evaluate  v[l] := spline(u[l], ...),	    l = 1,..,nu, i.e. 0:(nu-1)
//	* Nodes x[i], coef (y[i]; b[i],c[i],d[i]); i = 1,..,n , i.e. 0:(*n-1)
//	*/	
//	unsigned int i, j, k, l;
//	double ul, dx, tmp;
//
//	for( l = 0, i = 0; l < n2; l++ )
//	{
//		ul = u[l];
//		if( (ul < x[i]) | ( (i < (n1-1)) & (x[i+1] < ul) ) )
//		{
//			/* reset i  such that  x[i] <= ul <= x[i+1] : */
//			i = 0;
//			j = n1;
//			do 
//			{
//				k = ( i + j ) / 2;
//				if ( ul < x[k] ) j = k;
//				else i = k;
//			} 
//			while( j > ( i + 1 ) );
//		}
//		dx = ul - x[i];
//		/* for natural splines extrapolate linearly left */
//		tmp = ( ul < x[0] ) ? 0.0 : d[i];
//
//		v[l] = y[i] + dx * ( b[i] + dx * ( c[i] + dx*tmp ) );
//	}
//}
//
//void splineCubic ( 
//	std::vector<double> x
//	, std::vector<double> y
//	, std::vector<double> a1
//	, std::vector<double> a2
//	, std::vector<double> a3
//	, std::vector<double> a4
//	, std::vector<double> a5
//	, std::vector<double> b
//	, std::vector<double> ypp
//	, unsigned int N1
//	, std::vector<double> xi
//	, std::vector<double>& yi
//	, unsigned int N2
//	)
//{	
//	int i;
//	int j;
//	int jhi;
//	int k;
//	double yval;
//
//	splineCubicSet ( x, y, ypp, a1, a2, a3, a4, a5, b, N1 );
//
//	for ( i = 0; i < N2; i++ )
//	{
//		yi[i] = splineCubicVal ( x, y, ypp, xi[i], N1 );
//	}
//}
//
//void splineCubicSet ( 
//	std::vector<double> t
//	, std::vector<double> y
//	, std::vector<double>& ypp
//	, std::vector<double> a1
//	, std::vector<double> a2
//	, std::vector<double> a3
//	, std::vector<double> a4
//	, std::vector<double> a5
//	, std::vector<double> b
//	, int n
//	)
//{
//
//	int i;	
//	/*
//	Check.
//	*/
//	if ( n <= 1 )
//	{
//		fprintf ( stderr, "\n" );
//		fprintf ( stderr, "SPLINE_CUBIC_SET - Fatal error!\n" );
//		fprintf ( stderr, "  The number of data points N must be at least 2.\n" );
//		fprintf ( stderr, "  The input value is %d.\n", n );
//		exit ( 1 );
//	}
//
//	for ( i = 0; i < n - 1; i++ )
//	{
//		if ( t[i+1] <= t[i] )
//		{
//			fprintf ( stderr, "\n" );
//			fprintf ( stderr, "SPLINE_CUBIC_SET - Fatal error!\n" );
//			fprintf ( stderr, "  The knots must be strictly increasing, but\n" );
//			fprintf ( stderr, "  T(%d) = %g\n", i, t[i] );
//			fprintf ( stderr, "  T(%d) = %g\n", i+1, t[i+1] );
//			exit ( 1 );
//		}
//	}
//
//	for ( i = 0; i < n; i++ )
//	{
//		a1[i] = 0.0;
//		a2[i] = 0.0;
//		a3[i] = 0.0;
//		a4[i] = 0.0;
//		a5[i] = 0.0;
//	}
//	/*
//	Set up the first equation.
//	*/
//
//	b[0] = 0.0;
//	a3[0] = - ( t[2] - t[1] );
//	a4[0] =   ( t[2]        - t[0] );
//	a5[0] = - (        t[1] - t[0] );
//
//	/*
//	Set up the intermediate equations.
//	*/
//	for ( i = 1; i < n - 1; i++ )
//	{
//		b[i] = ( y[i+1] - y[i] ) / ( t[i+1] - t[i] )
//			- ( y[i] - y[i-1] ) / ( t[i] - t[i-1] );
//		a2[i] = ( t[i+1] - t[i]   ) / 6.0;
//		a3[i] = ( t[i+1] - t[i-1] ) / 3.0;
//		a4[i] = ( t[i]   - t[i-1] ) / 6.0;
//	}
//	/*
//	Set up the last equation.
//	*/
//	b[n-1] = 0.0;
//	a1[n-1] = - ( t[n-1] - t[n-2] );
//	a2[n-1] =   ( t[n-1]          - t[n-3] );
//	a3[n-1] = - (          t[n-2] - t[n-3] );
//
//	/*
//	Solve the linear system.
//	*/	
//	penta ( a1, a2, a3, a4, a5, b, ypp, n );
//}
//
//void penta (
//	std::vector<double> a1
//	, std::vector<double> a2
//	, std::vector<double> a3
//	, std::vector<double> a4
//	, std::vector<double> a5
//	, std::vector<double> b
//	, std::vector<double>& x
//	, int n
//	)
//{
//	int i;	
//	double xmult;
//
//	for ( i = 1; i < n - 1; i++ )
//	{
//		xmult = a2[i] / a3[i-1];
//		a3[i] = a3[i] - xmult * a4[i-1];
//		a4[i] = a4[i] - xmult * a5[i-1];
//		b[i] = b[i] - xmult * b[i-1];
//		xmult = a1[i+1] / a3[i-1];
//		a2[i+1] = a2[i+1] - xmult * a4[i-1];
//		a3[i+1] = a3[i+1] - xmult * a5[i-1];
//		b[i+1] = b[i+1] - xmult * b[i-1];
//	}
//
//	xmult = a2[n-1] / a3[n-2];
//	a3[n-1] = a3[n-1] - xmult * a4[n-2];
//	x[n-1] = ( b[n-1] - xmult * b[n-2] ) / a3[n-1];
//	x[n-2] = ( b[n-2] - a4[n-2] * x[n-1] ) / a3[n-2];
//	for ( i = n - 3; 0 <= i; i-- )
//	{
//		x[i] = ( b[i] - a4[i] * x[i+1] - a5[i] * x[i+2] ) / a3[i];
//	}
//}
//
//double splineCubicVal (
//	std::vector<double> t
//	, std::vector<double> y
//	, std::vector<double> ypp
//	, double tval
//	, int n
//	)
//{
//	double dt;
//	double h;
//	int i;
//	int ival;
//	double yval;
//	/*
//	Determine the interval [ T(I), T(I+1) ] that contains TVAL.
//	Values below T[0] or above T[N-1] use extrapolation.
//	*/
//	ival = n - 2;
//
//	for ( i = 0; i < n-1; i++ )
//	{
//		if ( tval < t[i+1] )
//		{
//			ival = i;
//			break;
//		}
//	}
//	/*
//	In the interval I, the polynomial is in terms of a normalized
//	coordinate between 0 and 1.
//	*/
//	dt = tval - t[ival];
//	h = t[ival+1] - t[ival];
//
//	yval = y[ival]
//	+ dt * ( ( y[ival+1] - y[ival] ) / h
//		- ( ypp[ival+1] / 6.0 + ypp[ival] / 3.0 ) * h
//		+ dt * ( 0.5 * ypp[ival]
//	+ dt * ( ( ypp[ival+1] - ypp[ival] ) / ( 6.0 * h ) ) ) );
//
//	return yval;
//}
//
//
///************************************************************/
//
//
//void splineCubicSpecial ( 
//std::vector<double> x
//, std::vector<double> y
//, std::vector<double> ypp
//, std::vector<double> a1
//, std::vector<double> a2
//, std::vector<double> a3
//, std::vector<double> a4
//, std::vector<double> a5
//, std::vector<double> b
//, unsigned int N1
//, std::vector<double> xi
//, std::vector<double>& yi
//, unsigned int N2
//)
//{	
//	int i;
//	int j;
//	int jhi;
//	int k;
//	double yval;
//	
//	/******************************************************************/
//	
//	//
//	//  Set up the data.
//	//
//
//	for ( i = 0; i < N1; i++ )
//	{
//		a1[i] = 0.0;
//		a2[i] = 0.0;
//		a3[i] = 0.0;
//		a4[i] = 0.0;
//		a5[i] = 0.0;
//	}
//	/*
//Set up the first equation.
//*/
//
//	b[0] = 0.0;
//	a3[0] = - ( x[2] - x[1] );
//	a4[0] =   ( x[2]        - x[0] );
//	a5[0] = - (        x[1] - x[0] );
//
//	/*
//Set up the last equation.
//*/
//	b[N1-1] = 0.0;
//	a1[N1-1] = - ( x[N1-1] - x[N1-2] );
//	a2[N1-1] =   ( x[N1-1]          - x[N1-3] );
//	a3[N1-1] = - (          x[N1-2] - x[N1-3] );
//	
//	/*
//Set up the intermediate equations.
//*/
//	for ( i = 1; i < N1 - 1; i++ )
//	{
//		b[i] = ( y[i+1] - y[i] ) / ( x[i+1] - x[i] )
//		- ( y[i] - y[i-1] ) / ( x[i] - x[i-1] );
//		a2[i] = ( x[i+1] - x[i]   ) / 6.0;
//		a3[i] = ( x[i+1] - x[i-1] ) / 3.0;
//		a4[i] = ( x[i]   - x[i-1] ) / 6.0;
//	}
//	
//
//	/*
//Solve the linear system.
//*/	
//	//penta ( a1, a2, a3, a4, a5, b, ypp, n );
//	
//	double xmult;
//
//	for ( i = 1; i < N1 - 1; i++ )
//	{
//		xmult = a2[i] / a3[i-1];
//		a3[i] = a3[i] - xmult * a4[i-1];
//		a4[i] = a4[i] - xmult * a5[i-1];
//		b[i] = b[i] - xmult * b[i-1];
//		
//		xmult = a1[i+1] / a3[i-1];
//		a2[i+1] = a2[i+1] - xmult * a4[i-1];
//		a3[i+1] = a3[i+1] - xmult * a5[i-1];
//		b[i+1] = b[i+1] - xmult * b[i-1];
//	}
//
//	xmult = a2[N1-1] / a3[N1-2];
//	a3[N1-1] = a3[N1-1] - xmult * a4[N1-2];
//	ypp[N1-1] = ( b[N1-1] - xmult * b[N1-2] ) / a3[N1-1];
//	ypp[N1-2] = ( b[N1-2] - a4[N1-2] * ypp[N1-1] ) / a3[N1-2];
//	
//	for ( i = N1 - 3; 0 <= i; i-- )
//	{
//		ypp[i] = ( b[i] - a4[i] * ypp[i+1] - a5[i] * ypp[i+2] ) / a3[i];
//	}
//	
//	/******************************************************************/
//	
//	double dt;
//	double h;	
//	int ival;	
//	
//	for ( i = 0; i < N2; i++ )
//	{
//		/*
//Determine the interval [ T(I), T(I+1) ] that contains TVAL.
//Values below T[0] or above T[N-1] use extrapolation.
//*/
//		ival = N1 - 2;
//
//		for ( j = 0; j < N1-1; j++ )
//		{
//			if ( xi[i] < x[j+1] )
//			{
//				ival = j;
//				break;
//			}
//		}
//		/*
//In the interval I, the polynomial is in terms of a normalized
//coordinate between 0 and 1.
//*/
//		dt = xi[i] - x[ival];
//		h = x[ival+1] - x[ival];
//
//		yi[i] = y[ival]
//		+ dt * ( ( y[ival+1] - y[ival] ) / h
//		- ( ypp[ival+1] / 6.0 + ypp[ival] / 3.0 ) * h
//		+ dt * ( 0.5 * ypp[ival]
//		+ dt * ( ( ypp[ival+1] - ypp[ival] ) / ( 6.0 * h ) ) ) );
//	}
//}
//
//void splineCubicSpecial ( 
//	std::vector<double> x
//	, std::vector<double> y
//	, std::vector<double> ypp
//	, std::vector<double> a1
//	, std::vector<double> a2
//	, std::vector<double> a3
//	, std::vector<double> a4
//	, std::vector<double> a5
//	, std::vector<double> b
//	, unsigned int N1
//	, std::vector<double> xi
//	, std::vector<double>& yi
//	, unsigned int N2
//	)
//{	
//	int i;
//	int j;
//	int jhi;
//	int k;
//	double yval;
//
//	/******************************************************************/
//
//	//
//	//  Set up the data.
//	//
//
//	for ( i = 0; i < N1; i++ )
//	{
//		a1[i] = 0.0;
//		a2[i] = 0.0;
//		a3[i] = 0.0;
//		a4[i] = 0.0;
//		a5[i] = 0.0;
//	}
//	/*
//	Set up the first equation.
//	*/
//
//	b[0] = 0.0;
//	a3[0] = - ( x[2] - x[1] );
//	a4[0] =   ( x[2]        - x[0] );
//	a5[0] = - (        x[1] - x[0] );
//
//	/*
//	Set up the last equation.
//	*/
//	b[N1-1] = 0.0;
//	a1[N1-1] = - ( x[N1-1] - x[N1-2] );
//	a2[N1-1] =   ( x[N1-1]          - x[N1-3] );
//	a3[N1-1] = - (          x[N1-2] - x[N1-3] );
//
//	/*
//	Set up the intermediate equations.
//	*/
//	for ( i = 1; i < N1 - 1; i++ )
//	{
//		b[i] = ( y[i+1] - y[i] ) / ( x[i+1] - x[i] )
//			- ( y[i] - y[i-1] ) / ( x[i] - x[i-1] );
//		a2[i] = ( x[i+1] - x[i]   ) / 6.0;
//		a3[i] = ( x[i+1] - x[i-1] ) / 3.0;
//		a4[i] = ( x[i]   - x[i-1] ) / 6.0;
//	}
//
//
//	/*
//	Solve the linear system.
//	*/	
//	//penta ( a1, a2, a3, a4, a5, b, ypp, n );
//
//	double xmult;
//
//	for ( i = 1; i < N1 - 1; i++ )
//	{
//		xmult = a2[i] / a3[i-1];
//		a3[i] = a3[i] - xmult * a4[i-1];
//		a4[i] = a4[i] - xmult * a5[i-1];
//		b[i] = b[i] - xmult * b[i-1];
//
//		xmult = a1[i+1] / a3[i-1];
//		a2[i+1] = a2[i+1] - xmult * a4[i-1];
//		a3[i+1] = a3[i+1] - xmult * a5[i-1];
//		b[i+1] = b[i+1] - xmult * b[i-1];
//	}
//
//	xmult = a2[N1-1] / a3[N1-2];
//	a3[N1-1] = a3[N1-1] - xmult * a4[N1-2];
//	ypp[N1-1] = ( b[N1-1] - xmult * b[N1-2] ) / a3[N1-1];
//	ypp[N1-2] = ( b[N1-2] - a4[N1-2] * ypp[N1-1] ) / a3[N1-2];
//
//	for ( i = N1 - 3; 0 <= i; i-- )
//	{
//		ypp[i] = ( b[i] - a4[i] * ypp[i+1] - a5[i] * ypp[i+2] ) / a3[i];
//	}
//
//	/******************************************************************/
//
//	double dt;
//	double h;	
//	int ival;	
//
//	for ( i = 0; i < N2; i++ )
//	{
//		/*
//		Determine the interval [ T(I), T(I+1) ] that contains TVAL.
//		Values below T[0] or above T[N-1] use extrapolation.
//		*/
//		ival = N1 - 2;
//
//		for ( j = 0; j < N1-1; j++ )
//		{
//			if ( xi[i] < x[j+1] )
//			{
//				ival = j;
//				break;
//			}
//		}
//		/*
//		In the interval I, the polynomial is in terms of a normalized
//		coordinate between 0 and 1.
//		*/
//		dt = xi[i] - x[ival];
//		h = x[ival+1] - x[ival];
//
//		yi[i] = y[ival]
//		+ dt * ( ( y[ival+1] - y[ival] ) / h
//			- ( ypp[ival+1] / 6.0 + ypp[ival] / 3.0 ) * h
//			+ dt * ( 0.5 * ypp[ival]
//		+ dt * ( ( ypp[ival+1] - ypp[ival] ) / ( 6.0 * h ) ) ) );
//	}
//}

void splineCubicSpecial ( 
	std::vector<double> x
	, std::vector<double> y
	, std::vector<double> ypp
	, std::vector<double> a1
	, std::vector<double> a2
	, std::vector<double> a3
	, std::vector<double> a4
	, std::vector<double> a5
	, std::vector<double> b
	, unsigned int N1
	, std::vector<double> xi
	, std::vector<double>& yi
	, unsigned int N2
	)
{	
	int i;
	int j;
	int jhi;
	int k;
	double yval;

	/******************************************************************/

	//
	//  Set up the data.
	//

	for ( i = 0; i < N2; i++ )
	{
		a1[i] = 0.0;
		a2[i] = 0.0;
		a3[i] = 0.0;
		a4[i] = 0.0;
		a5[i] = 0.0;
	}
	/*
	Set up the first equation.
	*/

	b[0] = 0.0;
	a3[0] = - ( x[2] - x[1] );
	a4[0] =   ( x[2]        - x[0] );
	a5[0] = - (        x[1] - x[0] );

	/*
	Set up the last equation.
	*/
	b[N1-1]  = 0.0;
	a1[N1-1] = - ( x[N1-1] - x[N1-2] );
	a2[N1-1] =   ( x[N1-1]          - x[N1-3] );
	a3[N1-1] = - (          x[N1-2] - x[N1-3] );

	/*
	Set up the intermediate equations.
	*/
	for ( i = 1; i < N1 - 1; i++ )
	{
		b[i]  = ( y[i+1] - y[i] ) / ( x[i+1] - x[i] )
			- ( y[i] - y[i-1] ) / ( x[i] - x[i-1] );
		a2[i] = ( x[i+1] - x[i]   ) / 6.0;
		a3[i] = ( x[i+1] - x[i-1] ) / 3.0;
		a4[i] = ( x[i]   - x[i-1] ) / 6.0;
	}

	/*
	Solve the linear system.
	*/	
	double xmult;

	for ( i = 1; i < N1 - 1; i++ )
	{
		xmult = a2[i] / a3[i-1];
		a3[i] = a3[i] - xmult * a4[i-1];
		a4[i] = a4[i] - xmult * a5[i-1];
		b[i]  =  b[i] - xmult *  b[i-1];

		xmult = a1[i+1] / a3[i-1];
		a2[i+1] = a2[i+1] - xmult * a4[i-1];
		a3[i+1] = a3[i+1] - xmult * a5[i-1];
		b[i+1]  =  b[i+1] - xmult *  b[i-1];
	}

	xmult = a2[N1-1] / a3[N1-2];
	a3[N1-1]  = a3[N1-1] - xmult * a4[N1-2];
	ypp[N1-1] = ( b[N1-1] - xmult * b[N1-2] ) / a3[N1-1];
	ypp[N1-2] = ( b[N1-2] - a4[N1-2] * ypp[N1-1] ) / a3[N1-2];

	for ( i = N1 - 3; 0 <= i; i-- )
	{
		ypp[i] = ( b[i] - a4[i] * ypp[i+1] - a5[i] * ypp[i+2] ) / a3[i];
	}

	/******************************************************************/

	double dt;
	double h;	
	int ival;	

	for ( i = 0; i < N2; i++ )
	{
		/*
		Determine the interval [ T(I), T(I+1) ] that contains TVAL.
		Values below T[0] or above T[N-1] use extrapolation.
		*/
		ival = N1 - 1;

		for ( j = 0; j < N1-1; j++ )
		{
			if ( xi[i] < x[j+1] )
			{
				ival = j;
				break;
			}
		}
		/*
		In the interval I, the polynomial is in terms of a normalized
		coordinate between 0 and 1.
		*/
		dt = xi[i] - x[ival];
		h  = x[ival+1] - x[ival];

		yi[i] = y[ival]
			+ dt * ( ( y[ival+1] - y[ival] ) / h
			- ( ypp[ival+1] / 6.0 + ypp[ival] / 3.0 ) * h
			+ dt * ( 0.5 * ypp[ival]
			+ dt * ( ( ypp[ival+1] - ypp[ival] ) / ( 6.0 * h ) ) ) );
	}
}

void penta( 
	std::vector<double> a1
	, std::vector<double> a2
	, std::vector<double> a3
	, std::vector<double> a4
	, std::vector<double> a5
	, std::vector<double> b	
	, std::vector<double>& ypp
	, unsigned int N2
	);

void splineCubicTaha( 
	std::vector<double> x
	, std::vector<double> y
	, std::vector<double> ypp
	, std::vector<double> a1
	, std::vector<double> a2
	, std::vector<double> a3
	, std::vector<double> a4
	, std::vector<double> a5
	, std::vector<double> b
	, unsigned int N1
	, std::vector<double> xi
	, std::vector<double>& yi
	, unsigned int N2
	)
{	
	for ( int i = 0; i < N2; i++ )
	{
		a1[i] = 0.0;
		a2[i] = 0.0;
		a3[i] = 0.0;
		a4[i] = 0.0;
		a5[i] = 0.0;
	}
	//
	//  Set up the first equation.
	//

	b[0] = 0.0;
	a3[0] = -( x[2] - x[1] );
	a4[0] =  ( x[2] - x[0] );
	a5[0] = -( x[1] - x[0] );

	//
	//  Set up the intermediate equations.
	//
	for ( int i = 1; i < N2 - 1; i++ )
	{
		b [i] = ( y[i + 1] - y[i] ) / ( x[i + 1] - x[i] ) - ( y[i] - y[i - 1] ) / ( x[i] - x[i - 1] );
		a2[i] = ( x[i + 1] - x[i] ) / 6.0;
		a3[i] = ( x[i + 1] - x[i - 1] ) / 3.0;
		a4[i] = ( x[i] - x[i - 1] ) / 6.0;
	}
	//
	//  Set up the last equation.
	//

	b [N2 - 1] = 0.0;
	a1[N2 - 1] = -( x[N2 - 1] - x[N2 - 2] );
	a2[N2 - 1] =  ( x[N2 - 1] - x[N2 - 3] );
	a3[N2 - 1] = -( x[N2 - 2] - x[N2 - 3] );

	//
	//  Solve the linear system.

	penta( a1, a2, a3, a4, a5, b, ypp, N2 );
	

	double dt;
	double h;
	int ival;

	//
	//  Determine the interval [ T(I), T(I+1) ] that contains TVAL.
	//  Values below T[0] or above T[N-1] use extrapolation.
	//

	ival = N2 - 2;
	for (int tval = 1; tval < N1 + 1; tval++)
	{
		for (int i = 0; i < N2 - 1; i++)
		{
			if (tval < x[i + 1])
			{
				ival = i;
				break;
			}
		}
		//
		//  In the interval I, the polynomial is in terms of a normalized
		//  coordinate between 0 and 1.
		//
		dt = tval - x[ival];
		h = x[ival + 1] - x[ival];

		yi[tval - 1] = y[ival] + dt * ((y[ival + 1] - y[ival]) / h
			- (ypp[ival + 1] / 6.0 + ypp[ival] / 3.0) * h
			+ dt * (0.5 * ypp[ival] + dt * ((ypp[ival + 1] - ypp[ival]) / (6.0 * h))));
	}
}

void penta( 
	std::vector<double> a1
	, std::vector<double> a2
	, std::vector<double> a3
	, std::vector<double> a4
	, std::vector<double> a5
	, std::vector<double> b	
	, std::vector<double>& ypp
	, unsigned int N2
	)
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    PENTA solves a pentadiagonal system of linear equations.
	//
	//  Discussion:
	//
	//    The matrix A is pentadiagonal.  It is entirely zero, except for
	//    the main diagaonal, and the two immediate sub- and super-diagonals.
	//
	//    The entries of Row I are stored as:
	//
	//      A(I,I-2) -> A1(I)
	//      A(I,I-1) -> A2(I)
	//      A(I,I)   -> A3(I)
	//      A(I,I+1) -> A4(I)
	//      A(I,I-2) -> A5(I)
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    07 June 2013
	//
	//  Author:
	//
	//    John Burkardt
	//
	//  Reference:
	//
	//    Cheney, Kincaid,
	//    Numerical Mathematics and Computing,
	//    1985, pages 233-236.
	//
	//  Parameters:
	//
	//    Input, int N, the order of the matrix.
	//
	//    Input, double A1[N], A2[N], A3[N], A4[N], A5[N], the nonzero
	//    elements of the matrix.  Note that the data in A2, A3 and A4
	//    is overwritten by this routine during the solution process.
	//
	//    Input, double B[N], the right hand side of the linear system.
	//
	//    Output, double PENTA[N], the solution of the linear system.
	//
{
	int i;	
	double xmult;

	for (i = 1; i < N2 - 1; i++)
	{
		xmult = a2[i] / a3[i - 1];
		a3[i] = a3[i] - xmult * a4[i - 1];
		a4[i] = a4[i] - xmult * a5[i - 1];
		b[i] = b[i] - xmult * b[i - 1];
		xmult = a1[i + 1] / a3[i - 1];
		a2[i + 1] = a2[i + 1] - xmult * a4[i - 1];
		a3[i + 1] = a3[i + 1] - xmult * a5[i - 1];
		b[i + 1] = b[i + 1] - xmult * b[i - 1];
	}

	xmult = a2[N2 - 1] / a3[N2 - 2];
	a3[N2 - 1] = a3[N2 - 1] - xmult * a4[N2 - 2];
	ypp[N2 - 1] = (b[N2 - 1] - xmult * b[N2 - 2]) / a3[N2 - 1];
	ypp[N2 - 2] = (b[N2 - 2] - a4[N2 - 2] * ypp[N2 - 1]) / a3[N2 - 2];
	for (i = N2 - 3; 0 <= i; i--)
	{
		ypp[i] = (b[i] - a4[i] * ypp[i + 1] - a5[i] * ypp[i + 2]) / a3[i];
	}
}