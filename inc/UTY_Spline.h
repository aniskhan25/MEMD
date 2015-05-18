#ifndef UTY_SPLINE_H
#define UTY_SPLINE_H

#include <vector>

void chermiteCatmullRom( 
	std::vector<double> x
	, std::vector<double> y
	, unsigned int n1
	, std::vector<double> xi
	, std::vector<double>& yi
	, unsigned int n2 
	);

void chermiteMonotone( 
	std::vector<double> x
	, std::vector<double> y
	, unsigned int n1
	, std::vector<double> xi
	, std::vector<double>& yi
	, unsigned int n2 
	);

void naturalSpline( 
	std::vector<double> x
	, std::vector<double> y
	, unsigned int n
	, std::vector<double>& b
	, std::vector<double>& c
	, std::vector<double>& d
	);

void splineEval( 
	std::vector<double> x
	, std::vector<double> y
	, std::vector<double> b
	, std::vector<double> c
	, std::vector<double> d
	, unsigned int n1
	, std::vector<double> u
	, std::vector<double>& v
	, unsigned int n2
	);

void splineCubic( 
	std::vector<double> x
	, std::vector<double> y
	, unsigned int N1
	, std::vector<double> xi
	, std::vector<double>& yi
	, unsigned int N2
	);

void splineCubicSet ( 
	std::vector<double> t
	, std::vector<double> y
	, std::vector<double>& ypp
	, int n
	);

double splineCubicVal ( 
	std::vector<double> t
	, std::vector<double> y
	, std::vector<double> ypp
	, double tval
	, int n
	);

void penta (
	std::vector<double> a1
	, std::vector<double> a2
	, std::vector<double> a3
	, std::vector<double> a4
	, std::vector<double> a5
	, std::vector<double> b
	, std::vector<double>& x
	, int n
	);

void splineCubic ( 
	std::vector<double> x
	, std::vector<double> y
	, std::vector<double> a1
	, std::vector<double> a2
	, std::vector<double> a3
	, std::vector<double> a4
	, std::vector<double> a5
	, std::vector<double> b
	, std::vector<double> ypp
	, unsigned int N1
	, std::vector<double> xi
	, std::vector<double>& yi
	, unsigned int N2
	);

void splineCubicSet ( 
	std::vector<double> t
	, std::vector<double> y
	, std::vector<double>& ypp
	, std::vector<double> a1
	, std::vector<double> a2
	, std::vector<double> a3
	, std::vector<double> a4
	, std::vector<double> a5
	, std::vector<double> b
	, int n
	);

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
	);

#endif