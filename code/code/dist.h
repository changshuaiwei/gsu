#ifndef _DIST_H_
#define _DIST_H_

#include <cmath>
#include <limits>
#include <iostream>
#include <fstream>
#include <vector>
#include "calculate.h"

using namespace std;

namespace Dist
{
	vector< vector<double> > euclid(vector< vector<double> > & dt);

	void meansd(vector<double> & dt, double & miu, double & sd);
	double euclid_AB(vector<double> & a, vector<double> & b);

	void projection(vector<double> Y, vector<double> & e, vector< vector<double> > & X);//require X is eigen vector

	void projectionSolu(vector<double> Y, vector<double> & e, 
		vector<double> & beta, vector<bool> & significant, vector< vector<double> > & X);//require X is eigen vector

	vector<int> projectionOnestep(vector<double> Y, vector<double> & e, vector< vector<double> > & X);//require X is eigen vector

	vector<int> projectionBackward(vector<double> Y, vector<double> & e, vector< vector<double> > & X);//require X is eigen vector



}


#endif