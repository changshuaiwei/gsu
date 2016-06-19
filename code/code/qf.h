#ifndef _QF_H_
#define _QF_H_


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <setjmp.h>
#include <vector>
#include <string>
#include <iostream>

using namespace std;

//#define UseDouble 0             /* all floating point double */
#define qfTRUE  1
#define qfFALSE 0
typedef int qfBOOL;
//#ifdef UseDouble
//typedef double real;
//#else
//typedef float real;
//#endif

#define qf_pi 3.14159265358979
#define qf_log28 .0866  /*  log(2.0) / 8.0  */

namespace QF{


	static double sigsq, lmax, lmin, mean, c;
	static double intl, ersm;
	static int count, r, lim;  static qfBOOL ndtsrt, fail;
	static int *n,*th; static double *lb,*nc;
	static jmp_buf env;

	static double exp1(double x);
	static void counter(void);
	static double square(double x);
	static double cube(double x);
	static double  log1(double x, qfBOOL first);
	static void order(void);
	static double   errbd(double u, double* cx);
	static double  ctff(double accx, double* upn);
	static double truncation(double u, double tausq);
	static void findu(double* utx, double accx);
	static void integrate(int nterm, double interv, double tausq, qfBOOL mainx);
	static double cfe(double x);
	double   qf(double* lb1, double* nc1, int* n1, int r1, double sigma, double c1,
		int lim1, double acc, double* trace, int* ifault);

	double qf(vector<double> lb1, vector<double> nc1, vector<int> n1,
		double sigma, double c1, int lim1, double acc, 
		vector<double> &trace, int &ifault);

	double qf(vector<double> lb1, vector<double> nc1, vector<int> n1,
		double sigma, double c1, int lim1, 
		vector<double> &trace, int &ifault);

	double qf(vector<double> lb1, double c1, vector<double> &trace, int &ifault);

	void status(vector<double> &trace, int & ifault);

}


#endif