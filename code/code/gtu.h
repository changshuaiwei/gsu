#ifndef _GTU_H_
#define _GTU_H_

#include "parset.h"
#include "snpdt.h"
#include "kinship.h"
#include <Eigen/Dense> //for matrix calculation
#include <boost/math/distributions/non_central_chi_squared.hpp> //for non-central chi squared
#include <boost/math/special_functions/fpclassify.hpp>
//using namespace Eigen;

class Dmat;
class Dvec;
class kinship;

///here use float, alternatively we can use MatrixXd and VectorXd;
//typedef Eigen::MatrixXd GTmat_;
//typedef Eigen::VectorXd GTvec_;
//typedef Eigen::ArrayXXd GTarr_;//two dims array
//typedef Eigen::MatrixXi GTimat_;

class GTU
{
public:
	GTU(){
		_datafile=0;
		_n_sub=0;
		_buf=0;
	}

	~GTU(){
		_datafile=0;
		_n_sub=0;
		delete _buf;
	}

	void initialize(snpdt * snp_data);
	void run();
	void wtResult(string outputfile);
protected:

	void normalizeY();
	virtual void transX();

	virtual void PJ_Tsim();
	virtual void PJ_Gsim();
	virtual void FamStruc();

	void getTraitSimSqr(GTmat_ & Y, GTmat_ & Tsim, GTvec_ & weight);
	void getTraitSimCov(GTmat_ & Y, GTmat_ & Tsim, GTvec_ & weight);

	void getTraitSimLap(GTmat_ & Y, GTmat_ & Tsim, GTvec_ & weight);

	void standardize(GTmat_ & a);
	void standardize(GTmat_ & a, GTmat_ & X);
	void standardizeX(GTmat_ & a);

	bool wIBS(vector<int> & idx, bool LK=false);

	void GTUcore();

	void EigenSmall(GTmat_ & a);
	void EigenLarge(GTmat_ & a);
	void EigenLargeCore(GTmat_ & a,GTmat_ & eigenvalues);
	double Liu(GTmat_ & a, double Q);
	double ncChiSqSurvival(double df, double ncp, double q);

	snpdt * _datafile;

	vector< vector<double> > _phe_ori;
	vector<double> _phe_wt;
	vector< vector<double> > _cov_ori;

	vector< vector<int> > _snp_set;

	vector<double> _pVec;
	vector<double> _QVec;

	vector<double> _Liu_df;
	vector<double> _Liu_q;
	vector<double> _Liu_ncp;

	int _n_sub;

	GTmat_ _Y;
	GTmat_ _X;

	GTmat_ _One;

	GTvec_ _Yweight;

	GTmat_ _Tsim;//final trait similarity
	GTmat_ _Gsim;//final genetic similarity

	///using tralan
	Dmat * _buf;
	double ** _bufPointer;

};

class FGTU : public GTU
{
public:
	FGTU(){

	}

	~FGTU(){

	}


protected:

	void FamStruc();
	void PJ_Tsim();
	void PJ_Gsim();
	void transX();
	GTmat_ _KSmat;

	GTmat_ _KSEigenValues;
	GTmat_ _KSEigenVectors;
	GTmat_ _weightKS;
};


#endif