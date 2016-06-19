#ifndef _HWU_H_
#define _HWU_H_



#include "dist.h"
#include "trl.h"
#include "qf.h"

#include "parset.h"
#include "snpdt.h"

#include <Eigen/Dense> //for matrix calculation
#include <boost/math/distributions/non_central_chi_squared.hpp> //for non-central chi squared
#include <boost/math/special_functions/fpclassify.hpp>

//#define _DEBUG_HWU_WEI_

//typedef Eigen::MatrixXd GTmat_;

//typedef Eigen::VectorXd GTvec_;

//typedef Eigen::ArrayXXd GTarr_;//two dims array

//typedef Eigen::MatrixXi GTimat_;

class Dmat;
class Dvec;

class HWU
{
public:
	HWU(){
		_datafile=0;
		//_buffer_cal=0;
		_outer_y=0;
		_weight=0;
		_weight_2=0;
		_outer_x=0;
		_Projection=0;
		_buffer=0;
		_weight_flag=false;
		//_buffer_on=false;

		_weight_corct=0;

	}
	~HWU(){
		_datafile=0;
		//delete [] _buffer_cal;
		//_outer_y->clear();
		//_weight->clear();
		//_outer_x->clear();
		delete _outer_y;
		delete _weight;
		delete _weight_2;
		delete _outer_x;
		delete _Projection;
		//delete [] _buffer;
	}

	void initialize(snpdt * snp_data);

	void weightFromSNP(int size=2);
	void weightFromCov();

	void weightFromGenomeWIBS();
	void weightFromGenomeIBS();
	void weightFastGenome();

	void updateWeight();

	void sgLocusAssoc();//using moment matching, get p-value for all snp, or replace some by davis method
	void sgLocusRank();//using davis method, get p-value for highest ranked U statistics

	double hwu_i(int i_snp);

	void wtResult(string outputfile);
	void wtWeight(string outputfile);
	void rdWeight(string inputfile);

protected:

	void getRankWU();
	void getRankP();
	void updateP();
	void updatePcut();

	void prepareY();
	void prepareX(int i_snp, vector<double> &x);
	vector<double> sry(vector<double> & y);
	vector<double> getrk(vector<double> & y);
	void projectY();
	void pjY();
	void projectY2();//only use the first PC;
	//void projectY3();//project Y on genetic and covariate; seems not good

	void projectY(vector< vector<double> > &P, vector<double> & newY, vector<double> & oldY);

	void rdAddCov(string filename);

	void projectionMat(vector< vector<double> > & X, vector< vector<double> > & P);
	void prepProject(vector< vector<double> > & P);

	//void outer_Y();
	void outer_YXW();

	void outer_X(vector<double> & x, vector< vector<double> > & ot_X);

	void outer_X(vector<double> & x);

	void update_X();

	double hwu_i_internal(vector<double> & x);
	double hwu_i_first_PC(vector<double> & x);
	double hwu_test(vector<double> & x, vector<double> & y);
	//double hwu_i_internalP(vector<double> & x); //this func is too slow!!

	double hwu_inter_U(vector<double> & x);
	double hwu_inter_p(vector<double> & x);
	double hwu_davis(vector<double> & x);

	void princ_comp(vector< vector<double> > & X);
	void princ_comp(vector< vector<double> > & X, vector< vector<double> > & mat);

	void stdW();
	
	snpdt * _datafile;
	int _n_snp;
	int _n_sub;

	/*
	vector<double> _y;//standardized ranked trait value
	vector<double> _org_y;

	vector< vector<double> > _outer_y;
	vector< vector<double> > _weight;
	vector<double> _vec_pvalue;

	//double * _buffer_cal;//work space for op
	*/

	vector<double> _y;//standardized ranked trait value
	vector<double> _org_y;

	Dmat * _outer_y;
	Dmat * _weight;
	Dmat * _weight_2;
	Dmat * _outer_x;
	Dmat * _Projection;
	vector< vector<float> > _PW;

	vector<double> _vec_pvalue;
	vector<double> _bk_vec_pvalue;
	vector<double> _weightU;
	vector<double> _Obs;
	vector<double> _expct;
	vector<double> _varian;
	vector<double> _bk_weightU;
	vector<int> _idx;

	vector<int> _x_int;
	vector<double> _x_double;
	vector<double> _first_PC;

	vector< vector<double> > _add_cov;
	vector<string> _cov_name;

	//bool _buffer_on;
	double ** _buffer;

	bool _weight_flag;

	double _weight_corct;

};




class egHWU : public HWU
{
public:
	void initialize(snpdt * snp_data);

	void initialize_min(snpdt * snp_data);

	void weightFromCov();
	void weightFromGenomeWIBS();
	void weightFastGenome();
	void updateWeight();
	void wtWeight(string outputfile);
	void rdWeight(string inputfile);

	void weightPC(int size, string outputfile);

	void sgLocusAssoc();
	void sgLocusRank();
	void wtResult(string outputfile);

protected:
	void rankY();
	void prepareY();
	void getRankZ();
	void getRankZfloat();

	void getRankZfltMp();//multiple thread

	void calZMp(int start, int length);
	//void calZMpp(int start, int length, GTmatF_ & Wt, GTmatF_ & Xf, GTmatF_ & Yf, GTmatF_ & XXinvf);
	void splitCal(int n_td, vector<int> & idx_start, vector<int> & idx_length);

	void D2F();

	void stdW();

	void wtZscore(string outputfile);

	void sortEigen(GTmat_ & egvec, GTmat_ & egval);//by abslute value

	void PCEigen(GTmat_ & mat, GTmat_ & PC, int size);

	inline double gsim(double g1, double g2);
	inline float gsim(float g1, float g2);
	
	double hwu_inter(vector<double> & x, GTmat_ & Wt);
	double hwu_inter_Z(vector<double> & x);//standardized score, used for rank
	double hwu_inter_p(vector<double> & x);//not accurate, use normal approx
	double hwu_liu(vector<double> & x);

	double hwu_inter_Zfloat(vector<float> & x);
	//double hwu_inter_ZfltMP(vector<float> & x, GTmatF_ & Wt, GTmatF_ & Xf, GTmatF_ & Yf, GTmatF_ & XXinvf);
	double hwu_interFloat(vector<float> & x, GTmatF_ & Wt);

	void transX();
	void EDweight(GTmat_ & cov, GTmat_ & Kappa, GTvec_ & weight);
	void EDweight(GTmat_ & cov, GTvec_ & weight);

	void StandMatByCol(GTmat_ &Mt);

	double Liu(GTmat_ & a, double Q);
	double Liu(GTmat_ & a, double Q, double & df_, double & ncp_, double & q_);
	double ncChiSqSurvival(double df, double ncp, double q);

	void EigenLargeCore(GTmat_ & a,GTmat_ & eigenvalues);
	void EigenSmallCore(GTmat_ & a,GTmat_ & eigenvalues);

	void EigenLargeCore(GTmat_ & a, GTmat_ & eigenvectors, GTmat_ & eigenvalues);
	void EigenSmallCore(GTmat_ & a, GTmat_ & eigenvectors, GTmat_ & eigenvalues);

	void prepareXfloat(int i_snp, vector<float> &x);

	vector<double> _q_vec;
	vector<double> _df_vec;
	vector<double> _ncp_vec;

	vector< vector<double> > _cov_ori;
	vector< vector<double> > _cov4wt_ori;
	vector<double> _wtcov4wt_ori;

	GTmat_ _Y;
	GTmat_ _X;

	GTmat_ _XXinv;

	GTmatF_ _Yflt;
	GTmatF_ _Xflt;

	GTmatF_ _XXinvflt;

	GTmat_ _One;

	GTmat_ _Kappa;

	GTmat_ _Kappaflt;

	int _SNPcounter;

	boost::shared_mutex _mutex;
	boost::shared_mutex _mutex2;
	//boost::shared_mutex _mutex3;
	//boost::shared_mutex _mutex4;

	//GTmat_ _Tsim;//final trait similarity
	//GTmat_ _Gsim;//final genetic similarity
};


#endif