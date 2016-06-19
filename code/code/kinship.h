#ifndef _KINSHIP_H_
#define _KINSHIP_H_

#include "parset.h"
#include "snpdt.h"

#include <Eigen/Dense> //for matrix calculation

//typedef Eigen::MatrixXd GTmat_;
//typedef Eigen::VectorXd GTvec_;
//typedef Eigen::ArrayXXd GTarr_;//two dims array
//typedef Eigen::MatrixXi GTimat_;

class Individual;

class kinship
{
public:
	kinship(){

	}

	~kinship(){

	}

	void initialize(snpdt *snp_data);
	
	void ksMAT(GTmat_ & ks);
	
	void find_sample();

	void construct_pedi();
	double cal_kinship(Individual* indiv1, Individual* indiv2);
	void common_ans(vector<Individual*> & ans, Individual* indiv1, Individual* indiv2);
	void ancestor(vector<Individual*> & ans, Individual* indiv);
	void cal_rout(vector< vector<Individual*> > & routs, Individual* indiv, Individual* ans);
	double cal_inbreeding(Individual* indiv);
	double rout_coeff(vector< vector<Individual*> > & routs1, vector< vector<Individual*> > & routs2);
	bool compare_rout(vector<Individual*> & rout1, vector<Individual*> & rout2);

private:
	void search_tree(Individual* tmp, Individual* aim, Individual* root);
	vector<Individual*> _ori_data;
	vector<Individual*> _sample;
	vector< vector<Individual*> > _tmp_routs;
	vector<int> _sample_idx;
};


#endif