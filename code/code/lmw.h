#ifndef _LMW_H_
#define _LMW_H_

#include <string>
#include <vector>
#include <iostream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <map>
#include <fstream>
#include <set>
#include <functional>
#include <new>
#include <sstream>

#include "parset.h"
#include "snpdt.h"

class snpdt;

class LMW
{
public:
	LMW();
	~LMW();

	void initialize(snpdt * snp_data);
	void iteration();
	void wtResult(string outputfile);
	void wtSubset(string ssfile);

	void give_LR(vector<int> & sel_genotype, vector<int> & sel_snp, vector< vector<int> > & no_exist, vector< vector<double> > & sel_LR_value);
	vector<double> apply_model(vector<int> & sel_genotype, vector<int> & sel_snp, vector< vector<int> > & no_exist, vector< vector<double> > & sel_LR_value);
	vector<double> apply_model(vector<int> & sel_genotype, vector<int> & sel_snp, vector< vector<int> > & no_exist, vector< vector<double> > & sel_LR_value, vector<double> & auc_var, int nc);
	void send_res(vector<int> & sel_genotye, vector<int> & sel_snp, vector<double> & vec_auc);
	void give_nom_LR(vector<double> & LR);

	void grow_tree(vector< vector<int> > & seq);
	vector<double> apply_test_LR(vector<int> & sel_genotype, vector<int> & sel_snp, vector< vector<int> > & no_exist,
		vector< vector<double> > & sel_LR_value, vector<double> & indi_LR);

protected:
	
	void step_auc();
	void patient_devide(bool apply=false);
	void col_auc(int col);
	void ranking_auc();
	void deviding_to_next_mis(bool apply =false );
	void out_nom_lr(string out);
	

	void apply_devide_mis(vector<double> & LR_value, vector<int> & no_exist);

	void rdm_step_auc();
	void ranking_auc2();

	virtual void temp_deviding_mis(int col, int genotype);
	virtual double auc_fm_tmpLR();
	virtual double auc_fm_LR();
	virtual double auc_fm_LR(double & var);
	virtual void auc_VC_fm_LR(vector<double> & pre_score, vector<double> & score, int step);

	snpdt * _datafile;

	vector<int> _sel_snp;
	vector<double> _vec_auc;
	vector<double> _var_auc;
	vector<double> _cov_auc;
	vector<int> _sel_genotype;

	int _marker_num;
	int _subj_num;
	double _n_health;
	double _n_disease;

	vector<bool> _aff;

	vector<double> _temp_sel_auc;
	vector<int> _temp_sel_genotype;

	vector< vector<bool> > _left_genotype;// for every col in the data, indicating whether it's selected
	vector < short > _left_genotype_num;// for every col in the data

	vector< vector<int> > _devided_patient;
	vector< vector<int> > _devided_health;
	vector< vector<int> > _pre_devided_patient;
	vector< vector<int> > _pre_devided_health;
	vector< vector<int> > _mis_patient;
	vector< vector<int> > _mis_health;
	vector<double> _mis_LR;

	vector<double> _cur_tmp_LR;
	vector<double> _pre_tmp_LR;

	vector<double> _LR_nom;
	vector<double> _tmp_LR_nom;

	vector<int> _temp_noexist_rcd; 
	vector<double> _temp_LR;
	vector<int> _temp_devided_np;
	vector<int> _temp_devided_nh;

	vector< vector<double> > _vec_LR_value;

	vector<int> _temp_col_indx;
	vector<int> _temp_col_rcd;

};




class FLMW : public LMW
{
public:
	FLMW ();
	~FLMW();
	void initialize(snpdt * snp_data);

protected:
	void gen_contrust();
	void temp_deviding_mis(int col, int genotype);
	double auc_fm_tmpLR();
	double auc_fm_LR();
	double auc_fm_LR(double & var);
	void auc_VC_fm_LR(vector<double> & pre_score, vector<double> & score, int step);

	vector< vector<int> > _compare_struct;
};

#endif

