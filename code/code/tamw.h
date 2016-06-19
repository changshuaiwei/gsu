#ifndef _TAMW_H_
#define _TAMW_H_

#include "lmw.h"

class TAMW
{
public:
	TAMW();
	~TAMW();
	void initialize(snpdt * data);
	void assembling();
	void wtResult(string file);
	void wtPredictorStat(string file);
	void wtConverge(string file);
	void apply_model(string model_file);

protected:
	void produce_sequence(int tn_col, int & seed);
	virtual void step_tree(int tree_idx);
	void set_n_classifier();
	void rt_frst(int tree_indx, string forest_file);
	virtual double cal_rf_auc();
	virtual double cal_rf_auc(double & var, double & P);
	bool rd_treemodel(ifstream & rd);
	void match_nameAllel();
	int burninDepth();
	virtual int step_tmp_tree();
	

	snpdt * _data;

	vector< vector<int> > _iter_seq;
	vector<int> _snp_map;

	vector<int> _tmp_btstping;
	vector<int> _tmp_bts;
	vector<int> _tmp_oob;

	vector<double> _temp_vec_auc;
	vector<string> _temp_sel_names;
	vector<string> _temp_sel_allels;
	vector<int> _temp_sel_genotype;
	vector<int> _temp_sel_snp;
	vector< vector<double> > _temp_LR_value;
	vector< vector<int> > _temp_noexist;

	vector< vector<double> > _snp_stat;
	vector< vector<double> > _indi_LRs;

	vector<double> _moniter_auc;
	vector<int> _n_tree;

	double _fin_auc;
	double _fin_auc_var;
	double _fin_auc_P;
};

class FTAMW : public TAMW
{
public:
	void initialize(snpdt * data);
	void apply_model(string model_file);

protected:
	void step_tree(int tree_idx);
	void getFamGrp(vector<int> subjIdx, vector< vector<int> > & grouped_indi_id, vector< vector<int> > & compare_struct);
	double cal_rf_auc();
	double cal_rf_auc(double & var, double & P);
	int step_tmp_tree();

	vector< vector<int> > _fam_group;
	vector< vector<int> >_compare_struct;

};

#endif