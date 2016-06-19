#ifndef _KFOLD_H_
#define _KFOLD_H_

#include "lmw.h"

using namespace std;

class LMW;
class FLMW;

class kfold{
public:
	kfold();
	~kfold();
	void initialize(snpdt * datafile, LMW * res);
	void crossVali(string lr_file);
	void noCvWtLr(string lr_file);
	void wtResult(string filename);
	void wtSubset(string);
	void applyModel(snpdt * datafile, string lr_file);

protected:
	void cal_average();
	void wt_lr(string lr_file);
	void rd_lr(string lr_file);
	void match_nameAllel();

	snpdt * _datafile;

	int _nc;

	vector<int> _scan_sel_genotype;
	vector<int> _scan_sel_snp;
	vector<double> _scan_vec_auc;

	vector<int> _finl_sel_genotype;
	vector<int> _finl_sel_snp;
	vector<double> _finl_vec_auc;

	vector< vector<int> > _Kf_sel_genotype;
	vector< vector<int> > _Kf_sel_snp;
	vector< vector<double> > _Kf_vec_auc;

	vector< vector<double> > _eval_vec_auc;

	vector<string> _sel_names;
	vector<string> _sel_allels;

	vector<int> _sel_genotype;
	vector<int> _sel_snp;
	vector<double> _vec_auc;
	vector<double> _eval_auc;

	//for missing data
	vector< vector<int> > _sel_noexist;
	vector< vector<double> > _sel_LR_value;
	vector<double> _nom_LR;

	vector<double> _everage_vec_auc;

	vector<double> _test_vec_auc;
	vector<double> _test_auc_p;
	vector<double> _test_auc_var;


};

class fkfold : public kfold
{
public:
	fkfold();
	~fkfold();
	void initialize(snpdt * datafile, FLMW * res);
	void crossVali(string lr_file);
	void noCvWtLr(string lr_file);
	void applyModel(snpdt * datafile, string lr_file);

protected:
	void getFamGrp();
	vector< vector<int> > _fam_group;
};

#endif

