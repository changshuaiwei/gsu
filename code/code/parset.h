#ifndef _PARSET_H_
#define _PARSET_H_

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

#include "calculate.h"
#include "pamu.h"

class CArg;

using namespace std;

class Pamu;

extern ofstream LOG;
extern Pamu * PP;

//globle parameter
class par
{
public:
//	static **
	static string output_file_name;//file name for output

	static bool silent;//whether cout;

	static double replace0_cor;// =0.5
	static double replace0_sml;// =0.00001

	static bool read_bitfile;
	static string fileroot;
	static string bitfilename;
	static string bitfilename_map;
	static string famfile;

	static bool qt;//quantitative trait
	static bool bt;//binary trait

	static bool coding01;//01 coding 0 for nd and 1 for d;(012 coding 0 for missing 1 for nd and 2 for d)
	static int missing_int;//missing code (if not specified otherwise)
	static string missing_str;//missing code (usually set as . )
	static string out_missing_phenotype;//"-9"
	static bool SNP_major;//
	static bool write_bitfile;//output the binary data format
	static bool out_SNP_major;//output binary file as snp-major

	static string recode_delimit;
	static string recode_indelimit;
	static string out_missing_genotype;

	static bool read_ped;
	static string pedfile;
	static string mapfile;
	static string genfile;
	static bool recode;

	static int most_nsnp;
	static double largest_auc;
	static bool out_nom_LR;
	static string out_nom_LR_f;
	static bool exclude_1;

	static bool lmw_run;
	static string lmw_scan_rst;
	static string lmw_cv_rst;

	static bool flmw_run;

	static int n_fold;
	static int seed;

	static bool choose_first_peak;

	static bool hwu_flt;

	static bool cross_vali;
	static string lmw_subset;
	static string lmw_lr;
	static bool match_name;

	static bool lmw_apl;
	static string lmw_apl_rst;

	static bool flmw_apl;

	static bool multi_core;
	static int n_core;
	static bool force_core;

	static bool gsim_add;
	static bool gsim_dist;

	static bool tamw_run;
	static int tree_depth;
	static int n_classifier;
	static bool clf_ln;
	static bool clf_sqrt;
	static string tamw_lr;
	static bool auto_tree_depth;
	static string tamw_rst;
	static string tamw_smr;

	static int between_ntree;
	static int max_ntree;

	static double thrh_Zscore;
	static int thrh_seltimes;

	static bool moniter;
	static string moniter_f;
	static bool vote;

	static bool tamw_apl;

	static bool ftamw_run;
	static bool ftamw_apl;
	static bool show_iteration;
	static bool show_crossvali;
	static bool show_ntree;
	static bool ckmissing;
	static double missing_rate;

	static bool burnin_tree_depth;
	static int ntree_burnin;

	static bool hwu_run;
	static string hwu_rst;
	static int trl_scheme;
	static int trl_ned;
	static double trl_init_g;
	static double trl_pres;
	static int trl_max;

	static bool IBSwt;
	static bool WIBSwt;
	static int IBS_N;
	static string cov_file;
	static bool cov_read;
	static bool cov2nd_read;
	static string cov2nd_file;
	static bool Wcov;

	static bool PC_Wmt;
	static int N_PC_Wmt;
	static string f_PC_Wmt;
	static bool write_Wmt;
	static bool read_Wmt;
	static string file_Wmt;
	static bool prjct_Y;
	static bool prjct_Y_New_Distr;
	static bool pj_1PC;
	static bool add_cov;
	//static int pj_cov_n;
	static string add_cov_f;
	static bool hwu_rk;

	static bool pjct_OneStep;
	static bool pjct_Backward;
	static bool pj_Y;
	static bool pj_Y_fwd;
	static double pj_threshold;

	static bool pj_rst;
	static string pj_f_Y;
	static string pj_f_PC;

	static bool appx_davis;
	static bool appx_davis_P;
	static bool mm_davis;
	static int n_davis_snp;
	static int n_davis_dim;

	static bool std_weight;

	static bool alt_pheno;
	static string alt_pheno_f;

	static bool ext_snp;
	static string ext_snp_f;
	static double cut_up;
	static bool cut_upB;

	static bool read_set;
	static string setfile;

	static bool read_setC;

	static string pheno_mat_f;
	static bool read_pheno_mat;
	static bool read_pheno_wt;
	static string pheno_wt_f;

	static bool gtu_run;
	static string gtu_rst;
	static int use_trl_threshold;
	static int gtu_trl_dim;

	static bool run_fmt;
	static bool run_csv2ped;
	static string fmt_csv;
	static string fmt_ped;
	
	static bool fgtu_run;
	static bool read_popu_info;
	static string f_popu_info;

	static bool gtu_rank;
	static bool lap_kernel;
	static bool gtu_Gproj;

	static bool kinsp;

};



//globle function
namespace gfun
{
	void NoMem();
	void shutdown();
	void error(string);
	void getOutFileName(CArg &a);
	void setInitialValue();
	void setPar(CArg &a);
	void checkFileExists(string);
	void printLOG(string);
	int split_string(const string &str, vector <string> &vec_str, string separator=" ,\t;\n");
};

//extended std
namespace std
{
	int str2int (const string &str);
	double str2double (const string &str);
	string int2str (int n);
	string double2str (double d);

};


class CArg
{
public:
	CArg(int,char**);

	int count()
	{ return _n; }

	bool any()
	{ return _n > 1 ? true : false; }

	bool find(string);

	void showArg();

	string value(string);

	vector<string> value(string,int);


	vector<string> a;

private:
	int _n;
	vector<bool> _parsed;
	vector<bool> _option;
	vector<string> _root_command;
	vector<string> _original;
	map<string,string> _optionLabel;

};

#endif