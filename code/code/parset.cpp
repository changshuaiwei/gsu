#include "parset.h"

/*****globle function*******/

void gfun::NoMem()
{
	cerr << "*****************************************************\n"
		<< "* FATAL ERROR    Exhausted system memory            *\n"
		<< "*                                                   *\n"
		<< "* You need a smaller dataset or a bigger computer...*\n"
		<< "*                                                   *\n"
		<< "* Forced exit now...                                *\n"
		<< "*****************************************************\n\n";
	exit(1);
}

void gfun::shutdown()
{
	time_t curr=time(0);
	string tdstamp = ctime(&curr);

	if (!par::silent) cout << "\nAnalysis finished: " + tdstamp +"\n";
	LOG << "\nAnalysis finished: " + tdstamp +"\n";

	LOG.close();

	PP->clear();


	//system("PAUSE");
	exit(0);
}

void gfun::error(string msg)
{
	cerr << "\nERROR: " << msg << "\n";
	LOG << "\nERROR: " << msg << "\n";  
	LOG.close();

	PP->clear();

	exit(1);
}

void gfun::getOutFileName(CArg &a)
{
	if (a.find("--out")) 
	{
		par::output_file_name = a.value("--out");
	}

	if (a.find("--silent")) 
	{
		par::silent = true;
	}

}

void gfun::setInitialValue()
{

}

void gfun::setPar(CArg &a)
{
	gfun::printLOG("Command: "); a.showArg();

	if (a.find("--bfile")) 
	{ 
		par::read_bitfile = true;
		par::fileroot = a.value("--bfile");
		par::bitfilename = par::fileroot + ".bed";
		par::bitfilename_map= par::fileroot + ".bim";
		par::famfile = par::fileroot + ".fam";
	}

	if (a.find("--file")) 
	{ 
		par::read_ped = true;
		par::fileroot = a.value("--file");
		par::pedfile = par::fileroot + ".ped";
		par::mapfile = par::fileroot + ".map";
		par::genfile = par::fileroot + ".gen";
	}

	if(a.find("--hwu"))
	{
		par::hwu_run=true;
		par::qt=true;
		par::hwu_rst=par::output_file_name + ".scan";
		par::use_trl_threshold=1600;
	}

	if(a.find("--lmw"))
	{
		par::lmw_run=true;
		par::lmw_scan_rst=par::output_file_name + ".scan";
		par::lmw_cv_rst=par::output_file_name + ".cv";
		par::lmw_subset=par::output_file_name + ".subset";
		par::lmw_lr=par::output_file_name + ".lmw.lr";
	}

	if(a.find("--lmw-apl")){
		par::lmw_apl=true;
		par::lmw_lr=a.value("--lmw-apl");
		par::lmw_apl_rst=par::output_file_name + ".apl";
		par::lmw_subset=par::output_file_name + ".subset";
	}

	if(a.find("--flmw"))
	{
		par::flmw_run=true;
		par::lmw_scan_rst=par::output_file_name + ".scan";
		par::lmw_cv_rst=par::output_file_name + ".cv";
		par::lmw_subset=par::output_file_name + ".subset";
		par::lmw_lr=par::output_file_name + ".flmw.lr";
	}

	if(a.find("--flmw-apl")){
		par::flmw_apl=true;
		par::lmw_lr=a.value("--flmw-apl");
		par::lmw_apl_rst=par::output_file_name + ".apl";
		par::lmw_subset=par::output_file_name + ".subset";
	}

	if(a.find("--tamw")){
		par::tamw_run=true;
		par::tamw_lr=par::output_file_name + ".tamw.lr";
		par::tamw_rst=par::output_file_name + ".tamw.rst";
		par::tamw_smr=par::output_file_name + ".tamw.smr";
		par::moniter_f=par::output_file_name + ".tamw.monitor";
	}

	if(a.find("--tamw-apl")){
		par::tamw_apl=true;
		par::tamw_lr=a.value("--tamw-apl");
		par::tamw_rst=par::output_file_name + ".tamw.rst";
		par::tamw_smr=par::output_file_name + ".tamw.smr";
		par::moniter_f=par::output_file_name + ".tamw.monitor";
	}

	if(a.find("--ftamw")){
		par::ftamw_run=true;
		par::tamw_lr=par::output_file_name + ".ftamw.lr";
		par::tamw_rst=par::output_file_name + ".ftamw.rst";
		par::tamw_smr=par::output_file_name + ".ftamw.smr";
		par::moniter_f=par::output_file_name + ".ftamw.monitor";
	}

	if(a.find("--ftamw-apl")){
		par::tamw_apl=true;
		par::tamw_lr=a.value("--ftamw-apl");
		par::tamw_rst=par::output_file_name + ".ftamw.rst";
		par::tamw_smr=par::output_file_name + ".ftamw.smr";
		par::moniter_f=par::output_file_name + ".ftamw.monitor";
	}

	if(a.find("--vote"))  par::vote=true;

	if(a.find("--converge")) par::moniter=true;

	if(a.find("--match-names")){
		par::match_name=true;
	}

	if(a.find("--maxauc")) par::largest_auc=str2double(a.value("--maxauc"));

	if(a.find("--no-cv")) par::cross_vali=false; else par::cross_vali=true;

	if(a.find("--lr-sub"))
	{
		par::out_nom_LR=true;
		par::out_nom_LR_f=par::output_file_name + ".nom.lr";
	}

	if(a.find("--seed")) par::seed=str2int(a.value("--seed"));

	if(a.find("--clsf")) {
		par::n_classifier=str2int(a.value("--clsf"));
		par::clf_ln=false; par::clf_sqrt=false;
	}

	if(a.find("--clsf-sqrt")){
		par::clf_sqrt=true; par::clf_ln=false;
	}

	if(a.find("--tree-depth")){
		par::tree_depth=str2int(a.value("--tree-depth"));
	}

	if(a.find("--ntree")){
		par::max_ntree=str2int(a.value("--ntree"));
	}

	if(a.find("--auto-depth")){
		par::auto_tree_depth=true;
	}

	if (a.find("--make-bed")) 
	{
		par::write_bitfile = true;
	}

	if (a.find("--recode")) 
	{
		par::recode = true;
	}

	if(a.find("--showcv"))
	{
		par::show_crossvali=true;
		par::show_iteration=true;
	}

	if(a.find("--showscan")){
		par::show_iteration=true;
	}

	if(a.find("--showntree")){
		par::show_ntree=true;
	}

	if(a.find("--maxcvauc")){
		par::choose_first_peak=false;
	}

	if(a.find("--hiorder")){
		par::most_nsnp=str2int(a.value("--hiorder"));
	}

	if(a.find("--ckmissing")){
		par::ckmissing=true;
		par::missing_rate=str2double(a.value("--ckmissing"));
	}

	if(a.find("--td-burnin")){
		par::burnin_tree_depth=true;
	}

	if(a.find("--nburnin")){
		par::ntree_burnin=str2int(a.value("--nburnin"));
	}

	if(a.find("--hz")){
		par::exclude_1=false;
	}

	if(a.find("--qt")){
		par::qt=true; //affect whether read as double
	}

	if(a.find("--schm")){
		par::trl_scheme=str2int(a.value("--schm"));
	}

	if(a.find("--maxop")){
		par::trl_max=str2int(a.value("--maxop"));
	}

	if(a.find("--pres")){
		par::trl_pres=str2double(a.value("--pres"));
	}

	if(a.find("--ned")){
		par::trl_ned=str2int(a.value("--ned"));
	}

	if(a.find("--trint")){
		par::trl_init_g=str2double(a.value("--trint"));
	}

	if(a.find("--cov")){
		par::cov_read=true;
		par::cov_file=a.value("--cov");

		//if(a.find("--cov-wt")) par::Wcov=true;
	}

	if(a.find("--cov-wt")){
		par::cov2nd_read=true;
		par::Wcov=true;
		par::cov2nd_file=a.value("--cov-wt");

		//if(a.find("--cov-wt")) par::Wcov=true;
	}

	if(a.find("--mp")){
		par::multi_core=true;
		par::n_core=str2int(a.value("--mp"));

		//The below will have memory leak
		//if(a.find("--force-core")){
		//	par::force_core=true;
		//}
	}

	if(a.find("--kp")){
		par::kinsp=true;
	}

	if(a.find("--IBS-wt")){
		par::IBSwt=true;
		par::IBS_N=str2int(a.value("--IBS-wt"));
	}
	if(a.find("--wIBS-wt")){
		par::WIBSwt=true;
		par::IBS_N=str2int(a.value("--wIBS-wt"));
	}

	if(a.find("--write-Wmt")) {
		par::write_Wmt=true;
		par::file_Wmt=par::output_file_name + ".wmt";
	}
	if(a.find("--read-Wmt")) {
		par::read_Wmt=true;
		par::file_Wmt=a.value("--read-Wmt");
	}

	if(a.find("--gsim-dist")){
		par::gsim_dist=true; par::gsim_add=false;
	}

	if(a.find("--pjct")){
		par::prjct_Y=true;
	}

	if(a.find("--pjsig")){
		par::prjct_Y=true;
		par::pjct_OneStep=true;
		par::pjct_Backward=false;
	}

	if(a.find("--pjbwd")){
		par::prjct_Y=true;
		par::pjct_OneStep=false;
		par::pjct_Backward=true;
	}

	if(a.find("--pj")){
		par::pj_Y=true;
	}

	if(a.find("--pj-rst")){
		par::pj_rst=true;
		par::pj_f_Y=par::output_file_name +".pj.Y.txt";
		par::pj_f_PC=par::output_file_name+".pj.PC.txt";
	}

	if(a.find("--PCwt")){
		par::PC_Wmt=true;
		par::f_PC_Wmt=par::output_file_name + ".PCwt.txt";
		par::N_PC_Wmt=str2int(a.value("--PCwt"));
	}

	if(a.find("--flt")){
		par::hwu_flt=true;
	}


	if(a.find("--pjf")){
		par::pj_Y=true;
		par::pj_Y_fwd=true;
	}

	if(a.find("--pj-cut")){
		par::pj_threshold=str2double(a.value("--pj-cut"));
	}

	//if(a.find("--pj-exact")){
	//	par::prjct_Y_New_Distr=true;
	//}

	if(a.find("--pjpc")){
		par::pj_1PC=true;
	}

	if(a.find("--pjG")){
		par::gtu_Gproj=true;
	}

	if(a.find("--add-cov")){
		par::add_cov=true;
		par::add_cov_f=a.value("--add-cov");
	}

	if(a.find("--rank")){
		par::hwu_rk=true;
		par::gtu_rank=true;
	}

	if(a.find("--dv")){
		par::appx_davis=true;
		par::n_davis_snp=str2int(a.value("--dv"));

		if(a.find("--dim")){
			par::n_davis_dim=str2int(a.value("--dim"));
		}
	}

	if(a.find("--dvp")){
		par::appx_davis_P=true;
		par::n_davis_snp=str2int(a.value("--dvp"));

		if(a.find("--dim")){
			par::n_davis_dim=str2int(a.value("--dim"));
		}
	}

	if(a.find("--stdw")){
		par::std_weight=true;
	}

	if(a.find("--pheno")){
		par::alt_pheno=true;
		par::alt_pheno_f=a.value("--pheno");
	}

	if(a.find("--ext-snp")){
		par::ext_snp=true;
		par::ext_snp_f=a.value("--ext-snp");
	}

	if(a.find("--no-md")){
		par::mm_davis=false;
	}

	if(a.find("--cutup")){
		par::cut_upB=true;
		par::cut_up=str2double(a.value("--cutup"));
	}

	if (a.find("--set")) 
	{
		par::read_set = true;
		par::setfile = a.value("--set");
	}

	if (a.find("--setc")) 
	{
		par::read_setC = true;
		par::setfile = a.value("--setc");
	}

	if(a.find("--phs"))
	{
		par::read_pheno_mat=true;
		par::pheno_mat_f=a.value("--phs");
	}

	if(a.find("--phs-wt")){
		par::read_pheno_wt=true;
		par::pheno_wt_f=a.value("--phs-wt");
	}

	if(a.find("--gtu"))
	{
		par::gtu_run=true;
		par::gtu_rst=par::output_file_name + ".scan";
	}

	if(a.find("--LK"))
	{
		par::lap_kernel=true;
	}

	if(a.find("--popu"))
	{
		par::read_popu_info=true;
		par::f_popu_info=a.value("--popu");
	}

	if(a.find("--fgtu"))
	{
		par::fgtu_run=true;
		par::gtu_rst=par::output_file_name + ".scan";
	}


	if(a.find("--fmt"))
	{
		par::run_fmt=true;
	}

	if(a.find("--c2p"))
	{
		par::run_csv2ped=true;
		par::fmt_csv=a.value("--c2p");
		par::fmt_ped=par::output_file_name + ".tped";
	}
}

void gfun::checkFileExists(string f)
{

	ifstream inp;

	inp.open(f.c_str(), ifstream::in);
	if(inp.fail())
	{
		inp.clear(ios::failbit);
		inp.close();
		string msg = "No file [ " + f + " ] exists.";
		gfun::error(msg);
	}
	inp.close();
	return;

}

void gfun::printLOG(string s)
{
	LOG << s;
	LOG.flush();

	if (!par::silent)
	{
		cout << s;
		cout.flush();
	}
}

int gfun::split_string(const string &str, vector <string> &vec_str, string separator)
{
	if(str.empty()) return 0;
	vec_str.clear();

	int i=0;
	bool look=false;
	string str_buf;
	string symbol_pool="`1234567890-=~!@#$%^&*()_+qwertyuiop[]\\asdfghjkl;'zxcvbnm,./QWERTYUIOP{}|ASDFGHJKL:\"ZXCVBNM<>? \t\n";
	string::size_type pos;

	for(i=0; i<separator.size(); i++) 
	{
		pos=symbol_pool.find(separator[i]);
		if( pos!=string::npos ) symbol_pool.erase(symbol_pool.begin()+pos);
	}

	for(i=0; i<str.size(); i++)
	{		
		if( symbol_pool.find(str[i])!=string::npos )
		{
			if(!look) look=true;
			str_buf += str[i];
		}
		else
		{
			if(look) 
			{
				look=false;
				vec_str.push_back(str_buf);
				str_buf.erase(str_buf.begin(), str_buf.end());
			}
		}
	}
	if(look) vec_str.push_back(str_buf);

	return vec_str.size();
}

/**********std**************/

int std::str2int (const string &str) {
	stringstream ss(str);
	int n;
	ss >> n;
	return n;
}

string std::int2str (int n) {
	stringstream ss;
	ss << n;
	return ss.str();
}

string std::double2str(double d){
	stringstream ss;
	ss << d;
	return ss.str();
}

double std::str2double(const string &str)
{
	stringstream ss(str);
	double n;
	ss >> n;
	return n;
}



/*****CArg class*******/

CArg::CArg(int n, char *argv[])
{
	_n = n;
	a.push_back(argv[0]);
	_parsed.resize(_n,false);
	_option.resize(_n,false);
	for (int i=1 ; i < _n ; i++ )
		a.push_back(argv[i]);
	_original = a;
}

void CArg::showArg()
{
	for(int i=0; i<a.size(); i++){
		gfun::printLOG(a[i]+" ");
	}
	gfun::printLOG("\n");
}

bool CArg::find(string s)
{  
	for (int i=0 ; i < _n ; i++ )
		if (a[i] == s) { _parsed[i]=true; return true; }
		return false;
}

string CArg::value(string s)
{
	for (int i=0 ; i < _n ; i++ )
		if (a[i] == s && (i+1 < _n) ) { _parsed[i+1]=true; return a[i+1]; }
		gfun::error("Missing an argument for "+s);
} 

vector<string> CArg::value(string s, int c)
{
	vector<string> r(0);

	for (int i=0 ; i < _n ; i++ )
		if (a[i] == s && (i+1 < _n) ) 
		{
			for (int j=1;j<=c;j++) 
			{
				if ( (i+j) < a.size() )
				{
					_parsed[i+j]=true;
					r.push_back(a[i+j]);
				}
				else
					gfun::error("Not enough arguments given for option: "+s+" ");
			} 
		}

		if (r.size() != c) gfun::error("Not enough arguments given for option: "+s+" ");
		return r;
} 



///////////////////////

string par::output_file_name="pamuRst";
bool par::silent=false;
double par::replace0_cor=0.5;
double par::replace0_sml=0.00001;
bool par::read_bitfile=false;
string par::fileroot="";
string par::bitfilename="";
string par::bitfilename_map="";
string par::famfile="";
bool par::qt=false;//quantitative trait
bool par::bt=true;
bool par::coding01=false;
int par::missing_int=-9;
string par::missing_str=".";
string par::out_missing_phenotype = "-9";
bool par::SNP_major=true;
bool par::write_bitfile=false;
bool par::out_SNP_major = true;
string par::recode_delimit = "\t";
string par::recode_indelimit = "";
string par::out_missing_genotype = "0";

bool par::read_ped = false;
string par::pedfile = "";
string par::mapfile = "";
string par::genfile = "";
bool par:: recode=false;

int par::most_nsnp=10;
double par::largest_auc=0.9;
bool par::out_nom_LR=false;
string par::out_nom_LR_f="pamuRst.nom.lr";
bool par::exclude_1=true;

bool par::lmw_run=false;
bool par::flmw_run=false;
string par::lmw_scan_rst="pamuRst.scan";
string par::lmw_cv_rst="pamuRst.cv";
int par::n_fold=10;
int par::seed=-2011;
bool par::choose_first_peak=true;
bool par::hwu_flt=false;
bool par::cross_vali=true;
string par::lmw_subset="pamuRst.subset";
string par::lmw_lr="pamuRst.lmw.lr";
bool par::match_name=false;
bool par::lmw_apl=false;
bool par::flmw_apl=false;
string par::lmw_apl_rst="pamuRst.apl";
int par::tree_depth=6;
int par::n_classifier=2;
bool par::clf_ln=true;
bool par::clf_sqrt=false;
string par::tamw_lr="pamuRst.tamw.lr";
bool par::auto_tree_depth=false;
int par::between_ntree=10;
int par::max_ntree=2000;
double par::thrh_Zscore=0;
int par::thrh_seltimes=5;
bool par::tamw_run=false;
string par::tamw_rst="pamuRst.tamw.rst";
string par::tamw_smr="pamuRst.tamw.smr";
bool par::moniter=false;
string par::moniter_f="pamuRst.tamw.monitor";
bool par::tamw_apl=false;
bool par::ftamw_run=false;
bool par::ftamw_apl=false;
bool par::show_iteration=false;
bool par::show_crossvali=false;
bool par::show_ntree=false;
bool par::ckmissing=false;
double par::missing_rate=0.1;
bool par::burnin_tree_depth=false;
int par::ntree_burnin=20;
bool par::vote=false;
bool par::hwu_run=false;
int par::trl_scheme=0;
int par::trl_ned=20;
double par::trl_init_g=1e+40;
double par::trl_pres=0.00001;
int par::trl_max=200;
string par::hwu_rst="hwu.rst";
bool par::IBSwt=false;
bool par::WIBSwt=false;
int par::IBS_N=-9;
string par::cov_file="cov.txt";
bool par::cov_read=false;
bool par::cov2nd_read=false;
string par::cov2nd_file="cov2.txt";
bool par::Wcov=false;
bool par::write_Wmt=false;
bool par::read_Wmt=false;
string par::file_Wmt="weight.wmt";
bool par::prjct_Y=false;
bool par::prjct_Y_New_Distr=false;
bool par::pj_1PC=false;
bool par::add_cov=false;
//int par::pj_cov_n=15;
string par::add_cov_f="add_cov_f.txt";
bool par::hwu_rk=false;
bool par::pjct_OneStep=false;
bool par::pjct_Backward=false;
bool par::pj_Y=false;
bool par::pj_rst=false;
string par::pj_f_Y="pj.Y.txt";
string par::pj_f_PC="pj.PC.txt";
bool par::PC_Wmt=false;
int par::N_PC_Wmt=20;
string par::f_PC_Wmt="";
bool par::pj_Y_fwd=false;
double par::pj_threshold=0.2;
bool par::appx_davis=false;
bool par::appx_davis_P=false;
bool par::mm_davis=true;
int par::n_davis_snp=1000;
int par::n_davis_dim=20;
bool par::multi_core=false;
int par::n_core=1;
bool par::force_core=false;
bool par::gsim_add=true;
bool par::gsim_dist=false;
bool par::std_weight=false;
bool par::alt_pheno=false;
string par::alt_pheno_f="";
bool par::ext_snp=false;
string par::ext_snp_f="";
double par::cut_up=1e-6;
bool par::cut_upB=false;
bool par::read_set=false;
string par::setfile="";
bool par::read_setC=false;
string par::pheno_mat_f="";
bool par::read_pheno_mat=false;
bool par::read_pheno_wt=false;
string par::pheno_wt_f="";
bool par::gtu_run=false;
string par::gtu_rst="";
int par::use_trl_threshold=5000;
int par::gtu_trl_dim=400;
bool par::run_fmt=false;

bool par::read_popu_info=false;
string par::f_popu_info="";
bool par::fgtu_run=false;

bool par::gtu_Gproj=false;
bool par::gtu_rank=false;

bool par::lap_kernel=false;

bool par::kinsp=false;

bool par::run_csv2ped=false;
string par::fmt_csv="csv.csv";
string par::fmt_ped="ped.ped";






