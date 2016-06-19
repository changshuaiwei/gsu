#include "pamu.h"

ofstream LOG;
Pamu * PP;

int main(int argc, char* argv[]) 
{

	set_new_handler(gfun::NoMem);

	snpdt SDT;//store all the data
	
	Pamu P; //perform all the function
	P.initialize(&SDT);
	PP=&P;



	//////set initial value for globe parameters
	gfun::setInitialValue();

	// Command line arguments

	CArg a(argc,argv);

	gfun::getOutFileName(a);
	LOG.open(string(par::output_file_name + ".log").c_str());
	P.printLOG("\n"
		"@----------------------------------------------------------@\n"
		"|        GSU!                                              |\n"
		"|----------------------------------------------------------|\n"
		"|        Generalized Similarity U Test                     |\n"
		"|----------------------------------------------------------|\n"
		"|        Copyright (C) 2016                                |\n"
		"|----------------------------------------------------------|\n"
		"|  Changshuai Wei, PhD.                                    |\n"
		"|        http://changshuaiwei.github.io/                   |\n"
		"@----------------------------------------------------------@\n"
		"\n");

	gfun::setPar(a);

	// Time stamp

	P.printLOG("Writing [ "+
		par::output_file_name + ".log ]\n");

	time_t curr=time(0);
	string tdstamp = (string)ctime(&curr);
	P.printLOG("Analysis started: " + tdstamp +"\n");

	///format
	if(par::run_fmt) P.fmtFMT();

	////input data
	if (par::read_bitfile) P.readBinData();
	else if(par::read_ped) P.readPedData();

	///get covariate data
	if(par::cov_read) P.readCovData();

	///get alternative phenotype
	if(par::alt_pheno) P.readPheData();

	///get multiple phenotypes
	if(par::read_pheno_mat) P.readPheMat();

	//get weight for multivariate phenotype
	if(par::read_pheno_wt) P.readPheWt();

	///read snp set info
	if(par::read_set || par::read_setC) P.readSnpSet();

	//read population infor
	if(par::read_popu_info) 
	{
		P.readPopuInfo();
	}


	///sample code
	//--gtu --bfile NIC_SUB --set set.txt --cov cov.txt --phs phe.txt --phs-wt wt.txt --rank --pjG --out NS

	if(par::gtu_run){
		P.assocGTU();
	}

	if(par::fgtu_run){
		P.assocFGTU();
	}

	//--bfile ch6 --kp --wIBS-wt --write-Wmt --PCwt 20 --out ccc
	if(par::kinsp){
		P.calKP();
	}
	
	//hwu --bfile ch6 --hwu --IBS-wt --out ccc
	//hwu --file NDs --cov cov.txt --hwu --dv 8 --wIBS-wt --add-cov cov2.txt --write-Wmt --out ccc
	//hwu --file NDs --cov cov.txt --hwu --dv 8 --gsim-dist --cov-wt cov2.txt --write-Wmt --out ccc
	//hwu --file NDs --hwu --dv 8 --read-Wmt ccc.wmt --PCwt 2 --flt --out ccc2
	//hwu --file NDs --hwu --dv 8 --read-Wmt ccc.wmt --PCwt 2 --mp 10 --force-core --out ccc2
	if(par::hwu_run){
		P.assocHWU();
	}

	////forward selection algorithm
	////sample command 1:pamu --bfile GVT2D_HPFS --lmw --out HPFS
	////sample command 2:pamu --bfile GVT2D_HPFS --lmw --no-cv --out HPFS
	if(par::lmw_run) {
		if(par::cross_vali) P.cvLMW();
		else P.assocLMW();
	}

	///Appling model
	///sample command:pamu --bfile ND_90 --lmw-apl ND.lmw.lr --out NDpdt
	if(par::lmw_apl) P.aplLMW();

	//forward for family data
	if(par::flmw_run) {
		if(par::cross_vali) P.cvFLMW();
		else P.assocFLMW();
	}

	if(par::flmw_apl) P.aplFLMW();

	//tree assembling method
	////sample command:pamu --bfile ND_90 --tamw --converge --tree-depth 8 --out NDta
	if(par::tamw_run) P.scanTAMW();
	////sample command:pamu --bfile ND_90 --tamw-apl NDta.tamw.lr --out NDta
	if(par::tamw_apl) P.aplTAMW();

	//family tree assembling method
	if(par::ftamw_run) P.scanFTAMW();
	if(par::ftamw_apl) P.aplFTAMW();

	/////output data
	if (par::write_bitfile) P.write_BITFILE();
	else if(par::recode) P.write_PEDFILE();

	gfun::shutdown();

}


void Pamu::initialize(snpdt * data)
{
	_data=data;
}

void Pamu::clear()
{
	_data->clear();
}

void Pamu::printLOG(string s)
{
	LOG << s;
	LOG.flush();

	if (!par::silent)
	{
		cout << s;
		cout.flush();
	}
}

void Pamu::readBinData()
{
	_data->readBinData();
}

void Pamu::readPedData()
{
	_data->readGenMapPed();
}

void Pamu::readCovData()
{
	_data->readCovFile(par::cov_file);

	if(par::cov2nd_read){
		_data->read2ndCovFile(par::cov2nd_file);
	}
}

void Pamu::readPheData()
{
	_data->readPheFile(par::alt_pheno_f);
}

void Pamu::readPopuInfo()
{
	_data->readPopuFile(par::f_popu_info);
}

void Pamu::readPheMat()
{
	_data->readPheMat(par::pheno_mat_f);
}

void Pamu::readPheWt()
{
	_data->readWtFile(par::pheno_wt_f);
}

void Pamu::readSnpSet()
{
	if(par::read_set){
		_data->readSnpSet();
	}else if(par::read_setC){
		_data->readSnpSet2Cols();
	}
}

void Pamu::write_BITFILE()
{
	vector<int> idx;
	if(par::ext_snp){
		idx=_data->readExtSnpF(par::ext_snp_f);
		_data->subset(idx);
		_data->write_BITFILE();
		_data->reset();
	}else{
		_data->write_BITFILE();
	}
	
}

void Pamu::write_PEDFILE()
{
	vector<int> idx;
	if(par::ext_snp){
		idx=_data->readExtSnpF(par::ext_snp_f);
		_data->subset(idx);
		_data->writeGenMapPed();
		_data->reset();
	}else{
		_data->writeGenMapPed();
	}
}

snpdt* Pamu::dataPointer()
{
	snpdt * tmp=_data;
	return tmp;
}

void Pamu::assocGTU()
{
	gfun::printLOG("\n*****GTU:Gene-Trait U*****\n");

	gfun::printLOG("Initializing...\n");
	GTU asso;
	asso.initialize(_data);

	time_t curr=time(0);
	string tdstamp = (string)ctime(&curr);
	gfun::printLOG("Scanning started at: " + tdstamp +"\n");

	asso.run();


	asso.wtResult(par::gtu_rst);

}

void Pamu::assocFGTU()
{
	gfun::printLOG("\n*****FGTU:Family Gene-Trait U*****\n");

	if(!par::read_popu_info){
		_data->readPopuFile(par::famfile);
	}

	gfun::printLOG("Initializing...\n");
	FGTU asso;
	asso.initialize(_data);

	time_t curr=time(0);
	string tdstamp = (string)ctime(&curr);
	gfun::printLOG("Scanning started at: " + tdstamp +"\n");

	asso.run();


	asso.wtResult(par::gtu_rst);

}

void Pamu::calKP()
{
	gfun::printLOG("\n*****Calculate Kinship related matrices*****\n");

	gfun::printLOG("Initializing...\n");

	//HWU asso;//this is the old class, will not use any more

	egHWU asso;
	asso.initialize_min(_data);

	time_t curr=time(0);
	string tdstamp = (string)ctime(&curr);
	gfun::printLOG("Scanning started at: " + tdstamp +"\n");

	gfun::printLOG("Get relateness Matrix...\n");
	
	if(par::WIBSwt) asso.weightFromGenomeWIBS();
	else if(par::read_Wmt) asso.rdWeight(par::file_Wmt);

	if(par::write_Wmt) asso.wtWeight(par::file_Wmt);
	if(par::PC_Wmt) asso.weightPC(par::N_PC_Wmt, par::f_PC_Wmt);

}

void Pamu::assocHWU()
{
	gfun::printLOG("\n*****HWU:Heterogeneity Weighted U for association*****\n");

	gfun::printLOG("Initializing...\n");

	//HWU asso;//this is the old class, will not use any more

	egHWU asso;
	asso.initialize(_data);

	gfun::printLOG("Get relateness Matrix...\n");
	//asso.weightFromSNP(4);
	//asso.weightFromGenomeIBS();
	if(par::Wcov) asso.weightFromCov();//this one uses the second sets of covariates
	else if(par::WIBSwt) asso.weightFromGenomeWIBS();
	else if(par::IBSwt) asso.weightFastGenome();
	else if(par::read_Wmt) asso.rdWeight(par::file_Wmt);

	if(par::add_cov) asso.updateWeight();

	if(par::write_Wmt) asso.wtWeight(par::file_Wmt);
	if(par::PC_Wmt) asso.weightPC(par::N_PC_Wmt, par::f_PC_Wmt);

	time_t curr=time(0);
	string tdstamp = (string)ctime(&curr);
	gfun::printLOG("Scanning started at: " + tdstamp +"\n");

	//gfun::printLOG("Scanning...\n");

	if(par::appx_davis || par::appx_davis_P){
		asso.sgLocusRank();
	}else{
		asso.sgLocusAssoc();
	}
	
	asso.wtResult(par::hwu_rst);

}

void Pamu::assocLMW()
{
	gfun::printLOG("\n*****LMW without Cross-validation*****\n");

	gfun::printLOG("Forward Scanning...\n");
	LMW asso;
	asso.initialize(_data);
	asso.iteration();
	asso.wtSubset(par::lmw_subset);
	asso.wtResult(par::lmw_scan_rst);

	kfold cv;
	cv.initialize(_data, &asso);
	cv.noCvWtLr(par::lmw_lr);
}

void Pamu::assocFLMW()
{

	gfun::printLOG("\n*****FLMW without Cross-validation*****\n");

	gfun::printLOG("Forward Scanning...\n");
	FLMW asso;
	asso.initialize(_data);
	asso.iteration();
	asso.wtSubset(par::lmw_subset);
	asso.wtResult(par::lmw_scan_rst);

	fkfold cv;
	cv.initialize(_data, &asso);
	cv.noCvWtLr(par::lmw_lr);
}

void Pamu::cvLMW()
{
	gfun::printLOG("\n*****LMW with Cross-validation*****\n");

	gfun::printLOG("Forward Scanning...\n");

	LMW asso;
	asso.initialize(_data);
	asso.iteration();
	asso.wtSubset(par::lmw_subset);
	asso.wtResult(par::lmw_scan_rst);

	gfun::printLOG("Cross validating...\n");
	kfold cv;
	cv.initialize(_data, &asso);
	cv.crossVali(par::lmw_lr);
	cv.wtResult(par::lmw_cv_rst);
}

void Pamu::cvFLMW()
{

	gfun::printLOG("\n*****FLMW with Cross-validation*****\n");

	gfun::printLOG("Forward Scanning...\n");
	FLMW asso;
	asso.initialize(_data);
	asso.iteration();
	asso.wtSubset(par::lmw_subset);
	asso.wtResult(par::lmw_scan_rst);

	gfun::printLOG("Cross validating...\n");
	fkfold cv;
	cv.initialize(_data, &asso);
	cv.crossVali(par::lmw_lr);
	cv.wtResult(par::lmw_cv_rst);
}

void Pamu::aplLMW()
{
	gfun::printLOG("\n*****Apply LMW*****\n");

	gfun::printLOG("Applying model...\n");
	kfold apl;
	apl.applyModel(_data,par::lmw_lr);
	apl.wtResult(par::lmw_apl_rst);
	apl.wtSubset(par::lmw_subset);
}

void Pamu::aplFLMW()
{
	gfun::printLOG("\n*****Apply FLMW*****\n");

	gfun::printLOG("Applying model...\n");
	fkfold apl;
	apl.applyModel(_data,par::lmw_lr);
	apl.wtResult(par::lmw_apl_rst);
	apl.wtSubset(par::lmw_subset);
}

void Pamu::scanTAMW()
{
	gfun::printLOG("\n*****TAMW*****\n");

	gfun::printLOG("Building assembling model...\n");
	TAMW ta;
	ta.initialize(_data);
	ta.assembling();
	ta.wtResult(par::tamw_rst);
	ta.wtPredictorStat(par::tamw_smr);
	if(par::moniter) ta.wtConverge(par::moniter_f);
}

void Pamu::aplTAMW()
{
	gfun::printLOG("\n*****Apply TAMW*****\n");

	gfun::printLOG("Applying assembling model...\n");
	TAMW ta;
	ta.initialize(_data);
	ta.apply_model(par::tamw_lr);
	ta.wtResult(par::tamw_rst);
	ta.wtPredictorStat(par::tamw_smr);
	if(par::moniter) ta.wtConverge(par::moniter_f);
}

void Pamu::scanFTAMW()
{
	gfun::printLOG("\n*****FTAMW*****\n");

	gfun::printLOG("Building assembling model...\n");

	FTAMW ta;
	ta.initialize(_data);
	ta.assembling();
	ta.wtResult(par::tamw_rst);
	ta.wtPredictorStat(par::tamw_smr);
	if(par::moniter) ta.wtConverge(par::moniter_f);
}

void Pamu::aplFTAMW()
{
	gfun::printLOG("\n*****Apply FTAMW*****\n");

	gfun::printLOG("Applying assembling model...\n");
	FTAMW ta;
	ta.initialize(_data);
	ta.apply_model(par::tamw_lr);
	ta.wtResult(par::tamw_rst);
	ta.wtPredictorStat(par::tamw_smr);
	if(par::moniter) ta.wtConverge(par::moniter_f);
}

void Pamu::fmtFMT()
{
	gfun::printLOG("\n*****Format Data*****\n");

	if(par::run_csv2ped){
		gfun::printLOG("\nFrom CSV to PED file\n");

		FMT::csv2tped(par::fmt_csv, par::fmt_ped);
	}

}









