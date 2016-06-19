#include "kfold.h"

kfold::kfold()
{
	_datafile=0;
}

kfold::~kfold()
{
	_datafile=0;

}

void kfold::initialize(snpdt * datafile, LMW * res)
{
	vector<int> sel_genotype;
	vector<int> sel_snp;
	vector<double> vec_auc;
	_datafile=datafile;
	res->send_res(sel_genotype,sel_snp,vec_auc);
	_scan_sel_genotype=sel_genotype;
	_scan_sel_snp=sel_snp;
	_scan_vec_auc=vec_auc;
}

void kfold::crossVali(string lr_file)
{
	//cout<<"\n\n\n\nBegin cross validation.....\n\n\n";
	int i=0, j=0;

	int total=_datafile->totalIndi(); vector<int> idx; vector< vector<int> > devide, devide_residual;
	for(i=0; i< total; i++) idx.push_back(i);
	Stat_fuc::ran_devide(idx,par::n_fold,devide,par::seed,devide_residual);



	for (i=0;i<par::n_fold;i++)
	{
		if(par::show_crossvali){
			cout<<"Performing cross validation "<<(i+1)<<"\n";
			cout.flush();
		}

		//cout<<"\n\n\nthe "<<(i+1)<<"th fold.....\n";
		_datafile->sampling(devide_residual[i]);

		LMW mdl_cons, mdl_LR;
		mdl_cons.initialize(_datafile);
		mdl_cons.iteration();
		mdl_cons.send_res(_sel_genotype,_sel_snp,_vec_auc);

		_Kf_sel_genotype.push_back(_sel_genotype);
		_Kf_sel_snp.push_back(_sel_snp);
		_Kf_vec_auc.push_back(_vec_auc);

		mdl_LR.initialize(_datafile);
		mdl_LR.give_LR(_sel_genotype, _sel_snp, _sel_noexist, _sel_LR_value);
		
		_datafile->reset();
		_datafile->sampling(devide[i]);


		LMW mdl_eval;

		mdl_eval.initialize(_datafile);
		_vec_auc=mdl_eval.apply_model(_sel_genotype, _sel_snp, _sel_noexist, _sel_LR_value);

		_eval_vec_auc.push_back(_vec_auc);

		_sel_LR_value.clear(); 
		_sel_noexist.clear();
		_sel_snp.clear();
		_sel_genotype.clear();
		_vec_auc.clear();

		_datafile->reset();

	}
	//cout<<"\ncalculating average and nc.....\n";
	cal_average();
	if(par::choose_first_peak) _nc=Stat_fuc::first_peak(_everage_vec_auc)+1;
	else _nc=Stat_fuc::max_index(_everage_vec_auc)+1;


	_sel_genotype=_scan_sel_genotype;
	_sel_snp=_scan_sel_snp;
	LMW mdl_LR;
	mdl_LR.initialize(_datafile);
	mdl_LR.give_LR(_sel_genotype,_sel_snp,_sel_noexist,_sel_LR_value);

	wt_lr(lr_file);

	for(i=0; i<_nc; i++){
		_finl_sel_genotype.push_back(_scan_sel_genotype[i]);
		_finl_sel_snp.push_back(_scan_sel_snp[i]);
		_finl_vec_auc.push_back(_scan_vec_auc[i]);
	}

	vector<double> auc_var;
	LMW mdl_eval;
	mdl_eval.initialize(_datafile);
	_vec_auc=mdl_eval.apply_model(_sel_genotype,_sel_snp,_sel_noexist,_sel_LR_value, auc_var, _nc);
	mdl_eval.give_nom_LR(_nom_LR);
	_test_vec_auc=_vec_auc;
	_test_auc_var=auc_var;
	double z=0, pauc=0;

	_test_auc_p.clear();
	for(i=0; i<_test_vec_auc.size(); i++){
		z=_test_vec_auc[i]-0.5;
		z/=sqrt(auc_var[i]);
		pauc=Stat_fuc::std_norm_p1(z);
		_test_auc_p.push_back(pauc);
	}
}

void kfold::noCvWtLr(string lr_file)
{
	_nc=_scan_sel_snp.size();
	_sel_genotype=_scan_sel_genotype;
	_sel_snp=_scan_sel_snp;
	LMW mdl_LR;
	mdl_LR.initialize(_datafile);
	mdl_LR.give_LR(_sel_genotype,_sel_snp,_sel_noexist,_sel_LR_value);
	mdl_LR.give_nom_LR(_nom_LR);

	wt_lr(lr_file);
}

void kfold::wt_lr(string lr_file)
{
	gfun::printLOG("Writing LR file at [ " + lr_file + " ]\n");

	int i=0, j=0;
	ofstream ff;
	ff.open(lr_file.c_str());

	ff<<"\n\n";

	ff<<"NC:\t"<<_nc<<"\n";

	ff<<"snpname:\t";
	for(i=0; i<_sel_snp.size(); i++){
		ff<<_datafile->snpName(_sel_snp[i])<<"\t";
	}
	ff<<"\n";

	ff<<"allels:\t";
	for(i=0; i<_sel_snp.size(); i++){
		ff<<_datafile->intToGenotype(_sel_snp[i],_sel_genotype[i])<<"\t";
	}
	ff<<"\n";

	ff<<"col_indx:\t";
	for(i=0; i<_sel_snp.size(); i++){
		ff<<_sel_snp[i]<<"\t";
	}
	ff<<"\n";

	ff<<"geno_type:\t";
	for(i=0; i<_sel_genotype.size(); i++){
		ff<<_sel_genotype[i]<<"\t";
	}
	ff<<"\n";


	ff<<"LR_value:\n";
	for(i=0; i<_sel_LR_value.size(); i++){
		ff<<i<<":\t";
		for(j=0; j<_sel_LR_value[i].size(); j++){
			ff<<_sel_LR_value[i][j]<<"\t";
		}
		ff<<"\n";
	}
	ff<<"\n";

	ff<<"no_exist:\n";
	for(i=0; i<_sel_noexist.size(); i++){
		ff<<i<<":\t";
		for(j=0; j<_sel_noexist[i].size(); j++){
			ff<<_sel_noexist[i][j]<<"\t";
		}
		ff<<"\n";
	}

}

void kfold::cal_average()
{
	int n_sel_snp=_scan_vec_auc.size();
	int i=0,j=0;
	_everage_vec_auc.resize(n_sel_snp);
	for(i=0; i<n_sel_snp; i++){
		int count=0;
		double ave=0.0;
		for(j=0; j<par::n_fold; j++){
			if(i<_eval_vec_auc[j].size()){
				if (_eval_vec_auc[j][i]>0.0 &&_eval_vec_auc[j][i]<=1)
				{
					ave+=_eval_vec_auc[j][i];
					count++;
				}
			}

		}

		if (count>0)
		{
			ave/=double(count);
		}else{
			ave=0.0;
		}

		_everage_vec_auc[i]=ave;
	}
}

void kfold::wtResult(string filename)
{
	gfun::printLOG("Writing Cross-Validation Result at [ " + filename + " ]\n");

	ofstream result(filename.c_str());
	if(!result) gfun::error("\ncould not open the result file\n");
	result<<"#result of associatio with testing (0=A1A1, 1=A1A2, 2=A2A2)\n";

	if(_nc>0 && _everage_vec_auc.size()>=_nc) {
		result<<"#Cross-validated average auc (at cut-off) is "<<_everage_vec_auc[_nc-1]<<"\n";
		result<<"#average cv auc:";
		for(int i=0; i<_everage_vec_auc.size(); i++){
			result<<_everage_vec_auc[i]<<"\t";
		}
		result<<"\n";
	}

	result<<"#order\t"<<"Chr\tName\tPos\tBp\tAllel1\tAllel2\t"<<"Geno\t"<<"AUC_tr\t"<<"AUC_ts\t"<<"AUC_var\t"<<"P-value\n";

	Locus * loc=0;
	for(int i=0; i<_nc; i++){
		loc=_datafile->getLocus(_finl_sel_snp[i]);

		result<<i+1<<"\t";//order

		result	<<loc->chr << "\t"
			<< loc->name << "\t"
			<< loc->pos << "\t"
			<< loc->bp << "\t"
			<< loc->allele1 <<"\t"
			<< loc->allele2 <<"\t";//snp information

		result	<<_finl_sel_genotype[i]<<"\t"
			<<_finl_vec_auc[i]<<"\t";//model information

		result<<_test_vec_auc[i]<<"\t"
			<<_test_auc_var[i]<<"\t"
			<<_test_auc_p[i]<<"\n";
	}

	result<<"\n\n#LR for each subject (" + par::missing_str + " for missing)\n";

	result<<"\n#FID\tIID\tPAT\tMAT\tGENDER\tPHENO\t"<<"LR"<<"\n";



	for (int i=0;i<_nom_LR.size();i++)
	{
		Individual * person = _datafile->getIndividual(i);
		result << person->fid << "\t"
			<< person->iid << "\t"
			<< person->pat << "\t"
			<< person->mat << "\t"
			<< person->sexcode << "\t"
			<< person->pheno_str << "\t";

		result<<_nom_LR[i]<<"\n";
	}

	result.close();
}

void kfold::rd_lr(string lr_file)
{

	//cout<<"\nreading the LR record file..\n";
	int i, j;
	string buf_str;
	vector<string> vec_str;
	int n_str;

	ifstream rd(lr_file.c_str());

	while (getline(rd,buf_str))
	{
		if(buf_str.size()==0) continue;
		else if(buf_str[0]=='#') continue;

		if(gfun::split_string(buf_str,vec_str)==0) continue;

		if(vec_str[0]=="NC:") _nc=atoi(vec_str[1].c_str());
		if(vec_str[0]=="snpname:"){
			_sel_names.clear();
			for(i=1; i<vec_str.size(); i++){
				_sel_names.push_back(vec_str[i]);
			}
		}
		if(vec_str[0]=="allels:"){
			_sel_allels.clear();
			for(i=1; i<vec_str.size(); i++){
				_sel_allels.push_back(vec_str[i]);
			}
		}
		if(vec_str[0]=="col_indx:"){
			_sel_snp.clear();
			for(i=1; i<vec_str.size(); i++){
				_sel_snp.push_back(atoi(vec_str[i].c_str()));
			}
		}
		if(vec_str[0]=="geno_type:"){
			_sel_genotype.clear();
			for(i=1; i<vec_str.size(); i++){
				_sel_genotype.push_back(atoi(vec_str[i].c_str()));
			}
		}
		if(vec_str[0]=="LR_value:"){
			_sel_LR_value.clear();
			_sel_LR_value.resize(_sel_snp.size());

			for(i=0; i<_sel_snp.size(); i++){
				getline(rd, buf_str);
				gfun::split_string(buf_str,vec_str);
				for(j=1; j<vec_str.size(); j++){
					_sel_LR_value[i].push_back(atof(vec_str[j].c_str()));
				}
			}
		}
		if(vec_str[0]=="no_exist:"){
			_sel_noexist.clear();
			_sel_noexist.resize(_sel_snp.size());

			for(i=0; i<_sel_snp.size(); i++){
				getline(rd, buf_str);
				gfun::split_string(buf_str,vec_str);
				for(j=1; j<vec_str.size(); j++){
					_sel_noexist[i].push_back(atoi(vec_str[j].c_str()));
				}
			}

		}
	}
}

void kfold::applyModel(snpdt * datafile, string lr_file)
{
	_datafile=datafile;
	rd_lr(lr_file);

	if(par::match_name) match_nameAllel();

	vector<double> auc_var;
	LMW mdl_eval;
	mdl_eval.initialize(_datafile);
	_vec_auc=mdl_eval.apply_model(_sel_genotype,_sel_snp,_sel_noexist,_sel_LR_value, auc_var, _nc);
	_test_vec_auc=_vec_auc;
	_test_auc_var=auc_var;

	int i=0;
	for(i=0; i<_nc; i++){
		_finl_sel_genotype.push_back(_sel_genotype[i]);
		_finl_sel_snp.push_back(_sel_snp[i]);
		_finl_vec_auc.push_back(_vec_auc[i]);
	}

	double z=0, pauc=0;

	_test_auc_p.clear();
	for(i=0; i<_test_vec_auc.size(); i++){
		z=_test_vec_auc[i]-0.5;
		z/=sqrt(auc_var[i]);
		pauc=Stat_fuc::std_norm_p1(z);
		_test_auc_p.push_back(pauc);
	}

}

void kfold::wtSubset(string ssfile)
{
	_datafile->subset(_finl_sel_snp);
	_datafile->writeGenMapPed(ssfile);
	_datafile->reset();
}

void kfold::match_nameAllel()
{
	_sel_snp.clear(); _sel_genotype.clear();
	_datafile->matchNameAllel(_sel_names,_sel_allels,_sel_snp,_sel_genotype);

	for(int i=0; i<_sel_names.size(); i++){
		if(_sel_snp[i]<0) gfun::error("Can't find snp:" + _sel_names[i] + " in the file");
		if(_sel_genotype[i]<0) gfun::error("Can't find corresponding allels at snp:"+ _sel_names[i] + " in the file");
	}

}