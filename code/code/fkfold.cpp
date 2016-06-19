#include "kfold.h"

fkfold::fkfold()
{

}
fkfold::~fkfold()
{

}
void fkfold::initialize(snpdt * datafile, FLMW * res)
{
	vector<int> sel_genotype;
	vector<int> sel_snp;
	vector<double> vec_auc;
	_datafile=datafile;
	res->send_res(sel_genotype,sel_snp,vec_auc);
	_scan_sel_genotype=sel_genotype;
	_scan_sel_snp=sel_snp;
	_scan_vec_auc=vec_auc;
	getFamGrp();
}

void fkfold::getFamGrp()
{
	_fam_group.clear();

	vector<int> family_id, family;
	vector< vector<int> > grouped_indi_id;
	int i=0, j=0, n_mk=0;
	bool ex_mk=false;

	int subj_num=_datafile->totalIndi();

	for(i=0; i<subj_num; i++){
		family.push_back(_datafile->fID(i));
	}


	for(i=0; i<subj_num; i++){
		ex_mk=false;
		for(j=0; j<family_id.size(); j++){
			if(family[i]==family_id[j]){
				ex_mk=true;
				n_mk=j;
				break;
			}
		}

		if(ex_mk){
			grouped_indi_id[n_mk].push_back(i);
		}else{
			family_id.push_back(family[i]);
			vector<int> tmp;
			tmp.push_back(i);
			grouped_indi_id.push_back(tmp);
		}
	}

	_fam_group=grouped_indi_id;
}

void fkfold::crossVali(string lr_file)
{
	//cout<<"\n\n\n\nBegin cross validation.....\n\n\n";
	int i=0, j=0, k=0;

	int total=_fam_group.size(); vector<int> idx; vector< vector<int> > devide, devide_residual;
	for(i=0; i< total; i++) idx.push_back(i);
	Stat_fuc::ran_devide(idx,par::n_fold,devide,par::seed,devide_residual);

	vector<int> tmp_indi;
	for(i=0; i<devide.size(); i++){
		tmp_indi.clear();
		for(j=0; j<devide[i].size(); j++){
			int fam_idx=devide[i][j];
			for(k=0; k<_fam_group[fam_idx].size(); k++){
				tmp_indi.push_back(_fam_group[fam_idx][k]);
			}
		}
		devide[i]=tmp_indi;
	}

	for(i=0; i<devide_residual.size(); i++){
		tmp_indi.clear();
		for(j=0; j<devide_residual[i].size(); j++){
			int fam_idx=devide_residual[i][j];
			for(k=0; k<_fam_group[fam_idx].size(); k++){
				tmp_indi.push_back(_fam_group[fam_idx][k]);
			}
		}
		devide_residual[i]=tmp_indi;
	}



	for (i=0;i<par::n_fold;i++)
	{

		if(par::show_crossvali){
			cout<<"Performing cross validation "<<(i+1)<<"\n";
			cout.flush();
		}

		//cout<<"\n\n\nthe "<<(i+1)<<"th fold.....\n";
		_datafile->sampling(devide_residual[i]);

		FLMW mdl_cons, mdl_LR;
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


		FLMW mdl_eval;

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
	FLMW mdl_LR;
	mdl_LR.initialize(_datafile);
	mdl_LR.give_LR(_sel_genotype,_sel_snp,_sel_noexist,_sel_LR_value);

	wt_lr(lr_file);

	for(i=0; i<_nc; i++){
		_finl_sel_genotype.push_back(_scan_sel_genotype[i]);
		_finl_sel_snp.push_back(_scan_sel_snp[i]);
		_finl_vec_auc.push_back(_scan_vec_auc[i]);
	}

	vector<double> auc_var;
	FLMW mdl_eval;
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


void fkfold::noCvWtLr(string lr_file)
{
	_nc=_scan_sel_snp.size();
	_sel_genotype=_scan_sel_genotype;
	_sel_snp=_scan_sel_snp;
	FLMW mdl_LR;
	mdl_LR.initialize(_datafile);
	mdl_LR.give_LR(_sel_genotype,_sel_snp,_sel_noexist,_sel_LR_value);
	mdl_LR.give_nom_LR(_nom_LR);

	wt_lr(lr_file);
}

void fkfold::applyModel(snpdt * datafile, string lr_file)
{
	_datafile=datafile;
	rd_lr(lr_file);

	if(par::match_name) match_nameAllel();

	vector<double> auc_var;
	FLMW mdl_eval;
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


