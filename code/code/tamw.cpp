#include "tamw.h"

TAMW::TAMW()
{

}

TAMW::~TAMW()
{

}

void TAMW::initialize(snpdt * data)
{
	_data=data;
	_indi_LRs.resize(_data->totalIndi());
	_snp_stat.resize(_data->totalLocus());

	set_n_classifier();
}

void TAMW::set_n_classifier()
{
	if(par::clf_sqrt){
		double temp=double(_data->totalLocus());
		temp=sqrt(temp);
		par::n_classifier=int(temp);
		if(par::n_classifier==0) par::n_classifier=1;
	}else{
		if(par::clf_ln){
			double temp=double(_data->totalLocus());
			temp=log(temp);
			par::n_classifier=int(temp);
			if(par::n_classifier==0) par::n_classifier=1;
		}
	}
}


void TAMW::produce_sequence(int tn_col, int & seed)
{

	int i=0, indx=0, j=0, in_idx=-1;
	_iter_seq.clear();
	_iter_seq.resize(par::tree_depth);
	double dtn_col=tn_col;

	vector<int> tmp_cols; tmp_cols.resize(tn_col,-9);
	_snp_map.clear();

	for(j=0; j<par::tree_depth; j++){
		i=0;
		while(i<par::n_classifier){

			indx=int(dtn_col*(Stat_fuc::ran1(seed)));

			if(indx>0 && indx<tn_col){

				if(tmp_cols[indx]==-9){
					in_idx++;
					tmp_cols[indx]=in_idx;
					_snp_map.push_back(indx);
					_iter_seq[j].push_back(in_idx);
				}else{
					_iter_seq[j].push_back(tmp_cols[indx]);
				}
				i++;
			}
		}
	}
}


void TAMW::step_tree(int tree_idx)
{
	int i=0, j=0;
	produce_sequence(_data->totalLocus(),par::seed);

	Stat_fuc::btstrp_sampling(_data->totalIndi(), _tmp_btstping, par::seed);

	_tmp_bts.clear(); _tmp_oob.clear();
	for(i=0; i<_tmp_btstping.size(); i++){
		if(_tmp_btstping[i]==0) _tmp_oob.push_back(i);
		else{
			for(j=0; j<_tmp_btstping[i]; j++)_tmp_bts.push_back(i);
		}
	}

	_data->subset(_tmp_bts,_snp_map);

	LMW cal_tree;
	cal_tree.initialize(_data);
	cal_tree.grow_tree(_iter_seq);
	cal_tree.send_res(_temp_sel_genotype, _temp_sel_snp, _temp_vec_auc);

	LMW cal_LR;
	cal_LR.initialize(_data);
	cal_LR.give_LR(_temp_sel_genotype, _temp_sel_snp,_temp_noexist, _temp_LR_value);

	_data->reset();
	_data->subset(_tmp_oob,_snp_map);

	if(par::auto_tree_depth){
		LMW apl_tree_tres;
		apl_tree_tres.initialize(_data);
		vector<double> tmp_auc;
		tmp_auc=apl_tree_tres.apply_model(_temp_sel_genotype, _temp_sel_snp, _temp_noexist, _temp_LR_value);

		double max=0;
		int maxidx=0;

		for(i=0; i<tmp_auc.size(); i++){
			if(tmp_auc[i]>max){
				max=tmp_auc[i]; maxidx=i;
			}else{
				if(par::choose_first_peak) break;
			}
		}
		int cut=maxidx+1;

		_temp_sel_genotype.resize(cut); _temp_sel_snp.resize(cut); _temp_noexist.resize(cut); _temp_LR_value.resize(cut);
	}

	rt_frst(tree_idx,par::tamw_lr);


	LMW apl_tree;
	apl_tree.initialize(_data);
	vector<double> temp_indi_LR;
	vector<double> temp_stat;

	temp_stat=apl_tree.apply_test_LR(_temp_sel_genotype, _temp_sel_snp, _temp_noexist, _temp_LR_value, temp_indi_LR);

	

	for(i=0; i<_temp_sel_snp.size(); i++){
		_snp_stat[_snp_map[_temp_sel_snp[i]]].push_back(temp_stat[i]);

	}

	if(par::vote){
		for(i=0; i<temp_indi_LR.size(); i++){
			int indx=_tmp_oob[i];
			double tmpscore=0.5; if(temp_indi_LR[i]>1) tmpscore=1; if(temp_indi_LR[i]<1) tmpscore=0;
			_indi_LRs[indx].push_back(tmpscore);
		}
	}else{
		for(i=0; i<temp_indi_LR.size(); i++){
			int indx=_tmp_oob[i];
			_indi_LRs[indx].push_back(temp_indi_LR[i]);
		}
	}
	

	_data->reset();

}

int TAMW::step_tmp_tree()
{
	int i=0, j=0;
	produce_sequence(_data->totalLocus(),par::seed);

	Stat_fuc::btstrp_sampling(_data->totalIndi(), _tmp_btstping, par::seed);

	_tmp_bts.clear(); _tmp_oob.clear();
	for(i=0; i<_tmp_btstping.size(); i++){
		if(_tmp_btstping[i]==0) _tmp_oob.push_back(i);
		else{
			for(j=0; j<_tmp_btstping[i]; j++)_tmp_bts.push_back(i);
		}
	}

	_data->subset(_tmp_bts,_snp_map);

	LMW cal_tree;
	cal_tree.initialize(_data);
	cal_tree.grow_tree(_iter_seq);
	cal_tree.send_res(_temp_sel_genotype, _temp_sel_snp, _temp_vec_auc);

	LMW cal_LR;
	cal_LR.initialize(_data);
	cal_LR.give_LR(_temp_sel_genotype, _temp_sel_snp,_temp_noexist, _temp_LR_value);

	_data->reset();
	_data->subset(_tmp_oob,_snp_map);


	LMW apl_tree_tres;
	apl_tree_tres.initialize(_data);
	vector<double> tmp_auc;
	tmp_auc=apl_tree_tres.apply_model(_temp_sel_genotype, _temp_sel_snp, _temp_noexist, _temp_LR_value);

	double max=0;
	int maxidx=0;

	for(i=0; i<tmp_auc.size(); i++){
		if(tmp_auc[i]>max){
			max=tmp_auc[i]; maxidx=i;
		}else{
			if(par::choose_first_peak) break;
		}
	}
	int cut=maxidx+1;
	_data->reset();

	return cut;
}

void TAMW::rt_frst(int tree_indx, string forest_file)
{
	//cout<<"\n\nnow output the "<<tree_indx<<"th tree to file...\n\n";

	long i=0, j=0;
	ofstream ff;
	if(tree_indx==0) ff.open(forest_file.c_str());
	else ff.open(forest_file.c_str(), ios_base::app);
	ff<<"\n\n";

	ff<<"tree_indx\t"<<tree_indx<<"\n";

	ff<<"snpname:\t";
	for(i=0; i<_temp_sel_snp.size(); i++){
		ff<<_data->snpName(_temp_sel_snp[i])<<"\t";
	}
	ff<<"\n";

	ff<<"allels:\t";
	for(i=0; i<_temp_sel_snp.size(); i++){
		ff<<_data->intToGenotype(_temp_sel_snp[i],_temp_sel_genotype[i])<<"\t";
	}
	ff<<"\n";

	ff<<"col_indx:\t";
	for(i=0; i<_temp_sel_snp.size(); i++){
		ff<<_snp_map[_temp_sel_snp[i]]<<"\t";
	}
	ff<<"\n";

	ff<<"geno_type:\t";
	for(i=0; i<_temp_sel_genotype.size(); i++){
		ff<<_temp_sel_genotype[i]<<"\t";
	}
	ff<<"\n";


		ff<<"LR_value:\n";
		for(i=0; i<_temp_LR_value.size(); i++){
			ff<<i<<":\t";
			for(j=0; j<_temp_LR_value[i].size(); j++){
				ff<<_temp_LR_value[i][j]<<"\t";
			}
			ff<<"\n";
		}
		ff<<"\n";

		ff<<"no_exist:\n";
		for(i=0; i<_temp_noexist.size(); i++){
			ff<<i<<":\t";
			for(j=0; j<_temp_noexist[i].size(); j++){
				ff<<_temp_noexist[i][j]<<"\t";
			}
			ff<<"\n";
		}


	ff<<"\ntree_end\n";
}

void TAMW::assembling()
{
	if(par::burnin_tree_depth){
		par::tree_depth=burninDepth();
	}

	int i=0;
	int distance=0;
	while(1){
		//cout<<"\n\nnow add the "<<i<<"th tree into forest...\n\n";

		if(par::show_ntree){
			cout<<"Construct tree "<<i+1;
			cout.flush();
		}

		step_tree(i);

		if(par::show_ntree){
			cout<<"\r";
			cout.flush();
		}

		if(distance%par::between_ntree==0){
			double auc=cal_rf_auc();
			_moniter_auc.push_back(auc);
			_n_tree.push_back(i);
		}
		distance++;

		if(i+1>=par::max_ntree){
			_fin_auc=cal_rf_auc(_fin_auc_var, _fin_auc_P);
			break;
		}
		i++;

		
		//cout<<i<<"\t";
	}

}

int TAMW::burninDepth()
{
	gfun::printLOG("Choosing Tree depth...");
	vector<double> depthrcd;
	int i=0;
	for(i=0; i<par::ntree_burnin; i++){
		depthrcd.push_back(double(step_tmp_tree()));
	}
	int depth=int(Stat_fuc::mean(depthrcd));
	if(depth==0) depth=1;
	if(depth>par::tree_depth) depth=par::tree_depth;

	gfun::printLOG(int2str(depth)+"\n");

	return depth;
}

double TAMW::cal_rf_auc()
{
	vector<double> ave_indi_LR; vector<bool> if_disease;
	double tmp_LR=0;

	for(int i=0; i<_indi_LRs.size(); i++){
		if(_indi_LRs[i].size()>0){
			tmp_LR=0;
			tmp_LR=Stat_fuc::mean(_indi_LRs[i]);
			ave_indi_LR.push_back(tmp_LR);
			if_disease.push_back(_data->bPheno(i));
		}
	}

	double auc=0;
	auc=Stat_fuc::auc_frm_LR(ave_indi_LR, if_disease);
	return auc;
}

double TAMW::cal_rf_auc(double & var, double & P)
{
	var=0; P=1;

	vector<double> ave_indi_LR; vector<bool> if_disease;
	double tmp_LR=0;

	for(int i=0; i<_indi_LRs.size(); i++){
		if(_indi_LRs[i].size()>0){
			tmp_LR=0;
			tmp_LR=Stat_fuc::mean(_indi_LRs[i]);
			ave_indi_LR.push_back(tmp_LR);
			if_disease.push_back(_data->bPheno(i));
		}
	}

	double auc=0;
	auc=Stat_fuc::auc_frm_LR(ave_indi_LR, if_disease, var);

	P=Stat_fuc::std_norm_p1((auc-0.5)/sqrt(var));
	return auc;
}

void TAMW::wtResult(string file)
{
	gfun::printLOG("Writing TAMW Result (by subjects) at [ " + file + " ]\n");

	int i, j;
	ofstream of;
	of.open(file.c_str());

	of<<"\n#AUC:\t"<<_fin_auc;
	of<<"\n#Var:\t"<<_fin_auc_var;
	of<<"\n#P-value:\t"<<_fin_auc_P;

	of<<"\n#LR for each subject (" + par::missing_str + " for missing)\n";

	vector<double> ave_LR; vector<double> emp95U, emp95L, empmed;
	vector<bool> ave_LR_mk;
	ave_LR.resize(_indi_LRs.size()); emp95U.resize(ave_LR.size()); emp95L.resize(ave_LR.size()); empmed.resize(ave_LR.size());
	ave_LR_mk.resize(ave_LR.size(),false);

	for(i=0; i<_indi_LRs.size(); i++){
		if(_indi_LRs[i].size()>0){
			double tmpU, tmpL;
			Stat_fuc::emprical_CI(_indi_LRs[i],tmpU,tmpL,0.95);
			ave_LR[i]=Stat_fuc::mean(_indi_LRs[i]);
			empmed[i]=Stat_fuc::median(_indi_LRs[i]);
			emp95U[i]=tmpU; emp95L[i]=tmpL;
			ave_LR_mk[i]=true;
		}
	}

//	of<<"\n#FID\tIID\tPAT\tMAT\tGENDER\tPHENO\t"<<"LRmean\tLRmedian\tLR95L\tLR95U"<<"\n";
	if(par::vote){
		of<<"\n#FID\tIID\tPAT\tMAT\tGENDER\tPHENO\t"<<"Votescore"<<"\n";
	}else{
		of<<"\n#FID\tIID\tPAT\tMAT\tGENDER\tPHENO\t"<<"LRmean"<<"\n";
	}
	

	for (i=0;i<ave_LR_mk.size();i++)
	{
		Individual * person = _data->getIndividual(i);
		of << person->fid << "\t"
			<< person->iid << "\t"
			<< person->pat << "\t"
			<< person->mat << "\t"
			<< person->sexcode << "\t"
			<< person->pheno_str << "\t";

		if(ave_LR_mk[i]){
			of<<ave_LR[i]<<"\n";
				//<<empmed[i]<<"\t"
				//<<emp95L[i]<<"\t"
				//<<emp95U[i]<<"\n";
		}else{
			of<<par::missing_str<<"\n";
				//<<par::missing_str<<"\t"
				//<<par::missing_str<<"\t"
				//<<par::missing_str<<"\n";
		}
	}


}


void TAMW::wtPredictorStat(string file)
{
	gfun::printLOG("Writing TAMW Result (by SNPs) at [ " + file + " ]\n");

	int i=0,j=0;

	vector<double> ave_Z, median_Z, Z95L, Z95U; vector<int> sel;
	ave_Z.resize(_snp_stat.size()); median_Z.resize(ave_Z.size()); Z95L.resize(ave_Z.size()); Z95U.resize(ave_Z.size());
	sel.resize(ave_Z.size(),0);

	for(i=0; i<_snp_stat.size(); i++){
		if(_snp_stat[i].size()!=0){
			double tmpU, tmpL;
			Stat_fuc::emprical_CI(_snp_stat[i],tmpU,tmpL,0.95);
			ave_Z[i]=Stat_fuc::mean(_snp_stat[i]);
			median_Z[i]=Stat_fuc::median(_snp_stat[i]);
			Z95L[i]=tmpL; Z95U[i]=tmpU;
			sel[i]=_snp_stat[i].size();
		}
	}


	vector<int> sel_idx;
	for(i=0; i<sel.size(); i++){
		if(sel[i]>par::thrh_seltimes && ave_Z[i]>par::thrh_Zscore){
			sel_idx.push_back(i);
		}
	}

	ofstream ff(file.c_str());

	ff<<"#Only show Predictors with selection times > "<<par::thrh_seltimes
		<<" and average Z score > "<< par::thrh_Zscore
		<<"\n";

	ff<<"#Chr\tName\tPos\tBp\tAllel1\tAllel2\t"<<"meanZ\tmedianZ\tZ95L\tZ95U\tnsel"<<"\n";

	Locus * loc=0;
	for(i=0; i<sel_idx.size(); i++){
		int snp=sel_idx[i];

		loc=_data->getLocus(sel_idx[i]);

		ff	<<loc->chr << "\t"
			<< loc->name << "\t"
			<< loc->pos << "\t"
			<< loc->bp << "\t"
			<< loc->allele1 <<"\t"
			<< loc->allele2 <<"\t";//snp information

		ff<< ave_Z[snp]<<"\t"
			<<median_Z[snp]<<"\t"
			<<Z95L[snp]<<"\t"
			<<Z95U[snp]<<"\t"
			<<sel[snp]<<"\n";

	}

	ff.close();
}

void TAMW::wtConverge(string file)
{
	ofstream ff(file.c_str());
	ff<<"#tree_idx\tAUC\n";
	for(int i=0; i<_n_tree.size(); i++){
		ff<<_n_tree[i]<<"\t"
			<<_moniter_auc[i]<<"\n";
	}
}

void TAMW::apply_model(string model_file)
{
	int i, j;
	ifstream rd(model_file.c_str());

	string buf_str;
	vector<string> vec_str;

	int tree_count=-1;
	int distance=0;
	bool get_tree=false;

	while (getline(rd,buf_str))
	{
		if(buf_str.size()==0) continue;
		else if(buf_str[0]=='#') continue;

		if(gfun::split_string(buf_str,vec_str)==0) continue;

		get_tree=false;

		if(vec_str[0]=="tree_indx"){
			get_tree=rd_treemodel(rd);
			if(par::match_name) match_nameAllel();
		}

		if(get_tree){

			LMW apl_tree;
			apl_tree.initialize(_data);
			vector<double> temp_indi_LR;
			vector<double> temp_stat;

			temp_stat=apl_tree.apply_test_LR(_temp_sel_genotype, _temp_sel_snp, _temp_noexist, _temp_LR_value, temp_indi_LR);

			for(i=0; i<_temp_sel_snp.size(); i++){
				_snp_stat[_temp_sel_snp[i]].push_back(temp_stat[i]);

			}
			for(i=0; i<temp_indi_LR.size(); i++){
				_indi_LRs[i].push_back(temp_indi_LR[i]);
			}


			tree_count++;
			if(distance%par::between_ntree==0){
				double auc=cal_rf_auc();
				_moniter_auc.push_back(auc);
				_n_tree.push_back(tree_count);
			}
			distance++;
			
			//cout<<tree_count<<"\t";
		}
	}

	_fin_auc=cal_rf_auc(_fin_auc_var, _fin_auc_P);
}

bool TAMW::rd_treemodel(ifstream & rd)
{
	int i, j;
	string buf_str;
	vector<string> vec_str;
	int n_str;
	bool finish_tree=false;

	while (getline(rd,buf_str))
	{
		if(buf_str.size()==0) continue;
		else if(buf_str[0]=='#') continue;

		if(gfun::split_string(buf_str,vec_str)==0) continue;

		if(vec_str[0]=="tree_end") { finish_tree=true; break; }

		if(vec_str[0]=="snpname:"){
			_temp_sel_names.clear();
			for(i=1; i<vec_str.size(); i++){
				_temp_sel_names.push_back(vec_str[i]);
			}
		}
		if(vec_str[0]=="allels:"){
			_temp_sel_allels.clear();
			for(i=1; i<vec_str.size(); i++){
				_temp_sel_allels.push_back(vec_str[i]);
			}
		}
		if(vec_str[0]=="col_indx:"){
			_temp_sel_snp.clear();
			for(i=1; i<vec_str.size(); i++){
				_temp_sel_snp.push_back(atoi(vec_str[i].c_str()));
			}
		}
		if(vec_str[0]=="geno_type:"){
			_temp_sel_genotype.clear();
			for(i=1; i<vec_str.size(); i++){
				_temp_sel_genotype.push_back(atoi(vec_str[i].c_str()));
			}
		}
		if(vec_str[0]=="LR_value:"){
			_temp_LR_value.clear();
			_temp_LR_value.resize(_temp_sel_snp.size());

			for(i=0; i<_temp_sel_snp.size(); i++){
				getline(rd, buf_str);
				gfun::split_string(buf_str,vec_str);
				for(j=1; j<vec_str.size(); j++){
					_temp_LR_value[i].push_back(atof(vec_str[j].c_str()));
				}
			}
		}
		if(vec_str[0]=="no_exist:"){
			_temp_noexist.clear();
			_temp_noexist.resize(_temp_sel_snp.size());

			for(i=0; i<_temp_sel_snp.size(); i++){
				getline(rd, buf_str);
				gfun::split_string(buf_str,vec_str);
				for(j=1; j<vec_str.size(); j++){
					_temp_noexist[i].push_back(atoi(vec_str[j].c_str()));
				}
			}

		}
	}

	return finish_tree;
}

void TAMW::match_nameAllel()
{
	_temp_sel_snp.clear(); _temp_sel_genotype.clear();
	_data->matchNameAllel(_temp_sel_names,_temp_sel_allels,_temp_sel_snp,_temp_sel_genotype);

	for(int i=0; i<_temp_sel_names.size(); i++){
		if(_temp_sel_snp[i]<0) gfun::error("Can't find snp:" + _temp_sel_names[i] + " in the file");
		if(_temp_sel_genotype[i]<0) gfun::error("Can't find corresponding allels at snp:"+ _temp_sel_names[i] + " in the file");
	}

}