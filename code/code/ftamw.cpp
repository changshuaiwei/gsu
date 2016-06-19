#include "tamw.h"

void FTAMW::initialize(snpdt * data)
{
	_data=data;
	_indi_LRs.resize(_data->totalIndi());
	_snp_stat.resize(_data->totalLocus());

	set_n_classifier();
	vector<int> subjIdx;
	for(int i=0; i<_data->totalIndi(); i++) subjIdx.push_back(i);
	getFamGrp(subjIdx, _fam_group,_compare_struct);

}

void FTAMW::getFamGrp(vector<int> subjIdx, vector< vector<int> > & grouped_indi_id, vector< vector<int> > & compare_struct)
{
	grouped_indi_id.clear();

	int subj_num=subjIdx.size();

	vector<bool> aff; aff.resize(subj_num);
	for(int i=0; i<aff.size(); i++) {
		int afft=_data->bPheno(subjIdx[i]);
		if(afft==1) aff[i]=true; else if(afft==0) aff[i]=false;
		else{
			gfun::error("do not allow missing phenotype now");
		}
	}

	vector<int> family_id, family;
	int i=0, j=0, n_mk=0;
	bool ex_mk=false;

	for(i=0; i<subj_num; i++){
		family.push_back(_data->fID(subjIdx[i]));
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

	compare_struct.clear();

	for(i=0; i<grouped_indi_id.size(); i++){
		vector<int> t1, t2;
		for(j=0; j<grouped_indi_id[i].size(); j++){
			if(aff[grouped_indi_id[i][j]]){
				t1.push_back(grouped_indi_id[i][j]);
			}else{
				t2.push_back(grouped_indi_id[i][j]);
			}
		}
		compare_struct.push_back(t1);
		compare_struct.push_back(t2);
	}

}

void FTAMW::step_tree(int tree_idx)
{
	int i=0, j=0, k=0;
	produce_sequence(_data->totalLocus(),par::seed);

	//Stat_fuc::btstrp_sampling(_data->totalIndi(), _tmp_btstping, par::seed);
	Stat_fuc::btstrp_sampling(_fam_group.size(), _tmp_btstping, par::seed);

	_tmp_bts.clear(); _tmp_oob.clear();
	for(i=0; i<_tmp_btstping.size(); i++){
		if(_tmp_btstping[i]==0) {
			for(j=0; j<_fam_group[i].size(); j++) _tmp_oob.push_back(_fam_group[i][j]);
		}
		else{
			int tmp_times=_tmp_btstping[i];
			for(j=0; j<tmp_times; j++){
				for(k=0; k<_fam_group[i].size(); k++) _tmp_bts.push_back(_fam_group[i][k]);
			}
		}
	}

	_data->subset(_tmp_bts,_snp_map);

	FLMW cal_tree;
	cal_tree.initialize(_data);
	cal_tree.grow_tree(_iter_seq);
	cal_tree.send_res(_temp_sel_genotype, _temp_sel_snp, _temp_vec_auc);

	FLMW cal_LR;
	cal_LR.initialize(_data);
	cal_LR.give_LR(_temp_sel_genotype, _temp_sel_snp,_temp_noexist, _temp_LR_value);

	_data->reset();
	_data->subset(_tmp_oob,_snp_map);

	if(par::auto_tree_depth){
		FLMW apl_tree_tres;
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


	FLMW apl_tree;
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

int FTAMW::step_tmp_tree()
{
	int i=0, j=0, k=0;
	produce_sequence(_data->totalLocus(),par::seed);

	//Stat_fuc::btstrp_sampling(_data->totalIndi(), _tmp_btstping, par::seed);
	Stat_fuc::btstrp_sampling(_fam_group.size(), _tmp_btstping, par::seed);

	_tmp_bts.clear(); _tmp_oob.clear();
	for(i=0; i<_tmp_btstping.size(); i++){
		if(_tmp_btstping[i]==0) {
			for(j=0; j<_fam_group[i].size(); j++) _tmp_oob.push_back(_fam_group[i][j]);
		}
		else{
			int tmp_times=_tmp_btstping[i];
			for(j=0; j<tmp_times; j++){
				for(k=0; k<_fam_group[i].size(); k++) _tmp_bts.push_back(_fam_group[i][k]);
			}
		}
	}

	_data->subset(_tmp_bts,_snp_map);

	FLMW cal_tree;
	cal_tree.initialize(_data);
	cal_tree.grow_tree(_iter_seq);
	cal_tree.send_res(_temp_sel_genotype, _temp_sel_snp, _temp_vec_auc);

	FLMW cal_LR;
	cal_LR.initialize(_data);
	cal_LR.give_LR(_temp_sel_genotype, _temp_sel_snp,_temp_noexist, _temp_LR_value);

	_data->reset();
	_data->subset(_tmp_oob,_snp_map);

	FLMW apl_tree_tres;
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

double FTAMW::cal_rf_auc()
{
	vector<double> ave_indi_LR; vector<int> subjIdx;
	vector< vector<int> > fam_grp, compare_struct;
	double tmp_LR=0;

	for(int i=0; i<_indi_LRs.size(); i++){
		if(_indi_LRs[i].size()>0){
			tmp_LR=0;
			tmp_LR=Stat_fuc::mean(_indi_LRs[i]);
			ave_indi_LR.push_back(tmp_LR);
			subjIdx.push_back(i);
		}
	}

	double auc=0;
	if(ave_indi_LR.size()==_data->totalIndi()){
		auc=Stat_fuc::auc_frm_LR(ave_indi_LR, _compare_struct);
	}else{
		getFamGrp(subjIdx, fam_grp, compare_struct);
		auc=Stat_fuc::auc_frm_LR(ave_indi_LR, compare_struct);
	}
	
	return auc;
}

double FTAMW::cal_rf_auc(double & var, double & P)
{
	var=0; P=1;

	vector<double> ave_indi_LR; vector<int> subjIdx;
	vector< vector<int> > fam_grp, compare_struct;
	double tmp_LR=0;

	for(int i=0; i<_indi_LRs.size(); i++){
		if(_indi_LRs[i].size()>0){
			tmp_LR=0;
			tmp_LR=Stat_fuc::mean(_indi_LRs[i]);
			ave_indi_LR.push_back(tmp_LR);
			subjIdx.push_back(i);
		}
	}

	double auc=0;
	if(ave_indi_LR.size()==_data->totalIndi()){
		auc=Stat_fuc::auc_frm_LR(ave_indi_LR, _compare_struct, var);
	}else{
		getFamGrp(subjIdx, fam_grp, compare_struct);
		auc=Stat_fuc::auc_frm_LR(ave_indi_LR, compare_struct, var);
	}

	P=Stat_fuc::std_norm_p1((auc-0.5)/sqrt(var));
	return auc;
}

void FTAMW::apply_model(string model_file)
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

			FLMW apl_tree;
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

			cout<<tree_count<<"\t";
		}
	}

	_fin_auc=cal_rf_auc(_fin_auc_var, _fin_auc_P);
}
