#include "lmw.h"

FLMW::FLMW()
{

}

FLMW::~FLMW()
{

}

void FLMW::initialize(snpdt * snp_data)
{
	LMW::initialize(snp_data);
	gen_contrust();
}

void FLMW::gen_contrust()
{
	vector<int> family_id, family;
	vector< vector<int> > grouped_indi_id;
	int i=0, j=0, n_mk=0;
	bool ex_mk=false;

	for(i=0; i<_subj_num; i++){
		family.push_back(_datafile->fID(i));
	}


	for(i=0; i<_subj_num; i++){
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

	_compare_struct.clear();

	for(i=0; i<grouped_indi_id.size(); i++){
		vector<int> t1, t2;
		for(j=0; j<grouped_indi_id[i].size(); j++){
			if(_aff[grouped_indi_id[i][j]]){
				t1.push_back(grouped_indi_id[i][j]);
			}else{
				t2.push_back(grouped_indi_id[i][j]);
			}
		}
		_compare_struct.push_back(t1);
		_compare_struct.push_back(t2);
	}

}

void FLMW::temp_deviding_mis(int col, int genotype)
{
	_tmp_LR_nom=_LR_nom;

	vector<int> p_temp_main, p_temp_other, h_temp_main, h_temp_other;

	int i=0,j=0,indx=0;
	for(i=0; i<_devided_patient.size(); i++){

		for(j=0; j<_devided_patient[i].size(); j++){

			indx=_devided_patient[i][j];
			int geno=_datafile->genotypeToInt(indx, col);
			if(geno==par::missing_int){
				continue;
			}else if (geno==genotype)
			{
				p_temp_main.push_back(indx);
			}else
			{
				p_temp_other.push_back(indx);
			}
		}

		for(j=0; j<_devided_health[i].size(); j++){
			indx=_devided_health[i][j];
			int geno=_datafile->genotypeToInt(indx, col);
			if(geno==-9){
				continue;
			}else if (geno==genotype)
			{
				h_temp_main.push_back(indx);
			}else
			{
				h_temp_other.push_back(indx);
			}

		}

		if(p_temp_main.size()!=0 || h_temp_main.size()!=0){

			double tmp_np=double(p_temp_main.size());
			double tmp_nh=double(h_temp_main.size());

			double temp;
			tmp_nh= (tmp_nh>0) ? tmp_nh : par::replace0_cor;
			tmp_np= (tmp_np>0) ? tmp_np : par::replace0_cor;
			temp=tmp_np/tmp_nh;

			temp=temp*_n_health/_n_disease;

			for(j=0; j<p_temp_main.size(); j++) _tmp_LR_nom[p_temp_main[j]]=temp;
			for(j=0; j<h_temp_main.size(); j++) _tmp_LR_nom[h_temp_main[j]]=temp;

		}

		if(p_temp_other.size()!=0 || h_temp_other.size()!=0){

			double tmp_np=double(p_temp_other.size());
			double tmp_nh=double(h_temp_other.size());

			double temp;
			tmp_nh= (tmp_nh>0) ? tmp_nh : par::replace0_cor;
			tmp_np= (tmp_np>0) ? tmp_np : par::replace0_cor;
			temp=tmp_np/tmp_nh;

			temp=temp*_n_health/_n_disease;
			for(j=0; j<p_temp_other.size(); j++) _tmp_LR_nom[p_temp_other[j]]=temp;
			for(j=0; j<h_temp_other.size(); j++) _tmp_LR_nom[h_temp_other[j]]=temp;
		}

		p_temp_main.clear();
		p_temp_other.clear();
		h_temp_main.clear();
		h_temp_other.clear();
	}

}




double FLMW::auc_fm_tmpLR()
{
	return Stat_fuc::auc_frm_LR(_tmp_LR_nom,_compare_struct);
}

double FLMW::auc_fm_LR()
{
	return Stat_fuc::auc_frm_LR(_LR_nom,_compare_struct);
}

double FLMW::auc_fm_LR(double & var)
{
	return Stat_fuc::auc_frm_LR(_LR_nom,_compare_struct,var);
}

void FLMW::auc_VC_fm_LR(vector<double> & pre_score, vector<double> & score, int step)
{
	double auc=0, var=0, cov=0;
	pre_score=score;
	Stat_fuc::score_frm_LR(_LR_nom, _compare_struct, score);	

	int i=0;
	for(i=0; i<score.size(); i++) auc+=score[i];
	auc/=double(score.size());

	for(i=0; i<score.size(); i++) var+=score[i]*score[i];
	var=var/double(score.size()) - auc*auc;

	if(step==0){
		cov=0;
	}else{
		for(i=0; i<score.size(); i++) cov+=score[i]*pre_score[i];
		cov= cov/double(score.size()) -auc*_vec_auc[step-1];
	}

	_vec_auc.push_back(auc); _var_auc.push_back(var); _cov_auc.push_back(cov);
}