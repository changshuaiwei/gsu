#include "lmw.h"

LMW::LMW()
{
	_datafile=0;
}

LMW::~LMW()
{
	_datafile=0;
}

void LMW::initialize(snpdt * snp_data)
{
	_datafile=snp_data;

	_subj_num=_datafile->totalIndi();
	_marker_num=snp_data->totalLocus();

	_aff.clear(); _aff.resize(_subj_num);
	for(int i=0; i<_aff.size(); i++) {
		int aff=_datafile->bPheno(i);
		if(aff==1) _aff[i]=true; else if(aff==0) _aff[i]=false;
		else{
			gfun::error("LMW do not allow missing phenotype now");
		}
	}

	vector <bool> temp; temp.push_back(true); temp.push_back(true); temp.push_back(true);
	for(int i=0; i<=_marker_num; i++){
		_left_genotype.push_back(temp);
		_left_genotype_num.push_back(3);
	}
	_n_health=double(_datafile->nHealth());
	_n_disease=double(_datafile->nDisease());

}

void LMW::iteration()
{

	int sel_auc_num;
	double last_auc_value=0.0;
	do 
	{
		sel_auc_num=_sel_snp.size();
		if(sel_auc_num!=0) last_auc_value=_vec_auc[sel_auc_num-1];

		if(sel_auc_num>=par::most_nsnp){
			break;
		}
		if ((par::largest_auc-last_auc_value)<1e-8 )
		{
			break;
		}
		if (last_auc_value<=0.0 && sel_auc_num!=0)
		{
			_vec_auc[sel_auc_num-1]=0.0;
			break;
		}
		
		if(par::show_iteration){
			cout<<"Searching order "<<sel_auc_num+1;
			cout.flush();
		}

		step_auc();

		if(par::show_iteration){
			cout<<"\r";
			cout.flush();
		}

		if(last_auc_value>=_vec_auc[_sel_snp.size()-1]){
			_vec_auc.pop_back();
			_sel_snp.pop_back();
			_sel_genotype.pop_back();
			break;
		}

		

	} while (1);
}

void LMW::step_auc()
{
	patient_devide();
	int i=0;
	for(i=0; i<_marker_num; i++){

		col_auc(i);
	}
	ranking_auc();// ranking auc and input the result to the "selected record"
	_temp_sel_auc.clear();
	_temp_sel_genotype.clear();
}

void LMW::patient_devide(bool apply)
{
	int sel_auc_num=_sel_snp.size();

	int n_devided=_devided_patient.size();
	if(n_devided!=_devided_health.size()) 
		gfun::error("\ndeviding error in the "+ int2str(sel_auc_num) + "step LMW::patient_devide!!\n");

	int i=0;
	if(sel_auc_num==0) {

		_devided_health.resize(1);
		_devided_patient.resize(1);


		for (i=0; i<_subj_num; i++)
		{
			int pheno=_datafile->bPheno(i);
			if(pheno==1){
				_devided_patient[0].push_back(i);
			}else if(pheno==0){
				_devided_health[0].push_back(i);
			}
		}
		
		double initial_LR=1;
		_cur_tmp_LR.push_back(initial_LR);

		_LR_nom.resize(_subj_num, initial_LR);



	}else {
		_pre_devided_health=_devided_health;
		_pre_devided_patient=_devided_patient;
		_devided_patient.clear();
		_devided_health.clear();

		deviding_to_next_mis(apply);

	}

}


void LMW::col_auc(int col)
{
	int i=0;

	vector<bool> left_geno=_left_genotype[col];

	if(par::exclude_1){
		left_geno[1]=false;
	}
	
	int aim_geno=0;
	int result_geno=0;
	double max_auc=0, temp=0;

	i=0;

	for(i=0; i<left_geno.size(); i++){
		if(left_geno[i]){
			aim_geno=i;
			temp_deviding_mis(col, aim_geno);
			temp=auc_fm_tmpLR();
			_temp_devided_np.clear();
			_temp_devided_nh.clear();
			_temp_LR.clear();

			if(max_auc<temp){
				max_auc=temp;
				result_geno=aim_geno;
			}
		}		
	}

	_temp_sel_auc.push_back(max_auc);
	_temp_sel_genotype.push_back(result_geno);

}

double LMW::auc_fm_tmpLR()
{
	return Stat_fuc::auc_frm_LR(_temp_LR, _temp_devided_np, _temp_devided_nh);
}

void LMW::temp_deviding_mis(int col, int genotype)
{
	bool disregard_missing=false;
	if(_devided_patient.size()==1) disregard_missing=true;

	_temp_devided_np.clear();
	_temp_devided_nh.clear();
	_temp_LR.clear();

	int p_temp_main=0;
	int p_temp_other=0;
	int h_temp_main=0;
	int h_temp_other=0;

	int p_temp_missing=0;
	int h_temp_missing=0;

	int i=0,j=0,indx=0;
	for(i=0; i<_devided_patient.size(); i++){

		for(j=0; j<_devided_patient[i].size(); j++){

			indx=_devided_patient[i][j];
			int geno=_datafile->genotypeToInt(indx, col);
			if(geno==par::missing_int){
				p_temp_missing++;
				continue;
			}else if (geno==genotype)
			{
				p_temp_main++;
			}else
			{
				p_temp_other++;
			}
		}

		for(j=0; j<_devided_health[i].size(); j++){
			indx=_devided_health[i][j];
			int geno=_datafile->genotypeToInt(indx, col);
			if(geno==-9){
				h_temp_missing++;
				continue;
			}else if (geno==genotype)
			{
				h_temp_main++;
			}else
			{
				h_temp_other++;
			}

		}

		if(p_temp_main!=0 || h_temp_main!=0){
			_temp_devided_np.push_back(p_temp_main);
			_temp_devided_nh.push_back(h_temp_main);

			double tmp_np=double(p_temp_main);
			double tmp_nh=double(h_temp_main);

			double temp;
			tmp_nh= (tmp_nh>0) ? tmp_nh : par::replace0_cor;
			tmp_np= (tmp_np>0) ? tmp_np : par::replace0_cor;
			temp=tmp_np/tmp_nh;
		

			temp=temp*_n_health/_n_disease;
			_temp_LR.push_back(temp);
		}

		if(p_temp_other!=0 || h_temp_other!=0){
			_temp_devided_np.push_back(p_temp_other);
			_temp_devided_nh.push_back(h_temp_other);

			double tmp_np=double(p_temp_other);
			double tmp_nh=double(h_temp_other);

			double temp;
			tmp_nh= (tmp_nh>0) ? tmp_nh : par::replace0_cor;
			tmp_np= (tmp_np>0) ? tmp_np : par::replace0_cor;
			temp=tmp_np/tmp_nh;
			
			temp=temp*_n_health/_n_disease;
			_temp_LR.push_back(temp);
		}

		if(!disregard_missing){
			if(p_temp_missing!=0 || h_temp_missing!=0){
				_temp_devided_np.push_back(p_temp_missing);
				_temp_devided_nh.push_back(h_temp_missing);
				double temp;
				temp=_cur_tmp_LR[i];
				_temp_LR.push_back(temp);
			}
		}

		p_temp_main=0;
		p_temp_other=0;
		h_temp_main=0;
		h_temp_other=0;

		p_temp_missing=0;
		h_temp_missing=0;
	}

	for(i=0; i<_mis_LR.size(); i++){
		_temp_devided_np.push_back(_mis_patient[i].size());
		_temp_devided_nh.push_back(_mis_health[i].size());
		_temp_LR.push_back(_mis_LR[i]);
	}
}



void LMW::ranking_auc()
{
	int i=0;
	int m_idx=0;
	
	m_idx=Stat_fuc::max_index(_temp_sel_auc);

	int geno=_temp_sel_genotype[m_idx];
	double auc=_temp_sel_auc[m_idx];

	_sel_snp.push_back(m_idx);
	_sel_genotype.push_back(geno);
	_vec_auc.push_back(auc);

	_left_genotype[m_idx][geno]=false;
	_left_genotype_num[m_idx]--;

}


void LMW::deviding_to_next_mis(bool apply/* =false */)
{

	bool disregard_missing=false;
	if(_pre_devided_patient.size()<=1) disregard_missing=true;

	_pre_tmp_LR=_cur_tmp_LR;
	_cur_tmp_LR.clear();

	int pre_size=_pre_devided_patient.size();
	if(pre_size!=_pre_devided_health.size()) gfun::error("\nprehealth should equal to prepatient!!\n");

	int sel_auc_num=_sel_snp.size();
	int aim_col=_sel_snp[sel_auc_num-1];
	int aim_genotype=_sel_genotype[sel_auc_num-1];

	vector<int> p_temp_main;
	vector<int> p_temp_other;
	vector<int> h_temp_main;
	vector<int> h_temp_other;

	vector<int> p_temp_missing;
	vector<int> h_temp_missing;


	if(apply) { _temp_noexist_rcd.clear(); _temp_LR.clear(); }

	int i=0,j=0, indx=0;
	for(i=0; i<pre_size; i++){

		for(j=0; j<_pre_devided_patient[i].size(); j++){
			indx=_pre_devided_patient[i][j];
			int geno=_datafile->genotypeToInt(indx,aim_col);
				if(geno==par::missing_int){
				p_temp_missing.push_back(indx);
				continue;
			}else if(geno==aim_genotype){
				p_temp_main.push_back(indx);
			}else{
				p_temp_other.push_back(indx);
			}
		}

		for(j=0; j<_pre_devided_health[i].size(); j++){
			indx=_pre_devided_health[i][j];
			int geno=_datafile->genotypeToInt(indx,aim_col);
			if(geno==par::missing_int){
				h_temp_missing.push_back(indx);
				continue;
			}else if(geno==aim_genotype){
				h_temp_main.push_back(indx);
			}else{
				h_temp_other.push_back(indx);
			}
		}

		if(apply){

			if(p_temp_main.size()!=0 || h_temp_main.size()!=0){
				_devided_patient.push_back(p_temp_main);
				_devided_health.push_back(h_temp_main);

				double tmp_np=double(p_temp_main.size());
				double tmp_nh=double(h_temp_main.size());

				double temp;
				
				tmp_nh= (tmp_nh>0) ? tmp_nh : par::replace0_cor;
				tmp_np= (tmp_np>0) ? tmp_np : par::replace0_cor;

				temp=tmp_np/tmp_nh;

				temp=temp*_n_health/_n_disease;
				_temp_LR.push_back(temp);

				for(j=0; j<p_temp_main.size(); j++) _LR_nom[p_temp_main[j]]=temp;
				for(j=0; j<h_temp_main.size(); j++) _LR_nom[h_temp_main[j]]=temp;

			}else{
				_temp_noexist_rcd.push_back(i*2);
			}

			if(p_temp_other.size()!=0 || h_temp_other.size()!=0){
				_devided_patient.push_back(p_temp_other);
				_devided_health.push_back(h_temp_other);

				double tmp_np=double(p_temp_other.size());
				double tmp_nh=double(h_temp_other.size());

				double temp;

				tmp_nh= (tmp_nh>0) ? tmp_nh : par::replace0_cor;
				tmp_np= (tmp_np>0) ? tmp_np : par::replace0_cor;
				temp=tmp_np/tmp_nh;

				temp=temp*_n_health/_n_disease;
				_temp_LR.push_back(temp);
				for(j=0; j<p_temp_other.size(); j++) _LR_nom[p_temp_other[j]]=temp;
				for(j=0; j<h_temp_other.size(); j++) _LR_nom[h_temp_other[j]]=temp;
			}else{
				_temp_noexist_rcd.push_back(i*2+1);
			}

		}else{
			if(p_temp_main.size()!=0 || h_temp_main.size()!=0){
				_devided_patient.push_back(p_temp_main);
				_devided_health.push_back(h_temp_main);

				double tmp_np=double(p_temp_main.size());
				double tmp_nh=double(h_temp_main.size());

				double temp=0;
				tmp_nh= (tmp_nh>0) ? tmp_nh : par::replace0_cor;
				tmp_np= (tmp_np>0) ? tmp_np : par::replace0_cor;
				temp=tmp_np/tmp_nh;

				temp=temp*_n_health/_n_disease;
				_cur_tmp_LR.push_back(temp);
				for(j=0; j<p_temp_main.size(); j++) _LR_nom[p_temp_main[j]]=temp;
				for(j=0; j<h_temp_main.size(); j++) _LR_nom[h_temp_main[j]]=temp;
			}

			if(p_temp_other.size()!=0 || h_temp_other.size()!=0){
				_devided_patient.push_back(p_temp_other);
				_devided_health.push_back(h_temp_other);

				double tmp_np=double(p_temp_other.size());
				double tmp_nh=double(h_temp_other.size());

				double temp=0;
				tmp_nh= (tmp_nh>0) ? tmp_nh : par::replace0_cor;
				tmp_np= (tmp_np>0) ? tmp_np : par::replace0_cor;
				temp=tmp_np/tmp_nh;

				temp=temp*_n_health/_n_disease;
				_cur_tmp_LR.push_back(temp);
				for(j=0; j<p_temp_other.size(); j++) _LR_nom[p_temp_other[j]]=temp;
				for(j=0; j<h_temp_other.size(); j++) _LR_nom[h_temp_other[j]]=temp;
			}

			if(!disregard_missing){
				if(p_temp_missing.size()!=0 || h_temp_missing.size()!=0){
					_mis_patient.push_back(p_temp_missing);
					_mis_health.push_back(h_temp_missing);
					double temp;
					temp=_pre_tmp_LR[i];
					
					_mis_LR.push_back(temp);
				}
			}

		}
		p_temp_main.clear();
		p_temp_other.clear();
		h_temp_main.clear();
		h_temp_other.clear();
		p_temp_missing.clear();
		h_temp_missing.clear();
	}
}

void LMW::out_nom_lr(string out)
{

	gfun::printLOG("Writing LR (by subjects) at [ " + out + " ]\n");

	int i, j;
	ofstream of;
	of.open(out.c_str());

	//of<<"\n#AUC:\t"<<_fin_auc;
	//of<<"\n#Var:\t"<<_fin_auc_var;
	//of<<"\n#P-value:\t"<<_fin_auc_P;

	of<<"\n#LR for each subject (" + par::missing_str + " for missing)\n";

	of<<"\n#FID\tIID\tPAT\tMAT\tGENDER\tPHENO\t"<<"LR"<<"\n";



	for (i=0;i<_LR_nom.size();i++)
	{
		Individual * person = _datafile->getIndividual(i);
		of << person->fid << "\t"
			<< person->iid << "\t"
			<< person->pat << "\t"
			<< person->mat << "\t"
			<< person->sexcode << "\t"
			<< person->pheno_str << "\t";

		of<<_LR_nom[i]<<"\n";
	}


}
void LMW::wtResult(string outputfile)
{
	gfun::printLOG("Writing Forward Selection Result at [ " + outputfile + " ]\n");

	ofstream result(outputfile.c_str());
	if(!result) gfun::error("\ncould not open the result file\n");
	result<<"#result of association scanning (0=A1A1, 1=A1A2, 2=A2A2)\n";
	result<<"#order\t"<<"Chr\tName\tPos\tBp\tAllel1\tAllel2\t"<<"Geno\t"<<"AUC\n";

	Locus * loc=0;
	for(int i=0; i<_vec_auc.size(); i++){
		loc=_datafile->getLocus(_sel_snp[i]);

		result<<i+1<<"\t";//order

		result	<<loc->chr << "\t"
				<< loc->name << "\t"
				<< loc->pos << "\t"
				<< loc->bp << "\t"
				<< loc->allele1 <<"\t"
				<< loc->allele2 <<"\t";//snp information

		result	<<_sel_genotype[i]<<"\t"
				<<_vec_auc[i]<<"\n";//model information
	}




	result<<"\n\n#LR for each subject (" + par::missing_str + " for missing)\n";

	result<<"\n#FID\tIID\tPAT\tMAT\tGENDER\tPHENO\t"<<"LR"<<"\n";



	for (int i=0;i<_LR_nom.size();i++)
	{
		Individual * person = _datafile->getIndividual(i);
		result << person->fid << "\t"
			<< person->iid << "\t"
			<< person->pat << "\t"
			<< person->mat << "\t"
			<< person->sexcode << "\t"
			<< person->pheno_str << "\t";

		result<<_LR_nom[i]<<"\n";
	}

	result.close();

}

void LMW::wtSubset(string ssfile)
{
	gfun::printLOG("Writing subset data for the selected SNP\n");

	_datafile->subset(_sel_snp);
	_datafile->writeGenMapPed(ssfile);
	_datafile->reset();
}


void LMW::give_LR(vector<int> & sel_genotype, vector<int> & sel_snp, vector< vector<int> > & no_exist, vector< vector<double> > & sel_LR_value)
{
	int i,n_step=sel_snp.size();

	//cout<<"\n\nExtract LR vectors from the training data set...\n\n\n";

	patient_devide();
	no_exist.clear();
	sel_LR_value.clear();

	for (i=0;i<n_step;i++)
	{
		//cout<<"\nextracting LR from the "<<(i+1)<<"th model....\n";

		{
			_sel_snp.push_back(sel_snp[i]);
			_sel_genotype.push_back(sel_genotype[i]);
			_vec_auc.push_back(0);

			bool apply=true;
			patient_devide(apply);
		}

		no_exist.push_back(_temp_noexist_rcd);
		sel_LR_value.push_back(_temp_LR);
	}

}

void LMW::give_nom_LR(vector<double> & LR)
{
	LR=_LR_nom;
}

vector<double> LMW::apply_model(vector<int> & sel_genotype, vector<int> & sel_snp, vector< vector<int> > & no_exist, vector< vector<double> > & sel_LR_value)
{
	_vec_LR_value=sel_LR_value;
	int i,n_step=sel_snp.size();
	//cout<<"\n\nApplying model to the evaluation data set...\n\n\n";

	patient_devide();


	for (i=0;i<n_step;i++)
	{

		//cout<<"\napplying the "<<(i+1)<<"th model....\n";
		//to be added
		_sel_snp.push_back(sel_snp[i]);
		_sel_genotype.push_back(sel_genotype[i]);
		apply_devide_mis(sel_LR_value[i], no_exist[i]);
		double auc=auc_fm_LR();
		_vec_auc.push_back(auc);
	}

	return _vec_auc;
}

double LMW::auc_fm_LR()
{
	return Stat_fuc::auc_frm_LR(_temp_LR, _temp_devided_np, _temp_devided_nh);
}

vector<double> LMW::apply_model(vector<int> & sel_genotype, vector<int> & sel_snp, vector< vector<int> > & no_exist, vector< vector<double> > & sel_LR_value, vector<double> & auc_var, int nc)
{
	auc_var.clear();
	_vec_LR_value=sel_LR_value;
	int i,n_step=sel_snp.size();
	double var;
	//cout<<"\n\nApplying model to the evaluation data set...\n\n\n";
	patient_devide();
	for (i=0;i<nc;i++)
	{

		//cout<<"\napplying the "<<(i+1)<<"th model....\n";

		_sel_snp.push_back(sel_snp[i]);
		_sel_genotype.push_back(sel_genotype[i]);
		apply_devide_mis(sel_LR_value[i], no_exist[i]);
		double auc=auc_fm_LR(var);
		_vec_auc.push_back(auc);

		auc_var.push_back(var);

	}

	return _vec_auc;
}

double LMW::auc_fm_LR(double & var)
{
	return Stat_fuc::auc_frm_LR(_temp_LR, _temp_devided_np, _temp_devided_nh, var);
}


void LMW::apply_devide_mis(vector<double> & LR_value, vector<int> & no_exist)
{
	_temp_devided_np.clear();
	_temp_devided_nh.clear();
	_temp_LR.clear();

	_pre_devided_health=_devided_health;
	_pre_devided_patient=_devided_patient;
	_devided_patient.clear();
	_devided_health.clear();

	bool disregard_missing=false;
	if(_pre_devided_patient.size()<=1) disregard_missing=true;

	int sel_auc_num=_sel_snp.size();
	int aim_col=_sel_snp[sel_auc_num-1];
	int aim_genotype=_sel_genotype[sel_auc_num-1];

	vector<int> p_temp_main;
	vector<int> p_temp_other;
	vector<int> h_temp_main;
	vector<int> h_temp_other;

	vector<int> p_temp_missing;
	vector<int> h_temp_missing;


	int pre_size=_pre_devided_patient.size();
	if(pre_size!=_pre_devided_health.size()) gfun::error("\nprehealth should equal to prepatient!!\n");

	int i=0,j=0, indx=0;
	int ne_idx=0;
	bool test=true;
	if(no_exist.size()==0) test=false;
	for(i=0; i<pre_size; i++){
		for(j=0; j<_pre_devided_patient[i].size(); j++){
			indx=_pre_devided_patient[i][j];
			char geno=_datafile->genotypeToInt(indx,aim_col);
			if(geno==par::missing_int){
				p_temp_missing.push_back(indx);
				continue;
			}else if(geno==aim_genotype){
				p_temp_main.push_back(indx);
			}else{
				p_temp_other.push_back(indx);
			}
		}

		for(j=0; j<_pre_devided_health[i].size(); j++){
			indx=_pre_devided_health[i][j];
			char geno=_datafile->genotypeToInt(indx,aim_col);
			if(geno==par::missing_int){
				h_temp_missing.push_back(indx);
				continue;
			}else if(geno==aim_genotype){
				h_temp_main.push_back(indx);
			}else{
				h_temp_other.push_back(indx);
			}
		}

		if(test && no_exist[ne_idx]==i*2){
			if(p_temp_main.size()!=0 || h_temp_main.size()!=0 ){
				for(j=0; j<p_temp_main.size(); j++) p_temp_missing.push_back(p_temp_main[j]);
				for(j=0; j<h_temp_main.size(); j++) h_temp_missing.push_back(h_temp_main[j]);
			}
			if(ne_idx+1<no_exist.size())ne_idx++;
			else test=false;
		}else{
			_devided_patient.push_back(p_temp_main);
			_devided_health.push_back(h_temp_main);
		}

		if(test && no_exist[ne_idx]==i*2+1){
			if(p_temp_other.size()!=0 || h_temp_other.size()!=0 ){
				for(j=0; j<p_temp_other.size(); j++) p_temp_missing.push_back(p_temp_other[j]);
				for(j=0; j<h_temp_other.size(); j++) h_temp_missing.push_back(h_temp_other[j]);
			}
			if(ne_idx+1<no_exist.size())ne_idx++;
			else test=false;
		}else{
			_devided_patient.push_back(p_temp_other);
			_devided_health.push_back(h_temp_other);
		}

		if(!disregard_missing){
			if(p_temp_missing.size()!=0 || h_temp_missing.size()!=0){
				_mis_patient.push_back(p_temp_missing);
				_mis_health.push_back(h_temp_missing);
				double temp;
				if(sel_auc_num==1){
					double tmp_nh, tmp_np; tmp_np=_pre_devided_patient[i].size(); tmp_nh=_pre_devided_health[i].size();

					tmp_nh= (tmp_nh>0) ? tmp_nh : par::replace0_cor;
					tmp_np= (tmp_np>0) ? tmp_np : par::replace0_cor;
					temp=tmp_np/tmp_nh;
				}else{
					temp=_vec_LR_value[sel_auc_num-2][i];
				}
				_mis_LR.push_back(temp);
			}
		}



		p_temp_main.clear();
		p_temp_other.clear();
		h_temp_main.clear();
		h_temp_other.clear();
		p_temp_missing.clear();
		h_temp_missing.clear();
	}

	//produce np,nh, and LR
	for(i=0; i<_devided_patient.size(); i++){
		_temp_devided_np.push_back(_devided_patient[i].size());
		_temp_devided_nh.push_back(_devided_health[i].size());

		_temp_LR.push_back(LR_value[i]);
	}

	for(i=0; i<_mis_patient.size(); i++){
		_temp_devided_np.push_back(_mis_patient[i].size());
		_temp_devided_nh.push_back(_mis_health[i].size());
		_temp_LR.push_back(_mis_LR[i]);
	}

	//produce LR
	for(i=0; i<_devided_patient.size(); i++){
		for(j=0; j<_devided_patient[i].size(); j++) _LR_nom[_devided_patient[i][j]]=LR_value[i];
		for(j=0; j<_devided_health[i].size(); j++) _LR_nom[_devided_health[i][j]]=LR_value[i];
	}

}

void LMW::send_res(vector<int> & sel_genotye, vector<int> & sel_snp, vector<double> & vec_auc)
{
	sel_genotye=_sel_genotype;
	sel_snp=_sel_snp;
	vec_auc=_vec_auc;
}
void LMW::grow_tree(vector< vector<int> > & seq)
{
	int n_size=_datafile->totalLocus();
	int depth_idx=0;
	int sel_auc_num;
	double last_auc_value=0.0;
	do 
	{
		sel_auc_num=_sel_snp.size();
		if(sel_auc_num!=0) last_auc_value=_vec_auc[sel_auc_num-1];

		if(sel_auc_num>=par::most_nsnp || sel_auc_num>=par::tree_depth){
			break;
		}
		if ((par::largest_auc-last_auc_value)<1e-8 )
		{
			break;
		}
		if (last_auc_value<=0.0 && sel_auc_num!=0)
		{
			_vec_auc[sel_auc_num-1]=0.0;
			break;
		}

		_temp_col_indx.clear(); _temp_col_indx.resize(n_size,0);
		for(int k=0; k<seq[depth_idx].size(); k++){
			_temp_col_indx[seq[depth_idx][k]]++;
		}

		rdm_step_auc();
		depth_idx++;

		if(last_auc_value>=_vec_auc[_sel_snp.size()-1]){
			_vec_auc.pop_back();
			_sel_snp.pop_back();
			_sel_genotype.pop_back();
			break;
		}

	} while (1);	
}

void LMW::rdm_step_auc()
{
	patient_devide();
	int i=0;
	for(i=0; i<_marker_num; i++){

		if(_temp_col_indx[i]==0) continue;

		_temp_col_rcd.push_back(i);


		col_auc(i);
	}
	ranking_auc2();// ranking auc and input the result to the "selected record"
	_temp_sel_auc.clear();
	_temp_sel_genotype.clear();
	_temp_col_rcd.clear();
}

void LMW::ranking_auc2()
{
	int i=0;
	int m_idx=0;

	m_idx=Stat_fuc::max_index(_temp_sel_auc);

	int col=_temp_col_rcd[m_idx];
	int geno=_temp_sel_genotype[m_idx];
	double auc=_temp_sel_auc[m_idx];

	_sel_snp.push_back(col);
	_sel_genotype.push_back(geno);
	_vec_auc.push_back(auc);

	_left_genotype[col][geno]=false;
	_left_genotype_num[col]--;

}

void LMW::auc_VC_fm_LR(vector<double> & pre_score, vector<double> & score, int step)
{
	double auc=0, var=0, cov=0;
	pre_score=score;
	Stat_fuc::score_frm_LR(_LR_nom,_aff,score);		

	if(step==0){
		auc=Stat_fuc::auc_frm_score(score,_aff,var);
		cov=0;
	}else{
		Stat_fuc::aucVarCov_frm_score(pre_score,_vec_auc[step-1],score,_aff,auc,var,cov);
	}
	_vec_auc.push_back(auc); _var_auc.push_back(var); _cov_auc.push_back(cov);
}

vector<double> LMW::apply_test_LR(vector<int> & sel_genotype, vector<int> & sel_snp, vector< vector<int> > & no_exist,
										 vector< vector<double> > & sel_LR_value, vector<double> & indi_LR)
{
	_vec_LR_value=sel_LR_value;
	int i,n_step=sel_snp.size();
	vector<double> score, pre_score;
	patient_devide();
	for (i=0;i<n_step;i++)
	{
		_sel_snp.push_back(sel_snp[i]);
		_sel_genotype.push_back(sel_genotype[i]);
		apply_devide_mis(sel_LR_value[i], no_exist[i]);

		auc_VC_fm_LR(pre_score,score,i);
	}


	vector<double> test_stat;
	double temp1=0, temp2=0;
	for(i=0; i<n_step; i++){
		if(i==0){
			temp1=_var_auc[i];
			temp2=_vec_auc[i]-0.5;
			if(temp1>0){
				temp1=sqrt(temp1);
				test_stat.push_back(temp2/temp1);
			}else{
				test_stat.push_back(0.0);
			}
		}else{
			temp1=(_var_auc[i]+_var_auc[i-1]-(2.0*_cov_auc[i]));
			temp2=_vec_auc[i]-_vec_auc[i-1];
			if(temp1>0){
				temp1=sqrt(temp1);
				test_stat.push_back(temp2/temp1);
			}else{
				test_stat.push_back(0.0);
			}
		}
	}

	indi_LR=_LR_nom;

	return test_stat;
}