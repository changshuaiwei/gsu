#include "kinship.h"

void kinship::initialize(snpdt *snp_data)
{
	snp_data->getPopuInfo(_ori_data);

	snp_data->getFamInfo(_sample);

	find_sample();

	for(int i=0; i<_ori_data.size(); i++){
		_ori_data[i]->inner_ID=i;
	}
}

void kinship::ksMAT(GTmat_ & ks)
{
	ks.resize(_sample_idx.size(), _sample_idx.size());

	for(int i=0; i<_sample_idx.size(); i++){
		double tmpks=1 + cal_inbreeding(_ori_data[_sample_idx[i]]);
		ks(i,i)=tmpks;
		for(int j=i+1; j<_sample_idx.size(); j++){
			double tmpks= 2 * cal_kinship(_ori_data[_sample_idx[i]], _ori_data[_sample_idx[j]]);
			ks(i,j)=tmpks;
			ks(j,i)=tmpks;
		}
	}
}

void kinship::find_sample()
{
	_sample_idx.clear();
	int i=0, j=0;

	for(i=0; i<_sample.size(); i++){

		bool is_here=false;
		
		for(j=0; j<_ori_data.size(); j++){
			if(_sample[i]->iid == _ori_data[j]->iid){
				_sample_idx.push_back(j);
				is_here=true;
				break;
			}
		}

		if(!is_here){
			_ori_data.push_back(_sample[i]);
			_sample_idx.push_back(j);
		}

	}
}

void kinship::construct_pedi()
{
	cout<<"\nnow we construct pedigree....\n";
	int i, j;
	for(i=0; i<_ori_data.size(); i++)
	{
		cout<<"now add the "<<i+1<<"th person....";
		if(_ori_data[i]->mat != "0" ){
			for(j=0; j<_ori_data.size(); j++){
				if(_ori_data[i]->mat == _ori_data[j]->iid ){
					_ori_data[i]->M_addr = _ori_data[j];
					(_ori_data[j]->children).push_back(_ori_data[i]);
					break;
				}
			}

			//this part may not be needed
			/*
			if(j==_ori_data.size()){
				Individual* tmp = new Individual();
				tmp->iid=_ori_data[i].fid;
				_ori_data.push_back(tmp);
				_ori_data[i]->M_addr=_ori_data[j];
				(_ori_data[j]->children).push_back(_ori_data[i]);
			}
			*/
		}

		if(_ori_data[i]->pat != "0"){
			for(j=0; j<_ori_data.size(); j++){
				if(_ori_data[i]->pat == _ori_data[j]->iid ){
					_ori_data[i]->F_addr=_ori_data[j];
					(_ori_data[j]->children).push_back(_ori_data[i]);
					break;
				}
			}
			//this part may not be needed
			/*
			if(j==_ori_data.size()){
				Individual* tmp = new Individual();
				tmp->iid=_ori_data[i]->pat;
				_ori_data.push_back(tmp);
				_ori_data[i]->F_addr=_ori_data[j];
				(_ori_data[j]->children).push_back(_ori_data[i]);
			}
			*/
		}

		cout<<"\r";
	}
}


double kinship::cal_kinship(Individual* indiv1, Individual* indiv2)
{
	vector<Individual*> ans;
	common_ans(ans, indiv1, indiv2);
	if(ans.size()==0) return 0;

	int i, j;
	double ksp=0;
	vector< vector<Individual*> > routs1, routs2;
	for(i=0; i<ans.size(); i++){
		//cout<<"\nget the rout for the first one";
		cal_rout(routs1, indiv1, ans[i]);
		//cout<<"\nget the rout for the second one";
		cal_rout(routs2, indiv2, ans[i]);
		//cout<<"\ncalculate inbreeding";
		double ibrd=cal_inbreeding(ans[i]);
		//cout<<"\ncalculate coeff";
		double coeff=rout_coeff(routs1, routs2);
		ksp+=(1.0+ibrd)*coeff;
	}
	return ksp;
}

void kinship::common_ans(vector<Individual*> & ans, Individual* indiv1, Individual* indiv2)
{
	vector<Individual*> ans1, ans2;

	if(indiv1->inner_ID==indiv2->inner_ID) return;

	ancestor(ans1, indiv1);
	ancestor(ans2, indiv2);

	int i, j;
	vector<int> tmp_inner_ID;
	tmp_inner_ID.resize(_ori_data.size()+1);


	for(i=0; i<ans1.size(); i++){
		int tmpid=ans1[i]->inner_ID;
		tmp_inner_ID[tmpid]=1;
	}

	for(i=0; i<ans2.size(); i++){
		int tmpid=ans2[i]->inner_ID;
		if(tmp_inner_ID[tmpid]==1){
			ans.push_back(ans2[i]);
		}
	}

	ans1.clear(); ans2.clear();
}

void kinship::ancestor(vector<Individual*> & ans, Individual* indiv)
{
	vector<Individual*> tmp;
	tmp.push_back(indiv);
	vector<Individual*> tmp1;
	ans.push_back(indiv);

	int i, j;
	while (tmp.size()>0)
	{
		for(i=0; i<tmp.size(); i++){
			if(tmp[i]->F_addr!=0){
				ans.push_back(tmp[i]->F_addr);
				tmp1.push_back(tmp[i]->F_addr);
			}
			if(tmp[i]->M_addr!=0){
				ans.push_back(tmp[i]->M_addr);
				tmp1.push_back(tmp[i]->M_addr);
			}
		}
		tmp.clear();
		tmp=tmp1;
		tmp1.clear();
	}
}

void kinship::cal_rout(vector< vector<Individual*> > & routs, Individual* indiv, Individual* ans)
{
	Individual* root=indiv;
	Individual* tmp=root;
	Individual* aim=ans;
	_tmp_routs.clear();

	search_tree(tmp, aim, root);

	routs=_tmp_routs;
	_tmp_routs.clear();
}

void kinship::search_tree(Individual* tmp, Individual* aim, Individual* root)
{
	if(tmp==aim){
		Individual* tmp0=tmp;
		vector<Individual*> tmp_rout;
		while(tmp0!=root){
			tmp0=tmp0->tmp_child;
			tmp_rout.push_back(tmp0);
		}
		_tmp_routs.push_back(tmp_rout);
		return;
	}else{ 
		if(tmp->F_addr!=0){
			tmp->F_addr->tmp_child=tmp;
			search_tree(tmp->F_addr, aim, root);
		}
		if(tmp->M_addr!=0){
			tmp->M_addr->tmp_child=tmp;
			search_tree(tmp->M_addr, aim, root);
		}
	}
}

double kinship::cal_inbreeding(Individual* indiv)
{
	if(indiv->F_addr==0 || indiv->M_addr==0){
		return 0;
	}else{
		double ibrd=cal_kinship(indiv->F_addr, indiv->M_addr);
		return ibrd;
	}
}

double kinship::rout_coeff(vector< vector<Individual*> > & routs1, vector< vector<Individual*> > & routs2)
{
	int i, j;
	double coeff=0;
	for(i=0; i<routs1.size(); i++){
		for(j=0; j<routs2.size(); j++){
			if(compare_rout(routs1[i], routs2[j])){
				double rout_lgth=routs1[i].size()+routs2[j].size()+1.0;
				coeff+=pow(0.5, rout_lgth);
			}
		}
	}

	return coeff;
}

bool kinship::compare_rout(vector<Individual*> & rout1, vector<Individual*> & rout2)//no common nodes on the two rout
{
	bool different=true;

	vector<int> indi_ID(_ori_data.size()+1, 0);
	int i, j;
	for(i=0; i<rout1.size(); i++){
		indi_ID[rout1[i]->inner_ID]=1;
	}
	for(i=0; i<rout2.size(); i++){
		if(indi_ID[rout2[i]->inner_ID]==1){
			different=false;
			break;
		}
	}
	if(rout1.size()==0 && rout2.size()==0) different=false;	

	return different;
}
