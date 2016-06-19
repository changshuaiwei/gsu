#include "hwu.h"

void HWU::initialize(snpdt * snp_data)
{
	_datafile=snp_data;
	_n_sub=_datafile->totalIndi();
	_n_snp=_datafile->totalLocus();

	vector<double> y;
	
	for(int i=0; i<_n_sub; i++){
		y.push_back(_datafile->qPheno(i));
	}

	_org_y=y;

	_y=sry(_org_y);
	//for(int i=0; i<_n_sub; i++) cout<<_y[i]<<"\t";

	outer_YXW();
	//_buffer_cal=new double [_n_sub*_n_sub];

}

vector<double> HWU::sry(vector<double> & y)
{
	vector<double> rk;
	rk.resize(y.size(),0);

	if(par::hwu_rk){
		vector<int> idx;
		Stat_fuc::indexx(y,idx);

		//get the rank
		int st=0, ed=-1;
		double tmp=y[idx[0]];
		for(int i_=0; i_<idx.size(); i_++){
			if(y[idx[i_]]==tmp){
				ed++;
			}else{
				double tmprk=(double(st+ed+2))/2;
				for(int j=st; j<=ed; j++){
					rk[idx[j]]=tmprk;
				}
				st=i_; ed=i_;
				tmp=y[idx[i_]];
			}

			if(i_==idx.size()-1){
				double tmprk=(double(st+ed+2))/2;
				for(int j=st; j<=ed; j++){
					rk[idx[j]]=tmprk;
				}
			}
		}
	}else{
		rk=y;
	}
	

	double miu, sd;

	Dist::meansd(rk,miu,sd);

	for(int i=0; i<rk.size(); i++){
		if(sd>0) rk[i]=(rk[i]-miu)/sd;
		else rk[i]=rk[i]-miu;
	}

	return rk;
}

vector<double> HWU::getrk(vector<double> & y)
{
	vector<double> rk;
	rk.resize(y.size(),0);

	vector<int> idx;
	Stat_fuc::indexx(y,idx);

	//get the rank
	int st=0, ed=-1;
	double tmp=y[idx[0]];
	for(int i_=0; i_<idx.size(); i_++){
		if(y[idx[i_]]==tmp){
			ed++;
		}else{
			double tmprk=(double(st+ed+2))/2;
			for(int j=st; j<=ed; j++){
				rk[idx[j]]=tmprk;
			}
			st=i_; ed=i_;
			tmp=y[idx[i_]];
		}

		if(i_==idx.size()-1){
			double tmprk=(double(st+ed+2))/2;
			for(int j=st; j<=ed; j++){
				rk[idx[j]]=tmprk;
			}
		}
	}

	return rk;
}

void HWU::outer_YXW()
{
	/*
	_outer_y.clear();
	vector<double> tmp(_y.size(),0);
	_outer_y.resize(_y.size(),tmp);

	for(int i=0; i<_y.size(); i++){
		for(int j=i; j<_y.size(); j++){
			_outer_y[i][j]=_y[i]*_y[j];
			_outer_y[j][i]=_outer_y[i][j];
		}
	}
	*/

	_outer_y=new Dmat(_y.size(),_y.size());
	_outer_x=new Dmat(_y.size(),_y.size());
	_weight=new Dmat(_y.size(),_y.size());
	_weight_2=new Dmat(_y.size(),_y.size());
	//_weight_2=new Dmat(_y.size(),_y.size());
	//int size=_y.size();
	//Dmat d1(size,size), d2(size,size), d3(size,size);
	//_outer_y=&d1; _outer_x=&d2; _weight=&d3;

	for(int i=0; i<_y.size(); i++){
		for(int j=i; j<_y.size(); j++){
			double tmp=_y[i]*_y[j];
			_outer_y->assign(i,j,tmp);
			_outer_y->assign(j,i,tmp);
		}
	}

}

void HWU::outer_X(vector<double> & x, vector< vector<double> > & ot_X)
{
	ot_X.clear();
	vector<double> tmp(x.size(),0);
	ot_X.resize(x.size(),tmp);

	for(int i=0; i<x.size(); i++){
		for(int j=i; j<x.size(); j++){
			ot_X[i][j]=x[i]*x[j];
			ot_X[j][i]=ot_X[i][j];
		}
	}
}

void HWU::outer_X(vector<double> & x)
{
	for(int i=0; i<x.size(); i++){
		for(int j=i; j<x.size(); j++){
			double tmp=x[i]*x[j];
			_outer_x->assign(i,j,tmp);
			_outer_x->assign(j,i,tmp);
		}
	}
}

void HWU::update_X()
{
	for(int i=0; i<_n_sub; i++){
		for(int j=i; j<_n_sub; j++){
			double tmp=(_outer_x->get(i,j))*(_weight->get(i,j));
			_outer_x->assign(i,j,tmp);
			_outer_x->assign(j,i,tmp);
		}
	}
}

void HWU::weightFromSNP(int size)
{
	vector< vector<double> > dt;
	for(int i=0; i<size; i++){
		int idx=-1;
		while(idx<0 || idx>=_n_snp){
			idx=int(double(_n_snp)*Stat_fuc::ran1(par::seed)+0.5);
		}
		
		vector<double> tmp;
		for(int j=0; j<_n_sub; j++){
			int tmptmp=_datafile->genotypeToInt(j,idx);
			if(tmptmp<0) tmptmp=int(3.0*Stat_fuc::ran1(par::seed));
			tmp.push_back(double(tmptmp));
		}
		dt.push_back(tmp);
	}


	_weight->copy(Dist::euclid(dt));
	_weight_flag=true;

	//for(int i=0; i<_n_sub; i++){
	//	for(int j=0; j<_n_sub; j++){
	//		double tmp=_weight->get(i,j);
	//		tmp*=tmp;
	//		_weight_2->assign(i,j,tmp);
	//	}
	//}

}

void HWU::weightFromCov()
{
	//cout<<"Relatedness From Covariates ..\n";
	gfun::printLOG("Relatedness From Covariates ..\n");

	vector< vector<double> > dt;
	_datafile->getCovMat(dt);

	_weight->copy(Dist::euclid(dt));
	_weight_flag=true;

#ifdef _DEBUG_HWU_WEI_
	cout<<"\n\n";
	for(int i=0; i<4; i++){
		for(int j=0; j<4; j++){
			cout<<_weight->get(i,j)<<"\t";
		}
		cout<<"\n";
	}
#endif

	//cout<<"\nRelateness Matrix Done.\n\n";
	gfun::printLOG("\nRelateness Matrix Done.\n\n");
}

//calculate kinship matrix: weighted IBS(average correlation)
//negative value is possible
void HWU::weightFromGenomeWIBS()
{
	//initialize tmp weight
	vector<float> tmp2, tmp(_n_sub,0);
	vector< vector<float> > weight(_n_sub,tmp);

	vector<int> mis_idx;
	float tmpgeno=0;
	float tmp_ave=0;
	float tmp_var=0;

	float count_snp=0;

	//cout<<"Relatedness From Whole Genome..\n";
	gfun::printLOG("Relatedness From Whole Genome..\n");

	vector<int> idx;
	double dtotal=_n_snp;

	for(int i=0; i<par::IBS_N; i++){
		int indx=(int)(Stat_fuc::ran1(par::seed)*(dtotal));
		while (indx<0 || indx>(_n_snp-1)) indx=(int)(Stat_fuc::ran1(par::seed)*(dtotal));
		idx.push_back(indx);
	}

	int totalSNP=0;
	if(par::IBS_N>0) totalSNP=par::IBS_N;
	else totalSNP=_n_snp;

	int i=0;
	for(int i_=0; i_<totalSNP; i_++){

		if (par::IBS_N>0) i=idx[i_];
		else i=i_;

		cout<<"include "<<i_+1<<"th SNP..";

		mis_idx.clear();
		tmp_ave=0;
		//tmp_var=0;
		for(int j=0; j<_n_sub; j++){
			tmpgeno=float(_datafile->genotypeToInt(j,i));
			if(tmpgeno<0){
				mis_idx.push_back(j);
			}else{
				tmp[j]=tmpgeno;
				tmp_ave+=tmpgeno;
				//tmp_var+=tmpgeno*tmpgeno;
			}
		}
		tmp_ave/=float(_n_sub-mis_idx.size());
		//tmp_var/=float(_n_sub-mis_idx.size());

		tmp_var=tmp_ave*(2.0-tmp_ave);
		//tmp_var-=tmp_ave*tmp_ave;
		tmp_var=sqrt(tmp_var);

		if(tmp_var==0) {
			cout<<"\r";
			cout.flush();
			continue;
		}

		count_snp+=1.0;

		for(int j=0; j<mis_idx.size(); j++){
			tmp[mis_idx[j]]=tmp_ave;
		}

		for(int j=0; j<_n_sub; j++){
			tmp[j]-=tmp_ave;
			tmp[j]/=tmp_var;
		}

		//now tmp is standardized, we can make outer product and sum to weight

		for(int j=0; j<weight.size(); j++){
			for(int k=j; k<weight[j].size(); k++){
				weight[j][k]+=(tmp[j]*tmp[k]);
				weight[k][j]=weight[j][k];
			}
		}

		cout<<"\r";
		cout.flush();

	}

	//average now
	vector<double> tmpd(_n_sub,0);
	vector< vector<double> > weightd(_n_sub,tmpd);

	for(int i=0; i<_n_sub; i++){
		for(int j=i; j<_n_sub; j++){
			weightd[i][j]=weight[i][j]/count_snp;
			weightd[j][i]=weightd[i][j];
		}
	}

	//trans to correlation (-1,+1), and then distance (0,1)
	vector<double> var_vec;
	for(int i=0; i<_n_sub; i++){
		var_vec.push_back(sqrt(weightd[i][i]));
	}
	for(int i=0; i<_n_sub; i++){
		weightd[i][i]=1.0;
		for(int j=i+1; j<_n_sub; j++){
			double tmpcov=var_vec[i]*var_vec[j];
			if(tmpcov>0){
				weightd[i][j]=(weightd[i][j]/tmpcov)*0.5 +0.5;
			}else{
				weightd[i][j]=0;
				//gfun::error("error in calculating distance matrix\n");
			}
			weightd[j][i]=weightd[i][j];
		}
	}



#ifdef _DEBUG_HWU_WEI_
	cout<<"\n\n";
	for(int i=0; i<4; i++){
		for(int j=0; j<4; j++){
			cout<<weightd[i][j]<<"\t";
		}
		cout<<"\n";
	}
	//system("pause");
#endif


	//assign to _weight
	_weight->copy(weightd);
	_weight_flag=true;

	//cout<<"\nRelateness Matrix Done.\n\n";
	gfun::printLOG("\nRelateness Matrix Done.\n\n");
	
}

void HWU::weightFromGenomeIBS()
{
	vector<float> tmp2, tmp(_n_sub,0);
	vector< vector<float> > weight(_n_sub,tmp);

	vector<int> mis_idx;
	float tmpgeno=0;
	float tmp_ave=0;
	float tmp_var=0;

	float count_snp=0;

	vector<int> idx;
	double dtotal=_n_snp;

	for(int i=0; i<par::IBS_N; i++){
		int indx=(int)(Stat_fuc::ran1(par::seed)*(dtotal));
		while (indx<0 || indx>(_n_snp-1)) indx=(int)(Stat_fuc::ran1(par::seed)*(dtotal));
		idx.push_back(indx);
	}

	int totalSNP=0;
	if(par::IBS_N>0) totalSNP=par::IBS_N;
	else totalSNP=_n_snp;

	int i;
	for(int i_=0; i_<totalSNP; i_++){

		if (par::IBS_N>0) i=idx[i_];
		else i=i_;

		mis_idx.clear();
		tmp_ave=0;
		for(int j=0; j<_n_sub; j++){
			tmpgeno=float(_datafile->genotypeToInt(j,i));
			if(tmpgeno<0){
				mis_idx.push_back(j);
			}else{
				tmp[j]=tmpgeno;
				tmp_ave+=tmpgeno;
			}
		}
		tmp_ave/=float(_n_sub-mis_idx.size());

		tmp_var=tmp_ave*(2.0-tmp_ave);
		tmp_var=sqrt(tmp_var);

		if(tmp_var==0) continue;

		count_snp+=1.0;

		for(int j=0; j<mis_idx.size(); j++){
			tmp[mis_idx[j]]=tmp_ave;
		}

		for(int j=0; j<_n_sub; j++){
			tmp[j]-=1.0;
		}

		for(int j=0; j<weight.size(); j++){
			for(int k=j; k<weight[j].size(); k++){
				weight[j][k]+=tmp[j]*tmp[k];
				weight[k][j]=weight[j][k];
			}
		}

	}

	//average now
	vector<double> tmpd(_n_sub,0);
	vector< vector<double> > weightd(_n_sub,tmpd);

	for(int i=0; i<_n_sub; i++){
		for(int j=i; j<_n_sub; j++){
			weightd[i][j]=weight[i][j]/count_snp/2.0+0.5;
			weightd[j][i]=weightd[i][j];
		}
	}

#ifdef _DEBUG_HWU_WEI_
	cout<<"\n\n";
	for(int i=0; i<4; i++){
		for(int j=0; j<4; j++){
			cout<<weightd[i][j]<<"\t";
		}
		cout<<"\n";
	}
#endif


	//assign to _weight
	_weight->copy(weightd);
	_weight_flag=true;
}


void HWU::weightFastGenome()
{
	vector<float> tmp(_n_sub,0);
	vector< vector<float> > weight(_n_sub,tmp);

	vector<int> tmp2(_n_sub,0);
	vector<bool> flag_9(_n_sub,0), flag_first(_n_sub,0), flag_second(_n_sub,0);
	int tmpgeno=0;
	const float one=1;

	//cout<<"Relatedness From Whole Genome..\n";
	gfun::printLOG("Relatedness From Whole Genome..\n");

	vector<int> idx;
	double dtotal=_n_snp;

	for(int i=0; i<par::IBS_N; i++){
		int indx=(int)(Stat_fuc::ran1(par::seed)*(dtotal));
		while (indx<0 || indx>(_n_snp-1)) indx=(int)(Stat_fuc::ran1(par::seed)*(dtotal));
		idx.push_back(indx);
	}

	int totalSNP=0;
	if(par::IBS_N>0) totalSNP=par::IBS_N;
	else totalSNP=_n_snp;

	int i;
	for(int i_=0; i_<totalSNP; i_++){

		if (par::IBS_N>0) i=idx[i_];
		else i=i_;

		cout<<"include "<<i_+1<<"th SNP..";

		
		for(int j=0; j<_n_sub; j++){
			tmp2[j]=_datafile->allel1ToInt(j,i);
		}

		for(int j=0; j<weight.size(); j++){
			for(int k=j; k<weight[j].size(); k++){
				if(tmp2[j]!=tmp2[k]) weight[j][k]+=one;
			}
		}

		for(int j=0; j<_n_sub; j++){
			tmp2[j]=_datafile->allel2ToInt(j,i);
		}

		for(int j=0; j<weight.size(); j++){
			for(int k=j; k<weight[j].size(); k++){
				if(tmp2[j]!=tmp2[k]) weight[j][k]+=one;
			}
		}
		
		/*
		for(int j=0; j<_n_sub; j++){
			tmpgeno=_datafile->genotypeToInt(j,i);
			if(tmpgeno==-9) flag_9[j]=true;
			else{
				flag_9[j]=false;
				if(tmpgeno<2) flag_first[j]=true; else flag_first[j]=false;
				if(tmpgeno>0) flag_second[j]=true; else flag_second[j]=false;
			}
		}

		for(int j=0; j<_n_sub; j++){
			if(!flag_9[j]){
				for(int k=j; k<_n_sub; k++){
					if(!flag_9[k]){
						float tmpadd=short((flag_first[j]!=flag_first[k]))+ short(flag_second[j]!=flag_second[k]);
						weight[j][k]+=tmpadd;
					}
					
				}
			}
			
		}
		*/

		cout<<"\r";
		cout.flush();

	}

	//average now
	vector<double> tmpd(_n_sub,0);
	vector< vector<double> > weightd(_n_sub,tmpd);

	for(int i=0; i<_n_sub; i++){
		for(int j=i; j<_n_sub; j++){
			weightd[i][j]=double(weight[i][j])/double(totalSNP)/2;
			weightd[i][j]=1-weightd[i][j];
			weightd[j][i]=weightd[i][j];
		}
	}

#ifdef _DEBUG_HWU_WEI_
	cout<<"\n\n";
	for(int i=0; i<4; i++){
		for(int j=0; j<4; j++){
			cout<<weightd[i][j]<<"\t";
		}
		cout<<"\n";
	}
#endif


	//assign to _weight
	_weight->copy(weightd);
	_weight_flag=true;

	//cout<<"\nRelateness Matrix Done.\n\n";
	gfun::printLOG("\nRelateness Matrix Done.\n\n");
}

void HWU::prepareY()
{
	if(par::prjct_Y && _weight_flag) projectY();
	if(par::pj_1PC && _weight_flag) projectY2();
	if(par::pj_Y && _weight_flag) pjY();

	//if(par::pj_cov && _weight_flag) projectY3();
	/*
	else if(par::prjct_Y_New_Distr){
	vector< vector<double> > X, P;
	if(_weight_flag) princ_comp(X);

	projectionMat(X,P);
	prepProject(P);


	}
	*/
}

void HWU::sgLocusRank()
{
	prepareY();
	stdW();

	if(par::appx_davis){
		getRankWU();
	}else if(par::appx_davis_P){
		getRankP();
	}

	vector<double> x;

	cout<<"\nP-value for ranked:\n";
	double pvalue;
	for(int i_=0; i_<_idx.size(); i_++){
		int i=_idx[i_];
		cout<<i_+1<<"th SNP";
		x.clear();
		//cout<<":"<<(_datafile->getLocus(i))->name;
		prepareX(i,x);
		//cout<<":Xprepared";

		//if(i_==230) system("pause");

		pvalue=hwu_davis(x);
		//cout<<":DavisDone";
		_vec_pvalue.push_back(pvalue);
		//wu=hwu_inter_U(x);
		//_weightU.push_back(wu);
		cout<<"\r";
		cout.flush();
	}

}

void HWU::getRankWU()
{
	cout<<"\nGet WU for\n";
	vector<double> x; double wu;
	_weightU.clear();
	for(int i=0; i<_n_snp; i++){
		cout<<i+1<<"th SNP";
		x.clear();
		prepareX(i,x);
		wu=hwu_inter_U(x);
		_weightU.push_back(wu);
		cout<<"\r";
		cout.flush();
	}

	cout<<"\nSorting WU\n";
	vector<int> tmpidx, idx;
	Stat_fuc::indexx(_weightU,tmpidx);

	int count=0;
	for(int i=tmpidx.size()-1; i>=0; i--){
		idx.push_back(tmpidx[i]);
		count++;

		if(count>=par::n_davis_snp) break;
	}
	_idx=idx;
	_bk_weightU=_weightU;
	_weightU.clear();
}

void HWU::getRankP()
{
	cout<<"\nGet pseudo-P for\n";
	vector<double> x; double p;
	_vec_pvalue.clear();
	for(int i=0; i<_n_snp; i++){
		cout<<i+1<<"th SNP";
		x.clear();
		prepareX(i,x);
		p=hwu_inter_p(x);
		_vec_pvalue.push_back(p);
		cout<<"\r";
		cout.flush();
	}

	cout<<"\nSorting pseudo-P\n";
	vector<int> tmpidx, idx;
	Stat_fuc::indexx(_vec_pvalue,tmpidx);

	int count=0;
	for(int i=0; i<tmpidx.size(); i++){
		idx.push_back(tmpidx[i]);
		count++;

		if(count>=par::n_davis_snp) break;
	}
	_idx=idx;
	_bk_vec_pvalue=_vec_pvalue;
	_vec_pvalue.clear();
}

void HWU::updateP()
{
	cout<<"\nSorting P-value\n";
	vector<int> tmpidx, idx;
	Stat_fuc::indexx(_vec_pvalue,tmpidx);

	int count=0;
	for(int i=0; i<tmpidx.size(); i++){
		idx.push_back(tmpidx[i]);
		count++;

		if(count>=par::n_davis_snp) break;
	}
	_idx=idx;
	_bk_vec_pvalue=_vec_pvalue;
	
	vector<double> x;
	cout<<"\nUpdate P-value for ranked:\n";
	double pvalue;
	for(int i_=0; i_<_idx.size(); i_++){
		int i=_idx[i_];
		cout<<i_+1<<"th SNP";
		x.clear();
		//cout<<":"<<(_datafile->getLocus(i))->name;
		prepareX(i,x);
		//cout<<":Xprepared";

		//if(i_==230) system("pause");

		pvalue=hwu_davis(x);
		//cout<<":DavisDone";
		_vec_pvalue[_idx[i_]]=pvalue;
		//wu=hwu_inter_U(x);
		//_weightU.push_back(wu);
		cout<<"\r";
		cout.flush();
	}

}

void HWU::updatePcut()
{
	
	vector<int> tmpidx, idx;

	for(int i=0; i<_vec_pvalue.size(); i++){
		if(_vec_pvalue[i]<par::cut_up) idx.push_back(i);
	}
	cout<<"\n"<<idx.size()<<" P-values <"<<par::cut_up<<"\n";
	_idx=idx;
	_bk_vec_pvalue=_vec_pvalue;

	vector<double> x;
	cout<<"\nUpdate P-value for:\n";
	double pvalue;
	for(int i_=0; i_<_idx.size(); i_++){
		int i=_idx[i_];
		cout<<i_+1<<"th SNP";
		x.clear();
		//cout<<":"<<(_datafile->getLocus(i))->name;
		prepareX(i,x);
		//cout<<":Xprepared";

		//if(i_==230) system("pause");

		pvalue=hwu_davis(x);
		//cout<<":DavisDone";
		_vec_pvalue[_idx[i_]]=pvalue;
		//wu=hwu_inter_U(x);
		//_weightU.push_back(wu);
		cout<<"\r";
		cout.flush();
	}
}
void HWU::sgLocusAssoc()
{
	prepareY();
	stdW();

	_vec_pvalue.clear();
	double pvalue=1;
	cout<<"\n";
	for(int i=0; i<_n_snp; i++){
		cout<<i+1<<"th SNP";
		pvalue=hwu_i(i);
		_vec_pvalue.push_back(pvalue);
		cout<<"\r";
		cout.flush();
	}

	if(par::mm_davis) {
		if(par::cut_upB){
			updatePcut();
		}else{
			updateP();
		}
	}
}

void HWU::wtResult(string outputfile)
{
	gfun::printLOG("Writing HWU Result at [ " + outputfile + " ]\n");

	ofstream result(outputfile.c_str());
	if(!result) gfun::error("\ncould not open the result file\n");
	result<<"#result of association scanning (0=A1A1, 1=A1A2, 2=A2A2)\n";
	
	if(par::appx_davis || par::appx_davis_P){ 
		result<<"#order\t"<<"Chr\tName\tPos\tBp\tAllel1\tAllel2\t"
		<<"U\t"
		<<"P_value\n";
	}else{
		result<<"#order\t"<<"Chr\tName\tPos\tBp\tAllel1\tAllel2\t"
			<<"U\tObs\tExp\tVar\t"
			<<"P_value";
		if(par::mm_davis){
			result<<"\tP_moment_match\n";
		}else{
			result<<"\n";
		}
	}

	Locus * loc=0;
	for(int i=0; i<_vec_pvalue.size(); i++){
		if(par::appx_davis || par::appx_davis_P){
			loc=_datafile->getLocus(_idx[i]);
		}else{
			loc=_datafile->getLocus(i);
		}
		

		result<<i+1<<"\t";//order

		result	<<loc->chr << "\t"
			<< loc->name << "\t"
			<< loc->pos << "\t"
			<< loc->bp << "\t"
			<< loc->allele1 <<"\t"
			<< loc->allele2 <<"\t";//snp information

		if(par::appx_davis || par::appx_davis_P){
			result	<<_weightU[i]<<"\t"
				<<_vec_pvalue[i]<<"\n";//pvalue
		}else{
			result	<<_weightU[i]<<"\t"
				<<_Obs[i]<<"\t"
				<<_expct[i]<<"\t"
				<<_varian[i]<<"\t"
				<<_vec_pvalue[i];//pvalue
			if(par::mm_davis){
				result<<"\t"<<_bk_vec_pvalue[i]<<"\n";
			}else{
				result<<"\n";
			}
		}
		
	}

	result.close();

}

void HWU::wtWeight(string outputfile)
{
	gfun::printLOG("Writing Weight Matrix at [ " + outputfile + " ]\n");

	ofstream result(outputfile.c_str());
	if(!result) gfun::error("\ncould not write the Weight Matrix file\n");
	//result<<"#result of association scanning (0=A1A1, 1=A1A2, 2=A2A2)\n";
	//result<<"#order\t"<<"Chr\tName\tPos\tBp\tAllel1\tAllel2\t"<<"P_value\n";

	//for(int i=0; i<_n_sub; i++){
	//	result<<"W"<<i+1<<"\t";
	//}
	//result<<"\n";

	result<<_n_sub<<"\n";

	for(int i=0; i<_n_sub; i++){
		for(int j=i; j<_n_sub; j++){
			result<<_weight->get(i,j)<<"\t";
		}
		result<<"\n";
	}

	result.close();

}

void HWU::rdWeight(string inputfile)
{
	gfun::printLOG("Reading Weight Matrix at [ " + inputfile + " ]\n");

	ifstream result(inputfile.c_str());
	if(!result) gfun::error("\ncould not read the Weight Matrix file\n");
	//result<<"#result of association scanning (0=A1A1, 1=A1A2, 2=A2A2)\n";
	//result<<"#order\t"<<"Chr\tName\tPos\tBp\tAllel1\tAllel2\t"<<"P_value\n";

	//vector<string> str_vec;
	//string buf_str;

	//getline(result,buf_str);
	//int dim=gfun::split_string(buf_str,str_vec);
	int dim=0;
	result>>dim;

	if(dim!=_n_sub) gfun::error("Dim of Matrix != sample size\n");

	vector<double> tmpvec(dim,0);
	vector< vector<double> > weight(dim,tmpvec);
	double tmp;
	for(int i=0; i<_n_sub; i++){
		for(int j=i; j<_n_sub; j++){
			result>>tmp;
			float tmptmp=tmp;
			weight[i][j]=tmptmp;
			weight[j][i]=weight[i][j];
			//_weight->assign(i,j,tmp);
			//_weight->assign(j,i,tmp);
		}
	}

	result.close();

	_weight->copy(weight);

	_weight_flag=true;

	cout<<"\nRelateness Matrix Done.\n\n";

}

void HWU::prepareX(int i_snp, vector<double> &x)
{
	//if(i_snp==157) system("pause");

	x.clear();
	vector<double> x_;
	vector<int> x_rcd;
	for(int i=0; i<_n_sub; i++){
		int geno=(_datafile->genotypeToInt(i,i_snp));
		x.push_back(double(geno));
		if(geno==-9){
			x_rcd.push_back(i);
		}else{
			x_.push_back(double(geno));
		}
	}

	//make sure minor allele is coded as 0
	double tst_minor=0;
	tst_minor=Stat_fuc::mean(x_);
	if(tst_minor>1.0){
		for(int i=0; i<x.size(); i++) x[i]=2-x[i];
	}

	//a very naive way of dealing missing value, by assigning average value
	for(int i=0; i<x_rcd.size(); i++){
		if(tst_minor>1.0) x[x_rcd[i]]=2-tst_minor;
		else x[x_rcd[i]]=tst_minor;
	}
}

double HWU::hwu_i(int i_snp)
{

	vector<double> x;
	prepareX(i_snp,x);

	//analysis
	
	/*
	if(par::prjct_Y_New_Distr){
		pvalue=hwu_i_internalP(x);
	}else{
		pvalue=hwu_i_internal(x);
	}
	*/

	double pvalue=0;
	if(par::pj_1PC && _weight_flag){
		pvalue=hwu_i_first_PC(x);
	}else{
		pvalue=hwu_i_internal(x);
	}
	return pvalue;

}


/*
double HWU::hwu_i_internalP(vector<double> & x)
{

	//using float to speed up
	float tmpsum=0;
	float tmpwt=0, tmpg=0;

	float diagsum=0, sumA2=0;

	float weightU=0;

	vector<float> tmp(_n_sub,0);
	vector< vector<float> > W(_n_sub,tmp);
	_PW.resize(_n_sub,tmp);
	
	//get weight matrix W
	for(int i=0; i<_n_sub; i++){
		if(x[i]!=0){
			for(int j=i; j<_n_sub; j++){
				if(x[j]!=0){
					if(_weight_flag) tmpwt=_weight->get(i,j);
					else tmpwt=1;

					if(tmpwt!=0){
						tmpg=x[i]*x[j];
						tmpwt*=tmpg;
						W[i][j]=tmpwt; W[j][i]=tmpwt;
						weightU+=(_outer_y->get(i,j)*tmpwt);
					}
				}
			}
		}	
	}

	//get matrix P*W
	//way too slow here!!!!!!!!!!!!!!!!
	for(int i=0; i<_n_sub; i++){
		for(int j=0; j<_n_sub; j++){
			tmpsum=0;
			for(int k=0; k<_n_sub; k++){
				if(W[i][k]!=0){
					tmpwt=_Projection->get(j,k);
					tmpsum+=W[i][k]*tmpwt;
				}
			}
			_PW[i][j]=tmpsum;
		}
		cout<<i<<"\t";
	}

	//get trace(PW) and trace(PWPW)
	sumA2=0;
	diagsum=0;
	for(int i=0; i<_n_sub; i++){
		diagsum+=_PW[i][i];
			for(int k=0; k<_n_sub; k++){
				sumA2+=_PW[i][k]*_PW[k][i];
			}
	}


	double expect=0, variance=0;

	expect=diagsum;

	variance=sumA2;
	
	variance*=2.0;

	

	if(variance<=0) return 1.0;

	double s_=variance/expect;
	double a_=expect/s_;


	//get observed value
	double obs=weightU + diagsum;

	cout<<"\t"<<obs<<"\t"<<variance<<"\t"<<expect<<"\n";

	if(obs<=0) return 1.0;

	double pvalue=1.0 - Stat_fuc::P_gamma(obs,a_,s_);

	return pvalue;
}


*/

double HWU::hwu_test(vector<double> & x,vector<double> & y)
{
	/*
	double tmpsum=0,tmp2sum=0,tmpsum2=0;
	double tmpwt=0, tmpwt2=0, tmpg=0;

	double diagsum=0, sumA=0, sumA2=0, sum2A=0;

	double weightU=0;
	*/

	//using float to speed up
	float tmpsum=0,tmp2sum=0,tmpsum2=0;
	float tmpwt=0, tmpwt2=0, tmpg=0;

	float diagsum=0, sumA=0, sumA2=0, sum2A=0;

	float weightU=0;

	for(int i=0; i<_n_sub; i++){
		tmpsum=0; tmp2sum=0;
		if(x[i]!=0){
			for(int j=0; j<_n_sub; j++){
				if(x[j]!=0){
					if(_weight_flag) tmpwt=_weight->get(i,j);
					else tmpwt=1;

					if(tmpwt!=0){
						tmpg=x[i]*x[j];
						tmpwt*=tmpg;
						tmpwt2=tmpwt*tmpwt;

						tmpsum+=tmpwt;
						tmp2sum+=tmpwt2;

						if(i==j) diagsum+=tmpwt;
						else{
							weightU+=(y[i]*y[j]*tmpwt);
						}

					}
				}
			}

			sumA+=tmpsum;
			sumA2+=tmp2sum;
			sum2A+=(tmpsum*tmpsum);

		}

		//cout<<weightU<<"\n";

	}

	double expect=0, variance=0;

	expect=diagsum-sumA/double(_n_sub);

	variance=sumA2 - 2*sum2A/double(_n_sub) + sumA*sumA/double(_n_sub)/double(_n_sub);


	variance*=2.0;

	//get observed value
	double obs=weightU + diagsum;

	//cout<<"\t"<<obs<<"\t"<<variance<<"\t"<<expect<<"\n";

	//_weightU.push_back(weightU);
	//_Obs.push_back(obs);
	//_expct.push_back(expect);
	//_varian.push_back(variance);

	if(variance<=0) return 1.0;
	if(obs<=0) return 1.0;

	double s_=variance/expect;
	double a_=expect/s_;

	double pvalue=1.0 - Stat_fuc::P_gamma(obs,a_,s_);

	return pvalue;
}

void HWU::stdW()
{
	if(!_weight_flag) return;
	if(!par::std_weight) return;

	cout<<"\nstandardizing weight matrix...\n";//to set for testing pure heterogeneity effect
	
	double sum=0;
	for(int i=0; i<_n_sub; i++){
		for(int j=0; j<_n_sub; j++){
			sum+=_weight->get(i,j);
		}
	}

	sum=sum/double(_n_sub)/double(_n_sub);

	for(int i=0; i<_n_sub; i++){
		for(int j=0; j<_n_sub; j++){
			double tmp=_weight->get(i,j)-sum;
			_weight->assign(i,j,tmp);
		}
	}

	cout<<"Done.\n";
}

double HWU::hwu_davis(vector<double> & x)
{
	double tmpwt, tmpg, weightU;

	weightU=0;
	for(int i=0; i<_n_sub; i++){
		if(x[i]!=0){
			for(int j=0; j<_n_sub; j++){
				if(x[j]!=0){
					if(_weight_flag) tmpwt=_weight->get(i,j);
					else tmpwt=1;

					if(tmpwt!=0){
						tmpg=x[i]*x[j];
						tmpwt*=tmpg;

						if(i!=j) weightU+=(_outer_y->get(i,j)*tmpwt);
					}
				}else{
					tmpwt=0;
				}
				
				_weight_2->assign(i,j,tmpwt);
				//_weight_2->assign(j,i,tmpwt);

			}
		}else{
			for(int j=0; j<_n_sub; j++){
				_weight_2->assign(i,j,0);
				//_weight_2->assign(j,i,0);
			}
		}
		//cout<<weightU<<"\n";
	}

	//weightU*=2.0;

	_buffer=_weight_2->getPointer();

	//get eigen value and eigen vectors
	vector<double> egval;
	vector< vector<double> > egvec;
	//TRL::eg_nuTRan(TRL::mt_op,ot_X,egval,egvec,10,0.05,1,30,60,10);

	//Dmat Dmt(ot_X);
	//TRL::eg_nuTRan(TRL::mt_op2,&Dmt,ot_X.size(),egval,egvec,10,0.05,1,30,60,10);
	//TRL::eg_nuTRan(TRL::mt_op2,_outer_x,_n_sub,egval,egvec,10,0.05,1,30,60,10);
	//TRL::eg_nuTRan(TRL::mt_op3,_buffer,_n_sub,egval,egvec,10,0.05,1,30,60,10);

	//TRL::eg_nuTRan(TRL::mt_op3,_buffer,_n_sub,egval,egvec,
	//	par::trl_ned,par::trl_pres,par::trl_scheme,100,200,par::trl_max);

	TRL::eg_nuTRan(TRL::mt_op3,_buffer,_n_sub,egval,egvec,
		par::n_davis_dim,par::trl_pres,par::trl_scheme,100,200,par::trl_max);

	_buffer=0;


	vector<double> c_vec;
	for(int i=0; i<egvec.size(); i++){
		c_vec.push_back(Stat_fuc::sum(egvec[i]));
	}

	vector< vector<double> > sig_mat(c_vec.size(),c_vec);
	for(int i=0; i<sig_mat.size(); i++){
		for(int j=i; j<sig_mat.size(); j++){
			if(i==j) sig_mat[i][i]=1.0-c_vec[i]*c_vec[i]/double(_n_sub);
			else{
				sig_mat[i][j]=0.0-c_vec[i]*c_vec[j]/double(_n_sub);
				sig_mat[j][i]=sig_mat[i][j];
			}
		}
	}

	vector< vector<double> > tmp_mat(sig_mat);
	//Mat_fuc::svdcmp(sig_mat,c_vec,tmp_mat);
	Mat_fuc::Symmeig smg(sig_mat);
	sig_mat=smg.z; c_vec=smg.d;

	//TRL::smg(tmp_mat,c_vec,sig_mat);

	//cout<<"C_VEC:\n";
	//for(int i=0; i<c_vec.size(); i++){
	//	cout<<c_vec[i]<<"\t";
	//}
	//cout<<"\n";
	
//**egval can be negative**//

	vector<double> lmbda0, lambda, minusplus(egval.size(),1);
	for(int i=0; i<egval.size(); i++){
		lambda.push_back(sqrt(abs(egval[i])));
		lmbda0.push_back(sqrt(c_vec[i]));
		if(egval[i]<0) minusplus[i]=-1;
	}

	for(int i=0; i<sig_mat.size(); i++){
		for(int j=0; j<sig_mat[i].size(); j++){
			sig_mat[i][j]*=lambda[i]*lmbda0[j];
		}
	}

	for(int i=0; i<sig_mat.size(); i++){
		for(int j=0; j<sig_mat.size(); j++){
			double tmpsum=0;
			for(int k=0; k<sig_mat.size(); k++){
				tmpsum+= sig_mat[k][j]*sig_mat[k][i]*minusplus[k];
			}
			tmp_mat[i][j]=tmpsum;
			//cout<<tmpsum<<"\t";
		}
		//cout<<"\n";
	}

	//Mat_fuc::svdcmp(tmp_mat,c_vec,sig_mat);
	Mat_fuc::Symmeig smg2(tmp_mat,false);
	c_vec=smg2.d;

	//TRL::smg(tmp_mat,c_vec,sig_mat);

	//cout<<"C_VEC:\n";
	//for(int i=0; i<c_vec.size(); i++){
	//	cout<<c_vec[i]<<"\t";
	//}
	//cout<<"\n";

	//system("pause");
/**/
	double adj=0;
	for(int i=0; i<egval.size(); i++){
		adj+=egval[i];
	}
	vector<double> trace; int fault;
	double pvalue = QF::qf(c_vec,weightU+adj,trace,fault);
	//double pvalue = QF::qf(egval,weightU+adj,trace,fault);

	_weightU.push_back(weightU);

	return pvalue;

}

double HWU::hwu_inter_U(vector<double> & x)
{
	/*
	double tmpsum=0,tmp2sum=0,tmpsum2=0;
	double tmpwt=0, tmpwt2=0, tmpg=0;

	double diagsum=0, sumA=0, sumA2=0, sum2A=0;

	double weightU=0;
	*/

	//using float to speed up
	float tmpsum=0,tmp2sum=0,tmpsum2=0;
	float tmpwt=0, tmpwt2=0, tmpg=0;

	float diagsum=0, sumA=0, sumA2=0, sum2A=0;

	float weightU=0;

	for(int i=0; i<_n_sub; i++){
		if(x[i]!=0){
			for(int j=i+1; j<_n_sub; j++){
				if(x[j]!=0){
					if(_weight_flag) tmpwt=_weight->get(i,j);
					else tmpwt=1;

					if(tmpwt!=0){
						tmpg=x[i]*x[j];
						tmpwt*=tmpg;

						weightU+=(_outer_y->get(i,j)*tmpwt);

					}
				}
			}
		}
		//cout<<weightU<<"\n";
	}

	weightU*=2.0;

	return weightU;
}

double HWU::hwu_inter_p(vector<double> & x)
{
	/*
	double tmpsum=0,tmp2sum=0,tmpsum2=0;
	double tmpwt=0, tmpwt2=0, tmpg=0;

	double diagsum=0, sumA=0, sumA2=0, sum2A=0;

	double weightU=0;
	*/

	//using float to speed up
	float tmpsum=0,tmp2sum=0,tmpsum2=0;
	float tmpwt=0, tmpwt2=0, tmpg=0;

	float diagsum=0, sumA=0, sumA2=0, sum2A=0;

	float weightU=0;

	for(int i=0; i<_n_sub; i++){
		tmpsum=0; tmp2sum=0;
		if(x[i]!=0){
			for(int j=0; j<_n_sub; j++){
				if(x[j]!=0){
					if(_weight_flag) tmpwt=_weight->get(i,j);
					else tmpwt=1;

					if(tmpwt!=0){
						tmpg=x[i]*x[j];
						tmpwt*=tmpg;
						tmpwt2=tmpwt*tmpwt;

						tmpsum+=tmpwt;
						tmp2sum+=tmpwt2;

						if(i==j) diagsum+=tmpwt;
						else{
							weightU+=(_outer_y->get(i,j)*tmpwt);
						}

					}
				}
			}

			sumA+=tmpsum;
			sumA2+=tmp2sum;
			sum2A+=(tmpsum*tmpsum);

		}

		//cout<<weightU<<"\n";

	}

	double expect=0, variance=0;

	expect=diagsum-sumA/double(_n_sub);

	variance=sumA2 - 2*sum2A/double(_n_sub) + sumA*sumA/double(_n_sub)/double(_n_sub);


	variance*=2.0;

	//get observed value
	double obs=weightU + diagsum;

	//cout<<"\t"<<obs<<"\t"<<variance<<"\t"<<expect<<"\n";

	//_weightU.push_back(weightU);
	//_Obs.push_back(obs);
	//_expct.push_back(expect);
	//_varian.push_back(variance);

	if(variance<=0) return 1.0;
	if(obs<=0) return 1.0;

	double s_=variance/expect;
	double a_=expect/s_;

	double pvalue=1.0 - Stat_fuc::P_gamma(obs,a_,s_);

	return pvalue;
}

double HWU::hwu_i_internal(vector<double> & x)
{
	/*
	double tmpsum=0,tmp2sum=0,tmpsum2=0;
	double tmpwt=0, tmpwt2=0, tmpg=0;

	double diagsum=0, sumA=0, sumA2=0, sum2A=0;

	double weightU=0;
	*/

	//using float to speed up
	float tmpsum=0,tmp2sum=0,tmpsum2=0;
	float tmpwt=0, tmpwt2=0, tmpg=0;

	float diagsum=0, sumA=0, sumA2=0, sum2A=0;

	float weightU=0;

	for(int i=0; i<_n_sub; i++){
		tmpsum=0; tmp2sum=0;
		if(x[i]!=0){
			for(int j=0; j<_n_sub; j++){
				if(x[j]!=0){
					if(_weight_flag) tmpwt=_weight->get(i,j);
					else tmpwt=1;

					if(tmpwt!=0){
						tmpg=x[i]*x[j];
						tmpwt*=tmpg;
						tmpwt2=tmpwt*tmpwt;

						tmpsum+=tmpwt;
						tmp2sum+=tmpwt2;

						if(i==j) diagsum+=tmpwt;
						else{
							weightU+=(_outer_y->get(i,j)*tmpwt);
						}

					}
				}
			}

			sumA+=tmpsum;
			sumA2+=tmp2sum;
			sum2A+=(tmpsum*tmpsum);

		}

		//cout<<weightU<<"\n";

	}

	double expect=0, variance=0;

	expect=diagsum-sumA/double(_n_sub);

	variance=sumA2 - 2*sum2A/double(_n_sub) + sumA*sumA/double(_n_sub)/double(_n_sub);


	variance*=2.0;

	//get observed value
	double obs=weightU + diagsum;

	//cout<<"\t"<<obs<<"\t"<<variance<<"\t"<<expect<<"\n";

	_weightU.push_back(weightU);
	_Obs.push_back(obs);
	_expct.push_back(expect);
	_varian.push_back(variance);

	if(variance<=0) return 1.0;
	if(obs<=0) return 1.0;

	double s_=variance/expect;
	double a_=expect/s_;

	double pvalue=1.0 - Stat_fuc::P_gamma(obs,a_,s_);

	return pvalue;
}

double HWU::hwu_i_first_PC(vector<double> & x)
{
	/*
	double tmpsum=0,tmp2sum=0,tmpsum2=0;
	double tmpwt=0, tmpwt2=0, tmpg=0;

	double diagsum=0, sumA=0, sumA2=0, sum2A=0;

	double weightU=0;
	*/

	//using float to speed up
	float tmpsum=0,tmp2sum=0,tmpsum2=0;
	float tmpwt=0, tmpwt2=0, tmpg=0;

	float diagsum=0, sumA=0, sumA2=0, sum2A=0;

	float weightU=0;

	for(int i=0; i<_n_sub; i++){
		tmpsum=0; tmp2sum=0;
		if(x[i]!=0){
			for(int j=0; j<_n_sub; j++){
				if(x[j]!=0){
					if(_weight_flag) tmpwt=_weight->get(i,j);
					else tmpwt=1;

					if(tmpwt!=0){
						tmpg=x[i]*x[j];
						tmpwt*=tmpg;
						tmpwt2=tmpwt*tmpwt;

						tmpsum+=tmpwt*_first_PC[j];//this is the difference
						tmp2sum+=tmpwt2;

						if(i==j) diagsum+=tmpwt;
						else{
							weightU+=(_outer_y->get(i,j)*tmpwt);
						}

					}
				}
			}

			sumA+=tmpsum*_first_PC[i];//this is the difference
			sumA2+=tmp2sum;
			sum2A+=(tmpsum*tmpsum);

		}

		//cout<<weightU<<"\n";

	}

	double expect=0, variance=0;

	expect=diagsum-sumA;

	variance=sumA2 - 2*sum2A + sumA*sumA;


	variance*=2.0;


	//get observed value
	double obs=weightU + diagsum;

	//cout<<"\t"<<obs<<"\t"<<variance<<"\t"<<expect<<"\n";

	_weightU.push_back(weightU);
	_Obs.push_back(obs);
	_expct.push_back(expect);
	_varian.push_back(variance);

	if(variance<=0) return 1.0;
	if(obs<=0) return 1.0;

	double s_=variance/expect;
	double a_=expect/s_;

	double pvalue=1.0 - Stat_fuc::P_gamma(obs,a_,s_);

	return pvalue;
}


void HWU::princ_comp(vector< vector<double> > & X)
{
	//first allocate, then delete, then assign. important!!!!!!

	/*
	double ** tmp_buf=new double* [_n_sub];
	for(int i=0; i<_n_sub; i++){
		double * tmptmp=new double[_n_sub];
		tmp_buf[i]=tmptmp;
	}

	for(int i=0; i<_n_sub; i++){
		for(int j=i; j<_n_sub; j++){
			_buffer[i*_n_sub+j]=_weight->get(i,j);
			_buffer[j*_n_sub+i]=_buffer[i*_n_sub+j];
		}
	}
	*/

	_buffer=_weight->getPointer();

	//get eigen value and eigen vectors
	vector<double> egval;
	vector< vector<double> > egvec;
	//TRL::eg_nuTRan(TRL::mt_op,ot_X,egval,egvec,10,0.05,1,30,60,10);

	//Dmat Dmt(ot_X);
	//TRL::eg_nuTRan(TRL::mt_op2,&Dmt,ot_X.size(),egval,egvec,10,0.05,1,30,60,10);
	//TRL::eg_nuTRan(TRL::mt_op2,_outer_x,_n_sub,egval,egvec,10,0.05,1,30,60,10);
	//TRL::eg_nuTRan(TRL::mt_op3,_buffer,_n_sub,egval,egvec,10,0.05,1,30,60,10);

	TRL::eg_nuTRan(TRL::mt_op3,_buffer,_n_sub,egval,egvec,
		par::trl_ned,par::trl_pres,par::trl_scheme,100,200,par::trl_max);

	_buffer=0;

	X=egvec;
}

void HWU::princ_comp(vector< vector<double> > & X, vector< vector<double> > & mat)
{
	//first allocate, then delete, then assign. important!!!!!!
	/*
	double * tmp_buf=new double[_n_sub*_n_sub];
	delete [] _buffer;
	_buffer=tmp_buf;

	for(int i=0; i<_n_sub; i++){
		for(int j=i; j<_n_sub; j++){
			_buffer[i*_n_sub+j]=mat[i][j];
			_buffer[j*_n_sub+i]=_buffer[i*_n_sub+j];
		}
	}
	*/

	_buffer=_weight->getPointer();


	//get eigen value and eigen vectors
	vector<double> egval;
	vector< vector<double> > egvec;
	//TRL::eg_nuTRan(TRL::mt_op,ot_X,egval,egvec,10,0.05,1,30,60,10);

	//Dmat Dmt(ot_X);
	//TRL::eg_nuTRan(TRL::mt_op2,&Dmt,ot_X.size(),egval,egvec,10,0.05,1,30,60,10);
	//TRL::eg_nuTRan(TRL::mt_op2,_outer_x,_n_sub,egval,egvec,10,0.05,1,30,60,10);
	//TRL::eg_nuTRan(TRL::mt_op3,_buffer,_n_sub,egval,egvec,10,0.05,1,30,60,10);
	TRL::eg_nuTRan(TRL::mt_op3,_buffer,_n_sub,egval,egvec,
		par::trl_ned,par::trl_pres,par::trl_scheme,100,200,par::trl_max);

	_buffer=0;

	X=egvec;
}

void HWU::pjY()
{
	//cout<<"PC\n"; system("pause");

	gfun::printLOG("\nProjection...");

	vector< vector<double> > X;
	princ_comp(X);
	vector<double> e;
	vector<double> one(_n_sub,1.0);
	
	vector< vector<double> > Xtmp;
	double pval=0;

	//_y=sry(e); do not rank, otherwise not good

	//can do standardize

	//cout<<"PJ\n"; system("pause");

	double miu, sd;


	if(par::pj_Y_fwd){

		vector<double>  beta; vector<bool>  significant;

#ifdef DEBUG_WEI_TRL
		cout<<"Adding the ";
#endif

		for(int i=0; i<X.size(); i++){
			Xtmp.push_back(X[i]);
			//Dist::projection(_y,e,Xtmp);
			Dist::projectionSolu(_y,e,beta,significant,Xtmp);

			if(!significant[significant.size()-1]){
				Xtmp.pop_back();
			}else{
				Dist::meansd(e,miu,sd);
				//cout<<"\nmiu:"<<miu<<"\tsd:"<<sd<<"\n";

				for(int j=0; j<e.size(); j++){
					if(sd>0) e[j]=(e[j]-miu)/sd;
					else e[j]=e[j]-miu;
				}

				pval=hwu_test(one,e);

#ifdef DEBUG_WEI_TRL
				cout<<i+1<<"\t";
				cout<<"(pvalue:"<<pval<<")\t";
#endif

				if(pval>par::pj_threshold){
					break;
				}

			}
			
		}

#ifdef DEBUG_WEI_TRL
		cout<<"th PCs\n";
#endif

	}else{

#ifdef DEBUG_WEI_TRL
		cout<<"Adding the ";
#endif

		for(int i=0; i<X.size(); i++){
			Xtmp.push_back(X[i]);
			Dist::projection(_y,e,Xtmp);

			Dist::meansd(e,miu,sd);
			//cout<<"\nmiu:"<<miu<<"\tsd:"<<sd<<"\n";

			for(int j=0; j<e.size(); j++){
				if(sd>0) e[j]=(e[j]-miu)/sd;
				else e[j]=e[j]-miu;
			}

			pval=hwu_test(one,e);

#ifdef DEBUG_WEI_TRL
			cout<<i+1<<"\t";
			cout<<"(pvalue:"<<pval<<")\t";
#endif

			if(pval>par::pj_threshold){
				break;
			}
		}

#ifdef DEBUG_WEI_TRL
		cout<<"th PCs\n";
#endif

	}

	gfun::printLOG("Done.(P="+ double2str(pval) +")\n");

	if(par::pj_rst){
		gfun::printLOG("Output the new phenotype and the PCs\n");
		
		Mat_fuc::printVec(par::pj_f_Y, e);
		Mat_fuc::printMat(par::pj_f_PC, Xtmp);

	}


	if(Xtmp.size()==0) e=_y;

	_y=e;

	for(int i=0; i<_y.size(); i++){
		for(int j=i; j<_y.size(); j++){
			double tmp=_y[i]*_y[j];
			_outer_y->assign(i,j,tmp);
			_outer_y->assign(j,i,tmp);
		}
	}

	//cout<<"PJ done\n"; system("pause");
}

void HWU::projectY()
{
	vector< vector<double> > X;
	princ_comp(X);
	vector<double> e;

	vector<int> idxxx;
	if(par::pjct_OneStep){
		idxxx=Dist::projectionOnestep(_y,e,X);
		cout<<"Choose "<<idxxx.size()<<" PCs\n";

#ifdef DEBUG_WEI_TRL
		cout<<"They are ";
		for(int i_=0; i_<idxxx.size(); i_++){
			cout<<idxxx[i_]+1<<"\t";
		}
		cout<<" th PCs\n";
#endif


	}else if(par::pjct_Backward){
		idxxx=Dist::projectionBackward(_y,e,X);
		cout<<"Choose "<<idxxx.size()<<" PCs\n";

#ifdef DEBUG_WEI_TRL
		cout<<"They are ";
		for(int i_=0; i_<idxxx.size(); i_++){
			cout<<idxxx[i_]+1<<"\t";
		}
		cout<<" th PCs\n";
#endif


	}else{
		Dist::projection(_y,e,X);
	}
	
	
	//_y=sry(e); do not rank, otherwise not good

	//can do standardize
	double miu, sd;

	Dist::meansd(e,miu,sd);
	//cout<<"\nmiu:"<<miu<<"\tsd:"<<sd<<"\n";

	for(int i=0; i<e.size(); i++){
		if(sd>0) e[i]=(e[i]-miu)/sd;
		else e[i]=e[i]-miu;
	}

	_y=e;

	for(int i=0; i<_y.size(); i++){
		for(int j=i; j<_y.size(); j++){
			double tmp=_y[i]*_y[j];
			_outer_y->assign(i,j,tmp);
			_outer_y->assign(j,i,tmp);
		}
	}
}

void HWU::projectY2()
{
	vector< vector<double> > X;
	princ_comp(X);
	vector<double> e;

	_first_PC=X[0];
	X.clear(); X.push_back(_first_PC);
	Dist::projection(_y,e,X);

	//_y=sry(e); do not rank, otherwise not good

	//do standardize, really?
	//double miu, sd;

	//Dist::meansd(e,miu,sd);

	//cout<<"\nmiu:"<<miu<<"\tsd:"<<sd<<"\n";

	//for(int i=0; i<e.size(); i++){
	//	if(sd>0) e[i]=(e[i]-miu)/sd;
	//	else e[i]=e[i]-miu;
	//}

	_y=e;

	for(int i=0; i<_y.size(); i++){
		for(int j=i; j<_y.size(); j++){
			double tmp=_y[i]*_y[j];
			_outer_y->assign(i,j,tmp);
			_outer_y->assign(j,i,tmp);
		}
	}
}

void HWU::updateWeight()
{
	cout<<"Update weight matrix\n";
	rdAddCov(par::add_cov_f);

	cout<<"Projecting Y\n";

	vector< vector<double> > mat;
	mat=Dist::euclid(_add_cov);

	if(_weight_flag){
		for(int i=0; i<_n_sub; i++){
			for(int j=i; j<_n_sub; j++){
				double tmp=mat[i][j]*_weight->get(i,j);
				_weight->assign(i,j,tmp);
				_weight->assign(j,i,tmp);
			}
		}
	}else{
		_weight->copy(mat);
		_weight_flag=true;
	}
	cout<<"Done Updating\n";
	
}
/*
void HWU::projectY3()
{
	rdAddCov(par::pj_cov_f);

	cout<<"Projecting Y\n";

	vector< vector<double> > mat;
	mat=Dist::euclid(_add_cov);

	vector< vector<double> > X0,X1;
	princ_comp(X0);
	princ_comp(X1,mat);

	for(int i=0; i<X1.size(); i++){
		X0.push_back(X1[i]);
	}

	projectionMat(X0,mat);

	vector<double> e;

	projectY(mat,e,_y);

	//do standardize, really?
	double miu, sd;

	Dist::meansd(e,miu,sd);

	//cout<<"\nmiu:"<<miu<<"\tsd:"<<sd<<"\n";

	for(int i=0; i<e.size(); i++){
		if(sd>0) e[i]=(e[i]-miu)/sd;
		else e[i]=e[i]-miu;
	}

	_y=e;

	for(int i=0; i<_y.size(); i++){
		for(int j=i; j<_y.size(); j++){
			double tmp=_y[i]*_y[j];
			_outer_y->assign(i,j,tmp);
			_outer_y->assign(j,i,tmp);
		}
	}

	cout<<"Projection Done\n";

}
*/

void HWU::rdAddCov(string filename)
{
	gfun::printLOG("Reading Additional Covariates in [ " 
		+ filename + " ] \n");

	gfun::checkFileExists( filename );

	ifstream COV;
	COV.open(filename.c_str());
	COV.clear();

	string str_buf;
	getline(COV,str_buf);
	int Ncov=gfun::split_string(str_buf,_cov_name);

	vector<double> tmp_cov(_n_sub,0);
	_add_cov.clear(); _add_cov.resize(Ncov,tmp_cov);

	int c=0;
	while(!COV.eof())
	{

		for(int i=0; i<Ncov; i++){
			COV >> _add_cov[i][c];
		}

		c++;

		//cout<<c<<"\n";

		if(c==_n_sub) break;
	}

	COV.clear();
	COV.close();

	gfun::printLOG(int2str(Ncov)+" Covariates in [ " 
		+ filename + " ] \n");
}

void HWU::projectionMat(vector< vector<double> > & X, vector< vector<double> > & P)
{
	//calculate I-X%*%solve(t(X)%*%X)%*%t(X)

	//use X1=[1,X];
	vector<double> tmp;
	if(X.size()==0) tmp.resize(_n_sub,1);
	else tmp.resize(X[0].size(),1);

	P.resize(tmp.size(), tmp);
	vector< vector<double> > X1;
	
	X1.push_back(tmp);
	for(int i=0; i<X.size(); i++){
		X1.push_back(X[i]);
	}

	//calculate XIns=solve(t(X1)%*%X1)
	tmp.resize(X1.size(),0);
	vector< vector<double> > XX(X1.size(),tmp), XIns(X1.size(),tmp);

	for(int i=0; i<X1.size(); i++){
		for(int j=i; j<X1.size(); j++){
			double tmpsum=0;
			for(int k=0; k<X1[0].size();k++){
				tmpsum+=X1[i][k]*X1[j][k];
			}
			XX[i][j]=tmpsum;
			XX[j][i]=tmpsum;
		}
	}

	//we can use TRLan maybe
	double rk=Mat_fuc::SVD_inverse1(XX,XIns);
	if(rk<=0 || rk>X1.size()) gfun::error("error in construct Projection Matrix\n");
	if(rk!=XX.size()){
		cout<<"Rank:"<<rk<<"\tNcol:"<<XX.size()<<"\n";
	}

	//use XX as slot
	XX=X1;
	for(int i=0; i<XX.size(); i++){
		for(int j=0; j<XX[0].size(); j++){
			double tmpsum=0;
			for(int k=0; k<XX.size(); k++){
				tmpsum+=X1[k][j]*XIns[k][i];
			}
			XX[i][j]=tmpsum;
		}
	}

	for(int i=0; i<XX[0].size(); i++){
		for(int j=i; j<XX[0].size(); j++){
			double tmpsum=0;
			for(int k=0; k<XX.size(); k++){
				tmpsum+=XX[k][i]*X1[k][j];
			}
			
			if(i==j){
				P[i][j]=1-tmpsum;
			}else{
				P[i][j]=0-tmpsum;
				P[j][i]=0-tmpsum;
			}
		}
	}
	
}

void HWU::projectY(vector< vector<double> > &P, vector<double> & newY, vector<double> & oldY)
{
	newY=oldY;
	for(int i=0; i<P.size(); i++){
		double tmpsum=0;
		for(int j=0; j<P.size(); j++){
			tmpsum+=P[i][j]*_y[i];
		}
		newY[i]=tmpsum;
	}

}

void HWU::prepProject(vector< vector<double> > & P)
{
	//already know _org_y
	_Projection=new Dmat(P);

	//get rank
	vector<double> rk;

	if(par::hwu_rk){
		rk=getrk(_org_y);
	}else{
		rk=_org_y;
	}

	_y.resize(rk.size());
	
	//project rank
	for(int i=0; i<P.size(); i++){
		double tmpsum=0;
		for(int j=0; j<P.size(); j++){
			tmpsum+=P[i][j]*rk[j];
		}
		_y[i]=tmpsum;
	}
	
	//check miu is 0, then devided by sd.
	double miu, sd;
	Dist::meansd(_y,miu,sd);

	if(miu>1e-2) gfun::error("error in Projecting trait values");
	else{
		if(sd!=0) sd=1.0/sd;

		for(int i=0; i<_y.size(); i++){
			_y[i]*=sd;
		}
	}

	//update _outer_y
	for(int i=0; i<_y.size(); i++){
		for(int j=i; j<_y.size(); j++){
			double tmp=_y[i]*_y[j];
			_outer_y->assign(i,j,tmp);
			_outer_y->assign(j,i,tmp);
		}
	}

	//initialze PW
	vector<float> tmpvec(_n_sub,0);
	_PW.resize(_n_sub,tmpvec);

}