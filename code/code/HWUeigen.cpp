#include "hwu.h"

void egHWU::initialize(snpdt * snp_data)
{
	_datafile=snp_data;
	_n_sub=_datafile->totalIndi();
	_n_snp=_datafile->totalLocus();
	
	_datafile->getPheVec(_org_y);
	_datafile->getCovMat(_cov_ori);

	rankY();
	_One=1*GTmat_::Constant(_Y.rows(),1,1);

	if(_n_sub>par::use_trl_threshold){
		_weight=new Dmat(_n_sub,_n_sub);
	}

	_Kappa.resize(_n_sub,_n_sub);
	_Kappa=_One * _One.transpose();
}

void egHWU::initialize_min(snpdt * snp_data)
{
	_datafile=snp_data;
	_n_sub=_datafile->totalIndi();
	_n_snp=_datafile->totalLocus();

	
	_One=1*GTmat_::Constant(_Y.rows(),1,1);

	_Kappa.resize(_n_sub,_n_sub);
	_Kappa=_One * _One.transpose();
}

void egHWU::rankY()
{
	vector<double> ytmp;
	ytmp=_org_y;
	Stat_fuc::crank(ytmp);

	_Y.resize(_n_sub,1);
	for(int i=0; i<_n_sub; i++){
		_Y(i,0)=ytmp[i];
	}

}

void egHWU::transX()
{
	int n=_n_sub;
	_X.resize(n,_cov_ori.size()+1);

	vector<double> tmp;
	tmp.resize(n,1);
	for(int j=0; j<tmp.size(); j++){
		_X(j,0)=tmp[j];
	}

	for(int i=0; i<_cov_ori.size(); i++){
		tmp=_cov_ori[i];
		for(int j=0; j<tmp.size(); j++){
			_X(j,i+1)=tmp[j];
		}
	}

	GTmat_ inv=_X.transpose()*_X;
	_XXinv=inv.inverse();
}

void egHWU::prepareY()
{
	//cout<<_Y<<"\n";
	_Y=_Y - (_X*_XXinv) * (_X.transpose() * _Y);

	//cout<<_Y<<"\n";
	double Yse= _Y.squaredNorm() / double (_n_sub - _X.cols()) ;

	//cout<<Yse<<"\n";
	if(Yse>1e-10){
		_Y = _Y / sqrt(Yse);
	}else{
		_Y = _Y * 0.0;
	}
	
	//cout<<_Y<<"\n";
}

void egHWU::sortEigen(GTmat_ & egvec, GTmat_ & egval)
{
	vector<double> val;
	vector<int> idx;
	for(int i=0; i<egval.rows(); i++){
		val.push_back(abs(egval(i,0)));
	}

	Stat_fuc::indexx(val,idx);

	GTmat_ Egvc, Egvl;
	Egvc=egvec; Egvl=egval;

	for(int i=0; i<val.size() ; i++){
		int tpidx=idx.size()-i-1;
		egval.row(i)=Egvl.row(idx[tpidx]);
		egvec.col(i)=Egvc.col(idx[tpidx]);
	}
}

void egHWU::stdW()
{
	if(!_weight_flag) return;
	if(!par::std_weight) return;

	cout<<"\nstandardizing weight matrix...\n";//to set for testing pure heterogeneity effect

	double inv=1.0/double(_Y.rows());

	_Kappa=_Kappa-(_One*inv)*(_One.transpose()*_Kappa);
	_Kappa=_Kappa-(_Kappa*_One)*(inv*_One.transpose());

	cout<<"Done.\n";
}


void egHWU::sgLocusAssoc()
{
	transX();
	prepareY();

	//cout<<_Y.sum()<<"\n";
	//cout<<_Y.squaredNorm()<<"\n";


	stdW();

	_vec_pvalue.clear();
	double pvalue=1;

	vector<double> x;
	cout<<"\n";
	for(int i=0; i<_n_snp; i++){
		cout<<i+1<<"th SNP";
		
		x.clear();
		//cout<<":"<<(_datafile->getLocus(i))->name;
		prepareX(i,x);

		pvalue=hwu_liu(x);

		_vec_pvalue.push_back(pvalue);
		cout<<"\r";
		cout.flush();
	}

}

void egHWU::weightPC(int size, string outputfile)
{
	GTmat_ PC;
	PCEigen(_Kappa,PC,size);

	gfun::printLOG("Writing PC for weight Result at [ " + outputfile + " ]\n");

	ofstream result(outputfile.c_str());
	if(!result) gfun::error("\ncould not open the result file\n");

	for(int i=0; i<PC.cols(); i++){
		result<<"PC"<<i+1<<"\t";
	}
	result<<"\n";

	for(int i=0; i<PC.rows(); i++){
		for(int j=0; j<PC.cols(); j++){
			result<<PC(i,j)<<"\t";
		}
		result<<"\n";
	}
	result<<"\n";

	//vector< vector<double> > cov_combine;
	vector<double> tmpvec;

	for(int i=0; i<PC.cols(); i++){
		for(int j=0; j<PC.rows(); j++){
			tmpvec.push_back(PC(j,i));
		}
		_cov_ori.push_back(tmpvec);
		tmpvec.clear();
	}

}

void egHWU::wtZscore(string outputfile)
{
	gfun::printLOG("Writing Zscore Result at [ " + outputfile + " ]\n");

	ofstream result(outputfile.c_str());
	if(!result) gfun::error("\ncould not open the result file\n");
	result<<"#result of Zscore scanning (0=A1A1, 1=A1A2, 2=A2A2)\n";

	result<<"#order\t"<<"Chr\tName\tPos\tBp\tAllel1\tAllel2\t"
		<<"Z_score\n";

	Locus * loc=0;
	for(int i=0; i<_bk_vec_pvalue.size(); i++){
			loc=_datafile->getLocus(i);


		result<<i+1<<"\t";//order

		result	<<loc->chr << "\t"
			<< loc->name << "\t"
			<< loc->pos << "\t"
			<< loc->bp << "\t"
			<< loc->allele1 <<"\t"
			<< loc->allele2 <<"\t";//snp information
		result
			<<(0.0-_bk_vec_pvalue[i])<<"\n";

	}

	result.close();
}

void egHWU::wtResult(string outputfile)
{
	gfun::printLOG("Writing HWU Result at [ " + outputfile + " ]\n");

	ofstream result(outputfile.c_str());
	if(!result) gfun::error("\ncould not open the result file\n");
	result<<"#result of association scanning (0=A1A1, 1=A1A2, 2=A2A2)\n";

	//if(par::appx_davis || par::appx_davis_P){ 
	//	result<<"#order\t"<<"Chr\tName\tPos\tBp\tAllel1\tAllel2\t"
			//<<"U\t"
	//		<<"P_value\n";
	//}else{
	//	result<<"#order\t"<<"Chr\tName\tPos\tBp\tAllel1\tAllel2\t"
			//<<"U\t"
	//		<<"P_value\n";
	//}

	result<<"#order\t"<<"Chr\tName\tPos\tBp\tAllel1\tAllel2\t"
		<<"Chi_df\t"<<"Chi_ncp\t"<<"Chi_q\t"
		<<"Z_score\t"
		<<"P_value\n";

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

		//if(par::appx_davis || par::appx_davis_P){
		//	result	//<<_weightU[i]<<"\t"
		//		<<_vec_pvalue[i]<<"\n";//pvalue
		//}else{
		//	result	//<<_weightU[i]<<"\t"
		//		<<_vec_pvalue[i]<<"\n";//pvalue
		//}
		result	<<_df_vec[i]<<"\t"<<_ncp_vec[i]<<"\t"<<_q_vec[i]<<"\t"
			<<_weightU[i]<<"\t"
					<<_vec_pvalue[i]<<"\n";

	}

	result.close();

	if(par::appx_davis || par::appx_davis_P){
		string outfile=outputfile+".zscore.txt";
		wtZscore(outfile);
	}
}


void egHWU::sgLocusRank()
{
	transX();
	prepareY();
	stdW();

	if(par::multi_core){
		D2F();
		getRankZfltMp();
	}else{
		if(!par::hwu_flt){
			getRankZ();
		}else{
			D2F();
			getRankZfloat();
		}
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

		pvalue=hwu_liu(x);
		_vec_pvalue.push_back(pvalue);
		cout<<"\r";
		cout.flush();
	}

}

void egHWU::getRankZ()
{
	cout<<"\nGet pseudo Z score for\n";
	vector<double> x; double p;
	_vec_pvalue.clear();
	for(int i=0; i<_n_snp; i++){
		cout<<i+1<<"th SNP";
		x.clear();
		prepareX(i,x);
		p=0.0 - hwu_inter_Z(x);
		_vec_pvalue.push_back(p);
		cout<<"\r";
		cout.flush();
	}

	cout<<"\nSorting pseudo Z\n";
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

void egHWU::D2F()
{
	_Yflt.resize(_Y.rows(),_Y.cols());
	for(int i=0; i<_Y.rows(); i++){
		for(int j=0;j<_Y.cols();j++){
			_Yflt(i,j)=float( _Y(i,j) );
		}
	}

	_Xflt.resize(_X.rows(),_X.cols());
	for(int i=0; i<_X.rows(); i++){
		for(int j=0;j<_X.cols();j++){
			_Xflt(i,j)=float( _X(i,j) );
		}
	}

	_XXinvflt.resize(_XXinv.rows(),_XXinv.cols());
	for(int i=0; i<_XXinv.rows(); i++){
		for(int j=0;j<_XXinv.cols();j++){
			_XXinvflt(i,j)=float( _XXinv(i,j) );
		}
	}

	_Kappaflt.resize(_Kappa.rows(),_Kappa.cols());
	for(int i=0; i<_Kappa.rows(); i++){
		for(int j=0;j<_Kappa.cols();j++){
			_Kappaflt(i,j)=float( _Kappa(i,j) );
		}
	}
}

void egHWU::getRankZfloat()
{
	cout<<"\nGet pseudo Z score for\n";
	vector<float> x; double p;
	_vec_pvalue.clear();
	for(int i=0; i<_n_snp; i++){
		cout<<i+1<<"th SNP";
		x.clear();
		prepareXfloat(i,x);
		p=0.0 - hwu_inter_Zfloat(x);
		_vec_pvalue.push_back(p);
		cout<<"\r";
		cout.flush();
	}

	cout<<"\nSorting pseudo Z\n";
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

void egHWU::getRankZfltMp()
{

	int n_td=boost::thread::hardware_concurrency()-1;

	gfun::printLOG("There are " + int2str(n_td + 1) + " cores, and we require " + int2str(par::n_core) + " cores\n");


	if(par::force_core){
		n_td=par::n_core;//might exceed memory when number is large.
	}else{
		if(n_td>par::n_core)	n_td=par::n_core;
		gfun::printLOG("Using " + int2str(n_td) + " cores instead of " + int2str(par::n_core) + " cores\n");
	}

	_vec_pvalue.resize(_n_snp);
	_SNPcounter=0;

	vector<int> idx_start, idx_length;
	splitCal(n_td,idx_start,idx_length);

	boost::thread_group threads;

	for(int i=0; i<idx_start.size(); i++){
		threads.add_thread(new boost::thread(&egHWU::calZMp, this, idx_start[i], idx_length[i]));

		//threads.create_thread( boost::bind (&egHWU::calZMp, boost::ref(this), boost::ref(idx_start[i]), boost::ref(idx_length[i]) ) );
	}

	threads.join_all();

	cout<<"\nSorting pseudo Z\n";
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
/*
void egHWU::calZMpp(int start, int length, GTmatF_ & Wt, GTmatF_ & Xf, GTmatF_ & Yf, GTmatF_ & XXinvf)
{


	//GTmatF_ Wt; GTmatF_ Yf; GTmatF_ Xf; GTmatF_ XXinvf;
	//Yf.resize(_n_sub,1); Xf.resize(_n_sub,_Xflt.cols()); XXinvf.resize(_Xflt.cols(),_Xflt.cols()); Wt.resize(_Kappaflt.rows(),_Kappaflt.cols());
	//Yf=_Yflt; Xf=_Xflt; XXinvf=_XXinvflt; Wt=_Kappaflt;

	vector<float> x; double p;
	int i_=0;
	for(int i=0; i<length; i++){
		i_=start + i;

		x.clear();

		{
			//boost::unique_lock<boost::shared_mutex> lock(_mutex);
			boost::shared_lock<boost::shared_mutex> lock(_mutex);
			prepareXfloat(i_,x);

			p=0.0 - hwu_inter_ZfltMP(x,Wt,Xf,Yf,XXinvf);

			_vec_pvalue[i_]=p;
			lock.unlock();
		}

		{
			boost::unique_lock<boost::shared_mutex> lock(_mutex2);
			_SNPcounter++;
			cout<<"\r"<<_SNPcounter<<"th SNP";
			cout.flush();
			lock.unlock();
		}
	}
}

*/
void egHWU::calZMp(int start, int length)
{

	vector<float> x; double p;
	int i_=0;
	for(int i=0; i<length; i++){
		i_=start + i;

		x.clear();

		{
			boost::shared_lock<boost::shared_mutex> lock(_mutex);
			prepareXfloat(i_,x);
			p=0.0 - hwu_inter_Zfloat(x);
			_vec_pvalue[i_]=p;
			lock.unlock();
		}
		
		{
			boost::unique_lock<boost::shared_mutex> lock(_mutex2);
			_SNPcounter++;
			cout<<"\r"<<_SNPcounter<<"th SNP";
			cout.flush();
			lock.unlock();
		}

	}
}

void egHWU::splitCal(int n_td, vector<int> & idx_start, vector<int> & idx_length)
{
	int length = int (double (_n_snp) / double (n_td)) ;

	idx_start.clear(); idx_length.clear();

	for(int i=0; i<n_td; i++){
		int tmps=i*length;
		int tmpl=length;

		if(i < n_td-1 ){
			idx_start.push_back(tmps);
			idx_length.push_back(tmpl);
		}else{
			idx_start.push_back(tmps);
			idx_length.push_back(_n_snp-tmps);
		}
	}
}

void egHWU::prepareXfloat(int i_snp, vector<float> &x)
{
	//if(i_snp==157) system("pause");

	x.clear();
	vector<double> x_;
	vector<int> x_rcd;
	for(int i=0; i<_n_sub; i++){
		int geno=(_datafile->genotypeToInt(i,i_snp));
		x.push_back(float(geno));
		if(geno==-9){
			x_rcd.push_back(i);
		}else{
			x_.push_back(float(geno));
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

double egHWU::hwu_inter(vector<double> & x, GTmat_ & Wt)
{
	Wt.resize(_n_sub,_n_sub);

	for(int i=0; i<_n_sub; i++){
		for(int j=i; j<_n_sub; j++){
			if(i==j){
				Wt(i,j)=0;
			}else{
				Wt(i,j)=_Kappa.coeff(i,j) * gsim(x[i] , x[j]);
				Wt(j,i)=Wt.coeff(i,j);
			}
		}
	}

	//cout<<_Y.sum()<<"\n";
	//cout<<_Y.squaredNorm()<<"\n";

	double U= (_Y.transpose() * (Wt * _Y) ).sum();

	Wt=Wt-(_X*_XXinv)*(_X.transpose()*Wt);
	Wt=Wt-(Wt*_X)*(_XXinv*_X.transpose());

	return U;
}

inline double egHWU::gsim(double g1, double g2)
{
	if(par::gsim_add){
		return g1*g2;
	}else{
		return exp(0.0-(g1-g2)*(g1-g2));
	}
}

inline float egHWU::gsim(float g1, float g2)
{
	if(par::gsim_add){
		return g1*g2;
	}else{
		return exp(0.0-(g1-g2)*(g1-g2));
	}
}

double egHWU::hwu_interFloat(vector<float> & x, GTmatF_ & Wt)
{
	Wt.resize(_n_sub,_n_sub);

	for(int i=0; i<_n_sub; i++){
		for(int j=i; j<_n_sub; j++){
			if(i==j){
				Wt(i,j)=0;
			}else{
				Wt(i,j)=_Kappaflt.coeff(i,j) * gsim ( x[i] , x[j] );
				Wt(j,i)=Wt.coeff(i,j);
			}
		}
	}

	//cout<<_Y.sum()<<"\n";
	//cout<<_Y.squaredNorm()<<"\n";

	double U= (_Yflt.transpose() * (Wt * _Yflt) ).sum();
	

	Wt=Wt-(_Xflt*_XXinvflt)*(_Xflt.transpose()*Wt);
	Wt=Wt-(Wt*_Xflt)*(_XXinvflt*_Xflt.transpose());

	return U;
}

double egHWU::hwu_inter_p(vector<double> & x)
{
	GTmat_ Wt;

	double U=hwu_inter(x,Wt);

	double mu=Wt.trace();
	double var=Wt.squaredNorm();

	double Z= (U-mu) /sqrt(var);

	double p_value=Stat_fuc::std_norm_p1(Z);
	return p_value;
}

double egHWU::hwu_inter_Z(vector<double> & x)
{
	GTmat_ Wt;

	//analyze time

	double U=hwu_inter(x,Wt);

	double mu=Wt.trace();
	double var=2 * Wt.squaredNorm();

	double Z= (U-mu) / sqrt(var);

	return Z;
}
/*
double egHWU::hwu_inter_ZfltMP(vector<float> & x, GTmatF_ & Wtf, GTmatF_ & Xf, GTmatF_ & Yf, GTmatF_ & XXinvf)
{
	GTmatF_ Wt;
	Wt.resize(_n_sub,_n_sub);
	for(int i=0; i<_n_sub; i++){
		for(int j=i; j<_n_sub; j++){
			if(i==j){
				Wt(i,j)=0;
			}else{
				Wt(i,j)=Wtf.coeff(i,j) * x[i] * x[j];
				Wt(j,i)=Wt.coeff(i,j);
			}
		}
	}

	//cout<<_Y.sum()<<"\n";
	//cout<<_Y.squaredNorm()<<"\n";

	double U= (Yf.transpose() * (Wt *Yf) ).sum();


	Wt=Wt-(Xf*XXinvf)*(Xf.transpose()*Wt);
	Wt=Wt-(Wt*Xf)*(XXinvf*Xf.transpose());


	double mu=Wt.trace();
	double var=2 * Wt.squaredNorm();

	double Z= (U-mu) / sqrt(var);

	return Z;

}
*/
double egHWU::hwu_inter_Zfloat(vector<float> & x)
{
	GTmatF_ Wt;

	double U=hwu_interFloat(x,Wt);

	double mu=Wt.trace();
	double var=2 * Wt.squaredNorm();

	double Z= (U-mu) / sqrt(var);

	return Z;
}

double egHWU::hwu_liu(vector<double> & x)
{
	GTmat_ Wt;

	double U=hwu_inter(x,Wt);

	GTmat_ Coef;

	if(_n_sub>par::use_trl_threshold){
		EigenLargeCore(Wt,Coef);

		//cout<<Coef.block(0,0,100,1)<<"\n";

		//cout<<Wt.trace()<<"\n";
		//cout<<Coef.sum()<<"\n";
		//cout<<Wt.squaredNorm()<<"\n";
		//cout<<Coef.squaredNorm()<<"\n";
	}else{
		EigenSmallCore(Wt,Coef);

		//cout<<Coef.block(0,0,300,1)<<"\n";

		//cout<<Wt.trace()<<"\n";
		//cout<<Coef.sum()<<"\n";
	}

	double df, ncp, q;
	double p_value=Liu(Coef,U,df,ncp,q);
	_df_vec.push_back(df);
	_ncp_vec.push_back(ncp);
	_q_vec.push_back(q);

	double mu=Wt.trace();
	double var=2 * Wt.squaredNorm();
	_weightU.push_back((U-mu)/sqrt(var));
	
	//_expct.push_back(mu);
	//_varian.push_back(var);

	return p_value;
}

double egHWU::Liu(GTmat_ & a, double Q, double & df_, double & ncp_, double & q_)
{
	double c1=0,c2=0,c3=0,c4=0;
	double tmp_c1=0,tmp_c2=0,tmp_c3=0,tmp_c4=0;

	for(int i=0; i<a.rows(); i++){
		for(int j=0; j<a.cols(); j++){
			tmp_c1=a(i,j);
			tmp_c2=tmp_c1*tmp_c1;
			tmp_c3=tmp_c2*tmp_c1;
			tmp_c4=tmp_c3*tmp_c1;

			c1+=tmp_c1;
			c2+=tmp_c2;
			c3+=tmp_c3;
			c4+=tmp_c4;
		}
	}

	double s1=0,s2=0,miuQ=0,sigmaQ=0,miuX=0,sigmaX=0,tstar=0,aa=0,delta=0,ll=0;

	s1 = c3/(pow(c2,1.5));
	s2 = c4/ (pow(c2,2.0));
	miuQ = c1;
	sigmaQ = sqrt(2*c2);
	tstar = (Q-miuQ)/sigmaQ;
	if(s1*s1 > s2){
		//cout<<"t1\t";
		aa=1.0/(s1 - sqrt( s1*s1 - s2) );
		delta= s1 * pow(aa, 3) - aa*aa;
		ll=aa*aa - 2 * delta;
	}else{
		//cout<<"t2\t";
		aa=1.0/s1;
		delta=0;
		ll=pow(c2,3)/pow(c3,2);
	}
	miuX=ll + delta;
	sigmaX= sqrt(2.0) * aa;

	double df, ncp, q;
	df=ll; ncp=delta; q=tstar*sigmaX + miuX;

	//cout<<df<<"\t"<<ncp<<"\t"<<q<<"\t";
	df_=df; ncp_=ncp; q_=q;

	double pp=ncChiSqSurvival(df,ncp,q);

	return pp;
}

double egHWU::Liu(GTmat_ & a, double Q)
{
	double c1=0,c2=0,c3=0,c4=0;
	double tmp_c1=0,tmp_c2=0,tmp_c3=0,tmp_c4=0;

	for(int i=0; i<a.rows(); i++){
		for(int j=0; j<a.cols(); j++){
			tmp_c1=a(i,j);
			tmp_c2=tmp_c1*tmp_c1;
			tmp_c3=tmp_c2*tmp_c1;
			tmp_c4=tmp_c3*tmp_c1;

			c1+=tmp_c1;
			c2+=tmp_c2;
			c3+=tmp_c3;
			c4+=tmp_c4;
		}
	}

	double s1=0,s2=0,miuQ=0,sigmaQ=0,miuX=0,sigmaX=0,tstar=0,aa=0,delta=0,ll=0;

	s1 = c3/(pow(c2,1.5));
	s2 = c4/ (pow(c2,2.0));
	miuQ = c1;
	sigmaQ = sqrt(2*c2);
	tstar = (Q-miuQ)/sigmaQ;
	if(s1*s1 > s2){
		//cout<<"t1\t";
		aa=1.0/(s1 - sqrt( s1*s1 - s2) );
		delta= s1 * pow(aa, 3) - aa*aa;
		ll=aa*aa - 2 * delta;
	}else{
		//cout<<"t2\t";
		aa=1.0/s1;
		delta=0;
		ll=pow(c2,3)/pow(c3,2);
	}
	miuX=ll + delta;
	sigmaX= sqrt(2.0) * aa;

	double df, ncp, q;
	df=ll; ncp=delta; q=tstar*sigmaX + miuX;

	//cout<<df<<"\t"<<ncp<<"\t"<<q<<"\t";

	double pp=ncChiSqSurvival(df,ncp,q);

	return pp;
}
double egHWU::ncChiSqSurvival(double df, double ncp, double q)
{
	if(q<0) q=0;

	double p_higher_tail=1.0 - cdf(boost::math::non_central_chi_squared(df, ncp), q);

	return p_higher_tail;
}


void egHWU::EigenSmallCore(GTmat_ & a,GTmat_ & eigenvalues)
{
	Eigen::SelfAdjointEigenSolver<GTmat_> es(_n_sub);

	es.compute(a,Eigen::EigenvaluesOnly);
	eigenvalues.resize(_n_sub,1);
	eigenvalues=es.eigenvalues();
}

void egHWU::EigenSmallCore(GTmat_ & a,GTmat_ & eigenvectors, GTmat_ & eigenvalues)
{
	Eigen::SelfAdjointEigenSolver<GTmat_> es(_n_sub);

	es.compute(a);
	eigenvalues.resize(_n_sub,1);
	eigenvalues=es.eigenvalues();
	eigenvectors.resize(a.rows(), a.rows());
	eigenvectors=es.eigenvectors();

	//cout<<eigenvectors.colwise().squaredNorm()<<"\n";
	//cout<<eigenvectors.transpose().colwise().squaredNorm()<<"\n";
	//cout<<eigenvectors * eigenvectors.transpose()<<"\n";
	//cout<<eigenvectors.colwise().mean()<<"\n";
	//cout<<eigenvectors.rowwise().mean()<<"\n";
}

void egHWU::PCEigen(GTmat_ & mat, GTmat_ & PC, int size)
{
	GTmat_ egvec, egval;

	//cout<<mat.rows()<<"\t"<<mat.cols()<<"\n";

	//cout<<mat.block(0,0,4,4);

	if(_n_sub>par::use_trl_threshold){
		EigenLargeCore(mat,egvec,egval);
	}else{
		EigenSmallCore(mat,egvec,egval);
	}

	//cout<<egval.block(0,0,100,1)<<"\n";

	//cout<<egvec.block(0,0,10,10)<<"\n";

	//cout<<mat.trace()<<"\t"<<egval.sum()<<"\t";

	sortEigen(egvec,egval);

	//cout<<egval<<"\n";
	//cout<<egvec.colwise().squaredNorm()<<"\n";


	if(size > egvec.cols()){
		cout<<"cumputed eigen-dim ("<<egvec.cols()<<") < demanded sizes (" << size <<"), use all the computed\n";
		PC=egvec;
	}else{
		cout<<"cumputed eigen-dim ("<<egvec.cols()<<") > demanded sizes (" << size <<"), use demanded\n";
		PC=egvec.block(0,0,mat.rows(),size);
	}

	//cout<<PC.colwise().mean()<<"\n";
	//cout<<PC.colwise().squaredNorm()<<"\n";

}

void egHWU::EigenLargeCore(GTmat_ & a, GTmat_ & eigenvectors, GTmat_ & eigenvalues)
{
	for(int i=0; i<a.rows(); i++){
		for(int j=0; j<a.cols(); j++){
			_weight->assign(i,j,a(i,j));
		}
	}

	_buffer=_weight->getPointer();

	//get eigen value and eigen vectors
	vector<double> egval;
	vector< vector<double> > egvec;

	TRL::eg_nuTRan(TRL::mt_op3,_buffer,_n_sub,egval,egvec,
		par::gtu_trl_dim,par::trl_pres,par::trl_scheme,500,1000,par::trl_max);

	//TRL::eg_nuTRan(TRL::mt_op3,_buffer,_n_sub,egval,egvec,
	//	par::gtu_trl_dim,par::trl_pres,par::trl_scheme,1000,2000,par::trl_max);

	_buffer=0;

	eigenvalues.resize(egval.size(),1);
	for(int i=0; i<egval.size(); i++){
		eigenvalues(i,0)=egval[i];
	}

	eigenvectors.resize(a.rows(), egval.size());
	for(int i=0; i<egval.size(); i++){
		for(int j=0; j<a.rows(); j++){
			eigenvectors(j,i)=egvec[i][j];
		}
	}
}

void egHWU::EigenLargeCore(GTmat_ & a,GTmat_ & eigenvalues)
{
	for(int i=0; i<a.rows(); i++){
		for(int j=0; j<a.cols(); j++){
			_weight->assign(i,j,a(i,j));
		}
	}

	_buffer=_weight->getPointer();

	//get eigen value and eigen vectors
	vector<double> egval;
	vector< vector<double> > egvec;

	TRL::eg_nuTRan(TRL::mt_op3,_buffer,_n_sub,egval,egvec,
		par::gtu_trl_dim,par::trl_pres,par::trl_scheme,500,1000,par::trl_max);

	//TRL::eg_nuTRan(TRL::mt_op3,_buffer,_n_sub,egval,egvec,
	//	par::gtu_trl_dim,par::trl_pres,par::trl_scheme,1000,2000,par::trl_max);

	_buffer=0;

	double sumeig=0, sumeigsq=0;
	for(int i=0; i<egval.size(); i++){
		sumeig+=egval[i];
		sumeigsq+=egval[i]*egval[i];
	}

	int nn=  (a.rows() - egval.size() );
	double mud=  ( a.trace() - sumeig );
	double vrd= ( a.squaredNorm() - sumeigsq ) ;

	vector<double> psedoEigen;
	Stat_fuc::RanNormal(0,1,nn,par::seed,psedoEigen);
	sumeig=0; sumeigsq=0;
	for(int i=0; i<psedoEigen.size(); i++){
		sumeig+=psedoEigen[i]; sumeigsq+=psedoEigen[i]*psedoEigen[i];
	}

	if( vrd - mud*mud/double(nn) > 0){
		for(int i=0; i<psedoEigen.size(); i++){
			psedoEigen[i] = ( psedoEigen[i] - sumeig/double(nn) ) / sqrt(sumeigsq - sumeig*sumeig/double(nn));
			psedoEigen[i] = psedoEigen[i] * sqrt( vrd - mud*mud/double(nn) ) + mud/double(nn);
			egval.push_back(psedoEigen[i]);
		}
	}else{
		for(int i=0; i<psedoEigen.size(); i++){
			egval.push_back(mud/double(nn));
		}
	}
	
	eigenvalues.resize(egval.size(),1);
	for(int i=0; i<egval.size(); i++){
			eigenvalues(i,0)=egval[i];
	}
	

	/*
	eigenvalues.resize(a.rows(),1);
	for(int i=0; i<a.rows(); i++){
		if(i<egval.size()){
			eigenvalues(i,0)=egval[i];
		}else{
			eigenvalues(i,0)=correction;
		}
	}
	*/

}

void egHWU::EDweight(GTmat_ & cov, GTvec_ & weight)
{
	_Kappa.resize(cov.rows(),cov.rows());
	weight=weight/(weight.sum());

	//cout<<weight<<"\n";

	for(int i=0; i<cov.rows(); i++){
		for(int j=i; j<cov.rows(); j++){
			double tmp_sum=0, tmp=0;
			for(int k=0; k<cov.cols(); k++){
				tmp=cov.coeff(i,k)-cov.coeff(j,k);
				tmp_sum-=weight(k)*tmp*tmp;
			}
			//cout<<tmp_sum<<"\n";
			_Kappa(i,j)=exp(tmp_sum);
			_Kappa(j,i)=_Kappa(i,j);
		}
	}
}


void egHWU::EDweight(GTmat_ & cov, GTmat_ & Kappa, GTvec_ & weight)
{
	Kappa.resize(cov.rows(),cov.rows());
	weight=weight/(weight.sum());

	for(int i=0; i<cov.rows(); i++){
		for(int j=i; j<cov.rows(); j++){
			double tmp_sum=0, tmp=0;
			for(int k=0; k<cov.cols(); k++){
				tmp=cov.coeff(i,k)-cov.coeff(j,k);
				tmp_sum-=weight(k)*tmp*tmp;
			}
			Kappa(i,j)=exp(tmp_sum);
			Kappa(j,i)=Kappa(i,j);
		}
	}
}

void egHWU::StandMatByCol(GTmat_ &Mt)
{
	Mt.rowwise() -= Mt.colwise().mean();

	GTmat_ sd(1, Mt.cols());
	sd=Mt.colwise().norm()/(sqrt(double(Mt.rows()-1)));

	for(int i=0; i<Mt.cols(); i++){
		if(sd(0,i) > 0){
			Mt.col(i) = Mt.col(i)/sd(0,i);
		}
	}
}

void egHWU::weightFromCov()
{
	//cout<<"Relatedness From Covariates ..\n";
	gfun::printLOG("Relatedness From Covariates ..\n");

	_datafile->get2ndCovMat(_cov4wt_ori);

	GTmat_ Cov4Wt; Cov4Wt.resize(_n_sub,_cov4wt_ori.size());

	for(int i=0; i<_cov4wt_ori.size(); i++){
		for(int j=0; j<_n_sub; j++){
			Cov4Wt(j,i)=_cov4wt_ori[i][j];
		}
	}

	GTvec_ WtCov4Wt; 
	WtCov4Wt.resize(Cov4Wt.cols());
	WtCov4Wt=1 * GTvec_::Constant(Cov4Wt.cols(),1);

	//cout<<WtCov4Wt<<"\n";

	//cout<<WtCov4Wt.sum()<<"\n";

	//cout<<Cov4Wt.colwise().mean()<<"\n";

	StandMatByCol(Cov4Wt);

	//cout<<Cov4Wt.colwise().mean()<<"\n";
	//cout<<Cov4Wt.colwise().squaredNorm()<<"\n";

	egHWU::EDweight(Cov4Wt,WtCov4Wt);

	//cout<<_Kappa.colwise().mean()<<"\n";
	//cout<<( _Kappa < 0.9 ).all()<<"\n";
	_weight_flag=true;

	gfun::printLOG("\nRelateness Matrix Done.\n\n");
}

void egHWU::rdWeight(string inputfile)
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
			double tmptmp=tmp;
			weight[i][j]=tmptmp;
			weight[j][i]=weight[i][j];
			//_weight->assign(i,j,tmp);
			//_weight->assign(j,i,tmp);
		}
	}

	result.close();

	//assign to _weight
	_Kappa.resize(_n_sub,_n_sub);
	for(int i=0; i<weight.size(); i++){
		for(int j=0; j<_n_sub; j++){
			_Kappa(i,j)=weight[i][j];
		}
	}

	_weight_flag=true;

	cout<<"\nRelateness Matrix Done.\n\n";

}

void egHWU::weightFastGenome()
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




	//assign to _weight
	_Kappa.resize(_n_sub,_n_sub);
	for(int i=0; i<weightd.size(); i++){
		for(int j=0; j<_n_sub; j++){
			_Kappa(i,j)=weightd[i][j];
		}
	}

	_weight_flag=true;

	//cout<<"\nRelateness Matrix Done.\n\n";
	gfun::printLOG("\nRelateness Matrix Done.\n\n");
}

//calculate kinship matrix: weighted IBS(average correlation)
//negative value is possible
void egHWU::weightFromGenomeWIBS()
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


	//assign to _weight
	_Kappa.resize(_n_sub,_n_sub);
	for(int i=0; i<weightd.size(); i++){
		for(int j=0; j<_n_sub; j++){
			_Kappa(i,j)=weightd[i][j];
		}
	}

	_weight_flag=true;


	//cout<<"\nRelateness Matrix Done.\n\n";
	gfun::printLOG("\nRelateness Matrix Done.\n\n");

}

void egHWU::updateWeight()
{
	cout<<"Update weight matrix\n";
	rdAddCov(par::add_cov_f);

	GTmat_ Cov4Wt; Cov4Wt.resize(_n_sub,_add_cov.size());

	for(int i=0; i<_add_cov.size(); i++){
		for(int j=0; j<_n_sub; j++){
			Cov4Wt(j,i)=_add_cov[i][j];
		}
	}

	GTvec_ WtCov4Wt; 
	WtCov4Wt.resize(_add_cov.size());
	WtCov4Wt=1 * GTvec_::Constant(_add_cov.size(),1);

	StandMatByCol(Cov4Wt);

	GTmat_ Kappa;
	EDweight(Cov4Wt,Kappa,WtCov4Wt);

	if(_weight_flag){
		for(int i=0; i<_n_sub; i++){
			for(int j=i; j<_n_sub; j++){
				_Kappa(i,j) = _Kappa(i,j) * Kappa(i,j);
				_Kappa(j,i) = _Kappa(i,j);
			}
		}
	}else{
		//assign to _weight
		_Kappa.resize(_n_sub,_n_sub);
		_Kappa=Kappa;

		_weight_flag=true;
	}
	cout<<"Done Updating\n";

}

void egHWU::wtWeight(string outputfile)
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
			result<<_Kappa(i,j)<<"\t";
		}
		result<<"\n";
	}

	result.close();

}