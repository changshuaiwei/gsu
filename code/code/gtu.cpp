#include "gtu.h"

void GTU::initialize(snpdt *snp_data)
{
	_datafile=snp_data;

	_datafile->getCovMat(_cov_ori);
	_datafile->getPheMat(_phe_ori);
	_datafile->getSnpSet(_snp_set);
	
	_n_sub=_phe_ori[0].size();

	FamStruc();

	if(_n_sub>par::use_trl_threshold){
		_buf=new Dmat(_n_sub,_n_sub);
	}

	normalizeY();
	
	if(par::read_pheno_wt){
		_datafile->getWtVec(_phe_wt);
		_Yweight.resize(_phe_wt.size());
		for(int i=0; i<_phe_wt.size(); i++){
			_Yweight(i)=_phe_wt[i];
		}
	}else{
		_Yweight=1*GTvec_::Constant(_Y.cols(),1);
	}
	

	
	_One=1*GTmat_::Constant(_Y.rows(),1,1);

	if(par::lap_kernel){
		getTraitSimLap(_Y,_Tsim,_Yweight);
	}else{
		getTraitSimSqr(_Y,_Tsim,_Yweight);
	}

//	cout<<"\n\nBefore:\n"<<_Tsim.block(0,0,5,5);

	PJ_Tsim();
}

void GTU::FamStruc()
{

}

void GTU::PJ_Tsim()
{
	//cout<<"Before:\n"<<_Tsim.block(0,0,5,5);
	standardize(_Tsim);

	//cout<<"\nAfter:\n"<<_Tsim.block(0,0,5,5);

	for(int i=0; i<_Tsim.rows(); i++) _Tsim(i,i)=0;

	transX();
	standardizeX(_Tsim);

	//get _X; covariate
}

void GTU::run()
{
	cout<<"\n";
	for(int i=0; i<_snp_set.size(); i++){
		
		cout<<i+1<<"th SNP: Genetic Similarity....";

		bool LK=par::lap_kernel;

		bool skip=wIBS(_snp_set[i], LK);

		cout<<"U statistic and p-value....";

		if(skip){
			_QVec.push_back(0);
			_pVec.push_back(1);

			_Liu_df.push_back(0);
			_Liu_ncp.push_back(0);
			_Liu_q.push_back(0);

		}else{
			GTUcore();
		}
		

		cout<<"\r";
		cout.flush();
	}
}

void GTU::wtResult(string outputfile)
{
	gfun::printLOG("Writing GTU Result at [ " + outputfile + " ]\n");

	ofstream result(outputfile.c_str());
	if(!result) gfun::error("\ncould not open the result file\n");
	result<<"#result of SNP-set association scanning\n";

	result<<"#tName\tnSNP\t"
		<<"Q\tP_value\t"
		<<"Liu_df\tLiu_ncp\tLiu_q\n";

	for(int i=0; i<_snp_set.size(); i++){
		result<<_datafile->getSnpSetNameI(i)<<"\t"
			<<_snp_set[i].size()<<"\t"
			<<_QVec[i]<<"\t"
			<<_pVec[i]<<"\t"
			<<_Liu_df[i]<<"\t"
			<<_Liu_ncp[i]<<"\t"
			<<_Liu_q[i]<<"\n";

	}

	result.close();
}

void GTU::normalizeY()
{

	if(par::gtu_rank){

		_Y.resize(_phe_ori[0].size(),_phe_ori.size());


		
		for(int i=0; i<_phe_ori.size(); i++){
			vector<double> tmp;
			Stat_fuc::unifQuantile(_phe_ori[i],tmp);
			//Stat_fuc::normalQuantile(_phe_ori[i],tmp);

			for(int j=0; j<tmp.size(); j++){
				_Y(j,i)=tmp[j];
			}

		}

	}else{

		GTmat_ Y;
		Y.resize(_phe_ori[0].size(),_phe_ori.size());

		for(int i=0; i<_phe_ori.size(); i++){
			for(int j=0; j<_phe_ori[i].size(); j++){
				Y(j,i)=_phe_ori[i][j];
			}
		}

		

		//cout<<Y;


		GTmat_ centered =Y.rowwise() -Y.colwise().mean();

		//cout<<centered;

		_Y=centered.colwise().normalized() * sqrt(_phe_ori[0].size()*1.0);
	}

	

	//cout<<_Y;
}

void GTU::transX()
{
	int n=_phe_ori[0].size();
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
}

void GTU::getTraitSimSqr(GTmat_ & Y, GTmat_ & Tsim, GTvec_ & weight)
{
	Tsim.resize(Y.rows(),Y.rows());
	weight=weight/(weight.sum());
	
	for(int i=0; i<Y.rows(); i++){
		for(int j=i; j<Y.rows(); j++){
			double tmp_sum=0, tmp=0;
			for(int k=0; k<Y.cols(); k++){
				tmp=Y.coeff(i,k)-Y.coeff(j,k);
				tmp_sum-=weight(k)*tmp*tmp;
			}
			Tsim(i,j)=exp(tmp_sum);
			Tsim(j,i)=Tsim(i,j);
		}
	}
}

void GTU::getTraitSimLap(GTmat_ & Y, GTmat_ & Tsim, GTvec_ & weight)
{
	Tsim.resize(Y.rows(),Y.rows());
	weight=weight/(weight.sum());

	for(int i=0; i<Y.rows(); i++){
		for(int j=i; j<Y.rows(); j++){
			double tmp_sum=0, tmp=0;
			for(int k=0; k<Y.cols(); k++){
				tmp=abs( Y.coeff(i,k)-Y.coeff(j,k) );
				tmp_sum-=weight(k)*tmp;
			}

			//cout<<exp(tmp_sum)<<"\t";

			Tsim(i,j)=exp(tmp_sum);
			Tsim(j,i)=Tsim(i,j);
		}
	}
	//cout<<"\n";
}

void GTU::getTraitSimCov(GTmat_ & Y, GTmat_ & Tsim, GTvec_ & weight)
{
	weight=weight/(weight.sum());

	GTmat_ tmp;
	tmp=Y.transpose();
	for(int i=0; i<tmp.rows(); i++){
		tmp.row(i)*=weight(i);
	}

	Tsim=Y*tmp;
}

void GTU::standardize(GTmat_ & a)
{
	double inv=1.0/double(_Y.rows());

	a=a-(_One*inv)*(_One.transpose()*a);
	a=a-(a*_One)*(inv*_One.transpose());
}

void GTU::standardize(GTmat_ & a, GTmat_ & X)
{
	GTmat_ inv=X.transpose()*X;
	inv=inv.inverse();

	a=a-(X*inv)*(X.transpose()*a);
	a=a-(a*X)*(inv*X.transpose());
}

void GTU::standardizeX(GTmat_ & a)
{
	GTmat_ inv=_X.transpose()*_X;
	inv=inv.inverse();

	a=a-(_X*inv)*(_X.transpose()*a);
	a=a-(a*_X)*(inv*_X.transpose());
}

bool GTU::wIBS(vector<int> & idx, bool LK)
{
	vector<double> p, weight;
	GTimat_ Z(_n_sub,idx.size());
	int i_, tmp_i, count;
	double tmp_d;
	vector<int> variantIdx;
	vector<int> tmp_vi;

	///attain the data and record the missing position
	for(int i=0; i<idx.size(); i++){
		i_=idx[i];
		tmp_vi.clear();
		tmp_d=0;
		for(int j=0; j<_n_sub; j++){
			tmp_i=_datafile->genotypeToInt(j,i_);

			//if(boost::math::isnan(tmp_i)) tmp_i=-9;

			if(tmp_i==-9){
				tmp_vi.push_back(j);
				Z(j,i)=0;
			}else{
				Z(j,i)=tmp_i;
			}
		}
		
		//calculate p
		if(_n_sub>tmp_vi.size()){
			tmp_d=double(Z.col(i).sum())/double(_n_sub-tmp_vi.size())/2.0;
			if(tmp_d>0.5){
				tmp_d=1-tmp_d;
				for(int j=0; j<_n_sub; j++) Z(j,i)=2-Z.coeff(j,i);
			}
		}
		

		if(tmp_d>0){
			///impute
			for(int j=0; j<tmp_vi.size(); j++){
				tmp_i=0;
				if(tmp_d>Stat_fuc::ran1(par::seed)){
					tmp_i++;
				}
				if(tmp_d>Stat_fuc::ran1(par::seed)){
					tmp_i++;
				}
				Z(tmp_vi[j],i)=tmp_i;
			}


			p.push_back(tmp_d);

			if(tmp_d>0) weight.push_back(sqrt(1.0/tmp_d));

			variantIdx.push_back(i);
		}

	}


	if(variantIdx.size()==0) return true; 

	double wtotal=0;
	for(int i=0; i<weight.size(); i++) wtotal+=weight[i];

	
	///get weighted IBS
	_Gsim.resize(_n_sub,_n_sub);
	for(int i=0; i<_n_sub; i++){
		for(int j=i+1; j<_n_sub; j++){
			tmp_d=0;

			if(variantIdx.size()>0){
				for(int k=0; k<variantIdx.size(); k++){
					tmp_d += weight[k] * abs ( Z.coeff(i,variantIdx[k]) - Z.coeff(j,variantIdx[k]) );
				}
				tmp_d=1.0-tmp_d/2.0/wtotal;
			}else{
				tmp_d=1.0;
			}

			
			if(LK){
				tmp_d= exp(2*tmp_d-2);
			}
			
			
			_Gsim(i,j)=tmp_d;
			_Gsim(j,i)=tmp_d;
		}

		_Gsim(i,i)=1.0;
	}

	///centered and corrected

	//cout<<"\n\nBefore:\n"<<_Gsim.block(0,0,5,5);

	//cout<<"\n\nBefore:\n"<<_Tsim.block(0,0,5,5);



	double tempa=_Gsim.sum();
	if(boost::math::isnan(tempa)){
		gfun::printLOG("\nThe Genetic similarity Matrix have NAN values!!! check the data!!!\n");
	}

	PJ_Gsim();

	return false;
}

void GTU::PJ_Gsim()
{
	standardize(_Gsim);

	//cout<<"\n\nBefore:\n"<<_Gsim.block(0,0,5,5);
	for(int i=0; i<_Gsim.rows(); i++) _Gsim(i,i)=0;
	
	if(par::gtu_Gproj){
		standardizeX(_Gsim);
	}else{
		standardize(_Gsim);
	}
	
	//standardize(_Gsim);

	//cout<<"\n\nBefore:\n"<<_Gsim.block(0,0,5,5);
}

void GTU::GTUcore()
{
	double Q=0;
	for(int i=0; i<_n_sub; i++){
		for(int j=0; j<_n_sub; j++){
			Q+=_Gsim.coeff(i,j)*_Tsim.coeff(i,j);
		}
	}

	GTmat_ Coef;

	if(_n_sub>par::use_trl_threshold){
		EigenLarge(Coef);
	}else{
		EigenSmall(Coef);
	}

	double p_value=Liu(Coef,Q);

	_QVec.push_back(Q);
	_pVec.push_back(p_value);

}

void GTU::EigenSmall(GTmat_ & a)
{
	GTmat_ Tcoef, Gcoef;

	double factor=1.0/double(_n_sub-_X.cols());

	Eigen::SelfAdjointEigenSolver<GTmat_> es(_n_sub);

	es.compute(_Tsim,Eigen::EigenvaluesOnly);
	Tcoef=es.eigenvalues();

	es.compute(_Gsim,Eigen::EigenvaluesOnly);
	Gcoef=es.eigenvalues();

	a=(Tcoef*factor)*Gcoef.transpose();
}

void GTU::EigenLarge(GTmat_ & a)
{
	GTmat_ Tcoef, Gcoef, Coef;

	double factor=1.0/double(_n_sub-_X.cols());

	EigenLargeCore(_Tsim,Tcoef);

	//cout<<"\n\n"<<Tcoef<<"\n\n";

	EigenLargeCore(_Gsim,Gcoef);

	//cout<<"\n\n"<<Gcoef<<"\n\n";

	a=(Tcoef*factor)*Gcoef.transpose();
}

void GTU::EigenLargeCore(GTmat_ & a,GTmat_ & eigenvalues)
{
	for(int i=0; i<a.rows(); i++){
		for(int j=0; j<a.cols(); j++){
			_buf->assign(i,j,a(i,j));
		}
	}

	_bufPointer=_buf->getPointer();

	//get eigen value and eigen vectors
	vector<double> egval;
	vector< vector<double> > egvec;

	TRL::eg_nuTRan(TRL::mt_op3,_bufPointer,_n_sub,egval,egvec,
		par::gtu_trl_dim,par::trl_pres,par::trl_scheme,1000,2000,par::trl_max);

	_bufPointer=0;

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
}


double GTU::Liu(GTmat_ & a, double Q)
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

	_Liu_q.push_back(q);
	_Liu_df.push_back(df);
	_Liu_ncp.push_back(ncp);

	double pp=ncChiSqSurvival(df,ncp,q);

	return pp;
}
double GTU::ncChiSqSurvival(double df, double ncp, double q)
{
	if(q<0) q=0;

	double p_higher_tail=1.0 - cdf(boost::math::non_central_chi_squared(df, ncp), q);

	return p_higher_tail;
}

