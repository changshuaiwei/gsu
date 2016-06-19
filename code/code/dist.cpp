#include "dist.h"

void Dist::meansd(vector<double> & dt, double & miu, double & sd)
{
	double tmp=0;
	for (int i=0; i<dt.size(); i++)
	{
		tmp+=dt[i];
	}
	tmp/=double(dt.size());
	miu=tmp;

	tmp=0;
	for(int i=0; i<dt.size(); i++){
		tmp+=(dt[i]-miu)*(dt[i]-miu);
	}
	tmp/=double(dt.size()-1);
	sd=sqrt(tmp);

}

double Dist::euclid_AB(vector<double> & a, vector<double> & b)
{
	double sum=0;
	for(int i=0; i<a.size(); i++){
		sum+=(a[i]-b[i])*(a[i]-b[i]);
	}
	sum=exp(-sum/double(a.size()));
	return sum;
}

vector< vector<double> > Dist::euclid(vector< vector<double> > & dt)
{
	int i, j, k;

	vector< vector<double> > dt0;
	vector<double> miu, sd, tmp;
	double tmp1,tmp2;

	//standardize and transpose
	tmp.resize(dt.size(),0); dt0.resize(dt[0].size(), tmp);

	for(i=0; i<dt.size(); i++){
		meansd(dt[i],tmp1,tmp2);
		miu.push_back(tmp1); sd.push_back(tmp2);
		for(j=0; j<dt[i].size(); j++){
			if(tmp2!=0) dt0[j][i]=(dt[i][j]-tmp1)/tmp2;
			else dt0[j][i]=(dt[i][j]-tmp1);
		}
	}

	vector< vector<double> > ds;
	vector<double> tmp3(dt0.size(),0);
	for(i=0; i<tmp3.size(); i++) ds.push_back(tmp3);

	for(i=0; i<dt0.size(); i++){
		for(j=i; j<dt0.size(); j++){
			ds[i][j]=euclid_AB(dt0[i],dt0[j]);
			ds[j][i]=ds[i][j];
		}
	}
	return ds;
}

void Dist::projection(vector<double> Y, vector<double> & e, vector< vector<double> > & X)//X is eigen vector
{
	e.resize(Y.size(), 0);
	for(int i=0; i<Y.size(); i++){
		double sum=0;
		for(int k=0; k<Y.size(); k++){
			double tmpsum=0;
			for(int j=0; j<X.size(); j++){
				tmpsum+=X[j][i]*X[j][k];
			}
			sum+=Y[k]*tmpsum;
		}
		e[i]=Y[i]-sum;
	}
}

void Dist::projectionSolu(vector<double> Y, vector<double> & e, 
					vector<double> & beta, vector<bool> & significant, vector< vector<double> > & X)//require X is eigen vector
{
	e.resize(Y.size(), 0);
	beta.resize(X.size(),0);

	for(int i=0; i<X.size(); i++){
		double tmpsum=0;
		for(int j=0; j<X[i].size(); j++){
			tmpsum+=X[i][j]*Y[j];
		}
		beta[i]=tmpsum;
	}

	for(int i=0; i<Y.size(); i++){
		double tmpsum=0;
		for(int j=0; j<X.size(); j++){
			tmpsum+=beta[j]*X[j][i];
		}
		e[i]=Y[i]-tmpsum;
	}

	double beta_sig=0;
	for(int i=0; i<e.size(); i++){
		beta_sig+=e[i]*e[i];
	}

	beta_sig=sqrt(beta_sig/double(Y.size()-X.size()));

	significant.resize(X.size(),false);

	for(int i=0; i<X.size(); i++){
		double pvalue=2.0 * Stat_fuc::std_norm_p1(abs(beta[i]/beta_sig));
		if(pvalue<0.05) significant[i]=true;
		else significant[i]=false;
	}
}


vector<int> Dist::projectionOnestep(vector<double> Y, vector<double> & e, vector< vector<double> > & X)
{
	e.resize(Y.size(), 0);
	vector< vector<double> > Xtmp, Xtmp2;
	Xtmp=X;
	vector<double> beta(Xtmp.size(),0);
	vector<bool> significant(Xtmp.size(),false);

	vector<int> idx, idx2;

	for(int i=0; i<X.size(); i++) idx.push_back(i);

	projectionSolu(Y,e,beta,significant,Xtmp);

	Xtmp2.clear(); idx2.clear();
	for(int i=0; i<significant.size(); i++){
		if(significant[i]) {
			Xtmp2.push_back(Xtmp[i]);
			idx2.push_back(idx[i]);
		}
	}

	if(Xtmp2.size()==0){
		e=Y;
	}else{
		projection(Y,e,Xtmp2);
	}
	

	return idx2;
}

vector<int> Dist::projectionBackward(vector<double> Y, vector<double> & e, vector< vector<double> > & X)
{
	e.resize(Y.size(), 0);
	vector< vector<double> > Xtmp, Xtmp2;
	Xtmp=X;
	vector<double> beta(Xtmp.size(),0);
	vector<bool> significant(Xtmp.size(),false);

	vector<int> idx, idx2;
	//idx.push_back(0);
	for(int i=0; i<X.size(); i++) idx.push_back(i);
/*
	while(1){
		projectionSolu(Y,e,beta,significant,Xtmp);

		if(!significant[significant.size()-1]) {
			if(significant.size()<2){
				break;
			}else{
				Xtmp.pop_back();
				idx.pop_back();
				projection(Y,e,Xtmp);
				break;
			}
			
		}

		if(Xtmp.size()==X.size()) break;
		
		int i_=Xtmp.size();
		Xtmp.push_back(X[i_]);
		idx.push_back(i_);
	}
*/	

	
	while(1){
		projectionSolu(Y,e,beta,significant,Xtmp);

		Xtmp2.clear(); idx2.clear();
		for(int i=0; i<significant.size(); i++){
			if(significant[i]) {
				Xtmp2.push_back(Xtmp[i]);
				idx2.push_back(idx[i]);
			}
		}

		if(idx.size()==idx2.size()){
			if(idx2.size()==0){
				e=Y;
			}
			break;
		}

		idx=idx2;
		Xtmp=Xtmp2;

	}

	return idx2;
}


