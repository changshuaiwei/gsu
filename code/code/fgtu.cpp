#include "gtu.h"

void FGTU::FamStruc()
{
	kinship ksp;
	ksp.initialize(_datafile);
	ksp.ksMAT(_KSmat);

	//cout<<_KSmat.block(520,520,10,10)<<"\n";

	//eigen dicomposition of the kinship matrix
	Eigen::SelfAdjointEigenSolver<GTmat_> es(_n_sub);

	es.compute(_KSmat);
	_KSEigenValues=es.eigenvalues();
	//cout<<_KSEigenValues - es.eigenvalues()<<"\n";

	//cout<<es.eigenvectors().col(1)<<"\n";
	//cout<<es.eigenvectors().cols()<<"\n";
	//cout<<es.eigenvectors().rows()<<"\n";
	//cout<<_n_sub<<"\n";

	_KSEigenVectors.resize(_n_sub,_n_sub);
	_KSEigenVectors=es.eigenvectors();
	
	//cout<<_KSEigenVectors.col(1)<<"\n";
	_weightKS=1*GTmat_::Constant(_n_sub,_n_sub,0);
	for(int i=0; i<_n_sub; i++){
		if(_KSEigenValues(i,0)>0){
			_weightKS(i,i)=1/sqrt(_KSEigenValues(i,0));
		}else{
			_weightKS(i,i)=0;
		}
	}

	_weightKS=_KSEigenVectors*_weightKS;


}

void FGTU::PJ_Tsim()
{
	_Tsim=_KSEigenVectors*_Tsim*_KSEigenVectors.transpose();
	GTmat_ newX;
	newX=_KSEigenVectors*_One;
	standardize(_Tsim, newX);

	for(int i=0; i<_Tsim.rows(); i++) _Tsim(i,i)=0;

	//newX=_KSEigenVectors*_X;
	standardize(_Tsim, _One);
}

void FGTU::PJ_Gsim()
{
	_Gsim=_weightKS* _Gsim * _weightKS.transpose();
	GTmat_ newX;
	newX=_KSEigenVectors*_One;

	standardize(_Gsim,newX);

	//cout<<"\n\nBefore:\n"<<_Gsim.block(0,0,5,5);
	for(int i=0; i<_Gsim.rows(); i++) _Gsim(i,i)=0;

	//
	standardize(_Gsim,_One);

	//cout<<"\n\nBefore:\n"<<_Gsim.block(0,0,5,5);
}

void FGTU::transX()
{
	int n=_phe_ori[0].size();
	_X.resize(n,_cov_ori.size()+1);

	vector<double> tmp;
	tmp.resize(n,1);
	for(int j=0; j<tmp.size(); j++){
		_X(j,0)=tmp[j];
	}

	if(_cov_ori.size()>0){
		GTmat_ X;
		X.resize(n,_cov_ori.size());
		for(int i=0; i<_cov_ori.size(); i++){
			tmp=_cov_ori[i];
			for(int j=0; j<tmp.size(); j++){
				X(j,i)=tmp[j];
			}
		}
		_X.block(0,1,n,_cov_ori.size()) = _KSEigenVectors * X;
	}
	
}