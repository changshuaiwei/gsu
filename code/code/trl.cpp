#include "trl.h"

void TRL::smg(vector< vector<double> > mat, vector<double> & egval, vector< vector<double> > & egvec)
{
	eg_nuTRan(TRL::mt_op,mat,egval,egvec,
		par::trl_pres,par::trl_scheme,100,200,par::trl_max);
}

void TRL::sort_EGrst(int comp, int ned, int nrow, double * eval, double *evec, 
				vector<double> & egval, vector< vector<double> > & egvec, bool allned)
{
	egval.clear(); egvec.clear();
	vector<double> tmp1, temp2;
	vector<int> idx;
	for (int i=0; i<comp; i++)
	{
		tmp1.push_back(abs(eval[i]));
	}

	idx.resize(tmp1.size());
	TRL::indexx(tmp1,idx);


	int tmpidx=0;
	for(int i_=0;i_<ned; i_++){

		int i=idx.size()-i_-1;

		if(i<0) break;
		

		tmpidx=idx[i];
		if( (!allned) && (abs(eval[tmpidx])<=1e-12) ) break;


		tmp1.clear(); tmp1.resize(nrow);
		for(int j=0; j<nrow; j++){

			tmp1[j]=evec[tmpidx*nrow+j];
		}

		double tmpsum=0;
		for(int j=0; j<tmp1.size(); j++){
			tmpsum+=tmp1[j]*tmp1[j];
		}

		if( (!allned) && (tmpsum<1e-2) ) break;

		if( (allned) && (tmpsum<1e-1) ){
			for(int j=0; j<tmp1.size(); j++){
				tmp1[j]=0;
			}
			egvec.push_back(tmp1);
			egval.push_back(0);
			continue;
		}

		egvec.push_back(tmp1);
		egval.push_back(eval[tmpidx]);
	}

}
void TRL::mt_op(const int nrow, const int ncol, double *xin, const int ldx,
		   double *yout, const int ldy, void* mvparam)
{
	vector< vector<double> > * mt;
	mt=(vector< vector<double> > *) mvparam;
	//
	// ..
	// .. local variables ..
	int i, j;
	//
	// ..
	// .. executable statements ..
	for( j=0; j<ncol; j++ ) {
		for( i=0; i<nrow; i++ ) {
			double tmp=0;
			for(int k=0; k<nrow; k++){
				tmp += (*mt)[i][k]*xin[j*ldx+k];
			}
			yout[j*ldy+i] = tmp;
		}
	}
}

void TRL::mt_op2(const int nrow, const int ncol, const double *xin, const int ldx,
				double *yout, const int ldy, void* mvparam)
{
	Dmat* mt=(Dmat*) mvparam;
	//
	// ..
	// .. local variables ..
	int i, j;
	//
	// ..
	// .. executable statements ..
	for( j=0; j<ncol; j++ ) {
		for( i=0; i<nrow; i++ ) {
			double tmp=0;
			for(int k=0; k<nrow; k++){
				tmp += (mt->get(i,k))*xin[j*ldx+k];
			}
			yout[j*ldy+i] = tmp;
		}
	}
}

void TRL::mt_op3(const int nrow, const int ncol, double *xin, const int ldx,
				 double *yout, const int ldy, void* mvparam)
{
	double ** mt=(double **) mvparam;
	double * xin_buffer=xin;
	//
	// ..
	// .. local variables ..
	int i, j;
	//
	// ..
	// .. executable statements ..

	//cout<<"\t"<<ncol<<"\t"<<nrow<<"\t"<<ldx<<"\t"<<ldy<<"\n";

	for( j=0; j<ncol; j++ ) {
		for( i=0; i<nrow; i++ ) {
			
			double tmp=0;
			for(int k=0; k<nrow; k++){
				//tmp += (mt[i*nrow+k])*xin[j*ldx+k];
				tmp += (mt[i][k])*xin[j*ldx+k];
			}
			
			yout[j*ldy+i] = tmp;
			
			//yout[j*ldy+i] = TRL::trl_ddot(nrow,&mt[i*nrow],1,&xin_buffer[j*ldx],1);
		}
	}
	
}

void TRL::mt_op4(const int nrow, const int ncol, double *xin, const int ldx,
				 double *yout, const int ldy, void* mvparam)
{
	double * mt=(double *) mvparam;
	double * xin_buffer=xin;
	//
	// ..
	// .. local variables ..
	int i, j;
	//
	// ..
	// .. executable statements ..

	
	integer_ nrow_=nrow, ncol_=ncol;
	char para1='N';	char para2='N';	double para3=1.0;	double para4=0.0;
	TRL::dgemm_(&para1,&para2,&nrow_,&ncol_,&nrow_,&para3,mt,&nrow_,xin_buffer,&nrow_,&para4,yout,&nrow_);
}

/*
void TRL::mt_op2(const int nrow, const int ncol, const double *xin, const int ldx,
				double *yout, const int ldy, void* mvparam)
{
	const double * mt=(double *) mvparam;
	//
	// ..
	// .. local variables ..
	int i, j;
	//
	// ..
	// .. executable statements ..
	cout<<mt[0]<<"\t"<<mt[nrow+1]<<"\t"<<mt[nrow*nrow-1]<<"\t"<<mt[20000]<<"\n";
	for( j=0; j<ncol; j++ ) {
		for( i=0; i<nrow; i++ ) {
			double tmp=0, tmp2=0;
			for(int k=0; k<nrow; k++){
				
				tmp2 = (mt)[i+nrow+k]*xin[j*ldx+k];
				tmp=tmp+tmp2;
			}		
			yout[j*ldy+i] = tmp;
		}
	}
}

*/

void TRL::eg_nuTRan(trl_matvec op, vector< vector<double> > & mt,
					vector<double> & egval, vector< vector<double> > & egvec, 
					double pres ,const int lohi, int mev, int maxlan, int maxmv)
{
	int ned=mt.size();
	int nrow=mt.size();
	if(nrow<maxlan){
		maxlan=nrow;
		mev=nrow;
	}
	if(nrow<ned) ned=nrow;

	int lwrk;
	// local variable declaration

	//double eval[mev], evec[mev*nrow], exact[mev];
	double *eval=0, *evec=0, *exact=0;
	eval = (double*)malloc(mev*sizeof(double));
	evec = (double*)malloc(mev*nrow*sizeof(double));
	exact = (double*)malloc(mev*sizeof(double));


	double *res=0, *wrk=0;
	trl_info info;
	int i, j, k, fp, check;
	//char name2[133], name[150];
	int tmp1, tmp2, nlen;


	lwrk=maxlan*(maxlan+10);

	if( lwrk > 0 ) {
		res = (double*)malloc(lwrk*sizeof(double));
		wrk = (double*)malloc(lwrk*sizeof(double));
	}


	TRL::trl_init_info( &info, nrow, maxlan, lohi, ned, pres, 7, maxmv, -1, &mt );

	TRL::trl_set_iguess( &info, 0, 1, 0, NULL );
	// the Lanczos recurrence is set to start with [1,1,...,1]^T
	memset(eval, 0, mev*sizeof(double) );


	clock_t t1, t2;
	t1 = clock();
	TRL::trlan(op, &info, nrow, mev, eval, evec, nrow, lwrk, res );
	t2 = clock();

	sort_EGrst(info.nec, ned, nrow, eval, evec,egval,egvec,true);

#ifdef DEBUG_WEI_TRL
	printf( "TRLan: %d secs\n",
		(int)((t2-t1)/CLOCKS_PER_SEC) );

	TRL::trl_print_info(&info, 3*nrow);

	cout<<"\n\n";
	for(i=0; i<info.nec;i++){
		cout<<eval[i]<<"\t";
	}
	cout<<"\n\n";

	cout<<"\n\n";
	for(i=0; i<egval.size();i++){
		cout<<egval[i]<<"\t";
	}
	cout<<"\n\n";

	cout<<"\nsum sqr of eigen vec\n";
	for(i=0; i<egvec.size();i++){
		double sum_egvec=0;
		for(j=0; j<egvec[i].size();j++){
			sum_egvec+=egvec[i][j]*egvec[i][j];
		}
		cout<<sqrt(sum_egvec)<<"\t";
	}
	cout<<"\n\n";
	//system("pause");


	/*
	int l=0;
	for(i=0; i<info.nec;i++){
	for(l=0; l<nrow; l++){
	cout<<evec[i*nrow + l]<<"\t";
	}
	cout<<"\n\n\n";
	}
	*/
#endif



	free(eval); free(evec); free(exact);

	if( lwrk > 0 ) {
		free(res);
		free(wrk);
	}

}




void TRL::eg_nuTRan(trl_matvec op, vector< vector<double> > & mt,
			   vector<double> & egval, vector< vector<double> > & egvec, 
			   int ned, double pres ,const int lohi, int mev, int maxlan, int maxmv)
{
	int nrow=mt.size();
	if(nrow<maxlan){
		maxlan=nrow;
		mev=nrow/2;
	}
	if(nrow<ned) ned=nrow;

	int lwrk;
	// local variable declaration

	//double eval[mev], evec[mev*nrow], exact[mev];
	double *eval=0, *evec=0, *exact=0;
	eval = (double*)malloc(mev*sizeof(double));
	evec = (double*)malloc(mev*nrow*sizeof(double));
	exact = (double*)malloc(mev*sizeof(double));


	double *res=0, *wrk=0;
	trl_info info;
	int i, j, k, fp, check;
	//char name2[133], name[150];
	int tmp1, tmp2, nlen;


	lwrk=maxlan*(maxlan+10);

	if( lwrk > 0 ) {
		res = (double*)malloc(lwrk*sizeof(double));
		wrk = (double*)malloc(lwrk*sizeof(double));
	}


	TRL::trl_init_info( &info, nrow, maxlan, lohi, ned, pres, 7, maxmv, -1, &mt );

	TRL::trl_set_iguess( &info, 0, 1, 0, NULL );
	// the Lanczos recurrence is set to start with [1,1,...,1]^T
	memset(eval, 0, mev*sizeof(double) );


	clock_t t1, t2;
	t1 = clock();
	TRL::trlan(op, &info, nrow, mev, eval, evec, nrow, lwrk, res );
	t2 = clock();

	sort_EGrst(info.nec, ned, nrow, eval, evec,egval,egvec);

#ifdef DEBUG_WEI_TRL
	printf( "TRLan: %d secs\n",
		(int)((t2-t1)/CLOCKS_PER_SEC) );

	TRL::trl_print_info(&info, 3*nrow);

	cout<<"\n\n";
	for(i=0; i<info.nec;i++){
		cout<<eval[i]<<"\t";
	}
	cout<<"\n\n";

	cout<<"\n\n";
	for(i=0; i<egval.size();i++){
		cout<<egval[i]<<"\t";
	}
	cout<<"\n\n";

	cout<<"\nsum sqr of eigen vec\n";
	for(i=0; i<egvec.size();i++){
		double sum_egvec=0;
		for(j=0; j<egvec[i].size();j++){
			sum_egvec+=egvec[i][j]*egvec[i][j];
		}
		cout<<sqrt(sum_egvec)<<"\t";
	}
	cout<<"\n\n";
	//system("pause");


	/*
	int l=0;
	for(i=0; i<info.nec;i++){
	for(l=0; l<nrow; l++){
	cout<<evec[i*nrow + l]<<"\t";
	}
	cout<<"\n\n\n";
	}
	*/
#endif



	free(eval); free(evec); free(exact);

	if( lwrk > 0 ) {
		free(res);
		free(wrk);
	}

}



void TRL::eg_nuTRan(trl_matvec op, void* mt, int nrow,
					vector<double> & egval, vector< vector<double> > & egvec, 
					int ned, double pres ,const int lohi, int mev, int maxlan, int maxmv)
{
	//int nrow=mt.size();
	if(nrow<maxlan){
		maxlan=nrow;
		mev=nrow/2;
	}
	if(nrow<ned) ned=nrow;

	int lwrk;
	// local variable declaration

	//double eval[mev], evec[mev*nrow], exact[mev];
	double *eval=0, *evec=0, *exact=0;
	eval = (double*)malloc(mev*sizeof(double));
	evec = (double*)malloc(mev*nrow*sizeof(double));
	exact = (double*)malloc(mev*sizeof(double));


	double *res=0, *wrk=0;
	trl_info info;
	int i, j, k, fp, check;
	//char name2[133], name[150];
	int tmp1, tmp2, nlen;


	lwrk=maxlan*(maxlan+10);

	if( lwrk > 0 ) {
		res = (double*)malloc(lwrk*sizeof(double));
		wrk = (double*)malloc(lwrk*sizeof(double));
	}


	TRL::trl_init_info( &info, nrow, maxlan, lohi, ned, pres, 7, maxmv, -1, mt );

	TRL::trl_set_iguess( &info, 0, 1, 0, NULL );
	// the Lanczos recurrence is set to start with [1,1,...,1]^T
	memset(eval, 0, mev*sizeof(double) );


	clock_t t1, t2;
	t1 = clock();
	TRL::trlan(op, &info, nrow, mev, eval, evec, nrow, lwrk, res );
	t2 = clock();

	sort_EGrst(info.nec, ned, nrow, eval, evec,egval,egvec);

#ifdef DEBUG_WEI_TRL
	printf( "TRLan: %d secs\n",
		(int)((t2-t1)/CLOCKS_PER_SEC) );

	TRL::trl_print_info(&info, 3*nrow);

	cout<<"\n\n";
	for(i=0; i<info.nec;i++){
		cout<<eval[i]<<"\t";
	}
	cout<<"\n\n";

	cout<<"\n\n";
	for(i=0; i<egval.size();i++){
		cout<<egval[i]<<"\t";
	}
	cout<<"\n\n";

	cout<<"\nsum sqr of eigen vec\n";
	for(i=0; i<egvec.size();i++){
		double sum_egvec=0;
		for(j=0; j<egvec[i].size();j++){
			sum_egvec+=egvec[i][j]*egvec[i][j];
		}
		cout<<sqrt(sum_egvec)<<"\t";
	}
	cout<<"\n\n";
	//system("pause");


	/*
	int l=0;
	for(i=0; i<info.nec;i++){
	for(l=0; l<nrow; l++){
	cout<<evec[i*nrow + l]<<"\t";
	}
	cout<<"\n\n\n";
	}
	*/
#endif



	free(eval); free(evec); free(exact);

	if( lwrk > 0 ) {
		free(res);
		free(wrk);
	}

}



/*
void TRL::eg_nuTRan3(trl_matvec op, vector< vector<double> > & mtt,
					vector<double> & egval, vector< vector<double> > & egvec, 
					int ned, double pres ,const int lohi, int mev, int maxlan, int maxmv)
{
	int nrow=mtt.size();
	if(nrow<maxlan){
		maxlan=nrow;
		mev=nrow/2;
	}
	if(nrow<ned) ned=nrow;


	int lwrk;
	// local variable declaration

	//double eval[mev], evec[mev*nrow], exact[mev];
	double *eval=0, *evec=0, *exact=0;
	eval = (double*)malloc(mev*sizeof(double));
	evec = (double*)malloc(mev*nrow*sizeof(double));
	exact = (double*)malloc(mev*sizeof(double));


	double *res=0, *wrk=0;
	trl_info info;
	int i, j, k, fp, check;
	//char name2[133], name[150];
	int tmp1, tmp2, nlen;

	//make a data pointer
	double *mt=new double[nrow*nrow];
	cout<<mt<<"\n";
	for(i=0; i<nrow; i++){
		for(j=0; j<nrow; j++){
			mt[i*nrow+j]=mtt[i][j];
		}
	}
	cout<<mt<<"\n";

	lwrk=maxlan*(maxlan+10);

	if( lwrk > 0 ) {
		res = (double*)malloc(lwrk*sizeof(double));
		wrk = (double*)malloc(lwrk*sizeof(double));
	}


	TRL::trl_init_info( &info, nrow, maxlan, lohi, ned, pres, 7, maxmv, -1, mt );

	cout<<info.mvparam<<"\n";

	TRL::trl_set_iguess( &info, 0, 1, 0, NULL );
	// the Lanczos recurrence is set to start with [1,1,...,1]^T
	memset(eval, 0, mev*sizeof(double) );


	clock_t t1, t2;
	t1 = clock();
	TRL::trlan(op, &info, nrow, mev, eval, evec, nrow, lwrk, res );
	t2 = clock();

	sort_EGrst(info.nec, ned, nrow, eval, evec,egval,egvec);

#ifdef DEBUG_WEI_TRL
	printf( "TRLan: %d secs\n",
		(int)((t2-t1)/CLOCKS_PER_SEC) );

	TRL::trl_print_info(&info, 3*nrow);

	cout<<"\n\n";
	for(i=0; i<info.nec;i++){
		cout<<eval[i]<<"\t";
	}
	cout<<"\n\n";

	cout<<"\n\n";
	for(i=0; i<egval.size();i++){
		cout<<egval[i]<<"\t";
	}
	cout<<"\n\n";

	cout<<"\nsum sqr of eigen vec\n";
	for(i=0; i<egvec.size();i++){
		double sum_egvec=0;
		for(j=0; j<egvec[i].size();j++){
			sum_egvec+=egvec[i][j]*egvec[i][j];
		}
		cout<<sqrt(sum_egvec)<<"\t";
	}
	cout<<"\n\n";
	system("pause");


	
#endif


	//delete data pointer
	delete [] mt;

	free(eval); free(evec); free(exact);

	if( lwrk > 0 ) {
		free(res);
		free(wrk);
	}

}


void TRL::eg_nuTRan2(trl_matvec op, double * mt, int nrow,
					vector<double> & egval, vector< vector<double> > & egvec, 
					int ned, double pres ,const int lohi, int mev, int maxlan, int maxmv)
{
	if(nrow<maxlan){
		maxlan=nrow;
		mev=nrow/2;
	}
	if(nrow<ned) ned=nrow;

	int lwrk;
	// local variable declaration

	//double eval[mev], evec[mev*nrow], exact[mev];
	double *eval=0, *evec=0, *exact=0;
	eval = (double*)malloc(mev*sizeof(double));
	evec = (double*)malloc(mev*nrow*sizeof(double));
	exact = (double*)malloc(mev*sizeof(double));


	double *res=0, *wrk=0;
	trl_info info;
	int i, j, k, fp, check;
	//char name2[133], name[150];
	int tmp1, tmp2, nlen;


	lwrk=maxlan*(maxlan+10);

	if( lwrk > 0 ) {
		res = (double*)malloc(lwrk*sizeof(double));
		wrk = (double*)malloc(lwrk*sizeof(double));
	}


	TRL::trl_init_info( &info, nrow, maxlan, lohi, ned, pres, 7, maxmv, -1, mt );

	TRL::trl_set_iguess( &info, 0, 1, 0, NULL );
	// the Lanczos recurrence is set to start with [1,1,...,1]^T
	memset(eval, 0, mev*sizeof(double) );


	clock_t t1, t2;
	t1 = clock();
	TRL::trlan(op, &info, nrow, mev, eval, evec, nrow, lwrk, res );
	t2 = clock();

	sort_EGrst(info.nec, ned, nrow, eval, evec,egval,egvec);

#ifdef DEBUG_WEI_TRL
	printf( "TRLan: %d secs\n",
		(int)((t2-t1)/CLOCKS_PER_SEC) );

	TRL::trl_print_info(&info, 3*nrow);

	cout<<"\n\n";
	for(i=0; i<info.nec;i++){
		cout<<eval[i]<<"\t";
	}
	cout<<"\n\n";

	cout<<"\n\n";
	for(i=0; i<egval.size();i++){
		cout<<egval[i]<<"\t";
	}
	cout<<"\n\n";

	cout<<"\nsum sqr of eigen vec\n";
	for(i=0; i<egvec.size();i++){
		double sum_egvec=0;
		for(j=0; j<egvec[i].size();j++){
			sum_egvec+=egvec[i][j]*egvec[i][j];
		}
		cout<<sqrt(sum_egvec)<<"\t";
	}
	cout<<"\n\n";


	

	
#endif



	free(eval); free(evec); free(exact);

	if( lwrk > 0 ) {
		free(res);
		free(wrk);
	}

}




*/
void TRL::trl_init_info(trl_info * info, int nrow, int mxlan, int lohi,
				   int ned, double tol, int restart, int maxmv,
				   int mpicom,void *mvparam)
{
	//
	// Purpose:
	// ========
	// Initializes a TRL_INFO variable. This function must be called before calling
	// any other user level routine in TRLAN package.
	//
	// Arguments:
	// ==========
	// info    (input) pointer to the structure trl_info_
	//          On entry, points to the data structure to store the information
	//          about the eigenvalue problem and the progress of TRLAN.
	//          On exit, points to the initialized data structure.
	//
	// nrow    (input) integer
	//          On entry, specifies the local dimension of the problem.
	//
	// mxlan   (input) integer
	//          On entry, specifies the maximum number of basis vectors to be used.
	//
	// lohi    (input)  integer
	//          On entry, specifies, which end of the spectrum to compute:
	//           lohi < 0, then lower end, the smallest eigenvalues
	//           lohi > 0, then high end, the largest eigenvalues
	//           lohi = 0, then either lower and or high end, whoever converges first
	//          are computed.
	//
	// ned      (input) integer
	//           On entry, specifies the number of wanted eigenvalues and eigenvectors.
	//
	// tol      (optional) double precision
	//           If provided, specifies the tolerance on residual norm. By default,
	//           tol is set to be sqrt(epsilon).
	//
	// trestart (optional) integer
	//           If provided, specifies the thick-restarting scheme, 1-4. Default is 1.
	//
	// mxmv     (optionial) integer
	//           If provided, specifies the maximum number of matrix-vector 
	//           multiplication allowed. By default, mxmv is set to be 
	//           (info->ntot)*(info->ned).
	//
	// mpicom   (optional) integer
	//           If provided, specifites the MPI communicator. By default, it is a 
	//           duplicate of MPI_COMM_WORLD. In sequential case, this is set to 0.
	//
	//
	// ..
	// .. executable statements ..
	//  va_list argptr;
	//  va_start( argptr,nopts );
	//  if( nopts > 0 ) {
	if (tol > 0) {
		//info->tol = va_arg( argptr, double );
		info->tol = tol;
		if (info->tol <= DBL_MIN) {
			info->tol = DBL_EPSILON;
		} else if (info->tol > 1.0) {
			info->tol = min(0.1, 1.0 / (info->tol));
		}
	} else {
		info->tol = sqrt(DBL_EPSILON);
	}
	//if( nopts > 1 ) {
	if (restart > 0) {
		//info->restart = va_arg( argptr,int );
		info->restart = restart;
	} else {
		info->restart = 0;
	}
	//if( nopts > 2 ) {
	if (maxmv > 0) {
		//info->maxmv = va_arg( argptr,int );
		info->maxmv = maxmv;
	} else {
		info->maxmv = min(max(info->ntot, 1000), 1000 * info->ned);
	}

	//should I change this to 0?
	info->mpicom = -INT_MAX;


	//va_end( argptr );
	// setup the rest of arguments
	info->maxlan = mxlan;
	if (mxlan <= ned) {
		info->maxlan = ned + max(ned, 6);
	}
	info->lohi = lohi;
	info->ned = ned;
	info->nloc = nrow;
	info->ntot = nrow;
	info->guess = 0;
	info->nec = 0;
	info->locked = info->nec;
	info->matvec = 0;
	info->nloop = 0;
	info->north = 0;
	info->nrand = 0;
	info->flop = 0;
	info->rflp = 0;
	info->flop_h = 0;
	info->rflp_h = 0;
	info->flop_r = 0;
	info->rflp_r = 0;
	info->clk_rate = CLOCKS_PER_SEC;
#ifdef __64INT
	info->clk_max = 9223372036854775807LL;
#else
	info->clk_max = (clock_t) (pow(2.0, (8.0 * sizeof(clock_t) - 1.0))-1.0);
	if( info->clk_max < 0 ) {
		if( sizeof(clock_t) == 8 ) {
			info->clk_max = 9223372036854775807LL;
		} else {
			printf( "error initializing clock.\n" );
		}
	}
#endif
	if( (double)(info->clk_max) <= 0 ) printf( "??\n" );
	//info->clk_max = -1;
	info->clk_tot = 0;
	info->clk_op = 0;
	info->clk_orth = 0;
	info->clk_res = 0;
	info->tick_t = 0;
	info->tick_o = 0;
	info->tick_h = 0;
	info->tick_r = 0;
	info->clk_in = 0;
	info->clk_out = 0;
	info->wrds_in = 0;
	info->wrds_out = 0;
	info->verbose = 0;
	info->stat = 0;
	info->anrm = 0;
	info->tmv = -1;
	info->trgt = -DBL_MAX;
	info->tres = -1.0;
	info->crat = 1.0;

	info->predicted_crate = 0.0;
	info->old_target = 0.0;
	info->target_id = 0;
	info->ref = 0.0;
	info->avgm = 0.0;
	info->k1 = 0;
	info->k2 = 0;
	info->k = 0;

	//add the pointer void *mvparam
	info->mvparam = mvparam;

	info->my_pe = 0;
	info->npes = 1;
	info->cpflag_ = 0;
	strcpy(info->oldcpf, "");
	// log file pointer
	info->log_io = 99;
	strcpy(info->log_file, "");
	// checkpoint file pointer
	info->cpio = 98;
	strcpy(info->cpfile, "");
	//
	// .. end of trl_init_info_ ..
	//
}


void TRL::trl_set_iguess(trl_info * info, int nec, int iguess, int nopts,
					char *cpf )
{
	/*
	// Purpose:
	// ========
	// Set up parameters related to initial guesses of the Lanczos iterations, i.e., the number of
	// eigenvector already converged (initially assumed to be zero) and whether the user has
	// provided initial guess vector (initially assumed no).  It can also tell TRLan to read check
	// point file that is different from the default name.
	//
	// Arguments:
	// ==========
	// info    (input/output) pointer to the structure trl_info_
	//          On entry, points to the data structure to store the current information about
	//          the eigenvalue problem and the progress of TRLAN.
	//
	// nec     (input) integer
	//          On entry, specifies the number of eigenvalues, that have been converged.
	//
	// iguess  (input) integer
	//          On entry, specifies one of the following options:
	//           iguess <= 0, user did not provide initial guess, use [1,1,..,1].
	//           iguess =  1, user has supplied initial guess, will only use the first one.
	//           iguess >  1, restart with previous check-point file.
	//
	// nopts   (optional)
	//          If provided, it specifies the name of checkpoints file.
	//
	// ..
	// .. executable statements ..
	*/
	/*
	char cpf[STRING_LEN];
	//printf( "nec=%d iguess=%d nopts=%d\n",nec,iguess,nopts );
	va_list argptr;
	va_start(argptr, nopts);
	if (nopts > 0) {
	strcpy(info->oldcpf, va_arg(argptr, char *));
	} else {
	strcpy(info->oldcpf, "");
	}
	va_end(argptr);
	*/
	/* assign nec and iguess flag_s to info */
	info->nec = nec;
	info->guess = iguess;
	if (strlen(info->oldcpf) > 0 && info->guess > 1) {
		/* check to make sure the files exist */
		trl_pe_filename(STRING_LEN, cpf, info->oldcpf, info->my_pe,
			info->npes);
		if ((info->cpt_fp = fopen(cpf, "r")) != NULL) {
			if (fclose(info->cpt_fp) != 0) {
				info->stat = -9;
			}
		} else {
			info->stat = -8;
		}
		info->stat = TRL::trl_sync_flag_(info->mpicom, info->stat);
	} else {
		info->stat = 0;
	}
	/*
	// .. end of trl_set_iguess_ ..
	*/
}


//
void TRL::trl_clear_counter(trl_info * info, int nrow, int mev, int lde)
{
	/*
	// Purpose:
	// ========
	// Clears the counters inside info and performs a minimal check on the input parameters.
	//
	// Arguments:
	// ==========
	// info    (input/output) pointer to the structure trl_info_
	//          On entry, points to the data structure to store the current information about
	//          the eigenvalue problem and the progress of TRLAN.
	//          On exit, the information is cleared.
	//
	// nrow    (input) integer
	//          On entry, specifies the number of rows that is on this proccesor.
	//
	// mev     (input) integer
	//          On entry, specifies the number of Ritz pairs, that can be stored in
	//          eval and evec.
	//
	// lde     (input) integer
	//          On entry, specifies the leading dimension of the array evec, i.e.,
	//          (lde >= nrow).
	//
	// ..
	// .. local scalars ..
	*/
	int ntmp;
	/*
	// ..
	// .. executable statements ..
	*/
	info->stat = 0;
	if (nrow != info->nloc || nrow > info->ntot) {
		printf("TRLAN: info not setup for this problem.\n");
		printf("       Please reset info before calling TRLAN.\n");
		info->stat = -1;
	}
	if (info->nec < 0)
		info->nec = 0;
	if (lde < nrow) {
		printf("TRLAN: leading dimension of EVEC to small.\n");
		info->stat = -2;
	}
	if (info->tol >= 1.0) {
		info->tol = sqrt(DBL_EPSILON);
	} else if (info->tol <= DBL_MIN) {
		info->tol = DBL_EPSILON;
	}
	if (info->ned + info->ned >= info->ntot) {
		printf
			("TRLAN: info->ned (%d) is large relative to the matrix dimension (%d)\n",
			info->ned, info->ntot);
		printf
			(" **    It is more appropriate to use LAPACK dsyev/ssyev.\n");
		if (info->ned > info->ntot) {
			info->ned = min(info->ntot - 1, info->maxlan - 3);
			printf("TRLAN: ** reduced ned to %d **\n", info->ned);
		}
	}
	if (mev < info->ned) {
		printf
			("TRLAN: array EVAL and EVEC can not hold wanted number of eigenpairs.\n");
		info->stat = -3;
	}
	if (info->ntot < 10) {
		printf
			("TRLAN is not designed to work with such a small matrix(%dx%d).  Please use LAPACK or EISPACK instead.\n",
			info->ntot, info->ntot);
		info->stat = -4;
	}
	info->nrand = info->stat;
	info->stat = TRL::trl_sync_flag_(info->mpicom, info->nrand);

	/* decide what is a good maximum basis size to use */
	if (info->maxlan < info->ned + 3) {
		info->maxlan =
			info->ned + min(info->ned,
			20) + (int) (log((double) info->ntot));
		info->maxlan = min(info->maxlan, info->ntot);
		printf("TRLAN: ** reset maxlan to %d! **\n", info->maxlan);
	}
	if (info->maxlan < mev) {
		ntmp = min(info->ntot, max(100 + info->ned, 10 * info->ned));
		if (mev < ntmp) {
			info->maxlan = mev;
		} else {
			info->maxlan = ntmp;
		}
	}
	if (info->maxlan < 5) {
		printf
			("TRLAN must have at least 5 basis vectors, it is currently %d.\n",
			info->maxlan);
		info->stat = -5;
	}

	/* clear regular counters */
	info->tmv = -1;
	info->klan = min(info->maxlan, info->ntot);
	if (info->restart >= 7) {
		info->klan =
			min(info->maxlan, max(100, min(info->klan, 2 * (info->ned))));
	}
	info->locked = info->nec;
	info->matvec = 0;
	info->nloop = 0;
	info->north = 0;
	info->nrand = 0;
	info->flop = 0;
	info->rflp = 0;
	info->flop_h = 0;
	info->rflp_h = 0;
	info->flop_r = 0;
	info->rflp_r = 0;
	info->tick_t = 0.0;
	info->clk_op = 0;
	info->tick_o = 0.0;
	info->clk_orth = 0;
	info->tick_h = 0.0;
	info->clk_res = 0;
	info->tick_r = 0.0;
	info->clk_in = 0;
	info->clk_out = 0;
	info->wrds_in = 0;
	info->wrds_out = 0;
	info->avgm = 0.0;
	return;
	/*
	// .. end of trl_clear_counter ..
	*/
}

void TRL::trl_g_sum(int mpicom, int nelm, double *x, double *y)
{
	//
	// Purpose:
	// ========
	// Performs global sum in the parallel environment, nothing is done here.
	//
	// Arguments:
	// ==========
	// mpicom    (input) integer
	//            On entry, specifies the MPI communicator.
	//
	// nelm      (input) integer
	//            On entry, specifies the number of elements in x and y.
	//
	// x         (input/output) double precision array of dimension nelm
	//            On entry, contains the resulting vector on this processor.
	//            On exit, contain the resulting vector of global sum.
	//
	// y         (workspace) double precision array of dimensioni nelm
	//
}

////
int TRL::trl_sync_flag_(int mpicom, int inflag_)
{
	//
	// Purpose:
	// ========
	// Given an integer_ value, returns the minimum value of all the PEs
	//
	// Arguments:
	// ==========
	// mpicom    (input) integer
	//            On entry, specifies the MPI communicator.
	//
	// inflag_    (inpuut) integer
	//            On entry, specifies the integer_ value from this processor.
	//
	return inflag_;
}

////
void TRL::trl_g_dot_(int mpicom, int nrow, double *v1, int ld1, int m1,
				double *v2, int ld2, int m2, double *rr, double *wrk)
{
	//
	// Purpose:
	// ========
	// Implements a distributed version of BLAS routine dgemv, which is used to compute
	// dot-products by TRLAN, i.e., wrk = [V1, V2]'*rr.
	//
	// Arguments:
	// ==========
	// mpicom     (input) integer
	//             On entry, specifies MPI communicator.
	//
	// nrow       (input) integer
	//             On entry, specifies, the number of rows on the local processor.
	//
	// v1         (input) double precision array of dimension (ld1,m1)
	//             On entry, contains the first part of the matrix.
	//
	// ld1        (input) integer
	//             On entry, specifies the leading dimension of v1.
	//
	// m1         (input) integer
	//             On entry, specifies the number of columns in v1.
	//
	// v2         (input) double precision array of dimension (ld2,m2)
	//             On entry, contains the second part of the matrix.
	//
	// ld2        (input) integer
	//             On entry, specifies the leading dimension of v2.
	//
	// m2         (input) integer
	//             On entry, specifies the number of columns in v2.
	//
	// rr         (input) double precision array of length (nrow)
	//             On entry, contains the vector to be multiplied.
	//
	// wrk        (output) double precision array of length (m1+m2)
	//             On exit, contains th results of this operation.  !! size not checked !!
	//
	// ..
	// .. local parameters ..
	char trans = 'T';
	double one = 1.0, zero = 0.0;
	integer_ c__1 = 1;
	//
	// ..
	// .. local scalars ..
	int i, nd;
	//
	// ..
	// .. executable statements ..
	nd = m1 + m2;
	// nothing to do if both m1 and m2 are zero
	if (nd <= 0)
		return;
	// make sure the array sizes are correct
	if (ld1 < nrow || ld2 < nrow) {
		fprintf(stderr, "trl_g_dot: incorrect array sizes\n");
		exit(0);
	}
	if (m1 > 2) {
		trl_dgemv(&trans, nrow, m1, one, v1, ld1, rr, c__1, zero, wrk,
			c__1);
	} else if (m1 == 2) {
		wrk[0] = zero;
		wrk[1] = zero;
		for (i = 0; i < nrow; i++) {
			wrk[0] += v1[i] * rr[i];
			wrk[1] += v1[ld1 + i] * rr[i];
		}
	} else if (m1 == 1) {
		wrk[0] = trl_ddot(nrow, v1, c__1, rr, c__1);
	}
	if (m2 > 2) {
		trl_dgemv(&trans, nrow, m2, one, v2, ld2, rr, c__1, zero, &wrk[m1],
			c__1);
	} else if (m2 == 2) {
		wrk[m1] = zero;
		wrk[nd - 1] = zero;
		for (i = 0; i < nrow; i++) {
			wrk[m1]     += v2[i]        * rr[i];
			wrk[nd - 1] += v2[ld2 + i] * rr[i];
		}
	} else if (m2 == 1) {
		wrk[m1] = trl_ddot(nrow, v2, c__1, rr, c__1);
	}
	//
	// .. end of trl_g_dot_ ..
	//
}

void TRL::trlan(trl_matvec op, trl_info * info, int nrow, int mev, double *eval,
		   double *evec, int lde, int lwrk, double *wrk )
{
	/*
	// Purpose: Top (user) level routines
	// ========
	// A thick-restart Lanczos routine for computing eigenvalues and
	// eigenvectors of a real symmetric operator/matrix (A).
	// -- only accept one input vector, the input vector is expected
	//    to be stored in the (nec+1)st column of EVEC.
	// -- it extends the Lanczos basis one vector at a time.
	// -- orthogonality among the Lanczos vectors are maintained using
	//    full re-orthogonalization when necessary.
	//
	// Requirements:
	// =============
	// 1) User supplies OP with the specified interface.
	// 2) If (info->nec>0), evec(1:nrow, 1:info->nec) must contain valid
	//    eigenvectors and eval(1:nec) must be the corresponding eigenvalues.
	//    These eigenpairs are assumed to have zero residual norms.
	// 3) lde >= nrow.
	// 4) The arrays evec and eval must be large enough to store the
	//    solutions, i.e., mev >= info->ned and mev >= info->nec.
	// 5) The array wrk may be of arbitrary size.  Internally, the workspace
	//    size is
	//        nrow*max(0,info->ned-size(evec,2))+maxlan*(maxlan+10)
	//
	//    If wrk is smaller than this, trlan routine will allocate additional
	//    workspace to accommodate.
	//
	// Arguments:
	// ==========
	// op      (input) function pointer
	//         On entry, points to a function that comptues op(X) == A*X,
	//         when given a set of vectors X.
	//         The operator that defines the eigenvalue problem is expected to
	//         have the following interface
	//          void op(nrow, ncol, xin, ldx, yout, ldy)
	//            nrow  (input) integer
	//                   On entry, specifies the number of rows in xin and yout.
	//            ncol  (input) integer
	//                   On entry, specifies, the number of columns in Xin and
	//                   Yout.
	//            xin   (input) double precision vector of length (ldx*ncol)
	//                   On entry, contatins the input vector to be multiplied.
	//                   It consists of ncol column vectors with each column
	//                   stored in consecutive order.
	//            ldx   (input) integer
	//                   On entry, specifies the leading dimension of the array
	//                   Xin, i.e., the i-th column vector starts with element
	//                   (i-1)*ldx+1 and ends with element (i-1)*ldx+nrow in Xin.
	//            yout  (output) double precision vector of length (ldy*ncol)
	//                   On exit, contains the result array, i.e., it stores the
	//                   result of matrix-vector multiplications.
	//            ldy   (input) integer
	//                   On entry, specifies the leading dimension of the array
	//                   yout, i.e., the i-th column vector starts with element
	//                   (i-1)*ldy+1 in Yout and ends with element (i-1)*ldy+nrow.
	//
	// info    (input) pointer to the structure trl_info_
	//          On entry, points to the data structure to store the information
	//          about the eigenvalue problem and the progress of TRLAN
	//
	// nrow    (input) integer
	//          On entry, specifies the number of rows that is on this processor.
	//
	// mev     (input) integer
	//          On entry, specifies the number of Ritz pairs, that can be stored in
	//          eval and evec.
	//
	// eval    (output) double precision vector of length (mev)
	//          On exist, stores the eigenvalues.
	//
	// evec    (output) double precision vector of lenvth (nrow*mev)
	//          On exit, stores the eigenvectors.
	//
	// lde     (input) integer
	//          On entry, specifies the leading dimension of the array evec, i.e.,
	//
	// lwrk    (optional) integer
	//          On entry, specifies, the size of WRK.  When both WRK and LWRK are
	//          present, then LWRK should correctly indicate the size of WRK. If WRK
	//          is present by not LWRK, the size of WRK is assumed to be MEV which is
	//          only enough to store the residual norms on exit.  If WRK is not
	//          present, LWRK is not used even if it is present.
	//          (lde >= nrow).
	//
	// wrk     (optional) workspace
	//          If it is provided and there is enough space, the residual norm of
	//          the converged eigenpairs will be stored at wrk(1:info->nec) on exit.
	//
	*/
	/*
	// ..
	// .. local scalars ..
	*/
	clock_t clk1;
	int ii, nbas, nmis, ibas, imis, ldb, lwrk0;
	double *base=0, *misc=0;
	/*
	// ..
	// .. executables statements ..
	*/
	imis = -1;			/* if this routine allocated misc, imis will be 0 */
	ibas = -1;			/* if this routine allocated base, ibas will be 0 */
	clk1 = clock();
	info->clk_tot = clk1;
	if (info->ned > mev) {
		printf
			("info->ned (%d) is larger than mev (%d) reducing info->ned to %d\n",
			info->ned, mev, mev);
		info->ned = mev;
	}
	/* there is nothing to do if there is no more eigenvalue to compute */
	if (info->ned <= info->nec || info->ned <= 0)
		goto end;
	/* determine the workspace size and leading dimensions
	if (nopts > 0) {
	va_list argptr;
	va_start(argptr, nopts);
	lwrk = va_arg(argptr, int);
	wrk = va_arg(argptr, double *);
	va_end(argptr);
	} else {
	lwrk = 0;
	wrk = NULL;
	}
	*/
	lwrk0 = lwrk;
	info->stat = 0;
	ldb = ((nrow + 3) / 4) * 4;
	if ((ldb % 4096) == 0)
		ldb = ldb + 4;
	trl_clear_counter(info, nrow, mev, lde);
	if (info->stat != 0)
		goto end;
	/*
	// Internally, the workspace is broken into two parts
	// one to store (maxlan+1) Lanczos vectors, and the other to
	// store all others (size maxlan*(maxlan+ned+14))
	// The next If-block decides how the two arrays are mapped.
	*/
	nbas = max(1, info->maxlan - mev + 1);
	ii = nbas * ldb;
	nmis = info->maxlan * (info->maxlan + 10);
	if (lwrk0 >= min(ii, nmis)) {
		/* use wrk either as base or misc or both depending its size    */
		if (lwrk0 >= ii + nmis) {
			/* WRK is large enough for both arrays */
			base = wrk;
			misc = &wrk[ii];
			nmis = lwrk0 - ii;
			//printf( "\n\n ******** Large enough workspace ********* \n\n" );
		} else if (lwrk0 >= max(ii, nmis)) {
			/* associate the larger one of base and misc to WRK */
			if (ii >= nmis) {
				base = wrk;
				misc = (double *) malloc(nmis * sizeof(double));
				if (misc == NULL)
					info->stat = -4;
				imis = 0;
			} else {
				misc = wrk;
				nmis = lwrk0;
				base = (double *) malloc(ii * sizeof(double));
				if (base == NULL)
					info->stat = -5;
				ibas = 0;
			}
		} else if (ii <= nmis) {
			/* base is smaller, associate base with WRK */
			base = wrk;
			misc = (double *) malloc(nmis * sizeof(double));
			if (misc == NULL)
				info->stat = -4;
			imis = 0;
		} else {
			/* misc is smaller, associate misc with WRK */
			misc = wrk;
			nmis = lwrk0;
			base = (double *) malloc(ii * sizeof(double));
			if (base == NULL)
				info->stat = -5;
			ibas = 0;
		}
	} else {
		/* have to allocate both base and misc */
		base = (double *) malloc(ii * sizeof(double));
		if (base == NULL)
			info->stat = -5;
		ibas = 0;
		misc = (double *) malloc(nmis * sizeof(double));
		if (misc == NULL)
			info->stat = -4;
		imis = 0;
	}
	memset(base, 0, ii * sizeof(double));
	memset(misc, 0, nmis * sizeof(double));
	/* make sure every process is successful so far */
	ii = trl_sync_flag_(info->mpicom, info->stat);
	info->stat = ii;
	if (ii != 0)
		goto end;
	/* open log and checkpoint files */
	trl_open_logfile(info);
	/* trl_open_cptfile(info); */
	if (info->verbose > 0) {
		trl_time_stamp(info->log_fp);
		trl_print_setup(info, nbas * ldb, nmis, lwrk0);
	}

	/* call trlanczos to do the real work  */
	//printf( "calling trlanczso (%d)\n",info->cpflag_ );
	trlanczos(op, info, nrow, mev, eval, evec, lde, base, ldb, nbas,
		misc, nmis);
#ifdef DEBUG
	printf( "DEBUG ** out of trlanczos (locked=%d) **\n",info->locked );
#endif

	/* close log and checkpoint files */
	trl_close_logfile(info);
	if (lwrk0 >= mev) {
		memcpy(wrk, misc, mev * sizeof(double));
	} else {
		memcpy(wrk, misc, lwrk0 * sizeof(double));
	}

	/* DONE, reclaim the space allocated */
end:
	if (imis == 0)
		free(misc);
	if (ibas == 0)
		free(base);
	clk1 = clock();
	if (clk1 < info->clk_tot) {
		info->tick_t +=
			(info->clk_max -
			info->clk_tot) / (double) (info->clk_rate);
		info->tick_t +=
			(info->clk_max + clk1) / (double) (info->clk_rate);
		info->clk_tot = 0;
	} else if (info->clk_tot < 0 && clk1 >= 0) {
		info->tick_t -= info->clk_tot / (double) (info->clk_rate);
		info->tick_t += clk1 / (double) (info->clk_rate);
		info->clk_tot = 0;
	} else {
		info->tick_t  += (clk1 - info->clk_tot) / (double) (info->clk_rate);
		info->clk_tot  = 0;
	}
	/*
	if (clk1 >= info->clk_tot) {
	//printf( "clk_tot = %d-%d = ",ii,info->clk_tot );
	info->clk_tot = ii - info->clk_tot;
	//printf( "%d\n",info->clk_tot );
	} else if (clk1 < 0) {
	//assuming ( info->clk_tot > 0 ), otherwise ?? wrap-around twice ??
	info->tick_t +=
	(info->clk_max - info->clk_tot) / (double) (info->clk_rate);
	info->clk_tot = info->clk_max + ii;
	} else {
	//assuming ( info->clk_tot < 0 ),
	info->tick_t +=
	(info->clk_max + info->clk_tot) / (double) (info->clk_rate);
	info->clk_tot = ii;
	}
	*/
	return;
	/*
	// .. end of trlan_ ..
	*/
}








void TRL::trl_set_restart(trl_info * info, double rfact)
{
	//
	// Purpose
	// =======
	// Set the (minimum) basis size for the restart 7 and 8 schemes, i.e., the (minimum) basis
	// size if rfact * (number of Ritz vectors kept)
	//
	// Arguments:
	// ==========
	// info    (input/output) pointer to the structure trl_info_
	//          On entry, points to the data structure to store the current information about
	//          the eigenvalue problem and the progress of TRLAN.
	//
	// rfact   (input) double precision
	//          On entry, specify the (minimum) basis size.
	//
	info->rfact = rfact;
}



void TRL::trl_set_debug(trl_info * info, int msglvl, char *filename)
{
	/*
	// Purpose:
	// ========
	// Set information related to debugging, the initialization routine trl_init_info sets the
	// parameters so that no debug information is printed.
	//
	// Arguments:
	// ==========
	// info    (input/output) pointer to the structure trl_info_
	//          On entry, points to the data structure to store the current information about
	//          the eigenvalue problem and the progress of TRLAN.
	//
	// msglvl  (input) integer
	//          On entry, specifies the new message/verbose level:
	//           msglvl <  0         : nothing is printed
	//           msglvl = 1, .., 10  : the higher the level, the more debug information is
	//                                 printed.
	//
	// file    (input) character string
	//          On entry, specifies the leading part of the debug file name. Each processor will
	//          generate its own debug file with the file name formed by appending the processor
	//          number to the string FILE. The actual name is composed by the routine
	//          TRL_PE_FILENAME in trlaux
	//
	// ..
	// .. executable statements ..
	*/
	info->verbose = msglvl;
	if (filename != NULL) {
		strcpy(info->log_file, filename);
		if (msglvl >= 0 && info->my_pe == 0) {
			printf("TRLan will write diagnostic messages to files with "
				"prefix %s.\n", info->log_file);
		}
	}
	/*
	// .. end of trl_set_debug_ ..
	*/
}

void TRL::trl_set_checkpoint(trl_info * info, int cpflag_, char *file)
{
	/*
	// Purpose:
	// ========
	// Set up the information related to check-pointing
	//
	// Arguments:
	// ==========
	// info    (input/output) pointer to the structure trl_info_
	//          On entry, points to the data structure to store the current information about
	//          the eigenvalue problem and the progress of TRLAN.
	//
	// cpflag_  (input) integer
	//          On entry, spcifies roughly how many titmes checkpoints are written.
	//
	// file    (input) character string
	//          On entry, specifies the leading part of the checkpoint file name. Each processor
	//          will generate its own debug file with the file name formed by appending the
	//          processor number to the string FILE. The actual name is composed by the routine
	//          TRL_PE_FILENAME in trlaux
	//
	// ..
	// .. executable statements ..
	*/
	info->cpflag_ = cpflag_;
	if (file != NULL) {
		strcpy(info->cpfile, file);
		if (info->verbose >= 0 && info->my_pe == 0) {
			printf("TRLan will write checkpoint to files with prefix %s.\n",
				info->cpfile);
		}
	}
	/*
	// .. end of trl_set_checkpoint_ ..
	*/
}




void TRL::trl_print_info(trl_info * info, int mvflop)
{
	/*
	// Purpose:
	// ========
	// Provides an uniform way of printing information stored in TRL_INFO_T.  It needs to be
	// called by all PEs. In parallel environment, when writing to standard outputd device, only
	// PE0 will actually write its local summary information. Note that this function must be
	// called on all PEs at the same time ***
	//
	// Arguments:
	// ==========
	// info    (input) pointer to the structure trl_info_
	//          On entry, points to the data structure to store the current information about
	//          the eigenvalue problem and the progress of TRLAN.
	//
	// mvflop  (input) integer
	//          On entry, specifies the number of floating-operations per MATVEC per PE.  This
	//          information has to be supplied by user, otherwise related entries are left blank
	//          in the print out.
	//
	// ..
	// .. local arrays ..
	*/
	double tmp1[12], tmp2[12];
	/*
	// ..
	// .. local variables ..
	*/
	int i;
	double t_tot, t_op, t_orth, t_res, t_in, t_out, rinv, r_tot, r_op,
		r_orth, r_res, r_in, r_out;
	/*
	// ..
	// .. executable statements ..
	*/
	if (info->clk_rate > 0) {
		rinv = 1.0 / (double) (info->clk_rate);
	} else {
		/* get clock rate */
		rinv = 1.0 / CLOCKS_PER_SEC;
	}
	t_op   = info->tick_o + info->clk_op   * rinv;
	t_tot  = info->tick_t + info->clk_tot  * rinv;
	t_res  = info->tick_r + info->clk_res  * rinv;
	t_orth = info->tick_h + info->clk_orth * rinv;
	t_in   = info->clk_in * rinv;
	t_out  = info->clk_out * rinv;
	//printf( "tick_t=%e clk_tot=%d\n",info->tick_t,info->clk_tot );
	if (t_op != 0 && mvflop != 0) {
		if (mvflop > 0) {
			r_op = ((double)mvflop) * ((double)info->matvec);
		} else {
			r_op = 0;
		}
	} else {
		r_op = 0;
	}
	if (t_orth != 0) {
		r_orth = info->rflp_h + info->flop_h;
	} else {
		r_orth = 0;
	}
	if (t_res != 0) {
		r_res = info->rflp_r + info->flop_r;
	} else {
		r_res = 0;
	}
	if (r_op > 0) {
		r_tot =
			((double) mvflop) * ((double) info->matvec) + info->rflp +
			((double) info->flop);
	} else {
		r_tot = 0;
	}
	if (info->clk_in > 0) {
		r_in = 8.0 * info->wrds_in;
	} else {
		r_in = 0;
	}
	if (info->clk_out > 0) {
		r_out = 8.0 * info->wrds_out;
	} else {
		r_out = 0;
	}
	tmp2[0] = t_tot;
	tmp2[1] = t_op;
	tmp2[2] = t_orth;
	tmp2[3] = t_res;
	tmp2[4] = t_in;
	tmp2[5] = t_out;
	tmp2[6] = r_tot;
	tmp2[7] = r_op;
	tmp2[8] = r_orth;
	tmp2[9] = r_res;
	tmp2[10] = r_in;
	tmp2[11] = r_out;
	//printf( "print info\n" );
	//printf( "calling g_sum\n" );
	trl_g_sum(info->mpicom, 12, tmp2, tmp1);
	if (info->log_fp == NULL) {
		trl_reopen_logfile(info);
	}
	trl_time_stamp(info->log_fp);
	//printf( "printing\n" );
	if (info->npes > 1) {
		fprintf(info->log_fp,
			"TRLAN execution summary (exit status = %d) on PE %d\n",
			info->stat, info->my_pe);
	} else {
		fprintf(info->log_fp,
			"TRLAN execution summary (exit status =%d)\n", info->stat);
	}
	if (info->lohi > 0) {
		fprintf(info->log_fp,
			"Number of LARGEST eigenpairs      %10d (computed) %11d (wanted)\n",
			info->nec, info->ned);
	} else if (info->lohi < 0) {
		fprintf(info->log_fp,
			"Number of SMALLEST eigenpairs    %10d (computed) %11d (wanted)\n",
			info->nec, info->ned);
	} else {
		fprintf(info->log_fp,
			"Number of EXTREME eigenpairs     %10d (computed) %11d (wanted)\n",
			info->nec, info->ned);
	}
	fprintf(info->log_fp,
		"Times the operator is applied:   %10d (MAX: %16d )\n",
		info->matvec, info->maxmv);
	fprintf(info->log_fp,
		"Problem size:                    %10d (PE: %4d) %11d (Global)\n",
		info->nloc, info->my_pe, info->ntot);
	fprintf(info->log_fp,
		"Convergence tolerance:           %10.3e (rel) %16.3e (abs)\n",
		info->tol, info->tol * info->anrm);
	fprintf(info->log_fp, "Maximum basis size:              %10d\n",
		info->maxlan);
	fprintf(info->log_fp, "Restarting scheme:               %10d\n",
		info->restart);
	fprintf(info->log_fp, "Number of re-orthogonalizations: %10d\n",
		info->north);
	fprintf(info->log_fp, "Number of (re)start loops:       %10d\n",
		info->nloop);
	if (info->nrand > 0) {
		fprintf(info->log_fp, "Number of random vectors used:   %10d\n",
			info->nrand);
	}
	if (info->npes > 1) {
		fprintf(info->log_fp, "Number of MPI processes:         %10d\n",
			info->npes);
	}
	fprintf(info->log_fp, "Number of eigenpairs locked:     %10d\n",
		info->locked);
	if (t_op > 0) {
		fprintf(info->log_fp,
			"OP(MATVEC):            %12.4e sec, %12.4e FLOP/S (%11.4e FLOP)\n",
			t_op, r_op / t_op, r_op);
	} else {
		fprintf(info->log_fp, "time in OP:            %12.4e sec\n", t_op);
	}
	if (t_orth > 0) {
		fprintf(info->log_fp,
			"Re-Orthogonalization:: %12.4e sec, %12.4e FLOP/S (%11.4e FLOP)\n",
			t_orth, r_orth / t_orth, r_orth);
	} else {
		fprintf(info->log_fp, "time in orth:          %12.4e sec\n",
			t_orth);
	}
	if (t_res > 0) {
		fprintf(info->log_fp,
			"Restarting::           %12.4e sec, %12.4e FLOP/S (%11.4e FLOP)\n",
			t_res, r_res / t_res, r_res);
	} else {
		fprintf(info->log_fp, "time in restarting:    %12.4e sec\n",
			t_res);
	}
	if (t_tot > 0) {
		fprintf(info->log_fp,
			"TRLAN on this PE:      %12.4e sec, %12.4e FLOP/S (%11.4e FLOP)\n",
			t_tot, r_tot / t_tot, r_tot);
	} else {
		fprintf(info->log_fp, "total time in TRLAN:   %12.4e sec\n",
			t_tot);
	}
	/*
	//if( info->verbose > 0 && info->log_fp != info->log_fp) {
	//  fprintf( info->log_fp, "Debug infomation written to files %s ####\n",info->log_file );
	//}
	*/
	if (info->guess > 1 && info->wrds_in > 0) {
		if (strlen(info->oldcpf) <= 0) {
			fprintf(info->log_fp,
				"TRLAN restarted with checkpoint files %s ####\n",
				info->oldcpf);
		} else {
			fprintf(info->log_fp,
				"TRLAN restarted with checkpoint files %s ####\n",
				info->cpfile);
		}
		fprintf(info->log_fp,
			"Bytes read   %12.5e, Time(sec): %12.5e, Rate(B/s): %12.5e\n",
			r_in, t_in, r_in / t_in);
	}
	if (info->clk_out > 0 && info->wrds_out > 0) {
		fprintf(info->log_fp, "Checkpoint files are %s ####\n",
			info->cpfile);
		fprintf(info->log_fp,
			"Bytes read   %12.5e, Time(sec): %12.5e, Rate(B/s): %12.5e\n",
			r_out, t_out, r_out / t_out);
	}
	if (info->npes > 1) {
		/* write global performance information */
		rinv = 1.0 / info->npes;
		for (i = 0; i < 12; i++) {
			tmp1[i] = tmp1[i] * rinv;
		}
		for (i = 0; i < 6; i++) {
			if (tmp1[i] > 0) {
				tmp1[i + 6] = tmp1[i + 6] / tmp1[i];
			} else {
				tmp1[i + 6] = 0.0;
			}
		}
		if (tmp1[4] == tmp1[5] && tmp1[4] == 0) {
			fprintf(info->log_fp,
				" -- Global summary -- \n" );
			fprintf(info->log_fp,
				"                       Overall,\t\t  MATVEC,\t  Re-orth,\t  Restart,\n");
			fprintf(info->log_fp,
				"Time(ave)             %11.4e,\t %11.4e,\t %11.4e,\t %11.4e\n",
				tmp1[0], tmp1[1], tmp1[2], tmp1[3]);
			fprintf(info->log_fp,
				"Rate(tot)             %11.4e,\t %11.4e,\t %11.4e,\t %11.4e\n",
				tmp1[6], tmp1[7], tmp1[8], tmp1[9]);
		} else {
			fprintf(info->log_fp,
				" -- Global summary -- \n" );
			fprintf(info->log_fp,
				"                       Overall,\t\t  MATVEC,\t  Re-orth,\t  Restart,\t  Read,\t  Write\n");
			fprintf(info->log_fp,
				"Time(ave)             %11.4e,\t %11.4e,\t %11.4e,\t %11.4e,\t %11.4e,\t %11.4e\n",
				tmp1[0], tmp1[1], tmp1[2], tmp1[3], tmp1[4], tmp1[5]);
			fprintf(info->log_fp,
				"Rate(tot)             %11.4e,\t %11.4e,\t %11.4e,\t %11.4e,\t %11.4e,\t %11.4e\n",
				tmp1[6], tmp1[7], tmp1[8], tmp1[9], tmp1[10],
				tmp1[11]);
		}
	}
	trl_close_logfile(info);
	return;
	/*
	// .. end of trl_print_info ..
	*/
}

void TRL::trl_terse_info(trl_info * info, FILE * iou)
{
	/*
	// Purpose:
	// ========
	// It is a more compact version of trl_print_info, i.e., this is a local routine, indivadual
	// PE can call it without regard of whether other PEs do the same and the output may be
	// written to a different I/O unit number than the log_io
	//
	// Arguments:
	// ==========
	// info    (input) pointer to the structure trl_info_
	//          On entry, points to the data structure to store the current information about
	//          the eigenvalue problem and the progress of TRLAN.
	//
	// iou     (input) pointer to a file
	//          On entry, points to the file, here the information is outputed.
	//
	// ..
	// .. local scalars ..
	*/
	int rate;
	double t_tot, t_op, t_orth, t_res;
	/*
	// ..
	// .. executable statements ..
	*/
	if (iou == NULL) {
		if (info->log_fp == NULL) {
			iou = stdout;
		} else {
			iou = info->log_fp;
		}
	}
	if (info->clk_rate > 0) {
		t_op = info->tick_o + info->clk_op / (double) (info->clk_rate);
		t_tot = info->tick_t + info->clk_tot / (double) (info->clk_rate);
		t_res = info->tick_r + info->clk_res / (double) (info->clk_rate);
		t_orth = info->tick_h + info->clk_orth / (double) (info->clk_rate);
	} else {
		/* get clock rate */
		rate = CLOCKS_PER_SEC;
		t_op = info->tick_o + info->clk_op / (double) (rate);
		t_tot = info->tick_t + info->clk_tot / (double) (rate);
		t_res = info->tick_r + info->clk_res / (double) (rate);
		t_orth = info->tick_h + info->clk_orth / (double) (rate);
	}
	if (info->lohi > 0) {
		fprintf(iou,
			"MAXLAN:%10d, Restart:%10d,   NED: + %7d,      NEC:%10d\n",
			info->maxlan, info->restart, info->ned, info->nec);
	} else if (info->lohi < 0) {
		fprintf(iou,
			"MAXLAN:%10d, Restart:%10d,   NED: - %7d,      NEC:%10d\n",
			info->maxlan, info->restart, info->ned, info->nec);
	} else {
		fprintf(iou,
			"MAXLAN:%10d, Restart:%10d,   NED: 0 %7d,      NEC:%10d\n",
			info->maxlan, info->restart, info->ned, info->nec);
	}
	fprintf(iou,
		"MATVEC:%10d,  Reorth:%10d, Nloop:   %7d,  Nlocked:%10d\n",
		info->matvec, info->north, info->nloop, info->locked);
	if (t_tot > 0.001 && max(t_tot, max(t_op, max(t_res, t_orth))) < 1000) {
		fprintf(iou,
			"Ttotal:%10.6f,    T_op:%10.6f, Torth:%10.6f,   Tstart:%10.6f\n",
			t_tot, t_op, t_orth, t_res);
	} else {
		fprintf(iou,
			"Ttotal:%10.3e,    T_op:%10.3e, Torth:%10.3e,   Tstart:%10.3e\n",
			t_tot, t_op, t_orth, t_res);
	}
	/*
	// .. end of trl_terse_info_ ..
	*/
}


void TRL::trl_print_setup(trl_info * info, int lbas, int lmis, int lwrk)
{
	/*
	// Purpose:
	// ========
	// Print the definition of the eigenvalue problme.
	//
	// Arguments:
	// ==========
	// info    (input) pointer to the structure trl_info_
	//          On entry, points to the data structure to store the current information about
	//          the eigenvalue problem and the progress of TRLAN.
	//
	// lbas    (input) integer
	//          On entry, specifies the size of workspace required to store Lanczos basis, i.e.,
	//          nrow*(maxlan-mev).
	//
	// lmis    (input) integer
	//          On entry, specifies the size of miscellenious workspace required to solve the
	//          eigen problem.
	//
	// lwrk    (input) integer
	//          On entry, specifies the size of workspace provided by the user.
	//
	// ..
	// .. executable statements ..
	// print the problem parameters
	*/
	if (info->lohi > 0) {
		fprintf(info->log_fp,
			"TRLAN is to compute %6d largest eigenpair(s).\n",
			info->ned);
	} else if (info->lohi < 0) {
		fprintf(info->log_fp,
			"TRLAN is to compute %6d smallest eigenpair(s).\n",
			info->ned);
	} else {
		fprintf(info->log_fp,
			"TRLAN is to compute %6d first converged eigenpair(s).\n",
			info->ned);
	}
	fprintf(info->log_fp,
		"Problem dimension: %9d (PE:%4d) %12d (Global)\n", info->nloc,
		info->my_pe, info->ntot);
	fprintf(info->log_fp, "Maximum basis size:                   %10d\n",
		info->maxlan);
	fprintf(info->log_fp, "Dynamic restarting scheme:            %10d\n",
		info->restart);
	fprintf(info->log_fp, "Maximum applications of the operator: %10d\n",
		info->maxmv);
	fprintf(info->log_fp, "Relative convergence tolerance: %10e\n",
		info->tol);
	/* initial guess */
	if (info->guess == 1) {
		fprintf(info->log_fp, "User provided the starting vector.\n");
	} else if (info->guess == 0) {
		fprintf(info->log_fp, "TRLAN uses [1,1,...] as starting vctor.\n");
	} else if (info->guess < 0) {
		fprintf(info->log_fp,
			"TRLAN generates a random starting vector.\n");
	} else if (info->oldcpf == 0 || strlen(info->oldcpf) == 0) {
		fprintf(info->log_fp,
			"Restarting with existing checkpoint files %s ####\n",
			info->oldcpf);
	} else {
		fprintf(info->log_fp,
			"Restarting with existing checkpoint files %s ####\n",
			info->cpfile);
	}
	if (info->cpflag_ > 0) {
		fprintf(info->log_fp,
			"TLRAN will write about %d sets of checkpointing files %s ####.\n",
			info->cpflag_, info->cpfile);
	}
	/* print the workspace size parameters */
	fprintf(info->log_fp, "(required) array BASE size is %d\n", lbas);
	fprintf(info->log_fp, "(required) array MISC size is %d\n", lmis);
	if (lwrk > 0) {
		fprintf(info->log_fp,
			"Caller has supplied a work array with %d elements.\n",
			lwrk);
	} else {
		fprintf(info->log_fp, "Caller did not supply work array.\n");
	}
}


void
TRL::trl_check_ritz(trl_matvec op, trl_info * info, int nrow, int ncol, double *rvec,
			   int ldrvec, double *alpha, int *check, double *beta,
			   double *eval, double *wrk, int lwrk)
{
	/*
	// Purpose:
	// ========
	// Performs a standard check on the computed Ritz pairs.
	//
	// Arguments:
	// ==========
	// op       (input) function pointer
	//           On entry, points to the matrix-vector multiplication routine.
	//
	// info     (input) pointer to the structure trl_info_
	//           On entry, points to the data structure to store the information
	//           about the eigenvalue problem and the progress of TRLAN
	//
	// nrow     (input) integer
	//           On entry, specifies the problem size.
	//
	// ncol     (input) integer
	//           On entry, specifies the number of Ritz values computed.
	//
	// rvec     (input) double precision array of dimension (nrow,ncol)
	//           On entry, specifies the array storing the Ritz vectors.
	//
	// alpha    (input) double precision array of dimension (ncol)
	//           On entry, contains he Ritz values computed.
	//
	// beta     (optional) double precision array of dimension (ncol)
	//           If provided, contaions the residual norms returned from a Lanczos routine.
	//
	// eval     (optional) double precision array of dimension (ncol)
	//           If provided, contains the actual eigenvalues.
	//
	// lwrk     (optional) integer
	//           If provided, specifies the size of workspace provided.
	//
	// wrk      (optional) double precision array of size(lwrk)
	//           If provided, double precidion workspace.
	//
	// ..
	// .. local parameters ..
	*/
	long c__1 = 1;
	int i__1 = 1;
	double d__1;
	/*
	// ..
	// .. local variables ..
	// aq -- store the result of matrix-vector multiplication, size nrow
	// rq -- store the Rayleigh-Quotient and the residual norms
	// gsumwrk -- workspace left over for trl_g_sum to use dimension of the input arrays
	*/
	double *aq, *rq, *gsumwrk, *res, *err;
	FILE *fp;
	int i, aqi, rqi, gsumwrki, icheck;
	double gapl, gapr;
	/*
	// ..
	// .. executable statements ..
	*/
	if (ncol <= 0)
		return;			/* there is nothing to do */

	/* figure out whether it is necessary to allocate new workspace */
	*check = 0;
	aqi = 0;
	rqi = 0;
	gsumwrki = 0;
	if (lwrk > nrow + (4 * ncol)) {
		aq = &wrk[0];
		rq = &wrk[nrow];
		gsumwrk = &wrk[nrow + (3 * ncol)];
	} else if (lwrk >= (nrow + ncol)) {
		aq = &wrk[0];
		gsumwrk = &wrk[nrow];
		rq = (double *) malloc(3 * ncol * sizeof(double));
		if (rq == NULL) {
			fprintf(stderr,
				"TRL_CHECK_RITZ: Failed to allocated workspace RQ.\n");
			exit(0);
		}
		rqi = 1;
	} else if (lwrk >= (4 * ncol)) {
		rq = &wrk[0];
		gsumwrk = &wrk[3 * ncol];
		aq = (double *) malloc(nrow * sizeof(double));
		if (aq == NULL) {
			fprintf(stderr,
				"TRL_CHECK_RITZ: Failed to allocated workspace AQ.\n");
			exit(0);
		}
		aqi = 1;
	} else if (lwrk >= ncol) {
		gsumwrk = wrk;
		aq = (double *) malloc(nrow * sizeof(double));
		if (aq == NULL) {
			fprintf(stderr,
				"TRL_CHECK_RITZ: Failed to allocated workspace AQ.\n");
			exit(0);
		}
		aqi = 1;
		rq = (double *) malloc(3 * ncol * sizeof(double));
		if (rq == NULL) {
			fprintf(stderr,
				"TRL_CHECK_RITZ: Failed to allocated workspace RQ.\n");
			free(aq);
			exit(0);
		}
		rqi = 1;
	} else {
		/* WRK not provided -- allocate space for AQ and RQ,  */
		/* gsumwrk points to the last third of RQ             */
		aq = (double *) malloc(nrow * sizeof(double));
		if (aq == NULL) {
			printf("TRL_CHECK_RITZ: Failed to allocated workspace AQ.\n");
			return;
		}
		aqi = 1;
		rq = (double *) malloc((3 * ncol) * sizeof(double));
		if (rq == NULL) {
			printf("TRL_CHECK_RITZ: Failed to allocated workspace RQ.\n");
			free(aq);
			return;
		}
		rqi = 1;
		gsumwrk = (double *) malloc(ncol * sizeof(double));
		if (gsumwrk == NULL) {
			printf
				("TRL_CHECK_RITZ: Failed to allocate workspace GSUMWRK.\n");
			free(aq);
			free(rq);
			return;
		}
		gsumwrki = 1;
	}
	memset(aq, 0, nrow * sizeof(double));
	memset(rq, 0, 2 * ncol * sizeof(double));
	memset(gsumwrk, 0, ncol * sizeof(double));
	/* go through each Ritz pair one at a time, compute Rayleigh  */
	/* quotient and the corresponding residual norm               */
	res = &rq[ncol];
	for (i = 0; i < ncol; i++) {
#ifdef TRL_FORTRAN_COMPATIBLE
		op(&nrow, &i__1, &rvec[i * ldrvec], &ldrvec, aq, &nrow);
#else
		op(nrow, i__1, rvec + i*ldrvec, ldrvec, aq, nrow, info->mvparam);
#endif
		/* Rayleigh quotient -- assuming rvec(:,i) has unit norm */
		rq[i] = trl_ddot(nrow, &rvec[i * ldrvec], c__1, aq, c__1);
		trl_g_sum(info->mpicom, 1, &rq[i], gsumwrk);
		d__1 = -rq[i]; /* indent separated =- into = - */
		trl_daxpy(nrow, d__1, &rvec[i * ldrvec], c__1, aq, c__1);
		res[i] = trl_ddot(nrow, aq, c__1, aq, c__1);
	}
	trl_g_sum(info->mpicom, ncol, res, gsumwrk);
	for (i = 0; i < ncol; i++) {
		res[i] = sqrt(res[i]);
	}
	/* compute the error estimate based on computed residual norms */
	/*  and the Ritz values                                        */
	err = &rq[2 * ncol];
	gapl = alpha[ncol - 1] - alpha[0];
	for (i = 0; i < ncol - 1; i++) {
		gapr = alpha[i + 1] - alpha[i];
		gapl = min(gapl, gapr);
		if (res[i] >= gapl) {
			err[i] = res[i];
		} else {
			err[i] = res[i] * res[i] / gapl;
		}
		gapl = gapr;
	}
	if (res[ncol - 1] >= gapl) {
		err[ncol - 1] = res[ncol - 1];
	} else {
		err[ncol - 1] = res[ncol - 1] * res[ncol - 1] / gapl;
	}
	/* if writing to stdout, only PE 0 does it */
	fp = info->log_fp;
	if (fp == NULL) {
		trl_reopen_logfile(info);
		fp = info->log_fp;
	}
	if (fp != stdout || info->my_pe <= 0) {
		if (info->stat != 0) {
#ifdef __DEBUG_OUT
			ferr = fopen("error.txt", "a");
			fprintf(ferr, " ** exit error (%d) ** \n", info->stat);
			fclose(ferr);
#endif
			*check = -4;
		}
		/* print out the information */
		fprintf(fp, "TRL_CHECK_RITZ: \n");
		fprintf(fp,
			"           Ritz value       res norm   res diff  est error  diff w rq  act. error\n");
#ifdef __DEBUG_OUT
		ferr = fopen("error.txt", "a");
#endif
		if (beta != NULL && eval != NULL) {
			for (i = 0; i < ncol; i++) {
				icheck = 0;
				fprintf(fp,
					"%21.14f    %11.3e%11.3e%11.3e%11.3e %11.3e%11.3e\n",
					alpha[i], res[i], beta[i] - res[i], err[i],
					rq[i] - alpha[i], eval[i] - alpha[i], eval[i]);
				/* check the accuracy of results.. */
				if (fabs(beta[i] - res[i]) > 0.00001) {
#ifdef __DEBUG_OUT
					fprintf(ferr,
						" res diff[%d] = |beta-res| = %5.3e - %5.3e = %5.3e > 0.00001\n",
						i, beta[i], res[i], fabs(beta[i] - res[i]));
#endif
					*check = *check - 1;
					icheck++;
				}

				if (fabs(rq[i] - alpha[i]) > nrow * nrow * info->tol) {
#ifdef __DEBUG_OUT
					fprintf(ferr,
						" diff rq[%d] = |rq-alpha| = %5.3e - %5.3e = %5.3e > nrow*nor*tau = %5.3e\n",
						i, rq[i], alpha[i], fabs(rq[i] - alpha[i]),
						nrow * nrow * info->tol);
#endif
					*check = *check - 1;
					icheck++;
				}

				if (fabs(eval[i] - alpha[i]) > 10 * nrow * nrow * info->tol ||
					fabs(eval[i] - alpha[i]) > 10 * err[i]) {
#ifdef __DEBUG_OUT
						fprintf(ferr,
							" act. error[%d] = |exact-alpha| = %5.3e - %5.3e = %5.3e > 10*nrow*nrow*tau =%5.3e or 10*est err = %5.3e\n",
							i, eval[i], alpha[i], fabs(eval[i] - alpha[i]),
							10 * nrow * nrow * info->tol, 10 * err[i]);
#endif
						*check = *check - 1;
						icheck++;
				}
			}

		} else if (beta != NULL) {
			for (i = 0; i < ncol; i++) {
				fprintf(fp, "%21.14f    %11.3e%11.3e%11.3e%11.3e\n",
					alpha[i], res[i], beta[i] - res[i], err[i],
					rq[i] - alpha[i]);
				/* check the accuracy of results.. */
				if (fabs(beta[i] - res[i]) > 0.00001) {
#ifdef __DEBUG_OUT
					fprintf(ferr,
						" res diff[%d] = |beta-res| = %5.3e - %5.3e = %5.3e > 0.00001\n",
						i, beta[i], res[i], fabs(beta[i] - res[i]));
#endif
					*check = *check - 1;
					icheck++;
				}

				if (fabs(rq[i] - alpha[i]) > nrow * nrow * info->tol) {
#ifdef __DEBUG_OUT
					fprintf(ferr,
						" diff rq[%d] = |rq-alpha| = %5.3e - %5.3e = %5.3e > nrow*nor*tau = %5.3e\n",
						i, rq[i], alpha[i], fabs(rq[i] - alpha[i]),
						nrow * nrow * info->tol);
#endif
					*check = *check - 1;
					icheck++;
				}
			}
		} else if (eval != NULL) {
			for (i = 0; i < ncol; i++) {
				fprintf(fp,
					"%21.14f     %11.3e           %11.3e%11.3e%11.3e%11.3e\n",
					alpha[i], res[i], err[i], rq[i] - alpha[i],
					eval[i] - alpha[i], eval[i]);
			}
		} else {
			for (i = 0; i < ncol; i++) {
				fprintf(fp, "%21.14f    %11.5e           %11.3e%11.3e\n",
					alpha[i], res[i], err[i], rq[i] - alpha[i]);
			}
		}
#ifdef __DEBUG_OUT
		fclose(ferr);
#endif
	}
	if (info->nec < info->ned)
		*check = 1;
	if (rqi > 0) {
		free(rq);
	}
	if (aqi > 0) {
		free(aq);
	}
	if (gsumwrki > 0) {
		free(gsumwrk);
	}
	trl_close_logfile(info);
	/*
	// .. end of trl_check_ritz_
	*/
}

void
TRL::trl_rayleigh_quotients(trl_matvec op, trl_info * info, int ncol, double *evec,
					   int lde, double *eres, double *base)
{
	/*
	// Purpose:
	// ========
	// Compute Rayleigh quotients, when it is given a set of Ritz vectors and Ritz values,
	// normalize the Ritz vectors, and compute their Rayleigh quotients to replace the Ritz values.
	//
	// Arguments:
	// ==========
	// op       (input) function pointer
	//           On entry, points to the matrix-vector multiplication routine.
	//
	// info     (input) pointer to the structure trl_info_
	//           On entry, points to the data structure to store the current information about
	//           the eigenvalue problem and the progress of TRLAN.
	//
	// evec     (input) double precision array of dimension (nloc,ncol)
	//           On entry, stores the portion of eigenvectors on this PE.
	//
	// base     (workspace)
	//           The workspace used to store results of MATVEC
	//
	// eres     (output) double precision array of dimension (ncol)
	//           On exist, store new Ritz values and new residual norms, i.e., if there are NEV
	//           Ritz pairs, eres(1:NEV) stores the new Rayleigh quotient and eres(nev+1:nev+nev)
	//           stores the new residual norms.
	//
	// base     (optional) double precision array od dimension (nloc)
	//           If provided, double precision workspace.
	//
	// ..
	// .. local parameters ..
	*/
	int c__1 = 1;
	int i__1 = 1;
	double d__1;
	/*
	// ..
	// .. local variables ..
	*/
	int i, nrow;
	double wrk[4], *avec;
	/*
	// ..
	// .. executable statements ..
	*/
	nrow = info->nloc;
	if (ncol <= 0)
		return;
	if (base != NULL) {
		avec = base;
	} else {
		avec = (double *) malloc(nrow * sizeof(double));
	}
	memset(avec, 0, nrow * sizeof(double));
	if (info->verbose >= 0) {
		FILE *fp = info->log_fp;
		if (fp == NULL) {
			trl_reopen_logfile(info);
		}
		fp = info->log_fp;
		fprintf(fp,
			"TRLAN computing Rayleigh Quotients for %d Ritz pairs\n",
			ncol);
	}
	/* loop through each vector to normalize the vector, compute Rayleigh  */
	/* quotient and compute residual norm of the new Ritz pairs            */
	for (i = 0; i < ncol; i++) {
		wrk[0] =
			trl_ddot(nrow, evec + i * lde, c__1, evec + i * lde, c__1);
#ifdef TRL_FORTRAN_COMPATIBLE
		op(&nrow, &i__1, evec + i * lde, &lde, avec, &nrow);
#else
		op(nrow, i__1, evec + i * lde, lde, avec, nrow, info->mvparam);
#endif
		wrk[1] = trl_ddot(nrow, evec+ i * lde, c__1, avec, c__1);
		trl_g_sum(info->mpicom, 2, wrk, &wrk[2]);
		info->matvec = info->matvec + 1;
		info->flop = info->flop + 4 * nrow;
		if (wrk[0] > 0.0) {
			eres[i] = wrk[1] / wrk[0];
			d__1 = -eres[i];
			trl_daxpy(nrow, d__1, evec + i * lde, c__1, avec, c__1);
			wrk[1] = trl_ddot(nrow, avec, c__1, avec, c__1);
			trl_g_sum(info->mpicom, 1, &wrk[1], &wrk[2]);
			wrk[0] = 1.0 / sqrt(wrk[0]);
			eres[ncol + i] = wrk[0] * sqrt(wrk[1]);
			d__1 = wrk[0];
			trl_dscal(nrow, d__1, evec + i * lde, c__1);
			info->flop = info->flop + 6 * nrow + 3;
		} else {
			eres[i] = -DBL_MAX;
			eres[ncol + i] = -DBL_MAX;
		}
	}
	if (base == NULL)
		free(avec);
	trl_close_logfile(info);
	/*
	// .. end of trl_rayleigh_quotients_ ..
	//
	*/
}


void
TRL::trl_ritz_projection(trl_matvec op, trl_info * info, int mev, double *evec,
					int lde, double *eres, double *wrk, int lwrk, double *base)
{
	/*
	// Purpose
	// =======
	// A separate Rayleigh-Ritz projection routine
	// Given a set of approximately orthonormal vectors (V), this routine
	// performs the following operations
	//  (1) V'*V ==> G
	//  (2) R'*R :=  G
	//  (3) V'*A*V => H1, inv(R')*H1*inv(R) => H
	//  (4) Y*D*Y' := H
	//  (5) V*inv(R)*Y => V, diag(D) => lambda,
	//      r(i) = ||A*V(:,i)-lambda(i)*V(:,i)||
	//
	// Arguments:
	// ==========
	// op       (input) function pointer
	//           On entry, points to the matrix-vector multiplication routine.
	//
	// info    (input) pointer to the structure trl_info_
	//          On entry, points to the data structure to store the information
	//          about the eigenvalue problem and the progress of TRLAN
	//
	// evec    (output) double precision vector of lenvth (nrow*mev)
	//          On exit, stores the eigenvectors.
	//
	// lde     (input) integer
	//          On entry, specifies the leading dimension of the array evec, i.e.,
	//
	// eres    (output) double precision vector of length (2*nev)
	//          the array to store new Ritz values and residual norms
	//
	// base    (workspace) double precision vector (nrow)
	//          Workspace to store the result of matrix-vector operation.
	//
	// wrk     (workspace) double precision vector of length (lwrk)
	//          Workspace to store projection matrix, etc.
	//
	// local variables
	*/
	//extern int dsyev_();
	//extern void dpotrf_();
	//extern void dtrtrs_();

	char trans = 'T', notrans = 'N', upper = 'U', job = 'V';
	double one = 1.0, zero = 0.0;
	int i__1 = 1;
	int i, j, ierr, nev, nsqr, nrow, iuau, irvv, lwrk2;
	double d__1;
	double *rvv, *uau, *wrk2, *avec;
	/*
	// ..
	// .. executable statements ..
	*/
	nrow = info->nloc;
	if (info->nec > 0) {
		nev = info->nec + 1;
	} else {
		nev = min(info->ned, mev - 1);
		if (info->lohi != 0)
			nev++;
	}
	nsqr = nev * nev;
	if (lwrk < 0) {
		lwrk = 0;
	}
	if (base != NULL) {
		avec = base;
	} else if (mev > nev) {
		avec = evec + (mev - 1) * lde;
	} else {
		avec = (double *) malloc(sizeof(double) * nrow);
	}
	if (info->verbose >= 0) {
		if (info->log_fp == NULL) {
			trl_reopen_logfile(info);
		}
		fprintf(info->log_fp,
			"TRLAN performing a Rayleigh-Ritz project for %d vectors.",
			nev);
	}
	/* memory allocation -- need 3*nev*nev elements, will allocate them     */
	/* in two consecutive blocks, uau(nev*nev), rvv(2*nev*nev)              */
	/* in actual use, rvv is further split in two until the last operation  */
	iuau = nsqr;
	irvv = nsqr + nsqr;
	if (lwrk >= iuau + irvv) {
		uau = wrk;
		rvv = &wrk[nsqr];
		wrk2 = &wrk[nsqr + nsqr];
		lwrk2 = lwrk - nsqr - nsqr;
	} else if (lwrk >= irvv) {
		rvv = wrk;
		wrk2 = &wrk[nsqr];
		lwrk2 = lwrk - nsqr;
		uau = (double *) malloc(nsqr * sizeof(double));
		if (uau == NULL) {
			info->stat = -231;
			goto end;
		}
	} else if (lwrk >= iuau) {
		uau = wrk;
		rvv = (double *) malloc((nsqr + nsqr) * sizeof(double));
		if (rvv == NULL) {
			info->stat = -232;
			goto end;
		}
		wrk2 = &rvv[nsqr];
		lwrk2 = nsqr;
	} else {
		uau = (double *) malloc(nsqr * sizeof(double));
		if (uau == NULL) {
			info->stat = -231;
			goto end;
		}
		rvv = (double *) malloc((nsqr + nsqr) * sizeof(double));
		if (rvv == NULL) {
			info->stat = -232;
			goto end;
		}
		wrk2 = &rvv[nsqr];
		lwrk2 = nsqr;
	}
	/* step (1) : V'*V ==> G */

	trl_dgemm(&trans, &notrans, nev, nev, nrow, one, evec, lde, evec, lde,
		zero, rvv, nev);
	trl_g_sum(info->mpicom, nsqr, rvv, wrk2);

	/* step (2) : Choleskey factorization of G */
	TRL::dpotrf_(&upper, &nev, rvv, &nev, &ierr);
	if (ierr != 0) {
		info->stat = -234;
		goto end;
	}
	/* step (3) : compute H_1 = V'*A*V                              */
	/* use the first nrow elements of avec to store the results of  */
	/* matrix-vector multiplication                                 */
	memset(wrk2, 0, lwrk2 * sizeof(double));
	for (i = 1; i <= nev; i++) {
#ifdef TRL_FORTRAN_COMPATIBLE
		op(&nrow, &i__1, evec + (i-1)*lde, &lde, avec, &nrow);
#else
		op(nrow, i__1, evec + (i-1)*lde, lde, avec, nrow, info->mvparam);
#endif
		trl_dgemv(&trans, nrow, i, one, evec, lde, avec, i__1, zero,
			&wrk2[(i - 1) * nev], i__1);
	}
	trl_g_sum(info->mpicom, nsqr, wrk2, uau);
	for (i = 1; i < nev; i++) {
		for (j = 0; j < i; j++) {
			wrk2[i + j * nev] = wrk2[(i - 1) * nev + j];
		}
	}
	/* compute solution of R^T H_2 = H_1 */
	TRL::dtrtrs_(&upper, &trans, &notrans, &nev, &nev, rvv, &nev, wrk2, &nev,
		&ierr);
	if (ierr != 0) {
		info->stat = -235;
		goto end;
	}
	/* compute solution of R^T H = H_2^T */
	for (i = 1; i < nev; i++) {
		for (j = 0; j < nev; j++) {
			uau[i + j * nev] = wrk2[(i - 1) * nev + j];
		}
	}
	TRL::dtrtrs_(&upper, &trans, &notrans, &nev, &nev, rvv, &nev, uau, &nev,
		&ierr);
	if (ierr != 0) {
		info->stat = -236;
		goto end;
	}
	/* solve the small symmetric eigenvalue problem */
	TRL::dsyev_(&job, &upper, &nev, uau, &nev, eres, wrk2, &nsqr, &ierr);
	if (ierr != 0) {
		info->stat = -237;
		goto end;
	}
	/* solve R Y = Y to prepare for multiplying with V */
	TRL::dtrtrs_(&upper, &notrans, &notrans, &nev, &nev, rvv, &nev, uau, &nev,
		&ierr);
	if (ierr != 0) {
		info->stat = -238;
		goto end;
	}
	/* call trl_ritz_vector to do the final multiplication */
	if (lwrk >= 3 * nsqr) {
		wrk2 = &wrk[nsqr];
	} else if (lwrk >= nsqr + nsqr) {
		wrk2 = wrk;
	} else {
		wrk2 = rvv;
	}
	i = lwrk2;
	trl_ritz_vectors(nrow, 0, nev, uau, nev, evec, lde, nev, avec, nrow,
		0, wrk2, i);
	/* compute the residual norms */
	for (i = 0; i < nev; i++) {
#ifdef TRL_FORTRAN_COMPATIBLE
		op(&nrow, &i__1, evec + i * lde, &lde, avec, &nrow);
#else
		op(nrow, i__1, evec + i * lde, lde, avec, nrow, info->mvparam);
#endif
		d__1 = eres[i];
		trl_daxpy(nrow, d__1, evec + i * lde, i__1, avec, i__1);
		eres[nev + i] = trl_ddot(nrow, avec, i__1, avec, i__1);
	}
	trl_g_sum(info->mpicom, nev, &eres[nev], avec);
	for (i = nev; i < nev + nev; i++) {
		if (eres[i] > 0.0) {
			eres[i] = sqrt(eres[i]);
		} else {
			eres[i] = -DBL_MAX;
		}
	}
	if (info->lohi < 0) {
		for (i = nev - 1; i < nev + nev - 2; i++) {
			eres[i] = eres[i + 1];
		}
	} else if (info->lohi > 0) {
		for (i = 0; i < nev - 1; i++) {
			eres[i] = eres[i + 1];
			memcpy(evec + i * lde, evec + (i + 1) * lde, nrow);
		}
		for (i = nev - 1; i < nev + nev - 2; i++) {
			eres[i] = eres[i + 2];
		}
	}
end:
	if (lwrk < iuau) {
		free(uau);
		free(rvv);
	} else if (lwrk < irvv) {
		free(rvv);
	} else if (lwrk < iuau + irvv) {
		free(uau);
	}
}



//
///////oxu func/////////
int TRL::indchar(char *a, char b)
{
	/*
	// Purpose:
	// ========
	// Used in trl_pe_filename to look for the first occurence of a character b in
	// a string a.
	//
	// Arguments:
	// ==========
	// a       (input) character string
	//          On entry, contains the string to search for a character b.
	//
	// b       (input) chararcter
	//          On entry, specifies the character to look for.
	//
	// ..
	// .. local scalars ..
	*/
	char *t = strchr(a, b);
	/*
	// ..
	// .. executable statements ..
	*/
	if (t != NULL) {
		return (1 + (t - a));
	} else {
		return strlen(a) + 1;
	}
}

void TRL::trl_pe_filename(int nlen, char *filename, char *base, int my_rank, int npe)
{
	/*
	// Purpose
	// =======
	// Generates file name from a base name and the PE number.
	//
	// Arguments:
	// ==========
	// nlen         (input) integer
	//               On entry, specifies the size of filiename.
	//
	// filename     (output) character string of length <= nlen.
	//               On exit, specifies the file name.
	//
	// base         (input) character string
	//               On entry, specifies the leading part of the file name.
	//
	// my_rank      (input) integer
	//               On entry, specifies the PE number.
	//
	// npe          (input) integer
	//               On entry, specifies the number of processors.
	//
	// ..
	// .. local variable ..
	*/
	int lead, ndig, len, off;
	char *format;
	/*
	// ..
	// .. executable statements ..
	*/
	ndig = 1;
	lead = npe;
	while (lead > 9) {
		lead /= 10;
		++ndig;
	}
	len = indchar(base, ' ') - 1;
	if (nlen < len + ndig + 1) {
		fprintf(stderr,
			"error: not enough space for filename (%d+%d chars).\n",
			len, ndig);
		exit(0);
	}
	memset(filename, 0, nlen * sizeof(char));
	strncpy(filename, base, len);
	off = 1 + ndig % 10;
	off = 5 + 2 * off;
	format = (char *) malloc(off * sizeof(char));
	sprintf(format, "%%s%%0%d.%dd", ndig, ndig);
	sprintf(filename, format, filename, my_rank);
	free(format);
	return;
	/*
	// .. end of trl_pe_filename ..
	*/
}


int TRL::close_file(FILE * fp, int err1, int err2)
{
	/*
	// Purpose
	// =======
	// Closes the file handler.
	//
	// Arguments
	// =========
	// fp      (input) pointer to file
	//          On entry, points to the file to be closed.
	//
	// err1    (input) integer
	//          On entry, specifies the return value when closing file failed.
	//
	// err2    (input) integer
	//          On entry, specifies the return value when closing file succeeded.
	//
	// ..
	// .. executable statements ..
	*/
	if (fclose(fp) != 0) {
		return err2;
	}
	return err1;
	/*
	// .. end of close_file ..
	*/
}

void TRL::trl_open_logfile(trl_info * info)
{
	/*
	// Purpose:
	// ========
	// Opens log file.
	//
	// Arguments:
	// ==========
	// info    (input) pointer to the structure trl_info_
	//          On entry, points to the data structure to store the information
	//          about the eigenvalue problem and the progress of TRLAN.
	//
	// ..
	// .. local arrays ..
	*/
	char filename[STRING_LEN];
	/*
	// ..
	// .. executable statements ..
	*/
	if (info->log_file != 0 && strlen(info->log_file) > 0) {
		trl_pe_filename(STRING_LEN, filename, info->log_file, info->my_pe,
			info->npes);
		info->log_fp = fopen(filename, "w");
	} else {
		info->log_fp = stdout;
	}
	/*
	// .. end of trl_open_logfilie_ ..
	*/
}

void TRL::trl_reopen_logfile(trl_info * info)
{
	/*
	// Purpose:
	// ========
	// Opens log file.
	//
	// Arguments:
	// ==========
	// info    (input) pointer to the structure trl_info_
	//          On entry, points to the data structure to store the information
	//          about the eigenvalue problem and the progress of TRLAN.
	//
	// ..
	// .. local arrays ..
	*/
	char filename[STRING_LEN];
	/*
	// ..
	// .. executable statements ..
	*/
	if (info->log_file != 0 && strlen(info->log_file) > 0) {
		trl_pe_filename(STRING_LEN, filename, info->log_file, info->my_pe,
			info->npes);
		info->log_fp = fopen(filename, "a");
	} else {
		info->log_fp = stdout;
	}
	/*
	// .. end of trl_open_logfilie_ ..
	*/
}

void TRL::trl_close_logfile(trl_info * info)
{
	/*
	// Purpose:
	// ========
	// Closes log file.
	//
	// Arguments:
	// ==========
	// info    (input) pointer to the structure trl_info_
	//          On entry, points to the data structure to store the information
	//          about the eigenvalue problem and the progress of TRLAN.
	//
	// ..
	// .. executable statements ..
	*/
	if (info->log_fp != NULL && info->log_fp != stdout) {
		fclose(info->log_fp);
	}
	info->log_fp = NULL;
	/*
	// .. end of trl_close_logfile ..
	*/
}



void TRL::trl_time_stamp(FILE * fp)
{
	/*
	// Purpose
	// =======
	// Print the current date and time.
	//
	// Arguments:
	// ==========
	// fp         (input) pointer to a file
	//             On entry, points to the file to output.
	//
	// ..
	// .. local variables ..
	*/
	time_t clck;
	/*
	// ..
	// .. executable statements ..
	*/
	clck = time(NULL);
	fprintf(fp, "                                                  %s",
		asctime(localtime(&clck)));
	/*
	// .. end of trl_time_stamp ..
	*/
}

int TRL::trl_read_checkpoint(char *filename, int nrow, double *evec, int lde,
						int mev, int *j1, double *base, int ldb, int nbas,
						int *j2, int nalpha, double *alpha, int nbeta,
						double *beta)
{
	/*
	// Purpose
	// =======
	// Read check-point file
	//
	// Arguments:
	// ==========
	// filename   (input) character string
	//             On entry, specifies the name of checkpoint file.
	//
	// nrow       (input) integer
	//             On entry, specifies the problem size.
	//
	// evec       (output) double precision array of dimensioni (lde,j1)
	//             On exit, contains the first part of basis vectors stored in checkpoint.
	//
	// lde         (input) integer
	//             On entry, specifies the leading dimension of evec.
	//
	// mev        (input) integer
	//             On entry, specifies the number of eigenvalues converged.
	//
	// j1         (input) integer
	//             On entry, specifies the last column index of evec, that contains a base vector.
	//
	// base       (output) double precision array of dimension (ldb,nbas)
	//             On exit, contains the second part of basis stored in the checkpoint.
	//
	// ldb        (input) integer
	//             On entry, specifies the leading dimension of base.
	//
	// nbas       (input) integer
	//             On entry, specifies the number of columns in base.
	//
	// j2         (input) integer
	//             On entry, specifies the last column index of base, that contains a base vector.
	//
	// nalpha     (input) integer
	//             On entry, specifies the size of alpha
	//
	// alpha      (output) double precision array of length (nalpha)
	//             On exit, contains the alpha values stored in checkpoint.
	//
	// nbeta      (input) integer
	//             On entry, specifies the size of beta.
	//
	// beta       (output) double precision array of length (nbeta)
	//             On exit, contains the beta values stored in checkpoint.
	//
	// ..
	// .. local variables ..
	*/
	int i, j;
	FILE *io_fp;
	/*
	// ..
	// .. executable statements ..
	*/
	if (lde < nrow || ldb < nrow) {
		printf("TRL_READ_CHECKPOINT: leading dimensions too small.\n");
		return -211;
	}
	/* open file */
	io_fp = fopen(filename, "r");
	if (io_fp == NULL) {
		printf
			("TRL_READ_CHECKPOINT: failed to open check-point file %s.\n",
			filename);
		return -212;
	}
	/* read size information */
	if (fread(j1, sizeof(*j1), 1, io_fp) <= 0) {
		return close_file(io_fp, -215, -216);
	}
	if (fread(j2, sizeof(*j2), 1, io_fp) <= 0) {
		return close_file(io_fp, -215, -216);
	}
	if (*j1 != nrow) {
		printf("TRL_READ_CHECKPOINT: Nrow mismatch.\n");
		return -213;
	}
	if (*j2 > min(nalpha, min(nbeta, mev + nbas - 1))) {
		printf("TRL_READ_CHECKPOINT: MAXLAN too small.");
		return -214;
	}
	/* can continue read all data */
	for (i = 0; i < *j2; i++) {
		if (fread(&alpha[i], sizeof(double), 1, io_fp) <= 0) {
			return close_file(io_fp, -215, -216);
		}
	}
	for (i = 0; i < *j2; i++) {
		if (fread(&beta[i], sizeof(double), 1, io_fp) <= 0) {
			return close_file(io_fp, -215, -216);
		}
	}
	*j1 = min(mev, *j2);
	*j2 = *j2 - *j1;
	if (*j1 < mev) {
		for (i = 0; i <= *j1; i++) {
			for (j = 0; j < nrow; j++) {
				if (fread
					(&evec[i * lde + j], sizeof(double), 1,
					io_fp) <= 0) {
						return close_file(io_fp, -215, -216);
				}
			}
		}
	} else {
		for (i = 0; i < *j1; i++) {
			for (j = 0; j < nrow; j++) {
				if (fread
					(&evec[i * nrow + j], sizeof(double), 1,
					io_fp) <= 0) {
						return close_file(io_fp, -215, -216);
				}
			}
		}
		for (i = 0; i < *j2; i++) {
			for (j = 0; j < nrow; j++) {
				if (fread
					(&base[i * ldb + j], sizeof(double), 1, io_fp) <= 0) {
						return close_file(io_fp, -215, -216);
				}
			}
		}
	}
	return close_file(io_fp, 0, -215);
	/*
	// .. end of trl_read_checkpoint ..
	*/
}


int TRL::trl_write_checkpoint(char *filename, int nrow, double *alpha,
						 double *beta, double *evec, int lde, int me,
						 double *base, int ldb, int nb)
{
	/*
	// Purpose
	// =======
	// Write a check-point file.
	//
	// Arguments
	// =========
	// filename   (input) character string
	//             On entry, specifies the name of checkpoint file.
	//
	// nrow       (input) integer
	//             On entry, specifies the problem size.
	//
	// alpha      (input) double precision array of length (me+nb-1)
	//             On entry, contains the alpha values computed so far.
	//
	// beta       (input) double precision array of length (me+ne-1)
	//             On entry, contains the beta values computed so far.
	//
	// evec       (input) double precision array of dimensioni (lde,me)
	//             On entry, contains the first part of basis vectors.
	//
	// lde        (input) integer
	//             On entry, specifies the leading dimension of evec.
	//
	// me         (input) integer
	//             On entry, specifies the last column index of evec, that contains a base vector.
	//
	// base       (input) double precision array of dimension (ldb,nb)
	//             On entry, contains the second part of basis.
	//
	// ldb        (input) integer
	//             On entry, specifies the leading dimension of base.
	//
	// nb         (input) integer
	//             On entry, specifies the last column index of base, that contains a base vector.
	//
	// ..
	// .. local variables ..
	*/
	int jnd, i, j;
	FILE *io_fp;
	/*
	// ..
	// .. executable statements ..
	*/
	jnd = me + nb - 1;
	io_fp = fopen(filename, "w");
	if (io_fp == NULL) {
		printf("TRL_WRITE_CHECKPOINT: failed to open file: %s.\n",
			filename);
		return -221;
	}
	if (fwrite(&nrow, sizeof(nrow), 1, io_fp) < 1) {
		return close_file(io_fp, -223, -222);
	}
	if (fwrite(&jnd, sizeof(jnd), 1, io_fp) < 1) {
		return close_file(io_fp, -223, -222);
	}

	for (i = 0; i < jnd; i++) {
		if (fwrite(&alpha[i], sizeof(double), 1, io_fp) < 1) {
			return close_file(io_fp, -223, -222);
		}
	}
	for (i = 0; i < jnd; i++) {
		if (fwrite(&beta[i], sizeof(double), 1, io_fp) < 1) {
			return close_file(io_fp, -223, -222);
		}
	}
	for (i = 0; i < me; i++) {
		for (j = 0; j < nrow; j++) {
			if (fwrite
				(&evec[i * lde + j], sizeof(double), 1,
				io_fp) < 1) {
					return close_file(io_fp, -223, -222);
			}
		}
	}
	for (i = 0; i < nb; i++) {
		for (j = 0; j < nrow; j++) {
			if (fwrite
				(&base[i * ldb + j], sizeof(double), 1,
				io_fp) < 1) {
					return close_file(io_fp, -223, -222);
			}
		}
	}
	return close_file(io_fp, 0, -223);
	/*
	// .. end of trl_write_checkpoint ..
	*/
}

void TRL::trl_print_int(trl_info * info, char *title, int size_array,
				   int *array, int inc)
{
	/*
	// Purpose:
	// ========
	// Print an integer_ array for debugging.
	//
	// Arguments:
	// ==========
	// info        (input) pointer to the structure trl_info_
	//              On entry, points to the data structure to store the information
	//              about the eigenvalue problem and the progress of TRLAN
	//
	// title       (input) character string
	//              On entry, specifies the title of the information to be printed.
	//
	// size_array  (input) integer
	//              On entry specifies, the number of integers to be printed.
	//
	// array       (input) integer_ array of length ((size_array-1)*inc+1)
	//              On entry, contains the integer_ to be printed.
	//
	// inc         (input) integer
	//              On entry, specifies how the index to array should be incremented.
	//
	// ..
	// .. local scalars ..
	*/
	int i;
	/*
	// ..
	// .. executable statements ..
	*/
	fprintf(info->log_fp, "PE %d : %s", info->my_pe, title);
	if (size_array > 2) {
		fprintf(info->log_fp, "\n");
	}
	for (i = 0; i < size_array; i += inc) {
		fprintf(info->log_fp, "%10d", array[i]);
		if ((i % 8) == 7)
			fprintf(info->log_fp, "\n");
	}
	if (((size_array - 1) % 8) != 7)
		fprintf(info->log_fp, "\n");
	/*
	// .. end of trl_print_int ..
	//
	*/
}

void TRL::trl_print_real(trl_info * info, char *title, int size_array,
					double *array, int inc)
{
	/*
	// Purpose
	// =======
	// Print a double precision array for debugging.
	//
	// Arguments:
	// ==========
	// info        (input) pointer to the structure trl_info_
	//              On entry, points to the data structure to store the information
	//              about the eigenvalue problem and the progress of TRLAN
	//
	// title       (input) character string
	//              On entry, specifies the title of the information to be printed.
	//
	// size_array  (input) integer
	//              On entry specifies, the number of doubles to be printed.
	//
	// array       (input) double array of length ((size_array-1)*inc+1)
	//              On entry, contains the doubles to be printed.
	//
	// inc         (input) integer
	//              On entry, specifies how the index to array should be incremented.
	//
	// ..
	// .. local scalars ..
	*/
	int i;
	/*
	// ..
	// .. executable statements ..
	*/
	fprintf(info->log_fp, "PE %d : %s", info->my_pe, title);
	if (size_array > 1) {
		fprintf(info->log_fp, "\n");
	}
	for (i = 0; i < size_array; i += inc) {
		fprintf(info->log_fp, " %10.7e", array[i]);
		if ((i % 8) == 7)
			fprintf(info->log_fp, "\n");
	}
	if (((size_array - 1) % 8) != 7)
		fprintf(info->log_fp, "\n");
	/*
	// .. end of trl_print_real ..
	*/
}

void TRL::trl_print_progress(trl_info * info)
{
	/*
	// Purpose
	// =======
	// Print the current progress of eigenvalue solution
	//
	// Arguments:
	// ==========
	// info        (input) pointer to the structure trl_info_
	//              On entry, points to the data structure to store the information
	//              about the eigenvalue problem and the progress of TRLAN
	//
	// ..
	// .. executable statements ..
	*/
	fprintf(info->log_fp, "MATVEC: %10d,    Nloop: %10d,     Nec: %10d\n",
		info->matvec, info->nloop, info->nec);
	fprintf(info->log_fp, "Reorth: %10d,    Nrand: %10d,    Ierr: %10d\n",
		info->north, info->nrand, info->stat);
	fprintf(info->log_fp,
		"Target: %10.3e,   ResNrm: %10.3e,    CFact: %10.3e\n",
		info->trgt, info->tres, info->crat);
	/*
	// .. end of trl_print_progress ..
	*/
}


void TRL::trl_check_orth(trl_info * info, int nrow, double *v1, int ld1,
					int j1, double *v2, int ld2, int j2, double *wrk,
					int lwrk)
{
	/*
	// Purpose:
	// ========
	// Check orthogonality of the basis.
	//
	// Arguments:
	// ==========
	// info     (input) pointer to the structure trl_info_
	//           On entry, points to the data structure to store the information 
	//           about the eigenvalue problem and the progress of TRLAN.
	//
	// nrow     (input) integer
	//           On entry, specifies the problem size, i.e., the number of rows in 
	//           v1 and v2.
	//
	// v1       (input) double precision array of diimension (ld1,j1)
	//           On entry, contains the first part of the basis.
	//
	// ld1      (input) integer
	//           On entry, specifies the leading diimension of v1.
	//
	// j1       (input) integer
	//           On entry, specifies the last column index of v1, containing the 
	//           basis.
	//
	// v2       (input) double precision array of dimension (ld2,j2)
	//           On entry, contains the second part of the basis.
	//
	// ld2      (input) integer
	//           On entry, specifies the leading dimension of v2.
	//
	// j2       (input) integer
	//           On entry, specifies the last column index of v2, containing the 
	//           basis.
	//
	// wrk      (workspace) double precision array of length (lwrk)
	//
	// lwrk     (input) integer
	//           On entry, specifies the size of workspace.
	//
	// ..
	// .. local parameters ..
	*/
	double one = 1.0, zero = 0.0;
	long c__1 = 1;
	/*
	// ..
	// .. local variables
	*/
	int i, j, k, jnd;
	double nrmfro, nrminf;
	/*
	// ..
	// .. executable statements ..
	*/
	jnd = j1 + j2;
	nrmfro = zero;
	nrminf = zero;
	if (jnd <= 0)
		return;
	if (lwrk < (jnd + jnd)) {
		fprintf(info->log_fp, "TRL_CHECK_ORTH: workspace too small.\n");
		return;
	}
	fprintf(info->log_fp,
		"TRL_CHECK_ORTH: check orthogonality of %d basis vectors.\n",
		jnd);
	/*
	// check orthognality of the basis vectors
	*/
	for (i = 0; i < j1; i++) {
		trl_g_dot_(info->mpicom, nrow, v1, ld1, i + 1, v2, ld2, 0,
			v1 + i * ld1, wrk);
		wrk[i] = wrk[i] - one;
		if (info->verbose > 7) {
			fprintf(info->log_fp, "Orthogonality level of v(%d) ..\n",
				i + 1);
			for (j = 0; j <= i; j++) {
				fprintf(info->log_fp, " %10.3e", wrk[j]);
				if ((j % 8) == 7)
					fprintf(info->log_fp, "\n");
			}
			if ((i % 8) != 7)
				fprintf(info->log_fp, "\n");
		}
		nrmfro =
			nrmfro + 2 * trl_ddot(i, wrk, c__1, wrk,
			c__1) + wrk[i] * wrk[i];
		if (i == 0) {
			wrk[i + 1] = fabs(wrk[i]);
		} else {
			wrk[i + 1] = max(wrk[i], wrk[i - 1]);
		}
		nrminf = max(nrminf, wrk[i + 1]);
	}
	for (i = 0; i < j2; i++) {
		j = j1 + i;
		trl_g_dot_(info->mpicom, nrow, v1, ld1, j1, v2, ld2, i + 1,
			v2 + i * ld2, wrk);
		wrk[j] = wrk[j] - one;
		if (info->verbose > 7) {
			fprintf(info->log_fp, "Orthogonality level of v(%d) ..\n",
				j + 1);
			for (k = 0; k <= j; k++) {
				fprintf(info->log_fp, " %10.3e", wrk[k]);
				if ((k % 8) == 7)
					fprintf(info->log_fp, "\n");
			}
			if ((j % 8) != 7)
				fprintf(info->log_fp, "\n");
		}
		nrmfro =
			nrmfro + 2 * trl_ddot(j, wrk, c__1, wrk,
			c__1) + wrk[j] * wrk[j];
		nrminf = max(nrminf, fabs(wrk[j]));
	}
	fprintf(info->log_fp,
		"Frobenius norm of orthogonality level %10i %4i  %14.5e\n",
		info->matvec, jnd, sqrt(nrmfro));
	fprintf(info->log_fp,
		"Maximum abs. value of orthogonality level is  %14.5e\n",
		nrminf);
	/*
	// .. end of trl_check_orth ..
	*/
}

void TRL::trl_check_recurrence(trl_matvec op,
					 trl_info * info, int nrow, int ncol, double *v1,
					 int ld1, int m1, double *v2, int ld2, int m2,
					 int kept, double *alpha, double *beta, double *wrk,
					 int lwrk)
{
	/*
	// Purpose
	// =======
	// Check Lanczos recurrence relation for debug purpose.
	//
	// Arguments:
	// ==========
	// op       (input) function pointer
	//           On entry, points to the matrix-vector multiplication routine.
	//
	// info     (input) pointer to the structure trl_info_
	//           On entry, points to the data structure to store the information about 
	//           the eigenvalue problem and the progress of TRLAN.
	//
	// nrow     (input) integer
	//           On entry, specifies the problem size, i.e., the number of ros in v1 
	//           and v2.
	//
	// ncol     (input) integer
	//           On entry, specifies the maximum number of eigenvalues that can be 
	//           stored.
	//
	// v1       (input) double precision array of dimension (ld1,m1)
	//           On entry, contains the first part of basis.
	//
	// ld1      (input) integer
	//           On entry, specifies the leading dimension of v1.
	//
	// m1       (input) integer
	//           On entry, specifies the last column index of v1 that contains a 
	//           base vector.
	//
	// v2       (input) double precision array of dimension (ld2,m2)
	//           On entry, contains the second part of basis.
	//
	// ld2      (input) integer
	//           On entry, specifies the leading dimension of v2.
	//
	// m2       (input) integer
	//           On entry, specifies the last column index of v2 that contains a 
	//           base vector.
	//
	// kept     (input) integer
	//           On entry, specifies the number of basis kept at the last restart.
	//
	// alpha    (input) integer
	//           On entry, contains the alpha values computed so far.
	//
	// beta     (input) integer
	//           On entry, contains the beta values computed so far.
	//
	// wrk      (workspace) double precision vector of length (lwrk)
	//
	// lwrk     (input) integer
	//           On entry, specifies the size of workspace.
	//
	// ..
	// .. local parameters ..
	*/
	long c__1 = 1;
	int i__1 = 1;
	double zero = 0.0, one = 1.0;
	/*
	// ..
	// .. local variables ..
	*/
	int i, ii, j, j1, j2, jnd, mv1;
	char title[STRING_LEN];
	double d__1;
	double *aq, *qkp1, *cs, *alf, *bet;
	/*
	// ..
	// .. executable statements ..
	*/
	mv1 = m1;
	if (m2 > 0) {
		j2 = m2 - 1;
		j1 = m1;
	} else {
		j2 = 0;
		j1 = m1 - 1;
	}
	jnd = j1 + j2;
	if (lwrk < jnd * 4 + max(jnd * 4, nrow)) {
		fprintf(info->log_fp,
			"TRL_CHECK_RECURRENCE: not enough workspace.\n");
		return;
	}
	if (lwrk >= jnd * 4 + nrow) {
		aq = wrk + (lwrk - nrow);
	} else {
		aq = (double *) malloc(nrow * sizeof(double));
		if (aq == NULL) {
			fprintf(info->log_fp,
				"TRL_CHECK_RECURRENCE: failed to allcoate workspace.\n");
			return;
		}
	}
	memset(wrk, 0, 4 * jnd * sizeof(double));
	cs = &wrk[jnd];
	alf = &wrk[2 * jnd];
	bet = &wrk[3 * jnd];
	/*
	// first type of relation
	// A q_i = Alpha_i q_i + Beta_i q_{k+1}
	*/
	if (kept < ncol) {
		qkp1 = v1 + kept * ld1;
	} else {
		qkp1 = v2 + (kept - j1) * ld2;
	}
	for (i = 0; i < min(j1, kept); i++) {
#ifdef TRL_FORTRAN_COMPATIBLE
		op(&nrow, &i__1, v1 + i * ld1, &ld1, aq, &nrow);
#else
		op(nrow, i__1, v1 + i * ld1, ld1, aq, nrow, info->mvparam);
#endif
		for (ii = 0; ii < nrow; ii++) {
			alf[i] += aq[ii] * v1[i * nrow + ii];
			aq[ii] -= alpha[i] * v1[i * nrow + ii];
			bet[i] += aq[ii] * aq[ii];
			cs[i] += aq[ii] * qkp1[ii];
			aq[ii] -= beta[i] * qkp1[ii];
			wrk[i] += aq[ii] * aq[ii];
		}
	}
	for (i = 0; i < (kept - j1); i++) {
		j = i + j1;
#ifdef TRL_FORTRAN_COMPATIBLE
		op(&nrow, &i__1, v2 + i * ld2, &ld2, aq, &nrow);
#else
		op(nrow, i__1, v2 + i * ld2, ld2, aq, nrow, info->mvparam);
#endif
		for (ii = 0; ii < nrow; ii++) {
			alf[j] += aq[ii] * v2[i * nrow + ii];
			aq[ii] -= alpha[j] * v2[i * nrow + ii];
			bet[j] += aq[ii] * aq[ii];
			cs[j] += aq[ii] * qkp1[ii];
			aq[ii] -= beta[j] * qkp1[ii];
			wrk[j] += aq[ii] * aq[ii];
		}
	}
	/*
	// the (k+1)st base vector need to orthogonalize against all previous
	// vectors
	*/
	if (jnd > kept) {
#ifdef TRL_FORTRAN_COMPATIBLE
		op(&nrow, &i__1, qkp1, &nrow, aq, &nrow);
#else
		op(nrow, i__1, qkp1, nrow, aq, nrow, info->mvparam);
#endif
		alf[kept] = trl_ddot(nrow, aq, c__1, qkp1, c__1);
		d__1 = -alpha[kept];
		trl_daxpy(nrow, d__1, qkp1, c__1, aq, c__1);
		for (i = 0; i < min(j1, kept); i++) {
			d__1 = -beta[i];
			trl_daxpy(nrow, d__1, v1 + i * ld1, c__1, aq, c__1);
		}
		for (i = 0; i < kept - j1; i++) {
			j = j1 + i;
			d__1 = -beta[j];
			trl_daxpy(nrow, d__1, v2+ i * ld2, c__1, aq, c__1);
		}
		bet[kept] = trl_ddot(nrow, aq, c__1, aq, c__1);
		if (kept + 2 <= j1) {
			cs[kept] =
				trl_ddot(nrow, aq, c__1, v1 + (kept + 1) * ld1, c__1);
			d__1 = -beta[kept];
			trl_daxpy(nrow, d__1, v1 + (kept + 1) * ld1, c__1, aq, c__1);
		} else {
			cs[kept] =
				trl_ddot(nrow, aq, c__1, v2 + (kept + 1 - j1) * ld2,
				c__1);
			d__1 = -beta[kept];
			trl_daxpy(nrow, d__1, v2 + (kept + 1 - j1) * ld2, c__1, aq,
				c__1);
		}
		wrk[kept] = trl_ddot(nrow, aq, c__1, aq, c__1);
	}
	/*
	// the third kind of relation -- normal three term recurrence
	// depending the fact that if the lower-bound of loop is less than
	// upper bound, the look should not be executed
	*/
	for (i = kept + 1; i < j1; i++) {
#ifdef TRL_FORTRAN_COMPATIBLE
		op(&nrow, &i__1, v1 + ld1 * i, &ld1, aq, &nrow);
#else
		op(nrow, i__1, v1 + ld1 * i, ld1, aq, nrow, info->mvparam);
#endif
		if (i < (mv1 - 1)) {
			for (ii = 0; ii < nrow; ii++) {
				alf[i] += aq[ii] * v1[i * ld1 + ii];
				aq[ii] -=
					(alpha[i] * v1[i * ld1 + ii] +
					beta[i - 1] * v1[(i - 1) * nrow + ii]);
				bet[i] += aq[ii] * aq[ii];
				cs[i] += aq[ii] * v1[(i + 1) * ld1 + ii];
				aq[ii] -= beta[i] * v1[(i + 1) * ld1 + ii];
				wrk[i] += aq[ii] * aq[ii];
			}
		} else {
			for (ii = 0; ii < nrow; ii++) {
				alf[i] += aq[ii] * v1[i * ld1 + ii];
				aq[ii] -=
					(alpha[i] * v1[i * ld1 + ii] +
					beta[i - 1] * v1[(i - 1) * ld1 + ii]);
				bet[i] += aq[ii] * aq[ii];
				cs[i] += aq[ii] * v2[ii];
				aq[ii] -= beta[i] * v2[ii];
				wrk[i] += aq[ii] * aq[ii];
			}
		}
	}
	for (i = max(0, kept - j1 + 1); i < j2; i++) {
		j = i + j1;
#ifdef TRL_FORTRAN_COMPATIBLE
		op(&nrow, &i__1, v2 + i * ld2, &ld2, aq, &nrow);
#else
		op(nrow, i__1, v2 + i * ld2, ld2, aq, nrow, info->mvparam);
#endif
		if (i > 0) {
			for (ii = 0; ii < nrow; ii++) {
				alf[j] += aq[ii] * v2[i * ld2 + ii];
				aq[ii] -=
					(beta[j - 1] * v2[(i - 1) * ld2 + ii] +
					alpha[j] * v2[i * ld2 + ii]);
				bet[j] += aq[ii] * aq[ii];
				cs[j] += aq[ii] * v2[(i + 1) * ld2 + ii];
				aq[ii] -= beta[j] * v2[(i + 1) * ld2 + ii];
				wrk[j] += aq[ii] * aq[ii];
			}
		} else {
			for (ii = 0; ii < nrow; ii++) {
				alf[j] += aq[ii] * v2[ii];
				aq[ii] -=
					(beta[j - 1] * v1[(j1 - 1) * ld1 + ii] +
					alpha[j] * v2[ii]);
				bet[j] += aq[ii] * aq[ii];
				cs[j] += aq[ii] * v2[ld2 + ii];
				aq[ii] -= beta[j] * v2[ld2 + ii];
				wrk[j] += aq[ii] * aq[ii];
			}
		}
	}

	trl_g_sum(info->mpicom, jnd * 4, wrk, &wrk[jnd * 4]);
	aq[0] = zero;
	for (ii = 0; ii < jnd; ii++) {
		aq[0] += wrk[ii];
	}
	aq[0] = sqrt(aq[0]);
	for (ii = 0; ii < jnd; ii++) {
		wrk[ii] = sqrt(wrk[ii]);
	}
	for (ii = 0; ii < jnd; ii++) {
		if (bet[ii] > zero) {
			if (beta[ii] < zero) {
				bet[ii] = -sqrt(bet[ii]);
			} else {
				bet[ii] = sqrt(bet[ii]);
			}
			cs[ii] = cs[ii] / bet[ii];
		} else {
			bet[ii] = zero;
		}
	}
	strcpy(title, "Alpha computed by TRLAN ..");
	trl_print_real(info, title, jnd, alpha, 1);
	strcpy(title, "Alpha computed explicitly in TRL_CHECK_RECURRENCE ..");
	trl_print_real(info, title, jnd, alf, 1);
	strcpy(title, "Differences in alpha ..");
	d__1 = -one;
	trl_daxpy(jnd, d__1, alpha, c__1, alf, c__1);
	trl_print_real(info, title, jnd, alf, 1);
	strcpy(title, "Beta computed by TRLAN ..");
	trl_print_real(info, title, jnd, beta, 1);
	strcpy(title, "Beta computed explicitly in TRL_CHECK_RECURRENCE ..");
	trl_print_real(info, title, jnd, bet, 1);
	strcpy(title, "Differences in beta ..");
	d__1 = -one;
	trl_daxpy(jnd, d__1, beta, c__1, bet, c__1);
	trl_print_real(info, title, jnd, bet, 1);
	strcpy(title, "Error in Lanczos recurrence (overall) =");
	trl_print_real(info, title, 1, aq, 1);
	if (info->verbose > 7) {
		strcpy(title,
			"|| A q_i - alpha_i q_i - beta_{i-1} q_{i-1} - beta_i q_{i+1} ||..");
		trl_print_real(info, title, jnd, wrk, 1);
		strcpy(title,
			"(A q_i - alpha_i q_i - beta_{i-1} q_{i-1})*q_{i+1}/beta_i ..");
		trl_print_real(info, title, jnd, cs, 1);
		strcpy(title, "Sine of the angles ..");
		for (ii = 0; ii < jnd; ii++) {
			cs[ii] = cs[ii] * cs[ii];
			if (cs[ii] < one) {
				cs[ii] = sqrt(one - cs[ii]);
			} else {
				cs[ii] = -one;
			}
		}
		trl_print_real(info, title, jnd, cs, 1);
	}
	if (lwrk < jnd * 4 + nrow)
		free(aq);
	/*
	// .. end of trl_check_recurrence ..
	*/
}

void TRL::trl_open_cptfile(trl_info * info)
{
	/*
	// Purpose:
	// ========
	// Opens check points file.
	//
	// Arguments:
	// ==========
	// info    (input) pointer to the structure trl_info_
	//          On entry, points to the data structure to store the information
	//          about the eigenvalue problem and the progress of TRLAN.
	//
	// ..
	// .. local arrays ..
	*/
	char filename[STRING_LEN];
	/*
	// ..
	// .. executable statements ..
	*/
	if (info->cpfile != 0 && strlen(info->cpfile) > 0) {
		trl_pe_filename(STRING_LEN, filename, info->cpfile, info->my_pe,
			info->npes);
		info->cpt_fp = fopen(filename, "w");
	} else {
		info->cpt_fp = stdout;
	}
	/*
	// .. end of trl_open_cptfile ..
	*/
}

void TRL::trl_close_cptfile(trl_info * info)
{
	/*
	// Purpose:
	// ========
	// Closes check point file.
	//
	// Arguments:
	// ==========
	// info    (input) pointer to the structure trl_info_
	//          On entry, points to the data structure to store the information
	//          about the eigenvalue problem and the progress of TRLAN.
	//
	// ..
	// .. executable statements ..
	*/
	if (info->cpt_fp != stdout) {
		fclose(info->cpt_fp);
	}
	info->cpt_fp = NULL;
	/*
	// .. end of trl_close_cptfile ..
	*/
}

///////////////
/////restart

void TRL::trl_shuffle_eig(int nd, int mnd, double *lambda, double *res,
					 trl_info * info, int *kept, int locked)
{
	// ..
	// .. local scalars ..
	int i, ncl, ncr, kl, kr, tind, minsep;
	double bnd;
	//
	// ..
	// .. executable statements ..
	// very small basis -- save the half with the smallest residual norms
	if (nd <= 5) {
		dsort2(nd, res, lambda);
		if (nd > 3) {
			*kept = 2;
		} else if (nd > 0) {
			*kept = 1;
		} else {
			*kept = 0;
		}
		if (*kept > 1)
			dsort2(*kept, lambda, res);
		return;
	}
	//
	// preparation for normal case, first determine what are converged.
	//    ncl are the index (base zero) of res converged from the left
	//    ncr are the index (base zero) of res converged from the right
	bnd = min(info->tol, DBL_EPSILON * info->anrm);
	ncr = 0;
	ncl = nd - 1;
	i = nd - 1;
	// determine how many has converged from the right
	while (i >= 0) {
		if (res[i] <= bnd) {
			i--;
		} else {
			ncr = i + 1;
			i = -1;
		}
	}
	i = 0;
	// determine how many has converged from the left
	while (i < nd) {
		if (res[i] <= bnd) {
			i++;
		} else {
			ncl = i - 1;
			i = nd;
		}
	}
	kl = ncl;
	kr = ncr;
	if (ncr > ncl) {
		// find the one that is closest to info->trgt
		tind = (kl + kr) / 2;
		while (lambda[tind] != info->trgt && kr > kl) {
			if (lambda[tind] < info->trgt) {
				kl = tind + 1;
				tind = (kl + kr) / 2;
			} else if (lambda[tind] > info->trgt) {
				kr = tind - 1;
				tind = (kl + kr) / 2;
			} else {
				kl = tind;
				kr = tind;
			}
		}
		// assign kl to the largest index of lambda that is smaller than
		// info->trgt
		if (lambda[tind] == info->trgt) {
			kl = tind - 1;
			while (kl >= 0 && lambda[kl] == info->trgt) {
				kl--;
			}
			// assign kr to the smallest index of lambda that is greater than
			// info->trgt
			kr = tind + 1;
			while (kr < nd && lambda[kr] == info->trgt) {
				kr++;
			}
		} else {
			kl = tind - 1;
			kr = tind + 1;
		}
		// initial assignment of kl and kr
		if (info->lohi > 0) {
			// large eigenvalues.
			kr = kl;
			kl = min(ncl, max(0, nd - info->ned) - 1);
		} else if (info->lohi < 0) {
			// small eigenvalues.
			kl = kr;
			kr = max(ncr, min(nd - info->nec, info->ned + 1) - 1);
		} else if (ncr - tind > tind - ncl) {
			// tind is closer to smallest lambda converged from left
			kl = kr;
			kr = max(ncr, min(nd - info->nec, info->ned + 1) - 1);
		} else {
			// tind is closer to largest lambda converged from right
			kr = kl;
			kl = min(ncl, max(0, nd - info->ned) - 1);
		}
	} else {
		// all have converged, keep all -- should not happen
		*kept = nd;
		return;
	}
	// We keep at least one from each end.
	/*
	if( kl < kr ) {
	kl ++;
	}
	if( kr > kl ) {
	kr --;
	}
	*/
	//
	// We also save Ritz values with the smallest residual
	/*
	if( kl < kr ) {
	if( res[kl+1] < res[kr-1] ) {
	kl ++;
	} else {
	kr --;
	}
	}
	*/
	//
	// normal cases, call subroutines to complete the tasks
	// [1 .. kl] and [kr .. nd] are saved for later
	// the initial values of kl and kr are simply ncl and ncr
	// they are further analyzed according to the restarting strategy
	// requested
	switch (info->restart) {
	case 0:
		if (info->restart <= -info->ned) {
			if (info->lohi >= 0) {
				kl = -1;
				kr = max(2, nd + info->restart);
			} else if (info->lohi < 0) {
				kl = min(-info->restart, nd - 3) - 1;
				kr = nd;
			} else {
				kl = min(nd - 3, -info->restart) / 2 - 1;
				kr = nd - kl - 2;
			}
		} else {
			trl_restart_fixed(nd, mnd, tind, lambda, res, info, &kl, &kr);
		}
		break;
	case 1:
		// fixed number beyond the currently converged ones
		trl_restart_fixed(nd, mnd, tind, lambda, res, info, &kl, &kr);
		break;
	case 2:
		// add the ones with smallest residual nroms
		trl_restart_small_res(nd, mnd, tind, lambda, res, info, &kl, &kr);
		break;
	case 3:
		if (info->nloop > 0) {
			// maximize the gap ratio
			trl_restart_max_gap_ratio(nd, tind, *kept, lambda, res, info,
				&kl, &kr);
		} else {
			// this is the first restart -- use trl_restart_fixed instead
			trl_restart_fixed(nd, mnd, tind, lambda, res, info, &kl, &kr);
		}
		break;
	case 4:
		if (info->nloop > 0) {
			// maximize [gap-ratio * (m-k)]
			trl_restart_max_progress(nd, tind, *kept, lambda, res, info,
				&kl, &kr);
		} else {
			// this is the first restart -- use trl_restart_fixed instead
			trl_restart_fixed(nd, mnd, tind, lambda, res, info, &kl, &kr);
		}
		break;
	case 5:
		if (info->nloop > 0) {
			// maximize [sqrt(gap tatio) * (m-k)]
			trl_restart_max_reduction(nd, tind, *kept, lambda, res, info,
				&kl, &kr);
		} else {
			// this is the first restart -- use trl_restart_fixed instead
			trl_restart_fixed(nd, mnd, tind, lambda, res, info, &kl, &kr);
		}
		break;
	case 6:
		// progressively vary the thickness
		trl_restart_scan(nd, res, info, *kept, &kl, &kr);
		break;
	default:
	case 7:
		if (info->nloop > 0) {
			trl_restart_max_gap_cost_ratio(nd, tind, info, lambda, res,
				&kl, &kr);
		} else {
			trl_restart_fixed(nd, mnd, tind, lambda, res, info, &kl, &kr);
		}
		break;
	case 8:
		if (info->nloop > 0) {
			trl_restart_max_gap_cost_ratio_static(nd, tind, info, lambda,
				res, &kl, &kr);
		} else {
			trl_restart_fixed(nd, mnd, tind, lambda, res, info, &kl, &kr);
		}
		break;
	}
	//
	// make sure kr > kl+minsep
	minsep = max(3, max(nd / 6, nd - 6 * info->ned));
	if (kr <= kl + minsep || kl + nd - kr + minsep > mnd) {
		if (ncl < kl && kl < kr && kr < ncr) {
			kl--;
			kr++;
		} else if (info->lohi > 0) {
			kr = max(minsep, min(nd / 3, ncr)) - 1;
			kl = -1;
		} else if (info->lohi < 0) {
			kl = min(nd - minsep, max((nd + nd) / 3, ncl + 2)) - 1;
			kr = nd;
		} else {
			kl = (nd - minsep) / 2 - 2;
			kr = (nd - minsep + 1) / 2;
		}
	}
	// copy the (kr:nd) elements to (kl+1:kl+nd-kr+2)
	// kr is temporarily decreased by 1 to make indices easier to compute
	kr--;
	for (i = 1; i < nd - kr; i++) {
		lambda[kl + i] = lambda[kr + i];
		res[kl + i] = res[kr + i];
	}
	*kept = kl + max(1, nd - kr);
	if (info->lohi < -1) {
		dsort2(*kept, lambda, res);
	}
	return;
}



static void
TRL::trl_restart_fixed(int nd, int mnd, int tind, double *lambda,
				  double *res, trl_info * info, int *kl, int *kr)
{
	//
	// Purpose
	// =======
	// Save fixed number of extra Ritz pairs.
	//
	// Arguments
	// =========
	// nd         (input) INTEGER
	//             On entry, specifies the number of Lanczos basis considered.
	//
	// mnd        (input) INTEGER
	//             On entry, specifies the number of maximum Lanczos basis.
	//             (lanczos basis size used)-(Ritz values locked).
	//
	// tind       (input) INTEGER
	//             On entry, specifies the index of lambda, that is closest to target.
	//
	// lambda     (input) DOUBLE PRECISION ARRAY of LENGTH nd
	//             On entry, contains the computed Ritz values.
	//
	// res        (input) DOUBLE PRECISION ARRAY of LENGTH nd
	//             On entry, contains the residual of the computed Ritz values.
	//
	// info       (input) POINTER to TRLANINFO structure
	//             On entry, points to the current TRLANINFO structure.
	//
	// kl         (input/output) POINTER to INTEGER
	//             On exit, points to the starting index of lambda to keep from left.
	//
	// kr         (input/output) POINTER to INTEGER
	//             On exit, points to the starting index of lambda to keep from right.
	//
	// January 99 -- added modification to ensure a minimal gap ratio is
	// maintained
	//
	// ..
	// .. local scalars ..
	int extra, i, kl0, kr0, minsep, kli, kri;
	double gmin;
	//
	// ..
	// .. Executable statements ..
	// the number of extra Ritz pairs to be saved
	kli = *kl;
	kri = *kr;
	kl0 = kli;
	kr0 = kri;
	extra =
		nint((mnd - info->nec) * (0.4 + 0.1 * info->ned / (double) mnd));
	if (extra > info->ned + info->ned && extra > 5) {
		gmin = ((double) mnd) / ((double) info->ned);
		extra =
			nint((extra + (log(gmin) * info->ned) * gmin) / (1.0 + gmin));
	}
	minsep = max(3, max(nd / 5, nd - 4 * info->ned));
	gmin = trl_min_gap_ratio_(info, nd, tind, res);
	if (info->lohi > 0) {
		kri = min(tind - 1, kri) - extra;
		kli = -1;
	} else if (info->lohi < 0) {
		kli = max(tind + 1, kli) + extra;
		kri = nd;
	} else {
		extra++;
		kli = kli + extra / 2;
		kri = kri - extra / 2;
		i = 0;
		while (kli > kl0 && kri < kr0 && i == 0) {
			if (10.0 * res[kli] < res[kri]) {
				// lambda converged much more from left, so shift to right
				if (res[kli + 1] < res[kri + 1]) {
					kli++;
					kri++;
				} else {
					i = -1;
				}
			} else if (10.0 * res[kri] < res[kli]) {
				// lambda converged much more from right, so shift to left
				if (res[kri - 1] < res[kli - 1]) {
					kri--;
					kli--;
				} else {
					i = -1;
				}
			} else {
				i = -1;
			}
		}
	}
	// adjust kl and kr until the minimal gap ratio is satisfied
	while (kli + minsep < kri &&
		gap_ratio_(max(0, kli), min(kri, nd - 1), tind,
		lambda) < gmin) {
			if (info->lohi < 0) {
				kli++;
			} else if (info->lohi > 0) {
				kri--;
			} else if (res[kli] < res[kri]) {
				kli++;
			} else {
				kri++;
			}
	}
	// make sure nearly degenerate/duplicated Ritz pairs are included
	// lambda(kl)+r(kl) > lambda(j) and
	//                lambda(kl)-r(kl) > lambda(j)-r(j)
	if (info->lohi > 0) {
		i = kri - 1;
		// (lambda[i],lambda[i]+res[i]) is
		//   in (lambda[kri]-res[kri],lambda[kri]+res[kri]).
		while (i > kli + minsep && lambda[kri] - res[kri] < lambda[i] &&
			lambda[kri] + res[kri] < lambda[i] + res[i]) {
				i--;
		}
		kri = i + 1;
	} else {
		kl0 = kli;
		i = kli + 1;
		while (i < kri - minsep && lambda[kli] + res[kli] > lambda[i] &&
			lambda[kli] - res[kli] > lambda[i] - res[i]) {
				i++;
		}
		kli = i - 1;
	}
	*kl = kli;
	*kr = kri;
}



static double TRL::trl_min_gap_ratio_(trl_info * info, int nd, int tind, double *res)
{
	//
	// Purpose
	// =======
	// Try to determine the minimal gap ratio need to compute all wanted
	// eigenvalues
	//
	// ..
	// .. local scalars ..
	double gamma;
	//
	// ..
	// .. executable statements ..
	gamma = info->maxmv * (info->nec + 1.0) / info->ned - info->matvec;
	if (gamma < info->klan) {
		gamma =
			max((info->maxmv -
			info->matvec) / ((double) (info->ned - info->nec)), 2.0);
	}
	return min(log(res[tind] / (info->tol * info->anrm)) / gamma, 0.5);
}


inline double TRL::gap_ratio_(int i, int j, int tind, double *lambda)
{
	// Internal functio used in trl_restart_fixed_
	return (lambda[i] - lambda[tind]) / (lambda[j] - lambda[tind]);
}

////
static void
TRL::trl_restart_small_res(int nd, int mnd, int tind, double *lambda,
					  double *res, trl_info * info, int *kl, int *kr)
{
	//
	// Purpose
	// =======
	// Save those that are close to converge (as close as the target).
	//
	// Arguments
	// =========
	// nd         (input) INTEGER
	//             On entry, specifies the number of Lanczos basis considered.
	//
	// mnd        (input) INTEGER
	//             On entry, specifies the number of maximum Lanczos basis.
	//             (lanczos basis size used)-(Ritz values locked).
	//
	// tind       (input) INTEGER
	//             On entry, specifies the index of lambda, that is closest to target.
	//
	// lambda     (input) DOUBLE PRECISION ARRAY of LENGTH nd
	//             On entry, contains the computed Ritz values.
	//
	// res        (input) DOUBLE PRECISION ARRAY of LENGTH nd
	//             On entry, contains the residual of the computed Ritz values.
	//
	// info       (input) POINTER to TRLANINFO structure
	//             On entry, points to the current TRLANINFO structure.
	//
	// kl         (input/output) POINTER to INTEGER
	//             On exit, points to the starting index of lambda to keep from left.
	//
	// kr         (input/output) POINTER to INTEGER
	//             On exit, points to the starting index of lambda to keep from right.
	//
	// ..
	// .. local variables ..
	int extra, i, j, ii, kl0, kr0, kli, kri, minsep, done;
	double acpt, resmax, gmin;
	//
	// ..
	// .. executable statements ..
	// the number of extra Ritz pairs to be saved
	minsep = max(3, max(nd / 5, nd - 4 * info->ned));
	extra =
		nint((mnd - info->nec) * (0.4 + 0.1 * info->ned / ((double) mnd)));
	if (extra > info->ned + info->ned && extra > 5) {
		gmin = ((double) mnd) / ((double) info->ned);
		extra =
			nint((extra + (log(gmin) * info->ned) * gmin) / (1.0 + gmin));
	}
	kli = *kl;
	kri = *kr;
	kl0 = kli;
	kr0 = kri;
	resmax = maxval(nd, res);
	acpt = resmax / res[tind];
	//
	// determine the number of Ritz pairs that has to be saved
	if (info->lohi > 0) {
		if (acpt < 0.999 && acpt >= 0.0) {
			ii = tind - 1;
			acpt = max(sqrt(acpt) * res[tind], res[ii] + res[ii]);
			acpt = min(acpt, resmax);
			kri = ii - 1;
			while (res[kri] < acpt && kri > kli + 3) {
				kri--;
			}
		} else {
			kri = kr0 - extra;
		}
		kri = max(2, kri);
		kli = min(kli, kri - 2);
	} else if (info->lohi < 0) {
		if (acpt < 0.999 && acpt >= 0.0) {
			ii = tind + 1;
			acpt = max(sqrt(acpt) * res[tind], res[ii] + res[ii]);
			acpt = min(acpt, resmax);
			kli = ii + 1;
			while (res[kli] < acpt && kli < kri - 3) {
				kli++;
			}
		} else {
			kli = kl0 + extra;
		}
		kli = min(nd - 4, kli);
		kri = max(kri, kli + 2);
	} else {
		// save whoever has smaller residual norms
		i = kli + 1;
		j = kri - 1;
		for (ii = 1; ii <= extra; ii++) {
			if (res[i] < res[j]) {
				kli = i;
				i++;
			} else if (res[i] > res[j]) {
				kri = j;
				j--;
			} else if (i < nd - j) {
				kli = i;
				i++;
			} else {
				kri = j;
				j--;
			}
		}
	}
	// adjust kl and kr until the minimal gap ratio is satisfied
	kl0 = kli;
	kr0 = kri;
	gmin = trl_min_gap_ratio_(info, nd, tind, res);
	done = 0;
	while (kli + minsep < kri &&
		gap_ratio_(max(0, kli), min(nd - 1, kri), tind, lambda) < gmin
		&& (done == 0)) {
			if (info->lohi < 0) {
				kli++;
			} else if (info->lohi > 0) {
				kri--;
			} else if (res[kli] < res[kri]) {
				kli++;
			} else if (kri < nd - 1) {
				kri++;
			} else {
				done = 1;
			}
	}
	// make sure nearly degenerate Ritz pairs are included
	// lambda(kl)+r(kl) > lambda(j) and
	//                lambda(kl)-r(kl) > lambda(j)-r(j)
	if (info->lohi > 0) {
		i = kr0 - 1;
		while (i > kli + minsep && lambda[kri] - res[kri] < lambda[i] &&
			lambda[kri] + res[kri] < lambda[i] + res[i]) {
				i--;
		}
		kri = min(kri, i + 1);
	} else {
		i = kl0 + 1;
		while (i < kri - minsep && lambda[kli] + res[kli] > lambda[i] &&
			lambda[kli] - res[kli] > lambda[i] - res[i]) {
				i++;
		}
		kli = max(kli, i - 1);
	}
	*kl = kli;
	*kr = kri;
}



static double TRL::maxval(int n, double *a)
{
	//
	// Purpose
	// =======
	// Returns maximum value in a
	//
	// Arguments
	// =========
	// n   (input) INTEGER
	//      On entry, specifies the size of the array.
	//
	// a   (input) DOUBLE PREICISION ARRAY of LENGTH n
	//      On entry, contains the values to search for the maximum value.
	//
	// ..
	// .. local scalars ..
	int i;
	double val;
	//
	// ..
	// .. executable statements ..
	if (n <= 0) {
		return 0.0;
	}
	val = a[0];
	for (i = 1; i < n; i++) {
		val = max(val, a[i]);
	}
	return val;
}

////
static void
TRL::trl_restart_max_gap_ratio(int nd, int tind, int kept, double *lambda,
						  double *res, trl_info * info, int *kl, int *kr)
{
	//
	// Purpose
	// =======
	// Search throught all pairs of (kl, kr) for the one with the maximum
	// gap ratio for the next Ritz pair(target).
	// This is an optimized version of the original version.  It only search
	// through nd number once. (Only single loop!)
	//
	// Arguments
	// =========
	// nd         (input) INTEGER
	//             On entry, specifies the number of Lanczos basis considered.
	//
	// mnd        (input) INTEGER
	//             On entry, specifies the number of maximum Lanczos basis.
	//             (lanczos basis size used)-(Ritz values locked).
	//
	// tind       (input) INTEGER
	//             On entry, specifies the index of lambda, that is closest to target.
	//
	// kept       (input) INTEGER
	//             On entry, specifies the number of Ritz values kept.
	//
	// lambda     (input) DOUBLE PRECISION ARRAY of LENGTH nd
	//             On entry, contains the computed Ritz values.
	//
	// res        (input) DOUBLE PRECISION ARRAY of LENGTH nd
	//             On entry, contains the residual of the computed Ritz values.
	//
	// info       (input) POINTER to TRLANINFO structure
	//             On entry, points to the current TRLANINFO structure.
	//
	// kl         (input/output) POINTER to INTEGER
	//             On exit, points to the starting index of lambda to keep from left.
	//
	// kr         (input/output) POINTER to INTEGER
	//             On exit, points to the starting index of lambda to keep from right.
	//
	//  statement function for computing gap ratio
#define gratio(i,j) ( lambda[(j)] == info->trgt ?  DBL_MAX : (lambda[(i)]-info->trgt)/(lambda[(j)]-info->trgt) )
	//
	// ..
	// .. local variables ..
	int i, j, lohi, klm, krm, kli, kri, igap;
	double bnd, tmp;
	//
	// ..
	// .. executable statements ..
	// determine the search range
	kli = *kl;
	kri = *kr;
	trl_restart_search_range_(nd, lambda, res, info, kli, kri, &lohi, tind,
		&klm, &krm);
	kli = klm;
	kri = krm;
	igap = max(min(nd - info->ned, nint((krm - klm) * 0.4)), 2);
	if (igap > 2 && info->matvec > info->maxlan) {
		if (info->clk_op + info->tick_o >
			10.0 * (info->clk_orth + info->tick_h + info->clk_res +
			info->tick_r)) {
				igap = max(2, nd - kept - 1);
		} else {
			bnd = trl_min_gap_ratio_(info, nd, tind, res);
			if (info->crat < bnd)
				igap = max(2, nd - kept - 1);
		}
	}
	if (kli + igap > kri) {
		*kl = kli;
		*kr = kri;
		return;
	}
	// two cases depending on lohi
	if (lohi > 0) {
		// target is at the high end of spectrum
		bnd = gratio(kri, kli);
		for (i = klm; i <= krm - igap; i++) {
			j = i + igap;
			tmp = gratio(j, i);
			if (tmp > bnd) {
				kli = i;
				kri = j;
				bnd = tmp;
			}
		}
	} else {
		bnd = gratio(kli, kri);
		for (i = klm; i <= krm - igap; i++) {
			j = i + igap;
			tmp = gratio(i, j);
			if (tmp > bnd) {
				kli = i;
				kri = j;
				bnd = tmp;
			}
		}
	}
	*kl = kli;
	*kr = kri;
}

static void
TRL::trl_restart_max_progress(int nd, int tind, int kept, double *lambda,
						 double *res, trl_info * info, int *kl, int *kr)
{
	//
	// Purpose
	// =======
	// Search for a pair (kl, kr) such that the reduction in residual norm
	// of the target (info->trgt) will be the largest before next restart
	// The merit function is [gap-ratio * (m-k)]
	//
	// Arguments
	// =========
	// nd         (input) INTEGER
	//             On entry, specifies the number of Lanczos basis considered.
	//
	// tind       (input) INTEGER
	//             On entry, specifies the index of lambda, that is closest to target.
	//
	// kept       (input) INTEGER
	//             On entry, specifies the number of Ritz values kept.
	//
	// lambda     (input) DOUBLE PRECISION ARRAY of LENGTH nd
	//             On entry, contains the computed Ritz values.
	//
	// res        (input) DOUBLE PRECISION ARRAY of LENGTH nd
	//             On entry, contains the residual of the computed Ritz values.
	//
	// info       (input) POINTER to TRLANINFO structure
	//             On entry, points to the current TRLANINFO structure.
	//
	// kl         (input/output) POINTER to INTEGER
	//             On exit, points to the starting index of lambda to keep from left.
	//
	// kr         (input/output) POINTER to INTEGER
	//             On exit, points to the starting index of lambda to keep from right.
	//
	// merit measure the factor of residual norm reduction
#define merit(i,j) ( (lambda[(i)]-info->trgt) * abs(j-i) / (lambda[(j)]-info->trgt) )
	//
	// ..
	// .. local variables ..
	int i, j, lohi, kli, kri, klm, krm, igap;
	double tmp, ss;
	//
	// ..
	// .. executable statements ..
	//
	// determine the search range
	trl_restart_search_range_(nd, lambda, res, info, *kl, *kr, &lohi, tind,
		&klm, &krm);
	//
	// perform the brute-force search
	kli = klm;
	kri = krm;
	igap = max(min(nd - info->ned, nint((kri - kli) * 0.4)), 2);
	if (igap > 2 && igap + kept > nd && info->crat > 0.0) {
		ss = trl_min_gap_ratio_(info, nd, tind, res);
		if (ss > info->crat)
			igap = max(2, nd - kept - 1);
	}
	if (lohi > 0) {
		ss = merit(kri, kli);
		for (i = klm; i <= krm - igap; i++) {
			for (j = i + igap; j <= krm; j++) {
				tmp = merit(j, i);
				if (tmp > ss) {
					ss = tmp;
					kli = i;
					kri = j;
				}
			}
		}
	} else {
		ss = merit(kli, kri);
		for (i = klm; i <= krm - igap; i++) {
			for (j = i + igap; j <= krm; j++) {
				tmp = merit(i, j);
				if (tmp > ss) {
					ss = tmp;
					kli = i;
					kri = j;
				}
			}
		}
	}
	*kl = kli;
	*kr = kri;
}


static void
TRL::trl_restart_max_reduction(int nd, int tind, int kept, double *lambda,
						  double *res, trl_info * info, int *kl, int *kr)
{
	//
	// Purpose
	// =======
	// Search for a pair (kl, kr) such that the reduction in residual norm
	// of the target (info->trgt) will be the largest before next restart
	// the merit function is [sqrt(gap ratio) * (m-k)]
	//
	// Arguments
	// =========
	// nd         (input) INTEGER
	//             On entry, specifies the number of Lanczos basis considered.
	//
	// tind       (input) INTEGER
	//             On entry, specifies the index of lambda, that is closest to target.
	//
	// kept       (input) INTEGER
	//             On entry, specifies the number of Ritz values kept.
	//
	// lambda     (input) DOUBLE PRECISION ARRAY of LENGTH nd
	//             On entry, contains the computed Ritz values.
	//
	// res        (input) DOUBLE PRECISION ARRAY of LENGTH nd
	//             On entry, contains the residual of the computed Ritz values.
	//
	// info       (input) POINTER to TRLANINFO structure
	//             On entry, points to the current TRLANINFO structure.
	//
	// kl         (input/output) POINTER to INTEGER
	//             On exit, points to the starting index of lambda to keep from left.
	//
	// kr         (input/output) POINTER to INTEGER
	//             On exit, points to the starting index of lambda to keep from right.
	//
	// merit measure the factor of residual norm reduction
	//#define merit_maxred(i,j) ( sqrt( (lambda[(i)]-info->trgt)/(lambda[(j)]-info->trgt)) * abs(j-i) );
#define merit_maxred(i,j) ( sqrt( (lambda[(i)]-info->trgt)/(lambda[(j)]-lambda[(i)])) * abs(j-i) );
	//
	// ..
	// .. Parameters ..
	const double PI = 3.141592654;
	//
	// ..
	// .. local variables ..
	int i, j, lohi, kli, kri, klm, krm, t, igap;
	double tmp, tmp2, tmp3, ss, up, dw, def1, def2, z1, z2;
	//
	// ..
	// .. executable statements ..
	// determine the search range
	kli = *kl;
	kri = *kr;
	trl_restart_search_range_(nd, lambda, res, info, kli, kri, &lohi, tind,
		&klm, &krm);
	// perform the brute-force search
	kli = klm;
	kri = krm;
	//
	// ** Static approach to decide the minimum gap ratio **
	//igap = max( min(nd-info->ned, nint((krm-klm)*info->rfact)), 2);
	//
	// ** Dynamic approach to decide the minimum gap ratio **
	def1 = 1.0;
	def2 = 0.7;
	up = 1.0;
	dw = 0.7;
	info->avgm =
		(info->avgm * ((double) info->nloop - 1) +
		(double) (info->klan - info->k + 1)) / ((double) info->nloop);
	z1 = (info->ptres) / (info->tol * info->anrm);
	z2 = 1.0 / info->cfac;
	if (z2 < 1.0) {
		tmp = def2;
		t = abs(krm - klm) * def2;
	} else if (z1 < 1.0) {
		tmp = def1;
		t = abs(krm - klm) * def1;
	} else {
		tmp3 =
			log(z1 +
			sqrt(z1 - 1.0) * sqrt(z1 + 1.0)) / (2.0 * (info->avgm));
		tmp2 =
			log(z2 + sqrt(z2 - 1.0) * sqrt(z2 + 1.0)) / (info->klan -
			info->k + 1);

		tmp = tmp2 / tmp3;
		tmp = pow(tmp, 2.0);
		tmp = atan(tmp) * (2.0 / PI);
		tmp = dw + (up - dw) * tmp;

		t = abs(krm - klm) * tmp;
	}
	igap = t;
	if (igap > 2 && igap + kept > nd && info->crat > 0.0) {
		ss = trl_min_gap_ratio_(info, nd, tind, res);
		if (ss > info->crat)
			igap = max(2, nd - kept - 1);
	}
	if (lohi > 0) {
		ss = merit_maxred(kri, kli);
		for (i = klm; i <= krm - igap; i++) {
			for (j = i + igap; j <= krm; j++) {
				tmp = merit_maxred(j, i);
				if (tmp > ss) {
					ss = tmp;
					kli = i;
					kri = j;
				}
			}
		}
	} else {
		ss = merit_maxred(kli, kri);
		for (i = klm; i <= krm - igap; i++) {
			for (j = i + igap; j <= krm; j++) {
				tmp = merit_maxred(i, j);
				if (tmp > ss) {
					ss = tmp;
					kli = i;
					kri = j;
				}
			}
		}
	}
	*kl = kli;
	*kr = kri;
}


static void
TRL::trl_restart_scan(int nd, double *res, trl_info * info, int kept,
				 int *kl, int *kr)
{
	//
	// Purpose
	// =======
	// This subroutine determines the number Ritz vectors to be saved by
	// adding some quantity to the current thickness.  If the thickness is
	// larger than nd-2, it is reset to something smaller.
	// The thickness is varied progressively to scan all possible values.
	//
	// Arguments
	// =========
	// nd         (input) INTEGER
	//             On entry, specifies the number of Lanczos basis considered.
	//
	// res        (input) DOUBLE PRECISION ARRAY of LENGTH nd
	//             On entry, contains the residual of the computed Ritz values.
	//
	// info       (input) POINTER to TRLANINFO structure
	//             On entry, points to the current TRLANINFO structure.
	//
	// kept       (input) INTEGER
	//             On entry, specifies the number of Ritz values kept.
	//
	// kl         (input/output) POINTER to INTEGER
	//             On exit, points to the starting index of lambda to keep from left.
	//
	// kr         (input/output) POINTER to INTEGER
	//             On exit, points to the starting index of lambda to keep from right.
	//
	// ..
	// .. local variables ..
	int extra, i, kl0, kr0, kli, kri;
	//
	// ..
	// .. executable statements ..
	kli = *kl;
	kri = *kr;
	// three cases have to be dealt with separately
	if (info->lohi < 0) {
		// smallest eigenvalues
		kri = nd;
		kli = kept + min(max(info->nec, 1), (nd - kept) / 2) - 1;
		if (kli <= 0) {
			if (nd > 6) {
				kli = nd / 2 - 1;
			} else if (nd > 2) {
				kli = 1;
			}
		} else if (kli + 3 > nd) {
			kli =
				info->nec + min(info->ned,
				min(10, (nd - info->ned) / 2)) - 1;
		}
	} else if (info->lohi > 0) {
		kli = -1;
		kri = kept + min(max(info->nec, 1), (nd - kept) / 2) - 1;
		if (kri <= 0) {
			if (nd > 6) {
				kri = nd / 2 - 1;
			} else if (nd > 2) {
				kri = 1;
			}
		} else if (kri + 4 > nd) {
			kri =
				info->nec + min(info->ned, min(10, (nd - info->ned) / 2));
		}
		kri = nd - kri - 1;
	} else {
		kl0 = kli;
		kr0 = kri;
		extra = kept + min(info->nec, (nd - kept) / 2) + 1;
		if (extra <= 0) {
			if (nd > 6) {
				extra = nd / 2;
			} else if (nd > 2) {
				extra = 2;
			}
		} else if (extra + 3 > nd) {
			extra =
				info->nec + min(info->ned, min(10, (nd - info->ned) / 2));
		}
		kli = max(kli, (extra / 2) - 1);
		kri = min(kri, (nd - extra / 2) - 1);
		i = 0;
		while (kli > kl0 && kri < kr0 && i == 0) {
			if (10.0 * res[kli] < res[kri]) {
				if (res[kli + 1] < res[kri + 1]) {
					kli++;
					kri++;
				} else {
					i = -1;
				}
			} else if (10.0 * res[kri] < res[kli]) {
				if (res[kri - 1] < res[kli - 1]) {
					kri--;
					kli--;
				} else {
					i = -1;
				}
			} else {
				i = -1;
			}
		}
	}
	*kl = kli;
	*kr = kri;
	//
	// .. end of trl_restart_scan_ ..
	//
}


static void
TRL::trl_restart_max_gap_cost_ratio(int n, int tind, trl_info * info,
							   double *lambda, double *res, int *kl, int *kr)
{
	//
	// Purpose
	// =======
	// Determine the size of the maximum Lanczos basis size and the number of Ritz
	// pairs to keep. The decision is base on the following criteria
	//   argmax( m, k ) exp( - gamma * (m-k) )/ (m-k)(m-k-1)
	// where gamma = (lambda(kl)-lambda(1))/(lambda(kr)-lambda(1)
	// the numerator of the criteria approximates the expected improvement in the
	// residual norm of the Ritz pairs at the next restart, and the denominator
	// approximates the cost of the reorthogonalization till the next restart.
	// The next basis size m is chosen between info-rfact * k and maxlan, where
	// k is the number of Ritz vectors kept.
	//
	// Arguments
	// =========
	// n         (input/output) INTEGER
	//            On entry, specifies the number of Lanczos basis considered.
	//            On exit, specifies the number of maximum Lanczos basis size
	//            for the next iterations.
	//
	// tind      (input) INTEGER
	//            On entry, specifies the index of the next target Ritz value.
	//
	//
	// lambda    (input) DOUBLE PRECISION ARRAY
	//            On entry, contains the computed Ritz values.
	//
	// res       (input) DOUBLE PRECISION ARRAY
	//            On entry, contains the residueal norm.
	//
	// kl        (input/output) Pointer to INTEGER
	//            On entry, specifies the largest Ritz value from left, that has
	//            been converged.
	//            On exit, specifies the largest Ritz value from left, that are kept.
	//
	// kr        (input/output) Pointer to INTEGER
	//            On entry, specifies the smallest Ritz value from right, that has
	//            been converged.
	//            On exit, specifies the smallest Ritz value from left, that are kept.
	//
	// ..
	// .. Parameters ..
	const double PI = 3.141592654;
	//
	// ..
	// .. Local scalars ..
	int i, j, k, nd, l1, l2, mn, mn1, mn2, t, k1, k2,
		min_k1, min_k2, min_m, min_l;
	double gamma0, gamma, val, min_val, min_gamma, min_gamma0, min_ratio,
		tmp, tmp2, tmp3, up, dw, def1, def2, z1, z2;
	//
	// ..
	// .. Executable statements ..
	nd = info->ned;
	mn = info->maxlan;
	trl_restart_search_range_(n, lambda, res, info, (*kl), (*kr),
		&(info->lohi), tind, &k1, &k2);
	// ** Static approach to decide the minimum gap ratio **
	// minimum gap beteween kl and kr for ok-conditioned systems,
	// i.e, diag(i) and diag(i^2)
	//t = nint(4.0*abs(k1-k2)/5.0);
	//t = nint((info->mgap)*abs(k1-k2));
	//t = nint(2.0*abs(k1-k2)/5.0);
	/*
	if( t > n-nd ) t=n-nd;
	if( t < 2 ) t = 2;
	if( info->lohi == -1 ) {
	mcnv0 = info->klan - (*kr);
	} else if( info->lohi == 1 ) {
	mcnv0 = (*kl)+1;
	}
	if( t > 2 && t+kept > nd && info->crat > 0.0 ) {
	min_val = trl_min_gap_ratio_(info, nd, tind, res);
	if( min_val > info->crat ) t = max(2, nd-kept-1);
	}
	*/
	//
	// ** Dynamic approach to decide the minimum gap ratio **
	def1 = 1.0;
	def2 = 0.7;
	up = 1.0;
	dw = 0.7;
	info->avgm =
		(info->avgm * ((double) info->nloop - 1) +
		(double) (info->klan - info->k + 1)) / ((double) info->nloop);
	z1 = (info->ptres) / (info->tol * info->anrm);
	z2 = 1.0 / info->cfac;
	if (z2 < 1.0) {
		tmp = def2;
		t = abs(k2 - k1) * def2;
	} else if (z1 < 1.0) {
		tmp = def1;
		t = abs(k2 - k1) * def1;
	} else {
		tmp3 =
			log(z1 +
			sqrt(z1 - 1.0) * sqrt(z1 + 1.0)) / (2.0 * (info->avgm));
		tmp2 =
			log(z2 + sqrt(z2 - 1.0) * sqrt(z2 + 1.0)) / (info->klan -
			info->k + 1);
		tmp = tmp2 / tmp3;
		tmp = pow(tmp, 2.0);
		tmp = atan(tmp) * (2.0 / PI);
		tmp = dw + (up - dw) * tmp;

		t = abs(k2 - k1) * tmp;
	}
	//
	// ** Static approach to decided the minimum gap ratio **
	//t = max (min (n - info->ned, nint ((k2 - k1) * info->rfact)), 2);
	//
	mn1 = nint(2.0 * (info->klan) / 5.0);
	min_val = 0.0;
	min_k1 = k1;
	min_k2 = k2;
	min_m = info->klan;
	//printf( "k1=%d k2=%d\n",k1,k2 );
	for (i = k1; i <= k2 - t; i++) {
		// considering how many Ritz vectors to keep from left.
		for (j = i + t; j <= k2; j++) {
			// considering how many Ritz vectors to keep from right.
			//
			/* note the two option, i.e., square-rooted, are hard-coded here !!! */
			gamma0 =
				sqrt((lambda[i + 1] - lambda[tind]) / (lambda[j - 1] -
				lambda[i + 1]));
			/* original object function */
			//gamma0 = sqrt((lambda[i+1]-lambda[tind])/(lambda[j-1]-lambda[tind]));
			/* Without square */
			//gamma0 = (lambda[i+1]-lambda[tind])/(lambda[j-1]-lambda[tind]);
			l1 = (i + 1) + (n - j);
			l2 = l1 + info->locked;

			//mn2 = mn1;
			//if( mn1 < (info->rfact*l2) ) mn2 = info->rfact*l2;
			//if( mn1 < (1.5*l2) ) mn2 = 1.5*l2;
			//if( mn1 < (l2+1) ) mn2 = l2+1;
			mn2 = l2 + 1;
			if (mn2 > mn)
				mn2 = mn;
			if (mn2 < nd)
				mn2 = nd;

			for (k = mn2; k <= mn; k++) {
				// considering k for the next Lanczos sizes
				gamma = ((k - l2) * gamma0);
				//* Approximate cost function *//
				//val = ( k * k - l2 * l2 + 2 * k * l1 ) / ( gamma );
				//* More exact cost function *//
				val = ((k - l2) * (k + l2 - 1) + k * l1) / (gamma);
				//* Including cosh function *//
				//val = ( (k-l2) * (k+l2-1) + k * l1 ) / ( cosh( 2.0 * gamma ) );
				//printf( "%d %e %e\n",( (k-l2) * (k+l2-1) + k * l1 ), gamma, cosh(2.0*gamma) );

				if (min_val == 0.0 || val < min_val) {
					min_k1 = i;
					min_k2 = j;
					min_m = k;
					min_val = val;
					min_gamma = gamma;
					min_gamma0 = gamma0;
					min_ratio =
						(lambda[i + 1] - lambda[tind]) / (lambda[j - 1] -
						lambda[tind]);
					min_l = l2;
				}
			}
		}
	}
	*kl = min_k1;
	*kr = min_k2;
	info->klan = min_m;
	info->k1 = min_k1;
	info->k2 = min_k2;
	info->k = min_l;
	//printf( "%e %d %d %d\n",tmp,min_k1,min_k2,min_m );
	info->predicted_crate = exp(-min_gamma);
	info->mgamma = min_gamma;
	info->gamma0 = min_gamma0;
	info->old_target = lambda[tind];
	info->target_id = tind;
	info->old_locked = info->locked;
}


static void
TRL::trl_restart_max_gap_cost_ratio_static(int n, int tind, trl_info * info,
									  double *lambda, double *res,
									  int *kl, int *kr)
{
	//
	// Purpose
	// =======
	// Determine the size of the maximum Lanczos basis size and the number of Ritz
	// pairs to keep. The decision is base on the following criteria
	//   argmax( m, k ) exp( - gamma * (m-k) )/ (m-k)(m-k-1)
	// where gamma = (lambda(kl)-lambda(1))/(lambda(kr)-lambda(1)
	// the numerator of the criteria approximates the expected improvement in the
	// residual norm of the Ritz pairs at the next restart, and the denominator
	// approximates the cost of the reorthogonalization till the next restart.
	// The next basis size m is fixed at info-rfact * k and maxlan, where k is the number
	// of Ritz vectors kept.
	//
	// Arguments
	// =========
	// n         (input/output) INTEGER
	//            On entry, specifies the number of Lanczos basis considered.
	//            On exit, specifies the number of maximum Lanczos basis size
	//            for the next iterations.
	//
	// tind      (input) INTEGER
	//            On entry, specifies the index of the next target Ritz value.
	//
	//
	// lambda    (input) DOUBLE PRECISION ARRAY
	//            On entry, contains the computed Ritz values.
	//
	// res       (input) DOUBLE PRECISION ARRAY
	//            On entry, contains the residueal norm.
	//
	// kl        (input/output) Pointer to INTEGER
	//            On entry, specifies the largest Ritz value from left, that has
	//            been converged.
	//            On exit, specifies the largest Ritz value from left, that are kept.
	//
	// kr        (input/output) Pointer to INTEGER
	//            On entry, specifies the smallest Ritz value from right, that has
	//            been converged.
	//            On exit, specifies the smallest Ritz value from left, that are kept.
	//
	// ..
	// .. Parameters ..
	const double PI = 3.141592654;
	//
	// ..
	// .. Local scalars ..
	int i, j, nd, l1, l2, mn, mn1, mn2, t, k1, k2,
		min_k1, min_k2, min_m, min_l;
	double gamma0, gamma, val, min_val, min_gamma, min_gamma0, min_ratio,
		tmp2, tmp3, tmp;
	//
	// ..
	// .. Executable statements ..
	nd = info->ned;
	mn = info->maxlan;
	trl_restart_search_range_(n, lambda, res, info, (*kl), (*kr),
		&(info->lohi), tind, &k1, &k2);
	//
	// ** Static approach to decide the minimum gap ratio **
	//t = nint(4.0*abs(k1-k2)/5.0);
	/*
	t = nint((info->mgap)*abs(k1-k2));
	if( t > n-nd ) t=n-nd;
	if( t < 2 ) t = 2;

	if( t > 2 && t+kept > nd && info->crat > 0.0 ) {
	min_val = trl_min_gap_ratio_(info, nd, tind, res);
	if( min_val > info->crat ) t = max(2, nd-kept-1);
	}
	*/
	//
	// ** Dynamic approach to decide the minimum gap ratio **
	if (info->crat > 0.0) {
		tmp2 = -log(info->crat);
		tmp3 = -(tmp2 * (info->maxlan)) / log(info->tol * info->anrm);
		if (tmp3 >= 0.5 || info->klan < info->maxlan) {
			tmp = 0.8;
			t = abs(k2 - k1) * 0.8;
		} else {
			tmp = trl_min_gap_ratio_(info, nd, tind, res);
			tmp = tmp2 / tmp;
			tmp = pow(tmp, 0.25);
			tmp = atan(tmp) * (4.0 / PI);
			t = abs(k2 - k1) * tmp;
		}
	} else {
		tmp = 0.8;
		t = abs(k2 - k1) * 0.8;
	}
	mn1 = nint(2.0 * (info->klan) / 5.0);
	min_val = 0.0;
	min_k1 = k1;
	min_k2 = k2;
	min_m = info->klan;
	for (i = k1; i <= k2 - t; i++) {
		// considering how many Ritz vectors to keep from left.
		for (j = i + t; j <= k2; j++) {
			// considering how many Ritz vectors to keep from right.
			gamma0 =
				sqrt((lambda[i + 1] - lambda[tind]) / (lambda[j - 1] -
				lambda[i + 1]));
			//gamma0 = sqrt((lambda[i+1]-lambda[tind])/(lambda[j-1]-lambda[tind]));
			//gamma0 = (lambda[i+1]-lambda[tind])/(lambda[j-1]-lambda[tind]);
			l1 = (i + 1) + (n - j);
			l2 = l1 + info->locked;
			mn2 = mn1;
			if (mn1 < ((info->rfact) * l2))
				mn2 = (info->rfact) * l2;
			if (mn2 > mn)
				mn2 = mn;
			if (mn2 < nd)
				mn2 = nd;
			gamma = (mn2 - l2 + 1) * gamma0;
			val = (mn2 * mn2 - l2 * l2 + 2 * mn2 * l1) / (gamma);
			if (min_val == 0.0 || val < min_val) {
				min_k1 = i;
				min_k2 = j;
				min_m = mn2;
				min_val = val;
				min_gamma = gamma;
				min_gamma0 = gamma0;
				min_ratio =
					(lambda[i + 1] - lambda[tind]) / (lambda[j - 1] -
					lambda[tind]);
				min_l = l2;
			}
		}
	}
	*kl = min_k1;
	*kr = min_k2;
	info->klan = min_m;
}


static void TRL::trl_restart_search_range_(int nd, double *lambda, double *res,
									  trl_info * info, int ncl, int ncr,
									  int *lohi, int tind, int *klm, int *krm)
{
	//
	// Purpose
	// =======
	// Determine the search range --
	// used by the schemes that performs brute-force search.
	//
	// ..
	// .. local variables ..
	int j, klmi, krmi;
	double bnd;
	//
	// ..
	// .. executables statements ..
	klmi = max(ncl, 0);
	krmi = min(ncr, nd - 1);
	bnd = info->tol * info->anrm;
	*lohi = info->lohi;
	// make sure there is some distance between the boundary and the
	// target Ritz value
	if (info->lohi > 0) {
		// save high end
		krmi =
			min(max
			(info->maxlan - info->ned,
			(info->maxlan + info->nec) / 2) - 1, min(krmi, tind - 1));
		while (krmi + krmi > ncl + ncr && res[krmi] < bnd) {
			krmi--;
		}
	} else if (info->lohi < 0) {
		// save lower end
		klmi = max(min(info->ned, (info->maxlan + info->nec) / 2) - 1,
			max(tind + 1, klmi));
		while (klmi + klmi < ncl + ncr && res[klmi] < bnd) {
			klmi++;
		}
	} else {
		// save both ends
		if (tind - klmi < krmi - tind) {
			// tind is closer to klmi
			*lohi = -1;
			klmi = tind + 1;
		} else {
			// tind is closer to krmi
			*lohi = 1;
			krmi = tind - 1;
		}
		j = info->locked + klmi + nd - krmi + 1;
		if (j > 0) {
			j = j / 2;
			// should be bounded by 0 and nd-1
			klmi = max(0, klmi - j);
			krmi = min(krmi + j, nd - 1);
		}
	}
	*klm = klmi;
	*krm = krmi;
}


///////////////

//

/////trl map////
void TRL::trl_dgemv(char *trans, int m, int n, double alpha, double *a, int lda,
			   double *x, int incx, double beta, double *y, int incy)
{

#ifdef BLAS_EXT
	extern int dgemv_();
#endif

	dgemv_(trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
}
double TRL::trl_ddot(integer_ n, double *dx, integer_ incx,
				double *dy, integer_ incy)
{
#ifdef BLAS_EXT
	extern double ddot_();
#endif

	return TRL::ddot_(&n, dx, &incx, dy, &incy);
}


void TRL::trl_daxpy(int n, double da, double *dx, int incx, double *dy,
			   int incy)
{

#ifdef BLAS_EXT
	extern int daxpy_();
#endif

	daxpy_(&n, &da, dx, &incx, dy, &incy);
}



void TRL::trl_dscal(int n, double da, double *dx, int incx)
{

#ifdef BLAS_EXT
	extern int dscal_();
#endif

	dscal_(&n, &da, dx, &incx);

}


void TRL::trl_dgemm(char *transa, char *transb, int m, int n, int k,
			   double alpha, double *a, int lda, double *b, int ldb,
			   double beta, double *c, int ldc)
{

#ifdef BLAS_EXT
	extern int dgemm_();
#endif

	dgemm_(transa, transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c,
		&ldc);
}

void TRL::trl_dcopy(int n, double *dx, int incx, double *dy, int incy)
{

#ifdef BLAS_EXT
	extern int dcopy_();
#endif

	dcopy_(&n, dx, &incx, dy, &incy);
}

///
////////////////


////************/////////
///////////core function////////////////

/*
// Purpose
// =======
// The actual work routine of restarted Lanczos program for real
// symmetric eigenvalue problems
//
// user may directly invoke this sunroutine but she/he is responsible
// for input correct arguments as needed
//
// 1) info needs to be initialized
// 2) if info->nec>0, the first nec elements of eval and first nec
//    columns of evec must contain valid eigenpairs
// 3) workspace of currect size
//    eval(mev)
//    evec(lde, mev) (lde >= nrow, mev >= ned)
//    base(ldb, info->maxlan-mev+1) (ldb>=nrow, not used if mev>maxlan)
//    wrk(lwrk) minimum amount of memory required by TRLANCZOS is
//    maxlan*(maxlan+10)
// 4) if log files are to be written, the user need to open files on IO
//    unit log_io so that the log gile may be written correctly.
// 5) the user must set the timing variable info->clk_tot and
//    info->clk_max using system_clock function call in order for this
//    subroutine to track timing results correctly
//
// Algorithm
// =========
//  0. initialize input vector
//  1. do while (more eigenvalues to compute .and. more MATVEC allowed)
//  2.    first step
//     o   alpha(k+1) = dot_product(q_{k+1}, Aq_{k+1})
//     o   rr = A*q_{k+1}-alpha(k+1)*q_{k+1}-\sum_{i=1}^{k} beta(i)q_i
//     o   (re-orthogonalization)
//  3.    do j = k+2, m
//     o     rr = Aq_j
//     o     alpha(j) = dot_product(q_j, rr)
//     o     rr = rr - alpha(j)*q_j - beta(j-1)*q_{j-1}
//     o     (re-orthogonalization)
//        end do j = k+2, m
//  4.    restarting
//     o   call dstqrb to decompose the tridiagonal matrix
//     o   perform convergence test
//     o   determine what and how many Ritz pairs to save
//     o   compute the Ritz pairs to be saved
//     end do while
//
// The re-orthogonalization procedure is implemented in trl_orth.  it
// produces a normalized vector rr that is guaranteed to be orthogonal
// to the current basis.  An error will be reported if it can not
// achieve its goal.
//
// Arguments
// =========
// ops     (input) Functin pointer.
//          On entry, points to the function that performs the matrix-vector 
//          operation. The operator that defines the eigenvalue problem is 
//          expected to have the following interface
//          void op(nrow, ncol, xin, ldx, yout, ldy)
//             nrow  (input) Integer
//                    On entry, specifies the number of rows in xin and xout.
//             ncol  (input) Integer
//                    On entry, specifies the number of columns in xin and xout.
//             xin   (input) double precision array of dimension (ldx,ncol)
//                    On entry, specifies the vector/vectors for which the 
//                    matrix-vector is performed.
//             ldx   (input) Integer
//                    On entry, specifies the leading dimension of xin
//             yout  (output) Double precision array of diimension (ldy,ncol)
//                    On exit, specifies the resulting vector/vectors.
//             ldy   (input) Integer
//                    On entry, specifies the leading dimension of yout.
//
// info    (input) Pointer to structure trl_info_
//          On entry, points to the current TRL_INFO.
//
// nrow    (input) Integer
//          On entry, specifies the number of rows in eigenvectors.
//
// mev     (input) Integer
//          On entry, specifies the number of columns allocated to store 
//          eigenvectors.
//
// eval    (output) Double array of dimension (mev)
//          On successful exit, contains the eigenvalues.
//
// evec    (output) Double array of dimension (lde,mev)
//          On successful exit, contains the eigenvectors.
//
// lde     (input) Integer
//          On entry, specifies the leading dimension of evec.
//
// base    (workspace) Double precision array of dimension (ldb,nbas)
//          Used to hold the lanczos vectors not fit in evec, i.e., 
//          nbas=info->maxlan-mev+1.
//
// ldb     (input) Integer
//          On entry, specifies the leading dimension of base.
//
// nbas    (input) Integer
//          On entry, specifies the number of columns in base.
//
// wrk     (workspace) Double precision array of dimension (lwrk)
//          Workspace for lanczos iterations.
//
// lwrk    (input) Integer
//          On entry, specifies the size of workspace provided.
//
*/
void
TRL::trlanczos(trl_matvec op, trl_info * info, int nrow, int mev, double *eval,
		  double *evec, int lde, double *base, int ldb, int nbas,
		  double *wrk, int lwrk)
{
	/*
	// ..
	// .. local parameters ..
	*/
	char notrans = 'N';
	int c__1 = 1;
	int i__1 = 1;
	double one = 1.0;
	/*
	// ..
	// .. local variables ..
	*/
	char title[STRING_LEN];
	int i, i1, i2, j1, j2, jnd, jml, j1n, j2n, kept, prek, ldqa, ldrr,
		count;
	int next_test, lwrk2, chkpnt, locked, degen;
	clock_t clk1;
	int *iwrk;
	double d__1;
	double *alpha, *beta, *rr, *rot, *alfrot, *betrot, *lambda, *res, *yy,
		*qa, *qb, *wrk2;
	/*
	// ..
	// .. executable statements ..
	//
	// initialize title and clock.
	*/
	strcpy(title, "");
	clk1 = 0;
	/*
	// alpha, beta, alfrot and betrot have fixed locations in wrk, i.e.,
	//   alpha: wrk(1:maxlan), beta: wrk(maxlan+1:2*maxlan),
	//   alfrot: wrk(2*maxlan+1:3*maxlan), betrot: wrk(3*maxlan+1:4*maxlan)
	*/
	alpha = &wrk[0];
	i1 = info->maxlan + 1;
	i2 = info->maxlan + info->maxlan;
	beta = &wrk[i1 - 1];
	i1 = i2 + 1;
	i2 = i2 + info->maxlan;
	alfrot = &wrk[i1 - 1];
	i1 = i2 + 1;
	i2 += info->maxlan;
	betrot = &wrk[i1 - 1];
	/*
	// allocate an integer_ workspace. iwrk holds...
	*/
	iwrk = (int *) malloc((4 * info->maxlan) * sizeof(int));
	memset(iwrk, 0, (4 * info->maxlan) * sizeof(int));
	/*
	// chkpnt specifies how often the check point should be written.
	*/
	if (info->cpflag_ <= 0) {
		/* check point is not written */
		chkpnt = info->maxmv + info->maxlan;
	} else {
		/* check point is written at every chpnt matrix-vector operations. */
		chkpnt = info->maxmv / info->cpflag_;
	}
	/*
	// locked specifies the number of eigenvalues converged.
	*/
	locked = info->nec;
	/*
	// assign values to alpha, beta
	// uses the assumption that the content of eval(1:nec) are eigenvalues
	// and their residual norms are zero
	*/
	if (locked > 0) {
		memcpy(alpha, eval, locked * sizeof(double));
		memset(beta, 0, locked * sizeof(double));
	}
	/*
	get valid initial guess for the Lanczos iterations
	wrk2 points to the end of available wrk 
	(first 4*maxlan hold alpha, beta, alfrot, and betrot)
	to the end of wrk.
	*** j1 and j2 are the number of colums, and not indices ***
	*/
	wrk2 = &wrk[i2];
	lwrk2 = lwrk - i2;
	/*
#ifdef DEBUG_WEI_TRL
	cout<<"Initial guess @ trlanzco";
#endif
	*/
	trl_initial_guess(nrow, evec, lde, mev, base, ldb, nbas, alpha,
		beta, &j1, &j2, info, wrk2, lwrk2);
	/*
	On return from trl_initial_guess, j1 is the last column index of
	evecsed, and j2 is the last column index of base used, i.e., jnd
	specifies the sizef the current lanczos basis.
	*/
	jnd = j1 + j2;
	kept = jnd;
	if (info->stat != 0) {
		if (info->stat < 0 && (info->verbose > 0 || info->my_pe == 0)) {
			qa = NULL;
			qb = NULL;
			rr = NULL;
			kept = 0;
			j1 = 0;
			j2 = 0;
			jnd = 0;
			TRL::log_error_state(info, kept, j1, j2, jnd, nrow, mev, eval,
				alpha, alfrot, beta, betrot, evec, base, qa,
				qb, rr, title, iwrk);
		}
		free(iwrk);
		return;
	}
	/*
	we will perform the first convergence test after next_test
	matrix-vector multiplications
	*/
	i1 = info->ned - jnd;
	next_test = i1 + min(i1, min(6, info->ned / 2));
	/*
	// *********************************************** //
	//            -- the TRLan outer loop --           //
	// restart if                                      //
	//    1. there is more eigenvalues to compute, and //
	//    2. more matrix-vector operations are allowed //
	// *********************************************** //
	*/
	//tick1 = clock();
	count = 0;
	degen = -1;
	while (info->matvec < info->maxmv
		&& (degen != -1 || info->nec < info->ned)) {
			/*
			// jnd is the size of the current lanczos basis, so increment it.
			*/
			//tick2 = clock();
			//tick1 = tick2;

			jnd++;
			/*
			// specify the workspace to hold the rotation matrix, that transform the 
			// tridiagonal matrix to a diagonal matrix, i.e., the size of the matrix is 
			// (jnd-locked) or the size of the current lanczos basis minus the number 
			// of the eigenvalues converged so far. The workspace of the rotation matrix 
			// is located at the end of wrk.
			*/
			i2 = lwrk - (jnd - locked) * (jnd - locked);
			rot = &wrk[i2];
			/*
			spcify the workspace required for the orthogonalization procedure.
			the workspace is after the space storing alpha, beta, alfrot, and
			betrot, but before the space storing the rotation matrx.
			*/
			i1 = 4 * info->maxlan + 1;
			wrk2 = &wrk[i1 - 1];
			lwrk2 = i2 - i1 + 1;
			/*
			check if there are enough workspace for the orthogonalization
			procedure.
			*/
			i1 = max(5 * info->maxlan - 3 * locked, 4 * info->maxlan);
			if (lwrk2 < i1) {
				info->stat = -11;
				return;
			}
			/*
			the first iteration of TRLan
			qa points to the last lanczos base vector.
			*/
			if (j1 < mev) {
				/* there is still enough space in evec, so use it. */
				j1++;
				qa = &evec[(j1 - 1) * lde];
				ldqa = lde;
			} else {
				/* no more space in evec, so use base. */
				j2++;
				qa = &base[(j2 - 1) * ldb];
				ldqa = nrow;
			}
			/*
			j1n and j2n specify the location of the next lanczos basis, and
			rr points the next lanczos base vector.
			*/
			if (j1 < mev) {
				/* there is still enough space in evec, so use it. */
				j1n = j1 + 1;
				j2n = 0;
				rr = &evec[(j1n - 1) * lde];
				ldrr = lde;
			} else {
				/* no more space in evec, so use base. */
				j1n = mev;
				j2n = j2 + 1;
				rr = &base[(j2n - 1) * ldb];
				ldrr = ldb;
			}
			/*
			perform matrix-vector multiplication, i.e., rr = A*qa
			record the total time and the time in MATVEC
			*/
			clk1 = clock();
#ifdef __CLK_RESTART
			if (clk1 <= info->clk_tot) {
				info->tick_t +=
					(info->clk_max -
					info->clk_tot) / (double) (info->clk_rate);
				info->tick_t += clk1 / (double) (info->clk_rate);
				info->clk_tot = clk1;
			}
#else
			if (clk1 < info->clk_tot) {
				info->tick_t +=
					(info->clk_max -
					info->clk_tot) / (double) (info->clk_rate);
				info->tick_t +=
					(info->clk_max + clk1) / (double) (info->clk_rate);
				info->clk_tot = clk1;
			} else if (info->clk_tot < 0 && clk1 >= 0) {
				info->tick_t -= info->clk_tot / (double) (info->clk_rate);
				info->tick_t += clk1 / (double) (info->clk_rate);
				info->clk_tot = clk1;
			}
#endif
#ifdef TRL_FORTRAN_COMPATIBLE
			op(&nrow, &i__1, qa, &ldqa, rr, &ldrr);
#else
			op(nrow, i__1, qa, ldqa, rr, ldrr, info->mvparam);
#endif
			add_clock_ticks(info, &(info->clk_op), &(info->tick_o), clk1);
			(info->matvec)++;
			/*
			// computed the next alpha = qa' * A * qa
			*/
			alpha[jnd - 1] = trl_ddot(nrow, qa, c__1, rr, c__1);
			trl_g_sum(info->mpicom, 1, &alpha[jnd - 1], wrk2);
			/*
			// Perform the Lanczos orthogonalization.
			// rr = rr - sum_{i=1,...j1} 
			//             beta(i)*evec(:,i) - sum_{1,...,j2} beta(j1+i)*base(:,i)
			// Just for a convenience beta[jnd-1]=alpha[jnd-1] just computed.
			*/
			beta[jnd - 1] = alpha[jnd - 1];
			info->flop = info->flop + nrow + nrow;
			/*
			// orthogonalize with lanczos vectors stored in evec, first.
			*/
			if (j1 > 2) {
				/*
				compute rr = rr - [evec(:,1),...,evec(:,i1)]*[beta(1),...,beta(i1)]'
				*/
				d__1 = -one;
				trl_dgemv(&notrans, nrow, j1, d__1, evec, lde, beta, c__1, one,
					rr, c__1);
				info->flop = info->flop + 2 * j1 * nrow;
			} else if (j1 == 1) {
				/*
				// there is no beta, so just compute
				//   rr = rr - alpha(1)*qa
				*/
				d__1 = -alpha[0];
				trl_daxpy(nrow, d__1, qa, c__1, rr, c__1);
				info->flop = info->flop + nrow + nrow;
			} else if (j1 == 2) {
				/*
				there is only one beta, so just do
				rr = rr - beta(1)*evec(1:nrow,1) - beta(2)*evec(1:nrow,2)
				*/
				d__1 = -beta[0];
				trl_daxpy(nrow, d__1, evec, c__1, rr, c__1);
				d__1 = -beta[1];
				trl_daxpy(nrow, d__1, &evec[lde], c__1, rr, c__1);
				info->flop = info->flop + 4 * nrow;
			}
			/*
			// orthogonalize with lanczos vectors stored in base, now.
			*/
			if (j2 > 2) {
				/*
				rr = rr - [evec(:,1),...,evec(:,i1)]*[beta(j1+1),...,beta(j1+j2)]'
				*/
				d__1 = -one;
				trl_dgemv(&notrans, nrow, j2, d__1, base, ldb, &beta[j1], c__1,
					one, rr, c__1);
				info->flop = info->flop + 2 * j2 * nrow;
			} else if (j2 == 1) {
				/*
				there is no beta, so just compute
				rr = rr - beta(jnd)*qa
				*/
				d__1 = -beta[jnd - 1];
				trl_daxpy(nrow, d__1, qa, c__1, rr, c__1);
				info->flop = info->flop + nrow + nrow;
			} else if (j2 == 2) {
				/*
				there is only one beta, so just do
				rr = rr - beta(j1+1)*base(1:nrow,1) - beta(jnd)*base(1:nrow,2)
				*/
				d__1 = -beta[j1];
				trl_daxpy(nrow, d__1, base, c__1, rr, c__1);
				d__1 = -beta[jnd - 1];
				trl_daxpy(nrow, d__1, &base[ldb], c__1, rr, c__1);
				info->flop = info->flop + 4 * nrow;
			}
			/*
			// perform re-orthogonalization (full-orthogonalization)
			*/
			info->flop_h = info->flop_h - info->flop;
			clk1 = clock();
			trl_orth(nrow, evec, lde, j1, base, ldb, j2, rr, kept, alpha,
				beta, wrk2, lwrk2, info);
			if (info->verbose > 8) {
				/* check orthogonality after the initilization step */
				trl_check_orth(info, nrow, evec, lde, j1n, base, ldb, j2n,
					wrk2, lwrk2);
			}
			add_clock_ticks(info, &(info->clk_orth), &(info->tick_h), clk1);
			info->flop_h = info->flop_h + info->flop;
			if (info->stat != 0)
				goto end;
			if (info->verbose > 5) {
				print_alpha_beta(info, title, jnd, alpha, beta);
			}
			/*
			// transform the matrix formed by alpha and beta into a
			// tridiagonal matrix, rot stores the transformation matrix
			*/
			/* the already-converged part is just diagonal. */
			memcpy(alfrot, alpha, locked * sizeof(double));
			memset(betrot, 0, locked * sizeof(double));
			/*
			// now, diagonalize the rest of matrix.
			*/
			i1 = jnd - locked;
			trl_tridiag(i1, &alpha[locked], &beta[locked], rot,
				&alfrot[locked], &betrot[locked], wrk2, lwrk2,
				&(info->stat));
			info->flop = info->flop + 8 * i1 * i1 * i1 / 3;	/* Golub:1996:MC, P415 */
			if (info->stat != 0)
				goto end;
			betrot[jnd - 1] = beta[jnd - 1];
			/*
			// **************************************************** //
			// regular iterations of Lanczos algorithm (inner loop) //
			// loop if                                              //
			//   1. there is space to store lanczos basis           //
			//   2. there is more eigenvalues to compute, and       //
			// **************************************************** //
			*/
			while (jnd < info->klan && (degen != -1 || info->nec < info->ned)) {
				/*
				// compute the kth lanczos vector.
				//  qb is (k-2)nd lanczos vector, and qa is (k-1)st lanczos vector.
				*/
				//printf( "   ** inner loop (%d) **\n",jnd  );
				qb = qa;
				qa = rr;
				/*
				// increment j1, j2, and jnd.
				*/
				j1 = j1n;
				j2 = j2n;
				jnd++;
				/*
				// find the next available space for the kth lanczos vector.
				*/
				if (j1n < mev) {
					/* there is still a space in evec. */
					j1n++;
					rr = &evec[(j1n - 1) * lde];
				} else {
					/* no more space in evec, so use base. */
					j2n++;
					if (j2n <= nbas) {
						rr = &base[(j2n - 1) * ldb];
					} else {
						info->stat = -1111;
						goto end;
					}
				}
				/*
				// perform the matrix-vector operation.
				*/
				clk1 = clock();
#ifdef TRL_FORTRAN_COMPATIBLE
				op(&nrow, &i__1, qa, &ldqa, rr, &ldrr);
#else
				op(nrow, i__1, qa, ldqa, rr, ldrr, info->mvparam);
#endif
				add_clock_ticks(info, &(info->clk_op), &(info->tick_o), clk1);
				info->matvec = info->matvec + 1;
				//
				/* compute alpha(jnd) = qa' * A * qa */
				alpha[jnd - 1] = trl_ddot(nrow, qa, c__1, rr, c__1);
				trl_g_sum(info->mpicom, 1, &alpha[jnd - 1], wrk2);
				/*
				// the Lanczos orthogonalization (three-term recurrence).
				//   rr = rr - alpha(jnd)*qa - beta(jnd-1)*qb
				*/
				d__1 = -alpha[jnd - 1];
				trl_daxpy(nrow, d__1, qa, c__1, rr, c__1);
				d__1 = -beta[jnd - 2];
				trl_daxpy(nrow, d__1, qb, c__1, rr, c__1);
				info->flop = info->flop + 6 * nrow;
				/*
				// re-orthogonalization, and compute beta(jnd)
				*/
				info->flop_h = info->flop_h - info->flop;
				clk1 = clock();
				trl_orth(nrow, evec, lde, j1, base, ldb, j2, rr, kept, alpha,
					beta, wrk2, lwrk2, info);
				add_clock_ticks(info, &(info->clk_orth), &(info->tick_h),
					clk1);
				info->flop_h = info->flop_h + info->flop;
				/*
				// copy alpha and beta into alfrot and betrot
				*/
				alfrot[jnd - 1] = alpha[jnd - 1];
				betrot[jnd - 1] = beta[jnd - 1];
				if (info->stat != 0)
					goto end;
				if (info->verbose > 4) {
					print_alpha_beta(info, title, jnd, alpha, beta);
				}
				/*
				perform convergence test once in a while
				*/
				if (info->matvec >= next_test) {
					if (info->verbose > 5) {
						print_all_alpha_beta(info, title, jnd, alfrot,
							betrot);
					}
					lambda = wrk2;
					res = &wrk2[jnd];
					/*
					At return of get_eval lambda are order in the ascending order
					*/
					trl_get_eval(jnd, locked, alfrot, betrot, lambda, res,
						&wrk2[jnd + jnd + 1], lwrk2 - jnd - jnd,
						&(info->stat));
					if (info->stat != 0)
						goto end;
					if (info->verbose > 2) {
						print_lambda_res(info, jnd, lambda, res);
					}
					i1 = min(mev, jnd);
					memcpy(eval, wrk2, i1 * sizeof(double));
					/*
					At return from convergence_test, lambda are order in the
					ascending order of the distance from ref if lohi < -1
					*/
					trl_convergence_test(jnd, lambda, res, info,
						&wrk2[jnd + jnd]);
					/*
					decide when to perform the next test
					*/
					degen = trl_check_dgen(info, jnd, lambda, res);
					if ((degen != -1 || info->nec < info->ned)
						&& info->nec > 0) {
							/*
							assuming a same number of matrix-vector product is
							required for each eigenvalues to converge.
							*/
							next_test =
								(double) (info->ned * info->matvec) /
								(double) (info->nec);
					} else if (info->nec == 0) {
						next_test = next_test + next_test;
						if (info->maxlan == info->ntot) {
							next_test =
								(int) ceil(0.5 * (info->maxlan + info->matvec));
						}
					}
					if (info->verbose > 0)
						trl_print_progress(info);
				}
			}
			/*
			// ************************************************************* //
			// end of inner (regular Lanczos three-term recurrence) loop     //
			// ************************************************************* //
			*/
			/*
			// error checking for debugging use
			*/
			//tick1 = clock();
			lambda = wrk2;
			res = &wrk2[jnd];
			if (info->verbose > 6) {
				wrk2 = &wrk2[jnd + jnd];
				i2 = lwrk2 - jnd - jnd;
				trl_check_orth(info, nrow, evec, lde, j1n, base, ldb, j2n,
					wrk2, i2);
				if (info->verbose > 7) {
					trl_check_recurrence(op, info, nrow, mev, evec, lde, j1n,
						base, ldb, j2n, kept, alpha, beta,
						wrk2, i2);
				}
			}
			/*
			// convert the integer_ counters to floating-point counters
			*/
			i2 = info->clk_max / 4;
			if (info->flop > i2) {
				info->rflp = info->rflp + info->flop;
				info->flop = 0;
			}
			if (info->flop_h > i2) {
				info->rflp_h = info->rflp_h + info->flop_h;
				info->flop_h = 0;
			}
			if (info->flop_r > i2) {
				info->rflp_r = info->rflp_r + info->flop_r;
				info->flop_r = 0;
			}
			if (info->clk_op > i2) {
				info->tick_o =
					info->tick_o + info->clk_op / (double) (info->clk_rate);
				info->clk_op = 0;
			}
			if (info->clk_orth > i2) {
				info->tick_h =
					info->tick_h + info->clk_orth / (double) (info->clk_rate);
				info->clk_orth = 0;
			}
			if (info->clk_res > i2) {
				info->tick_r =
					info->tick_r + info->clk_res / (double) (info->clk_rate);
				info->clk_res = 0;
			}
			info->flop_r = info->flop_r - info->flop;
			/*
			// *** Determine whether to restart ***
			// compute the Ritz values and Ritz vectors if they are not up to
			// date
			*/
			clk1 = clock();
			prek = kept;
			jml = jnd - locked;
			i2 = kept - locked + 1;
			if (degen != -1 || info->nec < info->ned) {
				/* need to compute the updated Ritz values and residual norms */
				wrk2 = &wrk[4 * info->maxlan + 2 * jnd];
				lwrk2 = lwrk - i2 * i2 - 4 * info->maxlan - 2 * jnd;
				if (lwrk2 < 3 * jnd) {
					info->stat = -12;
					goto end;
				}
				if (info->verbose > 5) {
					print_all_alpha_beta(info, title, jnd, alfrot, betrot);
				}
				/*
				Given tridiagonal matrix (diagonals stored in alfrot, and
				off-diagonals stored in betrot), computes Ritz value
				(approximate eigenvalues), using dstqrb, and retrned in
				lambda. the last components of the eigenvectors are also
				computed.  At return, lambda are stored in the ascending
				order.
				*/
				trl_get_eval(jnd, locked, alfrot, betrot, lambda, res, wrk2,
					lwrk2, &(info->stat));
				if (info->stat != 0)
					goto end;
				if (info->verbose > 2) {
					print_lambda_res(info, jnd, lambda, res);
				}
				/*
				At return, lambda are stored in the ascending order of the
				distance from ref if lohi < -1 otherwise, they are sorted in
				the ascending order of lambda.
				*/
				trl_convergence_test(jnd, lambda, res, info, wrk2);
				if (info->verbose > 0) {
					trl_print_progress(info);
				}
				degen = trl_check_dgen(info, jnd, lambda, res);
			}
			/*
			Given the tridiagonal matrix and Ritz values, compute the Ritz
			vectors (rotational matrix, used for tridiagonalization, is also
			applied).  Also, decide how many vectors to save if restart
			*/
			//tick2 = clock();
			//time4 += (tick2-tick1);
			if ((degen != -1 || info->nec < info->ned)
				&& info->matvec < info->maxmv) {
					/*
					prepare to restart, reorder the eigenvales based on the input
					parameter.  At return, lambda kept are ordered in the
					ascending order.
					*/
					trl_shuffle_eig(jml, info->klan - locked, &lambda[locked],
						&res[locked], info, &kept, locked);
					/*
					compute eigenvectors using dstein (inverse interations)
					*/
					if (kept * 3 < jml) {
						i1 = 4 * info->maxlan + jnd + kept * jml;
						yy = &wrk[4 * info->maxlan + jnd];
						wrk2 = &wrk[i1];
						lwrk2 = lwrk - i1 - i2 * i2;
						trl_get_tvec(jml, &alfrot[locked], &betrot[locked], 0, i2,
							rot, kept, &lambda[locked], yy, iwrk, wrk2,
							lwrk2, &(info->stat));
						if (info->stat == 0 && (locked + kept) > 0) {
							memcpy(alpha, lambda,
								(locked + kept) * sizeof(double));
						}
					}
					/*
					compute eigenvectors using dsyev (QR)
					*/
					if (kept * 3 >= jml || info->stat != 0) {
						if ((locked + kept) > 0)
							memcpy(alfrot, lambda,
							(locked + kept) * sizeof(double));
						i1 = 4 * info->maxlan + jml * jml;
						yy = &wrk[4 * info->maxlan];
						wrk2 = &wrk[i1];
						lwrk2 = lwrk - i1;

						trl_get_tvec_a(jml, prek - locked, &alpha[locked],
							&beta[locked], kept, &alfrot[locked], yy,
							wrk2, lwrk2, iwrk, &(info->stat));
					}
					if (info->stat != 0)
						goto end;
					for (i = 0; i < kept; i++) {
						beta[locked + i] = yy[(i + 1) * jml - 1] * betrot[jnd - 1];
					}
					if (jml > info->ned + (info->ned / 5 + 6)) {
						trl_set_locking(jml, kept, &alpha[locked], &beta[locked],
							yy, info->anrm, &i2);
					} else {
						i2 = 0;
					}
					/*
					generate Ritz vectos, reclaim the space pointed by ROT
					*/
					i1 = 4 * info->maxlan + kept * jml + jnd;
					wrk2 = &wrk[i1];
					lwrk2 = lwrk - i1;
					trl_ritz_vectors(nrow, locked, kept, yy, jml, evec, lde, j1,
						base, ldb, j2, wrk2, lwrk2);
					info->flop = info->flop + 2 * nrow * jml * kept;
					if (info->verbose > 0) {
						TRL::print_restart_state(info, title, nrow, mev, alpha, beta,
							betrot, evec, yy, kept, locked, iwrk,
							wrk2, i2, jml);
					}
					/*
					reset the counters and indices to the correct values for
					restarting
					*/
					kept += locked;
					locked += i2;
					info->locked = locked;
					jnd = kept;
					if (jnd <= mev) {
						j1 = jnd;
						j2 = 0;
					} else {
						j1 = mev;
						j2 = jnd - mev;
						if (j2 >= (nbas - 1)) {
							info->stat = -1111;
							goto end;
						}
					}
					if (info->nec > 0) {
						next_test =
							(int) (((double) (info->matvec * info->ned)) /
							((double) info->nec));
					} else {
						next_test = next_test + info->maxlan;
					}
					i1 = min(mev, jnd);
					if (i1 > 0)
						memcpy(eval, lambda, i1 * sizeof(double));
					/* copying the last Lanczos vector at the end of kept Ritz
					vectors */
					if (jnd < mev) {
						j1n = j1 + 1;
						j2n = 0;
						memcpy(&evec[(j1n - 1) * lde], rr, nrow * sizeof(double));
					} else {
						j1n = mev;
						j2n = j2 + 1;
						memcpy(&base[(j2n - 1) * ldb], rr, nrow * sizeof(double));
					}
					/*
					// write checkpoint files
					*/
					//printf( "%d %d:\n",info->matvec,chkpnt );
					if (info->matvec >= chkpnt) {
						TRL::write_checkpoint(info, title, nrow, alpha, beta, evec, lde,
							base, ldb, j1n, jnd, j2n);
						chkpnt = chkpnt + info->maxmv / info->cpflag_;
					}
			} else {
				/*
				all wanted eigenpairs converged or maximum MATVEC used sort
				the eigenvalues in final output order
				*/
				kept = min(info->nec, max(info->ned, mev - 1));
				info->nec = kept;
				if (kept == 0)
					kept = min(mev - 1, info->ned);
				trl_sort_eig(jnd, info->lohi, kept, info->ref, lambda, res);
				memcpy(eval, lambda, kept * sizeof(double));
				if (kept * 3 < jnd) {
					/*
					eigenvectors of the projection matrix (try inverse
					interations)
					*/
					i1 = kept * jnd + 4 * info->maxlan;
					yy = &wrk[4 * info->maxlan];
					wrk2 = &wrk[i1];
					lwrk2 = lwrk - i1 - i2 * i2;
					trl_get_tvec(jnd, alfrot, betrot, locked, i2, rot,
						kept, eval, yy, iwrk, wrk2, lwrk2,
						&(info->stat));
				}
				if (kept * 3 >= jnd || info->stat != 0) {
					/*
					// too many eigenvectors or inverse iterations have failed,
					// try QR
					*/
					i1 = 4 * info->maxlan + jnd * jnd;
					yy = &wrk[4 * info->maxlan];
					wrk2 = &wrk[i1];
					lwrk2 = lwrk - i1;
					trl_get_tvec_a(jnd, prek, alpha, beta, kept, eval, yy,
						wrk2, lwrk2, iwrk, &(info->stat));
					if (info->stat != 0)
						goto end;
				}
				if (kept > 0)
					memcpy(alpha, eval, kept * sizeof(double));
				for (i = 0; i < kept; i++) {
					beta[i] = betrot[jnd - 1] * yy[(1 + i) * jnd - 1];
				}
				/*
				// generate eigenvectos, reclaim the space pointed by ROT
				*/
				i1 = kept * jnd + 4 * info->maxlan;
				wrk2 = &wrk[i1];
				lwrk2 = lwrk - i1;

				trl_ritz_vectors(nrow, 0, kept, yy, jnd, evec, lde, j1, base,
					ldb, j2, wrk2, lwrk2);

				info->flop = info->flop + 2 * nrow * jml * kept;
				if (info->verbose > 1) {
					TRL::print_final_state(info, title, nrow, mev, eval, beta,
						evec, yy, kept, jml);
				}
				/*
				// reset the counters and indices to be used by check_orth and
				// check_recurrence
				*/
				jnd = kept;
				j1 = kept;
				j2 = 0;
				if (j1 < mev) {
					j1n = j1 + 1;
					j2n = 0;
					memcpy(&evec[(j1n - 1) * lde], rr, nrow * sizeof(double));
				} else {
					j1n = mev;
					j2n = 1;
					memcpy(base, rr, nrow * sizeof(double));
				}
				/*
				// write checkpoint files
				*/
				if (info->cpflag_ > 0) {
					TRL::write_checkpoint(info, title, nrow, alpha, beta, evec, lde,
						base, ldb, j1n, jnd, j2n);
				}
			}
			/*
			// check the orthogonality of the basis vectors before restarting
			*/
			if (info->verbose > 6) {
				trl_check_orth(info, nrow, evec, lde, j1n, base, ldb, j2n,
					wrk2, lwrk2);
				if (info->verbose > 7) {
					trl_check_recurrence(op, info, nrow, mev, evec, lde, j1n,
						base, ldb, j2n, kept, alpha, beta,
						wrk2, lwrk2);
				}
			}
			add_clock_ticks(info, &(info->clk_res), &(info->tick_r), clk1);
			info->flop_r = info->flop_r + info->flop;
			info->nloop = info->nloop + 1;
			//printf( "nloop: %d\n",info->nloop );
	}
	/*
	// ******************* //
	// end of restart_loop //
	// ******************* //
	*/
	/* write the estimated residual norms to the beginning of WRK */
	for (i = 0; i < j1; i++) {
		wrk[i] = fabs(beta[i]);
	}
end:
	if (info->stat < 0 && (info->verbose > 0 || info->my_pe == 0)) {
		TRL::log_error_state(info, kept, j1, j2, jnd, nrow, mev, eval, alpha,
			alfrot, beta, betrot, evec, base, qa, qb, rr,
			title, iwrk);
	}
	free(iwrk);
	return;
	/*
	// .. end of lanczos_ ..
	*/
}


void TRL::trl_initial_guess(int nrow, double *evec, int lde, int mev,
					   double *base, int ldb, int nbas, double *alpha,
					   double *beta, int *j1, int *j2, trl_info * info,
					   double *wrk, int lwrk)
{
	/*
	// Purpose
	// =======
	// check to make sure the initial guess vector contains valid nonzero numbers if not fill with
	// random numbers this routine will also read the checkpoint files to restore the previous state// of the Lancozs iterations
	//
	// Arguments
	// =========
	// nrow   (input) Integer
	//         On entry, specifies the number of rows in eigenvectors.
	//
	// evec   (input/output) Double array of dimension (lde,mev)
	//         On entry, the (nec+1)st column contains the initial guess.
	//
	// lde    (input) Integer
	//         On entry, specifies the leading dimention of evec.
	//
	// mev    (input) Integer
	//         On entry, specifies the number of Ritz vectors that can be stored in evec.
	//
	// base   (input/output) Double array of dimension (ldb,nbas)
	//         Stores the Ritz vectors, that cannot be stored in evec.
	//
	// ldb    (input) Integer
	//         On entry, specifies the leading dimention of base.
	//
	// nbas   (input) Integer
	//         On entry, specifies the number of Ritz vectors that can be stored in base
	//
	// alpha  (input/output) Double array of dimension (mev+nbas-1)
	//         On exit, stores alpha values if checkpoint files are provided.
	//
	// beta   (input/output) Double array of dimension (mev+nbas-1)
	//         On exit, stores beta values if checkpoint files are provided.
	//
	// j1     (output) Pointer to integer
	//         On exit, stores j1 (number of Ritz vectors in evec) if checkpoint files are
	//         provided.
	//
	// j2     (output) Pointer to integer
	//         On exit, stores j1 (number of Ritz vectors in base) if checkpoint files are
	//         provided.
	//
	// info   (input/output) Pointer to trl_info structure
	//         On entry, points to the data structure to store the current information about
	//         the eigenvalue problem and the progress of TRLAN.
	//
	// wrk    (workspace) Double array of dimension (lwrk)
	//
	// lwrk   (input) Integer
	//         Specifies the size of the workspace.
	//
	// Parameters
	*/
	long c__1 = 1;
	//
	// local variable
	//
	int i, j, k, nran, north;
	//double tmp, rnrm=1.09516e+292;
	//double tmp, rnrm=1.0e+10;
	//double tmp, rnrm=DBL_MAX;
	double tmp, rnrm=par::trl_init_g;
	clock_t ii, jj;
	char file[STRING_LEN];
	//
	// generate random seeds based on current clock ticks
	ii = clock();
	if (info->my_pe > 0) {
		ii = ii - (int) (info->my_pe * sqrt((double) ii));
	}
	srand48_new();
	//
	j = info->nec;
	if (info->guess > 1) {
		// retrieve a check-point file
		i = info->cpio;
		if (info->oldcpf != 0 && strlen(info->oldcpf) > 0) {
			trl_pe_filename(STRING_LEN, file, info->oldcpf, info->my_pe,
				info->npes);
		} else {
			trl_pe_filename(STRING_LEN, file, info->cpfile, info->my_pe,
				info->npes);
		}

		ii = clock();
		i = trl_read_checkpoint(file, nrow, &evec[j * lde], lde,
			mev - info->nec, j1, base, ldb, nbas, j2,
			(mev + nbas - 1 - j), &alpha[j],
			(mev + nbas - 1 - j), &beta[j]);
		info->stat = trl_sync_flag_(info->mpicom, i);
		jj = clock();
		if (jj > ii) {
			info->clk_in = jj - ii;
		} else {
			info->clk_in = (info->clk_max - ii) + jj;
		}
		info->wrds_in = (*j1 + *j2) * (nrow + nrow + 2) + nrow + 2;
		*j1 = *j1 + info->nec;
		if (info->stat != 0)
			return;
	} else {
		if (info->guess <= 0) {
			// generate an arbitrary initial starting vector
			// if (info->guess == 0), use the vector [1, 1, ...]^T
			// else perturb some random locations of the above vector
			for (k = 0; k < nrow; k++) {
				evec[j * lde + k] = 1.0;
			}
			nran = min(1 - info->guess, lwrk);
			nran = 2 * (nran / 2);
			if (nran > 0 && nran < nrow) {
				for (k = 0; k < nran; k++) {
					wrk[k] = drand48_new();
				}
				for (i = 0; i < nran - 1; i += 2) {
					ii = (int) (nrow * wrk[i]);
					evec[j * lde + ii] =
						evec[j * lde + ii] + wrk[i + 1] - 0.5;
				}
				info->flop = info->flop + nran + nran;
			} else if (nran >= nrow) {
				for (i = 0; i < nrow; i++) {
					evec[j * lde + i] = drand48_new();
				}
				TRL::trl_smooth_rr(nrow, &evec[(info->nec) * lde]);
				info->nrand++;
				info->flop += 4 * nrow;
			}
		}
		*j1 = info->nec;
		*j2 = 0;
	}
	tmp = 0.0;
	// make sure the norm of the next vector can be computed
	//wrk[0] = trl_ddot(nrow, &evec[j * lde], c__1, &evec[j * lde], c__1);
	double tmmmp = trl_ddot(nrow, &evec[j * lde], c__1, &evec[j * lde], c__1);
	wrk[0]=tmmmp;

	trl_g_sum(info->mpicom, 1, &wrk[0], &wrk[1]);
	info->flop = info->flop + nrow + nrow;
	if (wrk[0] >= DBL_MIN && wrk[0] <= DBL_MAX) {
		// set rnrm to let trl_CGS normalize evec(1:nrow, j)
		//cout<<"set rnrm to let trl_CGS normalize evec(1:nrow, j)\t"<<"wrk[0]"<<wrk[0]<<"\n";
		rnrm = sqrt(wrk[0]);
	} else {
		for (i = 0; i < nrow; i++) {
			evec[j * lde + i] = drand48_new();
			//cout<<evec[j * lde + i]<<"\t";
		}
		//cout<<"TRL::trl_smooth_rr(nrow, &evec[(info->nec) * lde]);";
		TRL::trl_smooth_rr(nrow, &evec[(info->nec) * lde]);
		info->nrand++;
		info->flop += 4 * nrow;
	}
	//
	// orthogonalize initial guess against all existing vectors
	//
	i = 0;
	tmp = 1.0;
	nran = info->nrand;
	north = info->north;
	//cout<<*j1<<"\t"<<mev<<"\n";
	//cout<<info->stat<<"\n";
	//-204
	
	if (*j1 < mev) {
		info->stat = trl_cgs(info, nrow, evec, lde, *j1, base, ldb, 0,
			&evec[(*j1) * lde], &rnrm, &tmp, &i, wrk);
	} else if (*j2 <= 0) {
		info->stat = trl_cgs(info, nrow, evec, lde, *j1, evec, lde, 0,
			base, &rnrm, &tmp, &i, wrk);

	} else {
		info->stat = trl_cgs(info, nrow, evec, lde, *j1, base, ldb, *j2,
			&base[(*j2) * ldb], &rnrm, &tmp, &i, wrk);
	}
	//cout<<info->stat<<"\n";

	info->flop =
		info->flop + 4 * nrow * ((info->north - north) * (*j1 + *j2) +
		info->nrand - nran)
		+ nrow;
	if (info->verbose > 6) {
		if (*j1 < mev) {
			i = *j1 + 1;
			ii = *j2;
		} else {
			i = *j1;
			ii = *j2 + 1;
		}
		trl_check_orth(info, nrow, evec, lde, *j1, base, ldb, ii, wrk,
			lwrk);
	}
	return;
	//
	// .. end of trl_initial_guess ..
	//
}



void TRL::trl_smooth_rr(int n, double *rr)
{
	/*
	// Purpose
	// =======
	// Smooth out a vector, i.e.,
	//  rr = rr + rr + Cshift(rr, 1) + Cshift(rr, -1) in Fortran.
	// Used in trl_initial_guess.
	//
	// Arguments
	// =========
	// n   (input) Integer
	//      On entry, specifies the dimension of rr.
	//
	// rr  (input/output) Double precision array of dimension (n)
	//      On entry, contains the initial state of the vector. On exit, contain the 
	//      vector after the smothing is applied.
	//
	// ..
	// .. local scalars ..
	*/
	int i;
	/*
	// ..
	// .. executable statements ..
	*/
	if (n <= 0)
		return;
	double rr1 = rr[0], rr2;
	rr2 = rr1;
	rr[0] = 2 * rr[0] + rr[2] + rr[n - 1];
	for (i = 1; i < n - 1; i++) {
		double tmp = rr[i];
		rr[i] = 2 * rr[i] + rr[i + 1] + rr2;
		rr2 = tmp;
	}
	rr[n - 1] = 2 * rr[n - 1] + rr[1] + rr2;
	/*
	// .. end of trl_smooth_rr ..
	*/
}


////
int TRL::trl_cgs(trl_info * info, int nrow, double *v1, int ld1, int m1,
			double *v2, int ld2, int m2, double *rr, double *rnrm,
			double *alpha, int *north, double *wrk)
{
	//
	// Purpose
	// =======
	// Perform full Gram-Schmidt routine -- orthogonalize a new vector against all existing vectors.//
	// Arguments
	// =========
	// info   (input) Pointer to structure trl_info_
	//         On entry, points to the current TRL_INFO.
	//
	// nrow   (input) Integer
	//         On entry, specifies the number of rows in eigenvectors.
	//
	// v1     (input) double precision array (ld1,m1)
	//         On entry, contains the first part of Lanczos basis computed.
	//
	// ld1    (input) Integer
	//         On entry, specifies the leading dimention of v1.
	//
	// m1     (input) Integer
	//         On entry, specifies the number of Lanczos basis in v1.
	//
	// v2     (input) double precision array (ld2,m2)
	//         On entry, contains the second part of Lanczos basis computed.
	//
	// ld2    (input) Integer
	//         On entry, specifies the leading dimention of v2.
	//
	// m2     (input) Integer
	//         On entry, specifies the number of Lanczos basis in v2.
	//
	// rr     (input/output) double precision array (nrow)
	//         On entry, contains the new Lanczos basis computed.
	//         On exit, contains the next Lanczos basis computed after the orthogonalization.
	//
	// rnrm   (output) double precision
	//         On entry, specifies the norm of the current Lanczos basis.
	//         On exit, specifies the norm of the new Lanczos basis.
	//
	// alpha  (input/output) double precision array (m1+m2)
	//         On entry, contains the alpha values, on exit, they are updated.
	//
	// north  (output)
	//         On exit, specifies the number of times the full-orthogonalization is applied.
	//
	// wrk    (workspace) complex array (m1+m2)
	//
	// ..
	// .. local parameters ..
	char notrans = 'N';
	double one = 1.0, zero = 0.0;
	int maxorth = 3;
	//
	// ..
	// .. local variables ..
	double d__1;
	int mpicom, myid, i, j, k, nold, irnd, cnt, ierr = 0;
	double tmp, old_rnrm;
	//
	// ..
	// .. executable statements ..
	mpicom = info->mpicom;
	myid = info->my_pe;
	nold = m1 + m2;
	//cout<<nold;
	if (ld1 < nrow || (ld2 < nrow && m2 > 0)) {
		return -201;
	}
	irnd = 0;
	ierr = 0;

	//cout<<"\n*rnrm: "<<*rnrm<<" DBL_MIN:"<<DBL_MIN<<"\n";

	if (nold > 0) {
		cnt = 0;
		while (cnt <= maxorth) {
			// compute [v1 v2]'*rr=wrk
			trl_g_dot_(mpicom, nrow, v1, ld1, m1, v2, ld2, m2, rr, wrk);
			if (m1 > 1) {
				d__1 = -one;
				trl_dgemv(&notrans, nrow, m1, d__1, v1, ld1, wrk, one,
					one, rr, one);
			} else if (m1 == 1) {
				d__1 = -wrk[0];
				trl_daxpy(nrow, d__1, v1, one, rr, one);
			}
			if (m2 > 1) {
				d__1 = -one;
				trl_dgemv(&notrans, nrow, m2, d__1, v2, ld2, &wrk[m1],
					one, one, rr, one);
			} else if (m2 == 1) {
				d__1 = -wrk[nold - 1];
				trl_daxpy(nrow, d__1, v2, one, rr, one);
			}
			/*
			if( irnd == 0) {
			tmp = (ld1*nold)*DBL_EPSILON*max(fabs(*alpha), *rnrm);
			if( fabs(wrk[nold-1]) > tmp && tmp > zero) {
			return -202;
			}
			*alpha = *alpha + wrk[nold-1];
			}
			*/
			(*north)++;
			cnt = cnt + 1;
			tmp = trl_ddot(nold, wrk, one, wrk, one);
			wrk[0] = trl_ddot(nrow, rr, one, rr, one);
			trl_g_sum(mpicom, 1, wrk, &wrk[1]);
			*rnrm = sqrt(wrk[0]);
			old_rnrm = sqrt(wrk[0] + tmp);
			//
			// decisions about whether to re-orthogonalize is based on
			// relative size of tmp and wrk(1) (R norm square)
			//printf( "wrk[0]=%e tmp=%e\r\n",wrk[0],tmp );
			//if( wrk[0] > tmp) {
			if (DBL_EPSILON * wrk[0] > tmp) {
				// no need for more orthogonalization
				cnt = maxorth + 1;
				//            } else if( ((wrk[0] <= DBL_EPSILON*tmp && cnt > 1) ||
				//                      !( wrk[0] > DBL_MIN)) && irnd < maxorth ) {
				//
				//            } else if( ((wrk[0] <= DBL_EPSILON*tmp && cnt > 1) ||
				//                      !( wrk[0] > 100.0*DBL_EPSILON*DBL_EPSILON*info->ntot)) &&
				//                      irnd < maxorth ) {
			} else if (((cnt > 1 && !
				(tmp <=
				info->ntot * DBL_EPSILON * DBL_EPSILON *
				(wrk[0] + tmp))) || !(wrk[0] > DBL_MIN))
				&& irnd < maxorth) {
					// the input vector is so nearly linear dependent on the
					// existing vectors, we need to perturb it in order to
					// generate a new vector that is orthogonal to the existing
					// ones
					// the perturbation is done in two stages:
					// -- perturbing only one number
					// -- call random_number to generate a whole random vector
					cnt = 0;
					irnd++;
					info->nrand++;
					if (irnd == 1 && *rnrm > 0.0
						&& *rnrm > DBL_EPSILON * old_rnrm) {
							// old version:  modify only one element
							// new version:  modify (nrow * epsilon * rnrm / old_rnrm ) elements
							tmp = drand48_new();
							i = (int) (nrow * tmp);
							k = i +
								(int) (max
								(1.0,
								(nrow * DBL_EPSILON * (*rnrm)/old_rnrm  )));
								//(nrow * DBL_EPSILON * old_rnrm / *rnrm)));
							for (j = i; j < k && j<nrow ; j++) {
								tmp = drand48_new();
								while (fabs(tmp - 0.5) <= DBL_EPSILON) {
									//tmp = drand48_new();
									tmp = drand48_new();
								}
								rr[j] += (*rnrm) * (tmp - 0.5);
							}
					} else {
						// fill with random numbers produced by intrinsic function
						for (i = 0; i <= myid; i++) {
							tmp = drand48_new();
						}
						i = (int) (nrow * tmp);
						tmp = drand48_new();
						j = (int) (nrow * tmp);
						if (i < j) {
							for (k = i; k <= j; k++) {
								rr[k] = drand48_new();
							}
						} else if (i > j) {
							for (k = j; k <= i; k++) {
								rr[k] = drand48_new();
							}
						} else {
							for (k = 0; k < nrow; k++) {
								rr[k] = drand48_new();
							}
						}
					}
					//rr = rr + rr + Cshift(rr, 1) + Cshift(rr, -1)
					trl_smooth_rr(nrow, rr);
			}
		}
		// failed to reduce the dot-products between the new vector
		// and the old vectors to small number -- usually an indication of
		// problem with orthogonality of the old vectors
		if (!(wrk[0] >= tmp))
			ierr = - 203;
	}
	//
	// normalization
	//
	if (ierr == 0) {
		
		if ( *rnrm > DBL_MIN ) {
			//cout<<"\n*rnrm: "<<*rnrm<<" > DBL_MIN:"<<DBL_MIN<<"\n";
			tmp = one / *rnrm;
			trl_dscal(nrow, tmp, rr, one);
		} else{
			//cout<<"\n*rnrm: "<<*rnrm<<" < DBL_MIN:"<<DBL_MIN<<"\n";
			return -204;
		}
		//tmp = one / *rnrm;
		//trl_dscal(nrow, tmp, rr, one);
	}
	if (irnd > 0)
		*rnrm = zero;
	return ierr;
	//
	// .. end of trl_cgs ..
	//
}


void TRL::log_error_state(trl_info * info, int kept, int j1, int j2, int jnd,
					 int nrow, int mev, double *eval, double *alpha,
					 double *alfrot, double *beta, double *betrot,
					 double *evec, double *base, double *qa, double *qb,
					 double *rr, char *title, int *iwrk)
{
	/*
	// Purpose
	// =======
	// Dump important variables when return due to error
	//
	// Arguments
	// =========
	// info       (input) Pointer to structure trl_info_
	//             On entry, points to the current TRL_INFO.
	//
	// kept       (input) Integer
	//             On entry, specifies the number of lanczos vector kept at the last 
	//             restart.
	//
	// j1         (input) Integer
	//             On entry, specifies the last column of evec, that contains a lanczos 
	//             vector.
	//
	// j2         (input) Integer
	//             On entry, specifies the last column of base, that contains a lanczos 
	//             vector.
	//
	// jnd        (input) Integer
	//             On entry, specifies the number of lanczos vectors computed.
	//
	// nrow       (input) Integer
	//             On entry, specifies the number of rows in the eigenvectors.
	//
	// mev        (input) Integer
	//             On entry, specifies the maximum number of eigenvalues allowed.
	//
	// eval       (input) Double precision array of dimension (mev)
	//             On entry, contains the eigenvalues computed.
	//
	// alpha      (input) Double precision array of dimension (info->maxlan)
	//             On entry, contains the values of alpha computed.
	//
	// alfrot     (input) Double precision array of dimension (info->maxlan)
	//             On entry, contains the values of alpha after rotation 
	//             (diagonalization).
	//
	// beta       (input) Double precision array of dimension (info->maxlan)
	//             On entry, contains the values of beta computed.
	//
	// betrot     (input) Double precisino array of dimension (info->maxlan)
	//             On entry, contains the values of beta after rotation 
	//             (diagonalization).
	//
	// evec       (input) Double precision array of dimension (nrow,mev)
	//             On entry, contains the eigevectors computed.
	//
	// base       (input) Double precision array of dimension (nrow,nbas)
	//             On entry, contains the lanczos vectors, that did not fit in evec.
	//
	// qa         (input) Double precision array of dimension (nrow)
	//             On entry, contains the lanczos vector from the last iteration.
	//
	// qb         (input) Double precision array of dimension (nrow)
	//             On entry, contains the lanczos vector from the two iterations ago.
	//
	// rr         (input) Double precision array of dimension (nrow)
	//             On entry, contains the current lanczos vector being computed.
	//
	// title      (workspace) String length of (STRING_LEN)
	//             On entry, provides a space to store the title of the information 
	//             printed out.
	//
	// iwrk       (workspace) Integer array of dimension ()
	//
	// ..
	// .. local variables ..
	*/
	FILE *fp = info->log_fp;
	/*
	// ..
	// .. executable statements ..
	*/
	trl_time_stamp(fp);
	strcpy(title, "Dumping the content of the variables on error..");
	iwrk[0] = info->stat;
	trl_print_int(info, title, 1, iwrk, 1);
	trl_terse_info(info, fp);
	fprintf(fp, "This Lanczos iteration started with %d vectors.\n", kept);
	fprintf(fp, "There are %d (%d, %d) Lanczos vectors currently.\n", jnd,
		j1, j2);
	if (jnd != j1 + j2)
		jnd = j1 + j2;
	if (jnd < 0 || jnd > info->klan)
		jnd = 0;
	strcpy(title, "Content of eval ..");
	trl_print_real(info, title, mev, eval, 1);
	if (jnd > 0) {
		sprintf(title, "Alpha(1:%d) .. ", jnd);
		trl_print_real(info, title, jnd, alpha, 1);
		sprintf(title, " Beta(1:%d) .. ", jnd);
		trl_print_real(info, title, jnd, beta, 1);
		sprintf(title, "Alfrot(1:%d) .. ", jnd);
		trl_print_real(info, title, jnd, alfrot, 1);
		sprintf(title, "betrot(1:%d) .. ", jnd);
		trl_print_real(info, title, jnd, betrot, 1);
	}
	if (j1 > 0) {
		strcpy(title, "the First row of evec ..");
		trl_print_real(info, title, j1, evec, nrow);
		sprintf(title, "row %d of evec ..", nrow);
		trl_print_real(info, title, j1, &evec[nrow - 1], nrow);
	}
	if (j2 > 0) {
		strcpy(title, "the First row of base ..");
		trl_print_real(info, title, j2, base, nrow);
		sprintf(title, "row %d of base ..", nrow);
		trl_print_real(info, title, j2, &base[nrow - 1], nrow);
	}
	if (qb != NULL) {
		sprintf(title, "Content of qb (q_%d) ..", jnd - 1);
		trl_print_real(info, title, nrow, qb, 1);
	}
	if (qa != NULL) {
		sprintf(title, "Content of qa (q_%d) ..", jnd);
		trl_print_real(info, title, nrow, qa, 1);
	}
	if (rr != NULL) {
		sprintf(title, "Content of rr (residual == q_%d) ..", jnd + 1);
		trl_print_real(info, title, nrow, rr, 1);
	}
	if (info->my_pe == 0 && info->log_fp != stdout) {
		printf("TRLanczos returned with error\n");
		printf("Contents of most variables are dumped to log file %s.\n",
			info->log_file);
	}
	/*
	//  .. end of print_error_state_ ..
	*/
}


void TRL::add_clock_ticks(trl_info * info, clock_t *time, double *rtime,
					 clock_t clk1)
{
	/*
	// ..
	// .. local variables ..
	*/
	clock_t clk2, clk3;
	/*
	// ..
	// .. executable statements ..
	*/
	clk2 = clock();
#ifdef __CLK_RESTART
	if (clk2 <= clk1) {
		clk3  = (info->clk_max - clk1);
		clk3 += clk2;
	} else {
		clk3 = clk2 - clk1;
	}
#else
	if (clk2 < clk1) {
#ifdef DEBUG
		printf( "DEBUG -- 1: clk1: %e clk2: %e clk3: %e\n",
			(double)clk1,(double)clk2,(double)(info->clk_max) );
#endif
		clk3  = (info->clk_max - clk1);
		clk3 += (info->clk_max + clk2);
	} else if (clk1 < 0 && clk2 >= 0) {
#ifdef DEBUG
		printf( "1: clk1: %e clk2: %e clk3: %e\n",
			(double)clk1,(double)clk2,(double)(info->clk_max) );
#endif
		clk3  = -clk1;
		clk3 +=  clk2;
	} else {
		clk3 = clk2 - clk1;
	}
#endif
	if (clk3 + (*time) >= (*time)) {
		*time = clk3 + *time;
	} else {
		*rtime = (*rtime) + ((*time) + clk3) / (double) (info->clk_rate);
		*time = 0;
	}
	/*
	//  .. end of add_clock_ticks ..
	*/
}


void TRL::trl_orth(int nrow, double *v1, int ld1, int m1, double *v2, int ld2,
			  int m2, double *rr, int kept, double *alpha, double *beta,
			  double *wrk, int lwrk, trl_info * info)
{
	//
	// Purpose
	// =======
	// Applies full re-orthogonalization;
	//  1. if (global re-orthogonalization is needed)
	//      call trl_cgs
	//    else
	//      perform extended local re-reorthogonalization
	//    endif
	//  2. perform normalization
	//
	// Arguments:
	// ==========
	// nrow   (input) Integer
	//         On entry, specifies the number of rows in eigenvectors.
	//
	// v1     (input) double precision array (ld1,m1)
	//         On entry, contains the first part of Lanczos basis computed.
	//
	// ld1    (input) Integer
	//         On entry, specifies the leading dimention of v1.
	//
	// m1     (input) Integer
	//         On entry, specifies the number of Lanczos basis in v1.
	//
	// v2     (input) double precision array (ld2,m2)
	//         On entry, contains the second part of Lanczos basis computed.
	//
	// ld2    (input) Integer
	//         On entry, specifies the leading dimention of v2.
	//
	// m2     (input) Integer
	//         On entry, specifies the number of Lanczos basis in v2.
	//
	// rr     (input/output) double precision array (nrow)
	//         On entry, contains the new Lanczos basis computed.
	//         On exit, contains the next Lanczos basis computed after the orthogonalization.
	//
	// kept   (input) Integer
	//         On etnry, specifies the number of Ritz vectors kept.
	//
	// alpha  (input/output) double precision array (m1+m2)
	//         On entry, contains the alpha values, on exit, they are updated.
	//
	// beta   (input/output) double precision array (m1+m2)
	//         On entry, contains the beta values, on exit, they are updated if necessary,
	//         (full orthogonalization).
	//
	// wrk    (workspace) complex array (lwrk)
	//
	// lwrk   (input) Integer
	//         Specifies the size of workspace.
	//
	// info   (input) Pointer to structure trl_info_
	//         On entry, points to the current TRL_INFO.
	//
	// ..
	// .. local parameters ..
	double zero = 0.0, one = 1.0;
	long c__1 = 1;
	//
	// ..
	// .. local variables ..
	double d__1;
	int i, usecgs, jnd, jm1, no, nr;
	double tmp;
	double *qa, *qb;
	//
	// ..
	// .. executable statements ..
	//
	// check for workspace size
	jnd = m1 + m2;
	jm1 = jnd - 1;
	tmp = zero;
	if (ld1 >= nrow && ld2 >= nrow && lwrk >= max(4, jnd + jnd)) {
		info->stat = 0;
	} else {
		info->stat = -101;
		return;
	}
	//
	// compute the norm of the vector RR
	//
	wrk[0] = trl_ddot(nrow, rr, c__1, rr, c__1);
	//cout<<"\nnorm of vector:"<<wrk[0]<<"\n";
	trl_g_sum(info->mpicom, 1, &wrk[0], &wrk[1]);
	if (!(wrk[0] >= zero) || !(wrk[0] <= DBL_MAX)) {
		info->stat = -102;
		
		return;
	}
	beta[jnd - 1] = sqrt(wrk[0]);
	tmp = alpha[jnd - 1] * alpha[jnd - 1];
	if (jm1 > kept) {
		tmp += (beta[jm1 - 1] * beta[jm1 - 1]);
		info->flop += (2 * nrow + 4);
	} else if (kept > 0) {
		tmp += trl_ddot(jm1, beta, c__1, beta, c__1);
		info->flop += (2 * (nrow + kept + 2));
	}

	if (jm1 == kept) {
		usecgs = 1;
	} else if (jnd >= info->ntot) {
		usecgs = 0;
	} else if (DBL_EPSILON * wrk[0] >= tmp) {
		double anorm = 0.0;
		for (i = 0; i < jnd; ++i) {
			d__1 = fabs(alpha[i]) + fabs(beta[i]);
			if (d__1 > anorm)
				anorm = d__1;
		}
		usecgs = (beta[jm1] < DBL_EPSILON * anorm * info->ntot);
	} else {
		usecgs = 1;
	}
	//
	// whether to perform full re-orthogonalization or extended local
	// re-orthogonalization
	if (usecgs != 0) {
		// perform global re-orthogonalization
		nr = info->nrand;
		no = info->north;
		info->stat = trl_cgs(info, nrow, v1, ld1, m1, v2, ld2, m2, rr,
			&beta[jnd - 1], &alpha[jnd - 1],
			&(info->north), wrk);
		info->flop =
			info->flop + 4 * nrow * ((info->north - no) * jnd +
			info->nrand - nr) + nrow;
	} else if (jnd > 1) {
		// perform local re-orthogonalization against two previous vectors
		if (m2 > 1) {
			qa = &v2[(m2 - 1) * ld2];
			qb = &v2[(m2 - 2) * ld2];
		} else if (m2 == 1) {
			qa = v2;
			qb = &v1[(m1 - 1) * ld1];
		} else {
			qa = &v1[(m1 - 1) * ld1];
			qb = &v1[(jm1 - 1) * ld1];
		}
		wrk[0] = zero;
		wrk[1] = zero;
		for (i = 0; i < nrow; i++) {
			wrk[0] = wrk[0] + qa[i] * rr[i];
			wrk[1] = wrk[1] + qb[i] * rr[i];
		}
		trl_g_sum(info->mpicom, 2, &wrk[0], &wrk[2]);
		alpha[jnd - 1] = alpha[jnd - 1] + wrk[0];
		d__1 = -wrk[0];
		trl_daxpy(nrow, d__1, qa, c__1, rr, c__1);
		d__1 = -wrk[1];
		trl_daxpy(nrow, d__1, qb, c__1, rr, c__1);
		tmp = one / beta[jnd - 1];
		trl_dscal(nrow, tmp, rr, c__1);
		info->flop = info->flop + 9 * nrow;
	} else {
		// perform local re-orthogonalization against the only vector
		if (m1 == 1) {
			qa = v1;
		} else {
			qa = v2;
		}
		wrk[0] = trl_ddot(nrow, qa, c__1, rr, c__1);
		trl_g_sum(info->mpicom, 1, &wrk[0], &wrk[1]);
		alpha[jnd - 1] = alpha[jnd - 1] + wrk[0];
		d__1 = -wrk[0];
		trl_daxpy(nrow, d__1, qa, c__1, rr, c__1);
		tmp = one / beta[jnd - 1];
		trl_dscal(nrow, tmp, rr, c__1);
		info->flop = info->flop + 5 * nrow;
	}
	// when beta(jnd) is exceedingly small, it should be treated as zero
	if (info->stat == 0) {
		if (beta[jnd - 1] <= DBL_EPSILON * fabs(alpha[jnd - 1])) {
			beta[jnd - 1] = zero;
		} else if (jnd >= info->ntot) {
			beta[jnd - 1] = zero;
		}
	}
	//
	// .. end of trl_orth ..
	//
}



void TRL::print_alpha_beta(trl_info * info, char *title, int i,
					  double *alpha, double *beta)
{
	/*
	// Purpose
	// =======
	// Print the Ith alpha and beta value to the log file. Function trl_print_real is 
	// used.
	//
	// Arguments
	// =========
	// info    (input) Pointer to structure trl_info_
	//          On entry, points the current TRL_INFO. The information is printed out 
	//          to the log file specified in trl_info.
	//
	// title   (workspace) String of length (STRING_LEN)
	//          On entry, provides the space to store the title to print out, i.e., 
	//          "alpha(jnd) =" and "beta(jnd) =".
	//
	// i       (input) Integer
	//          On entry, specifies the index of alpha and beta to print out.
	//
	// alpha   (input) Double array of dimension (info->maxlan)
	//          On entry, contains the alpha values.
	//
	// beta    (input) Double array of dimension (info->maxlan)
	//          On entry, contains the beta values.
	//
	// ..
	// .. executable statements ..
	*/
	sprintf(title, " alpha(%d) =", i);
	trl_print_real(info, title, 1, &alpha[i - 1], 1);
	sprintf(title, "  beta(%d) =", i);
	trl_print_real(info, title, 1, &beta[i - 1], 1);
	/*
	// .. end of print_alpha_beta ..
	*/
}


void TRL::trl_tridiag(int nd, double *alpha, double *beta, double *rot,
				 double *alfrot, double *betrot, double *wrk, int lwrk,
				 int *ierr)
{
	//
	// Purpose
	// =======
	// transforms an real symmetric arrow matrix into a symmetric tridiagonal matrix.
	//
	// Arguments
	// =========
	// nd       (input) integer
	//           On entry, specifies the dimention of the arrow matrix.
	//
	// alpha    (input) double precision array (nd)
	//           On entry, contains the alpha values.
	//
	// beta     (input) double precision array (nd)
	//           On entry, contains the beta values
	//
	// rot      (workspace) double precision array (nd*nd)
	//           Used to store the arrow matrix.
	//
	// alfrot   (output) double precision array (nd)
	//           On exit, contains alpha values after rotation.
	//
	// betrot   (output) double precision array (nd)
	//           On exit, contains beta values after rotation.
	//
	// wrk      (workspace) double precision array (lwrk)
	//
	// lwrk     (input) integer
	//           Specifies the size of workspace.
	//
	// ierr     (output) integer
	//           Returns the error from LAPACK calls.
	//
	// ..
	// .. CLAPACK subroutines..
	//extern int dsytrd_();
	//extern int dorgtr_();
	//
	// ..
	// .. local parameters ..
	char upper = 'U';
	//
	// ..
	// .. local variables ..
	int i, lwrk2;
	//
	// special case, nd == 1;
	if (nd == 0) {
		return;
	} else if (nd <= 1) {
		rot[0] = 1.0;
		alfrot[0] = alpha[0];
		betrot[0] = beta[0];
		*ierr = 0;
		return;
	}
	if (lwrk < nd + nd) {
		*ierr = -11;
		return;
	} else {
		*ierr = 0;
	}
	//
	// first form the array matrix as a full matrix in rot
	// alpha on diagonal, and beta on the last off-diagonal column and row
	memset(rot, 0, (nd * nd) * sizeof(double));
	for (i = 0; i < nd; i++) {
		rot[i * nd + i] = alpha[i];
	}
	for (i = 0; i < nd - 1; i++) {
		rot[(nd - 1) * nd + i] = beta[i];
		rot[(i + 1) * nd - 1] = beta[i];
	}
	lwrk2 = lwrk - nd;
	//
	// call LAPACK routines to reduce the matrix into tridiagonal form
	// and generate the rotation matrix
	TRL::dsytrd_(&upper, &nd, rot, &nd, alfrot, betrot, wrk, &wrk[nd], &lwrk2,
		ierr);
	if (*ierr != 0) {
		*ierr = -112;
		return;
	}
	betrot[nd - 1] = beta[nd - 1];
	TRL::dorgtr_(&upper, &nd, rot, &nd, wrk, &wrk[nd], &lwrk2, ierr);
	if (*ierr != 0) {
		*ierr = -113;
		return;
	}
	//
	// .. end of trl_tridiag ..
	//
}


void TRL::print_all_alpha_beta(trl_info * info, char *title, int jnd,
						  double *alfrot, double *betrot)
{
	/*
	// Purpose
	// =======
	// Print all computed alpha and beta to the log file. Function trl_print_real is 
	// used.
	//
	// Arguments
	// =========
	// info     (input) Pointer to structure trl_info_
	//           On entry, points to the current TRL_INFO. The information is printed 
	//           out to the log file specified in info.
	//
	// title    (workspace) String of length (12+digit of jnd)
	//           On entry, provides the space to store the title to print out, i.e., 
	//           "alfrot(1:jnd)..", and "beta(1:jnd).."
	//
	// jnd      (input) Integer
	//           On entry, specifies the number of alpha and beta computed so far.
	//
	// alfrot   (input) Double precision array of dimension (info->maxlan)
	//           On entry, contains the alpha computed so far.
	//
	// betrot   (input) Double precision array of dimension (info->maxlan)
	//           On entry, contains the beta computed so far.
	//
	// ..
	// .. executable statements ..
	*/
	sprintf(title, "alfrot(1:%d)..", jnd);
	trl_print_real(info, title, jnd, alfrot, 1);
	sprintf(title, "betrot(1:%d)..", jnd);
	trl_print_real(info, title, jnd, betrot, 1);
	/*
	// .. end of print_all_alpha_beta ..
	*/
}

void TRL::print_lambda_res(trl_info * info, int jnd, double *lambda,
					  double *res)
{
	/*
	// Purpose
	// =======
	// Print the lambda and its residual norm computed so far. Function trl_print_real 
	// is used.
	//
	// Arguments
	// =========
	// info      (input) Pointer to sructure trl_info
	//            On entry, points to the current TRL_INFO. The information is printed 
	//            out to the log file specified in info.
	//
	// jnd       (input) Integer
	//            On entry, specifies the number of lambda computed so far.
	//
	// lambda    (input) Double precision array of dimension (info->maxlan)
	//            On entry, contains the lambda computed so far.
	//
	// res       (input) Double precision array of dimension (info->maxlen)
	//            On entry, contains the residual norm of lambda computed so far.
	//
	// ..
	// .. executable statements ..
	*/
	trl_print_real(info, "Current eigenvalues..", jnd, lambda, 1);
	trl_print_real(info, "Current residual norms..", jnd, res, 1);
	/*
	// .. end of print_lambda_res
	*/
}

void TRL::trl_get_eval(int nd, int locked, double *alpha, double *beta,
				  double *lambda, double *res, double *wrk, int lwrk,
				  int *ierr)
{
	//
	// Purpose:
	// =======
	// Evaluates the eigenvalues and their corresponding residual norms of a
	// real symmetric tridiagonal matrix
	//
	// it returns eigenvalues in two sections
	//  1) the first section is the locked eigenvalues, their residual norms are zero
	//  2) the second section contains the new Ritz values in their ascending order.
	// res will contain corresponding residual norms
	//
	// Arguments;
	// ==========
	// nd         (input) integer
	//             On entry specifies, the size of alpha and beta.
	//
	// locked     (input) integer
	//             On entry, specifies the number of Ritz values locked.
	//
	// alpha      (input) double preicsion array (nd)
	//             On entry, contains the alpha values.
	//
	// beta       (input) double precision array (nd)
	//             On entry, contains the beta values.
	//
	// lambea     (output) double precision array (nd)
	//             On exit, contains the Ritz values.
	//
	// res        (output) double precision array (nd)
	//             On exit, contains the residual norm of the Ritz values.
	//
	// wrk        (workspace) double precision array (lwrk)
	//
	// lwrk       (input) integer
	//             Specifies the size of workspace.
	//
	// ..
	// .. local variables ..
	int i;
	integer_ d__1, d__2;
	//
	// ..
	// .. executable statements ..
	if (lwrk > 3 * nd) {
		*ierr = 0;
	} else {
		*ierr = -121;
		return;
	}
	memcpy(lambda, alpha, nd * sizeof(double));
	memcpy(wrk, &beta[locked], (nd - locked) * sizeof(double));
	d__1 = (long) (nd - locked);
	dstqrb_(&d__1, &lambda[locked], wrk, &res[locked], &wrk[nd], &d__2);
	*ierr = (int) d__2;
	if (*ierr == 0) {
		memset(res, 0, locked * sizeof(double));
		for (i = locked; i < nd; i++) {
			res[i] = beta[nd - 1] * fabs(res[i]);
		}
	} else {
		*ierr = -122;
	}
	//
	// .. end of trl_get_eval ..
	//
}



void TRL::trl_convergence_test(int nd, double *lambda, double *res,
						  trl_info * info, double *wrk)
{
	//
	// Purpose:
	// ========
	// count the numer of wanted eigenvalues that have small residual norms
	// eigenvalues are assumed to be order from small to large
	//
	// Arguments;
	// ==========
	// nd       (input) integer
	//           On entry, specifies the size of lambda.
	//
	// lambda   (input) double precision array (nd)
	//           On entry, contains the Ritz values.
	//
	// res      (input) double precision array (nd)
	//           On entry, contains the residual norm of the Ritz values.
	//
	// info     (input) Pointer to trl_info structure
	//           On entry, points to the data structure to store the current information about
	//           the eigenvalue problem and the progress of TRLAN.
	//
	// wrk      (workspace) double precision array (2*nd)
	//
	// ..
	// .. local scalars ..
	double bnd;
	int i, j, ncl, ncr;
	//
	// copy lambda and res to wrk, sort them in ascending order of lambda
	if (info->lohi < -1) {
		dsort2s(nd, info->ref, lambda, res);
	}
	memcpy(&wrk[nd], lambda, nd * sizeof(double));
	for (i = 0; i < nd; i++) {
		//wrk[nd+i] = lambda[i];
		wrk[i] = fabs(res[i]);
	}
	// sort in the ascending order of wrk[nd:2*nd]=lambda
	if (info->lohi == -2) {
		// around ref
		dsort2s(nd, info->ref, &wrk[nd], wrk);
	} else if (info->lohi == -3) {
		// larger than ref
		dsort2su_(nd, info->ref, &wrk[nd], wrk);
	} else if (info->lohi == -4) {
		// smaller than ref
		dsort2sd(nd, info->ref, &wrk[nd], wrk);
	} else {
		dsort2(nd, &wrk[nd], wrk);
	}
	//
	// determine the convergence rate of the previous target
	if (info->tmv > 0 && info->matvec > info->tmv) {
		j = 0;
		bnd = fabs(lambda[j] - info->trgt);
		for (i = 0; i < nd; i++) {
			if (fabs(lambda[i] - info->trgt) < bnd) {
				bnd = fabs(lambda[i] - info->trgt);
				j = i;
			}
		}
		if (info->tres > res[j]) {
			bnd = res[j] / info->tres;
			if (bnd > 0.0) {
				info->crat =
					exp(log(bnd) / (double) (info->matvec - info->tmv));
				info->cfac = bnd;
			} else {
				info->crat = 1.0;
				info->cfac = bnd;
			}
		} else {
			info->crat = 1.0;
			info->cfac = bnd;
		}
	}
	//
	// find out who has converged at the lower end of the spectrum
	info->anrm =
		max(info->anrm, max(fabs(wrk[nd + 1]), fabs(wrk[nd + nd - 1])));
	bnd = DBL_MIN + info->tol * info->anrm;
	ncl = 0;
	ncr = nd;
	if (info->lohi <= 0) {
		ncl = nd - 1;
		i = 0;
		while (i < nd) {
			if (wrk[i] < bnd) {
				if (info->lohi == -3 && wrk[i + nd] < info->ref) {
					ncl = i - 1;
					i = nd;
				} else if (info->lohi == -4 && wrk[i + nd] > info->ref) {
					ncl = i - 1;
					i = nd;
				} else {
					i++;
				}
			} else {
				ncl = i - 1;
				i = nd;
			}
		}
	}
	// find out who has converged at the high end of the spectrum
	if (info->lohi >= 0) {
		ncr = 0;
		i = nd - 1;
		while (i >= 0) {
			if (wrk[i] < bnd) {
				i--;
			} else {
				ncr = i + 1;
				i = -1;
			}
		}
	}
	// determine the number of wanted eigenvalues that have converged
	// compute the next target
	// ncl = index of wrk corresponding to the smallest eig converged.
	// ncr = index of wrk corresponding to the laragest eig converged.
	info->tmv = info->matvec;
	info->ptres = info->trgt;
	if (info->lohi < 0) {
		info->nec = ncl + 1;
		info->trgt = wrk[nd + min(nd - 1, ncl + 1)];
		info->tres = wrk[min(nd - 1, ncl + 1)];
	} else if (info->lohi > 0) {
		info->nec = nd - ncr;
		info->trgt = wrk[nd + max(0, ncr - 1)];
		info->tres = wrk[max(0, ncr - 1)];
	} else {
		if (ncr <= ncl) {
			ncl = nd / 2;
			ncr = ncl + 1;
			info->trgt = wrk[nd + (nd + 1) / 2 - 1];
			info->tres = wrk[(nd + 1) / 2 - 1];
		} else if (wrk[ncl + 1] <= wrk[ncr - 1]) {
			info->trgt = wrk[nd + ncl + 1];
			info->tres = wrk[ncl + 1];
		} else {
			info->trgt = wrk[nd + ncr - 1];
			info->tres = wrk[ncr - 1];
		}
		info->nec = ncl + nd - ncr + 1;
		//for( i=ncl; i<ncr-1; i++ ) {
		for (i = ncl + 1; i < ncr; i++) {
			if (wrk[i] < bnd)
				info->nec = info->nec + 1;
		}
	}
	//
	// .. end of trl_convergence_test ..
	//
}


int TRL::trl_check_dgen(trl_info * info, int jnd, double *lambda,
				   double *res)
{
	/*
	if( info->nec >= info->ned ) {
	if( (info->lohi == -1 || info->lohi == 0 ) &&
	(lambda[info->nec]  +res[info->nec]   > lambda[info->nec+1] &&
	lambda[info->nec+1]-res[info->nec+1] < lambda[info->nec]) ) {
	return 1;
	}
	if( (info->lohi == 1 || info->lohi == 0 ) &&
	(lambda[jnd-info->nec+1]-res[jnd-info->nec+1] < lambda[jnd-info->nec] &&
	lambda[jnd-info->nec]  +res[jnd-info->nec]   > lambda[jnd-info->nec+1]) ) {
	return 1;
	}
	}
	*/
	return -1;
}




void TRL::trl_get_tvec(int nd, double *alpha, double *beta, int irot, int nrot,
				  double *rot, int nlam, double *lambda, double *yy,
				  int *iwrk, double *wrk, int lwrk, int *ierr)
{
	//
	// Purpose:
	// ========
	// generating eigenvectors of the projected eigenvalue problem acorrding to the given 
	// Ritz values using LAPACK routine DSTEIN (inverse iterations).
	//
	// Arguments;
	// ==========
	// nd        (input) integer
	//            On entry, specifies the size of alpha and beta.
	//
	// alpha     (input) doubel precision array (nd)
	//            On entry, contains the alpha values.
	//
	// beta      (input) double precision array (nd)
	//            On entry, contains the beta values.
	//
	// irot      (input) integer
	//            On entry, specifies the starting column index of yy to apply the rotation.
	//
	// nrot      (input) integer
	//            On entry, specifies the ending column index of yy to apply the rotation.
	//
	// rot       (input) double precision array (nrot, nrot)
	//            On entry, contains the rotation matrix.
	//
	// nlam      (input) integer
	//            On entry, specifies the size of lambda.
	//
	// lambda    (input) double precision array (nlam)
	//            On entry, contains the Ritz values.
	//
	// yy        (output) double precision array (nd,nlam)
	//            On exit, contains the eigenvector of the tri-diagonal matrix.
	//
	// iwrk      (workspace) integer_ array (4nd)
	// wrk       (workspace) double precision (lwrk>=5nd)
	// lwrk      (input) integer
	//            specifies the size of workspace.
	//
	// ierr      (output) integer
	//            Error from Lapack subroutines.
	//
	//
	// local variables
	char notrans = 'N';
	int c__1 = 1;
	double zero = 0.0, one = 1.0;
	int i, j, k, ncol, ii, ioff;
	//
	// conventional external subprograms
	//extern int dstein_();
	//
	if (nlam <= 0) {
		*ierr = 0;
		return;
	}
	if (lwrk > 5 * nd) {
		*ierr = 0;
	} else {
		*ierr = -131;
		return;
	}
	//
	// set up IBLOCK and ISPLIT for calling dstein
	for (i = 0; i < nd; i++) {
		iwrk[i] = 1;
		iwrk[nd + i] = nd;
	}
	TRL::dstein_(&nd, alpha, beta, &nlam, lambda, iwrk, &iwrk[nd], yy, &nd, wrk,
		&iwrk[2 * nd], &iwrk[3 * nd], ierr);

	if (*ierr != 0) {
		printf("TRL_GET_TVEC: dstein failed with error code %d\n", *ierr);
		*ierr = -132;
		return;
	}
	//
	// apply the rotations to the IROT+1:IROT+NROT rows of YY
	// generates results 'NCOL' columns at a time
	if (nrot > 1) {
		ncol = lwrk / nrot;
		for (i = 1; i <= nlam; i += ncol) {
			j = min(nlam, i + ncol - 1);
			k = j - i + 1;
			if (k > 1) {
				trl_dgemm(&notrans, &notrans, nrot, k, nrot, one, rot,
					nrot, &yy[(i - 1) * nd + irot], nd, zero, wrk,
					nrot);
				for (ii = i - 1; ii < j; ii++) {
					ioff = (ii - i + 1) * nrot;
					memcpy(&yy[ii * nd + irot], &wrk[ioff],
						nrot * sizeof(double));
				}
			} else {
				trl_dgemv(&notrans, nrot, nrot, one, rot, nrot,
					&yy[(i - 1) * nd + irot], c__1, zero, wrk, c__1);
				memcpy(&yy[(i - 1) * nd + irot], wrk,
					nrot * sizeof(double));
			}
		}
	}
	//
	// .. end of trl_get_tvec ..
	//
}

////
void TRL::trl_get_tvec_a(int nd, int kept, double *alpha, double *beta,
					int nlam, double *lambda, double *yy, double *wrk,
					int lwrk, int *iwrk, int *ierr)
{
	//
	// Purpose
	// =======
	// compute all eigenvalues and eigenvectors of the projected matrix
	// use LAPACK routine DSYEV
	// The eigenvectors corresponding to lambda(1:nlam) are placed at the
	// first nlam*nd locations of yy on exit.
	//
	// Arguments;
	// ==========
	// nd        (input) integer
	//            On entry, specifies the size of alpha and beta.
	//
	// kept      (input) integer
	//            On entry, specifies the number of Ritz values kept.
	//
	// alpha     (input) doubel precision array (nd)
	//            On entry, contains the alpha values.
	//
	// beta      (input) double precision array (nd)
	//            On entry, contains the beta values.
	//
	// nlam      (input) integer
	//            On entry, specifies the size of lambda.
	//
	// lambda    (input) double precision array (nlam)
	//            On entry, contains the Ritz values.
	//
	// yy        (output) double precision array (nd,nlam)
	//            On exit, contains the eigenvector of the arrow-head matrix.
	//
	// iwrk      (workspace) integer_ array (nd)
	// wrk       (workspace) double precision (lwrk)
	// lwrk      (input) integer
	//            specifies the size of workspace.
	//
	// ierr      (output) integer
	//            Error from Lapack subroutines.
	//
	// ..
	// .. CLAPACK subroutines ..
	//extern void dsyev_();
	//
	// ..
	// .. local parameters ..
	char job = 'V';
	char upl = 'U';
	//
	// local variables
	int i, j, i2, j2, ii;
	double tmp;
	//
	// ..
	// .. executables statements ..
	//
	// fill yy with the projection matrix, then call DSYEV

	if (nlam <= 0) {
		*ierr = 0;
		return;
	}
	if (lwrk >= nd + nd + nd) {
		*ierr = 0;
	} else {
		*ierr = -141;
		return;
	}
	memset(yy, 0, (nd * nd) * sizeof(double));
	j = 0;
	for (i = 0; i < nd; i++) {
		yy[j] = alpha[i];
		j += (nd + 1);
	}
	if (kept > 0) {
		memcpy(&yy[kept * nd], beta, kept * sizeof(double));
	}
	for (i = kept; i < nd - 1; i++) {
		yy[(i + 1) * (nd + 1) - 1] = beta[i];
	}

	TRL::dsyev_(&job, &upl, &nd, yy, &nd, alpha, wrk, &lwrk, ierr);

	if (*ierr != 0) {
		printf("Error from dsyev: %d.\n", *ierr);
		*ierr = -142;
		return;
	}
	if (nlam >= nd)
		return;
	//
	// reorder the eigenvectors
	// both lambda(1:kept) and alpha are in ascending order
	//
	tmp =
		max(alpha[nd - 1] - alpha[0],
		max(fabs(alpha[nd - 1]), fabs(alpha[0])));
	tmp = DBL_EPSILON * tmp * nd;
	j = 0;
	i = 0;
	while (i < nlam) {
		// move j so that alpha(j) is within tmp distance away
		ii = j;
		j = nd - 1;
		while (ii < nd) {
			if (alpha[ii] < lambda[i] - tmp) {
				ii++;
			} else {
				j = ii;
				ii = nd;
			}
		}
		if (alpha[j] > lambda[i] + tmp) {
			*ierr = -143;
			return;
		}
		// identify the group size in lambda
		ii = i + 1;
		i2 = nlam - 1;
		while (ii < nlam) {
			if (lambda[ii] <= lambda[i] + tmp) {
				ii++;
			} else {
				i2 = ii - 1;
				ii = nd;
			}
		}
		// identify the group size in alpha
		ii = j + 1;
		j2 = nd - 1;
		while (ii < nd) {
			if (alpha[ii] <= lambda[i] + tmp) {
				ii++;
			} else {
				j2 = ii - 1;
				ii = nd;
			}
		}
		// assign the index values
		if (j2 == j && i2 == i) {
			iwrk[i] = j;
		} else if (j2 - j == i2 - i) {
			//iwrk(i:i2) = (/(ii, ii=j, j2)/)
			for (ii = i; ii < i2; ii++) {
				iwrk[ii] = j + ii - i;
			}
		} else if (j2 - j > i2 - i) {
			//j2 = j + i2 - i;
			//iwrk(i:i2) = (/(ii, ii=j, j2)/)
			for (ii = i; ii < i2; ii++) {
				iwrk[ii] = j + ii - i;
			}
		} else if (j2 < nd) {
			i2 = i + j2 - j;
			//iwrk(i:i2) = (/(ii, ii=j, j2)/)
			for (ii = i; ii < i2; ii++) {
				iwrk[ii] = j + ii - i;
			}
		} else {
			*ierr = -144;
			return;
		}
		i = i2 + 1;
		j = j2 + 1;
	}
	// perform the actual copying
	for (i = 0; i < nlam; i++) {
		// for( i=1; i<=nlam; i++ ) {
		j = iwrk[i];
		if (j > i) {
			alpha[i] = alpha[j];
			memcpy(&yy[i * nd], &yy[j * nd], nd * sizeof(double));
		}
	}
	//
	// .. end of trl_get_tvec_a ..
	//
}

////
void TRL::trl_set_locking(int jnd, int nlam, double *lambda, double *res,
					 double *yy, int anrm, int *locked)
{
	//
	// Purpose
	// =======
	// Move the Ritz pairs with extremely small residual norms to the front of the 
	// arrays so that locking can be performed cleanly.
	//
	// Arguments
	// =========

	//  double: lambda(nlam), res(nlam), yy(jnd*nlam)
#define small(tmp,eps) (fabs(tmp) >= (eps) ? (eps)*fabs(tmp) : (eps)*(eps)*(anrm))
	//
	// ..
	// .. local parameters ..
	double zero = 0.0;
	//
	// ..
	// .. local variables ..
	int i, j, ii, ioff, ti, tj;
	double tmp, eps;
	//
	// ..
	// .. Executable statements ..
	eps = DBL_EPSILON;
	i = 0;
	j = nlam - 1;

	//ti = (fabs(1) <= small(2, 1));
	ti = (fabs(res[i]) <= small(lambda[i], eps));
	tj = (fabs(res[j]) <= small(lambda[j], eps));
	while (i < j) {
		if (ti != 0) {
			// res[i] is very small, so lock ith lambda
			res[i] = zero;
			i = i + 1;
			if (i <= j) {
				ti = (fabs(res[i]) <= small(lambda[i], eps));
			} else {
				ti = 0;
			}
		} else {
			if (tj != 0) {
				// res[i] (small ones) is still large,
				// but res[j] (big ones) is very small,
				// so swap res[j] with res[i], and lock res[i].
				// swap the eigenvectors accordingly.
				tmp = lambda[i];
				lambda[i] = lambda[j];
				lambda[j] = tmp;
				res[j] = res[i];
				res[i] = zero;
				ioff = (j - i) * jnd;
				for (ii = (i + 1) * jnd - jnd; ii < (i + 1) * jnd; ii++) {
					tmp = yy[ii];
					yy[ii] = yy[ii + ioff];
					yy[ii + ioff] = tmp;
				}
				i++;
				if (i <= j) {
					ti = (fabs(res[i]) <= small(lambda[i], eps));
				} else {
					ti = 0;
				}
			}
			j--;
			if (j > i) {
				tj = (fabs(res[j]) <= small(lambda[j], eps));
			} else {
				tj = 0;
			}
		}
	}
	if (ti != 0) {
		*locked = i + 1;
	} else {
		*locked = i;
	}
	//
	// .. end of trl_set_locking ..
	//
}

////
void TRL::trl_ritz_vectors(int nrow, int lck, int ny, double *yy, int ldy,
					  double *vec1, int ld1, int m1, double *vec2,
					  int ld2, int m2, double *wrk, int lwrk)
{
	//
	// Purpose
	// =======
	// compute the Ritz vectors from the basis vectors and the eigenvectors of the projected system
	// the basis vectors may be stored in two separete arrays the result need to be stored back in them
	// lwrk should be no less than ny (lwrk>=ny) ***NOT checked inside***
	//
	// Arguments
	// =========
	// nrow   (input) Integer
	//         On entry, specifies the number of rows in eigenvectors.
	//
	// lck    (input) Integer
	//         On entry, specifies the number of Ritz values locked.
	//
	// ny     (input) Integer
	//         On entry, specifies the number of columns in yy.
	//
	// yy     (input) double precision array (ldy,ny)
	//         On entry, contains the eigenvector of the "tri-diagonal" matrix.
	//
	// ldy    (input) Integer
	//         On entry. specify the leading dimention of yy.
	//
	// vec1   (input) double precision array (ld1,m1)
	//         On entry, contains the first part of Lanczos basis.
	//
	// m1     (input) Integer
	//         On entry, specifies the number of Lanczos basis stored in vec1.
	//
	// ld1    (input) Integer
	//         On entry, specifies the leading dimention of vec1.
	//
	// vec2   (input) double precision array (ld2,m2)
	//         On entry, contains the second part of Lanczos basis.
	//
	// m2     (input) Integer
	//         On entry, specifies the number of Lanczos basis stored in vec2.
	//
	// ld2    (input) Integer
	//         On entry, specifies the leading dimention of vec2.
	//
	// wrk    (workspace) double precision array (lwrk)
	// yy2    (workspace) double precision array (ldy,ny)
	// lwrk   (input)
	//         Specifies the size of the workspace.
	//
	// ..
	// .. local parameters ..
	char notrans = 'N';
	double zero = 0.0, one = 1.0;
	int c__1 = 1;
	//
	// .. local variables ..
	int i, j, k, stride, ii, jl1, jl2, il1, il2, kv1;
	//
	// ..
	// .. executable statements ..
	// vec1*yy and vec2*yy where vec1 and vec2 are kept-locked
	// m1 number of kept in vec1
	// m2 number of kept in vec2
	if (lck <= m1) {
		// all locked are in vec1
		il1 = lck + 1;
		jl1 = m1 - lck;
		il2 = 1;
		jl2 = m2;
	} else {
		// all kept in vec1 are locked
		il1 = m1 + 1;
		jl1 = 0;
		il2 = lck - m1 + 1;
		jl2 = m1 + m2 - lck;
	}

	if (jl1 == 0 && jl2 == 0)
		return;

	kv1 = min(m1 - il1 + 1, ny);
	memset(wrk, 0, lwrk * sizeof(double));
	if (ny > 1) {
		stride = lwrk / ny;
		for (i = 0; i < nrow; i += stride) {
			j = min(nrow - 1, i + stride - 1);
			k = j - i + 1;

			if (jl1 > 0) {
				// compute wrk = vec1(i:j,:)*yy
				// (Note the leading dimension of vec1 is ld1. This effectively shift 
				// the vec1(i:j) to the top of the matrix.
				trl_dgemm(&notrans, &notrans, k, ny, jl1, one,
					&vec1[(il1 - 1) * ld1 + i], ld1, yy, ldy, zero,
					wrk, k);
			} else {
				memset(wrk, 0, lwrk * sizeof(double));
			}

			if (jl2 > 0) {
				trl_dgemm(&notrans, &notrans, k, ny, jl2, one,
					&vec2[(il2 - 1) * ld2 + i], ld2, &yy[jl1], ldy,
					one, wrk, k);
			}

			for (ii = 0; ii <= kv1 - 1; ii++) {
				memcpy(&vec1[(ii + il1 - 1) * ld1 + i], &wrk[ii * k],
					k * sizeof(double));
			}

			for (ii = 0; ii <= (ny - kv1 - 1); ii++) {
				memcpy(&vec2[(ii + il2 - 1) * ld2 + i],
					&wrk[(kv1 + ii) * k], k * sizeof(double));
			}

		}
	} else if (ny == 1) {
		stride = lwrk;
		for (i = 0; i < nrow; i += stride) {
			j = min(nrow - 1, i + stride - 1);
			k = j - i + 1;
			if (jl1 > 0) {
				trl_dgemv(&notrans, k, jl1, one,
					&vec1[(il1 - 1) * ld1 + i], ld1, yy, c__1, zero,
					wrk, c__1);
				if (jl2 > 0) {
					trl_dgemv(&notrans, k, jl2, one,
						&vec2[(il2 - 1) * ld2 + i], ld2, &yy[jl1],
						c__1, one, wrk, c__1);
				}
			} else {
				trl_dgemv(&notrans, k, jl2, one,
					&vec2[(il2 - 1) * ld2 + i], ld2, &yy[jl1], c__1,
					zero, wrk, c__1);
			}
			if (kv1 > 0) {
				memcpy(&vec1[(il1 - 1) * ld1 + i], wrk,
					k * sizeof(double));
			} else {
				memcpy(&vec2[(il2 - 1) * ld2 + i], wrk,
					k * sizeof(double));
			}
		}
	}
	//
	// .. end of trl_ritzs_vectors_ ..
	//
}


void TRL::trl_sort_eig(int nd, int lohi, int nec, double ref, double *lambda,
				  double *res)
{
	//
	// Purpose:
	// ========
	// sort the eigenvalues so that the wanted eigenvalues are ouputed to the user in
	// front of the arrays. the final Ritz values are in ascending order so that DSTEIN 
	// can be used to compute the eigenvectors
	//
	// Arguments;
	// ==========
	// nd       (input) integer
	//           On entry, specifies the size of lambda.
	//
	// lohi     (input) integer
	//           On entry, specifies which eigenvalues are desired.
	//
	// nec      (input) integer
	//           On entry, specifies how many Ritz values have been converged.
	//
	// lambda   (input) double precision array (nd)
	//           On entry, contains the Ritz values.
	//
	// res      (input) double precision array (nd)
	//           On entry, contains the residual norm of the Ritz values.
	//
	// ..
	// .. local scalars ..
	int i, j;
	//
	// ..
	// .. executable statements ..
	if (lohi == 0) {
		// sort the eigenvalues according to their absolute residual values
		// to get those converged first
		dsort2a(nd, res, lambda);
		// sort the first nec eigenvalue in the order of lambda
		dsort2(nec, lambda, res);
	} else {
		// sort the eigenvalues and residual norms in ascending order of the
		// eigenvalues
		if (lohi == -2) {
			// around ref
			dsort2s(nd, ref, lambda, res);
			dsort2(nec, lambda, res);
		} else if (lohi == -3) {
			// larger than ref
			dsort2su_(nd, ref, lambda, res);
			dsort2(nec, lambda, res);
		} else if (lohi == -4) {
			// smaller than ref
			dsort2sd(nd, ref, lambda, res);
			dsort2(nec, lambda, res);
		} else {
			dsort2(nd, lambda, res);
			if (lohi > 0) {
				// move the largest ones to the front (still ascending order)
				j = nd - nec;
				for (i = 0; i < nec; i++) {
					res[i] = res[j];
					lambda[i] = lambda[j];
					j++;
				}
			}
		}
	}
	//
	// .. end of trl_sort_eig ..
	//
}


void TRL::print_final_state(trl_info * info, char *title, int nrow, int mev,
					   double *eval, double *beta, double *evec,
					   double *yy, int kept, int jml)
{
	/*
	// Purpose
	// =======
	// print the final state
	//
	// Arguments
	// =========
	// info    (input) Pointer to structure trl_info_
	//          On entry, points to the current TRL_INFO.
	//
	// title   (workspace) String of length (STRING_LEN)
	//          On entry, provides space to store the title of the information to 
	//          print out.
	//
	// nrow    (input) Integer
	//          On entry, specifies the number of rows in the eigenvectors.
	//
	// mev     (input) Integer
	//          On entry, specifies the maximum number of eigenvalues allowed.
	//
	// eval    (input) Double precision array of dimension (mev)
	//          On entry, contains the eigenvalues computed.
	//
	// beta    (input) Double precision array of dimension (info->maxlan)
	//          On entry, contains the value of beta computed.
	//
	// evec    (input) Double precision array of dimension (nrow,mev)
	//          On entry, contains the eigenvectors computed.
	//
	// yy      (input) Double precision array of dimension (nrow,jml)
	//          On entry, contains the litz vectors computed at the last restart.
	//
	// kept    (input) Integer
	//          On entry, specifies the number of lanczos vectors kept at the last 
	//          restart.
	//
	// jml     (input) Integer
	//          On entry, specifies the number of new lanczos vectors computed at 
	//          the last restart.
	//
	// ..
	// .. local scalars ..
	*/
	int j1;
	/*
	// ..
	// .. executable statements ..
	*/
	strcpy(title, "Final eigenvalues  (in ascending order)..");
	trl_print_real(info, title, kept, eval, 1);
	if (info->verbose > 4) {
		strcpy(title, "Final residual norms..");
		trl_print_real(info, title, kept, beta, 1);
	}
	if (info->verbose > 8) {
		for (j1 = 0; j1 < min(kept, info->verbose); j1++) {
			sprintf(title, "Eigenvector %d of Q''AQ ..", j1);
			trl_print_real(info, title, jml, &yy[j1 * jml], 1);
		}
	}
	if (info->verbose > 10) {
		int j1n = min(nrow, info->verbose);
		for (j1 = 0; j1 < min(kept, mev); j1++) {
			sprintf(title, "Ritz vector %d (1:%d) ..", j1, j1n);
			trl_print_real(info, title, j1n, &evec[j1 * nrow], 1);
		}
	}
	/*
	// .. end of print_final_state ..
	*/
}

/*
// Purpose
// =======
// Output a check point.
//
// Arguments
// =========
// info     (input) Pointer to structure trl_info_
//           On entry, points to the current TRL_INFO.
//
// title    (input) String of length TITLE_LEN.
//           On entry, provides space to store the title of the information to 
//           print out.
//
// nrow     (input) Integer
//           On entry, specifies the number of rows in the eigenvectors.
//
// alpha    (input) Double precision array of length (info->maxlan)
//           On entry, contains the values of alpha computed.
//
// beta     (input) Double precision array of length (info->maxlan)
//           On entry, contains the value of beta computed.
//
// evec     (input) Double precision array of length (nrow,mev)
//           On entry, contains the eigenvectors computed.
//
// base     (input) Double precision array of length (ldb,nbase)
//           On entry, contains the lanczos vectors, that did not fit in evec.
//
// lde      (input) Integer
//           On entry, specifies the leading dimension of evec.
//
// j1n      (input) Integer
//           On entry, specifies the column index of evec, that stores the current 
//           lanczos vector.
//
// jnd      (input) Integer
//           On entry, specifies the number of lanczos vector computed so far.
//
// ldb      (input) Integer
//           On entry, specifies the leading dimension of base.
//
// j2n      (input) Integer
//           On entry, specifies the column index of base, the stores the current 
//           lanczos vector.
//
*/
void TRL::write_checkpoint(trl_info * info, char *title, int nrow,
					  double *alpha, double *beta, double *evec, int lde,
					  double *base, int ldb, int j1n, int jnd, int j2n)
{
	/*
	// ..
	// .. local variables ..
	*/
	int ii, c1, c2;
	/*
	// ..
	// .. executable statements ..
	*/
	trl_pe_filename(138, title, info->cpfile, info->my_pe, info->npes);
	c1 = clock();
	ii = trl_write_checkpoint(title, nrow, alpha, beta, evec, lde, j1n,
		base, ldb, j2n);
	c2 = clock();
	if (c2 > c1) {
		info->clk_out = info->clk_out + (c2 - c1);
	} else {
		info->clk_out = info->clk_out + ((info->clk_max - c1) + c2);
	}
	info->wrds_out = info->wrds_out + jnd * (nrow + nrow + 2) + nrow + 2;
	info->stat = trl_sync_flag_(info->mpicom, ii);
	/*
	//  .. end of print final_state_ ..
	*/
}


void TRL::print_restart_state(trl_info * info, char *title, int nrow,
						 int mev, double *alpha, double *beta,
						 double *betrot, double *evec, double *yy,
						 int kept, int locked, int *iwrk, double *wrk2,
						 int i2, int jml)
{
	/*
	// Purpose
	// =======
	// Print the current solution status to the log file.
	//
	// Arguments
	// =========
	// info     (input) Pointer to structure trl_info_
	//           On entry, points to the current TRL_INFO.
	//
	// title    (workspace) String of length (STRING_LEN)
	//           On entry, provides space to store title to print out.
	//
	// nrow     (input) Integer
	//           On entry, specifies the number of rows in the lanczos vectors.
	//
	// mev      (input) Integer
	//           On entry, specifies the maximum number of eigenvalues allowed.
	//
	// alpha    (input) Double precision array of dimension (info->maxlan)
	//           On entry, contains the value of alphas computed so far.
	//
	// beta     (input) Double precision array of dimension (info->maxlan)
	//           On entry, contains the value of beta computed so far.
	//
	// betrot   (input) Double precision array of dimension (info->maxlan)
	//           On entry, contains the value of beta rotated.
	//
	// evec     (input) Double precision array of dimension (nrow,mev)
	//           On entry, contains the eigenvectors computed.
	//
	// yy       (input) Double precision array of dimension (nrow,jml)
	//           On entry, contains the litz vectors of the tridiagonal matrix
	//           computed after the previous restart.
	//
	// kept     (input) Integer
	//           On entry, specifies the number of lanczos vector kept at the restart.
	//
	// locked   (input) Integer
	//           On entry, specifies the number of eigenvalues converged so far.
	//
	// iwrk     (workspace) Integer array of dimension (4*maxlan)
	//           Integer workspace used to..
	//
	// wrk2     (workspace) Double precision array of dimension
	//           Double precision workspace used to
	//
	// i2       (input) Integer
	//           On entry, specifies
	//
	// jml      (input) Integer
	//           On entry, specifies the number of litz vectors computed for
	//           the current restart.
	//
	// ..
	// .. local parameters ..
	*/
	long c__1 = 1;
	/*
	// ..
	// .. local scalars ..
	*/
	int i, j1, j2;
	/*
	// ..
	// .. executable statements ..
	*/
	iwrk[0] = kept + locked;
	iwrk[1] = locked + i2;
	strcpy(title, "Number of saved and locked Ritz pairs ..");
	trl_print_int(info, title, 2, iwrk, 1);
	if (info->verbose > 2) {
		if (iwrk[1] == 0) {
			strcpy(title, "Ritz values saved (ascending ordered) ..");
		} else {
			strcpy(title, "Ritz values saved (may not be ordered) ..");
		}
		trl_print_real(info, title, kept + locked, alpha, 1);
		strcpy(title, "Residual norms of the saved Ritz pairs ..");
		for (i = 0; i < (kept + locked); i++) {
			betrot[i] = fabs(beta[i]);
		}
		trl_print_real(info, title, kept + locked, betrot, 1);
	}
	if (info->verbose > 7) {
		for (j1 = 0; j1 < min(kept, info->verbose); j1++) {
			for (j2 = 0; j2 <= j1; j2++) {
				wrk2[j2] =
					trl_ddot(jml, &yy[j2 * jml], c__1, &yy[j1 * jml],
					c__1);
			}
			wrk2[j1] = wrk2[j1] - 1;
			sprintf(title, "Orthogonality level of y(%d) ..", j1 + 1);
			trl_print_real(info, title, j1 + 1, wrk2, 1);
		}
	}
	if (info->verbose > 10) {
		for (j1 = 0; min(kept, info->verbose); j1++) {
			sprintf(title, "eigenvector %d of Q'AQ ..", j1);
			trl_print_real(info, title, jml, &yy[(j1 - 1) * jml], 1);
		}
	}
	if (info->verbose > 10) {
		int j1n = min(nrow, info->verbose);
		for (j1 = 0; j1 < min(kept + locked, mev); j1++) {
			sprintf(title, "Ritz vector %d (1:%d) ..", j1, j1n);
			trl_print_real(info, title, j1n, &evec[j1 * nrow], 1);
		}
	}
	/*
	//  .. end of print_restart_state ..
	*/
}




