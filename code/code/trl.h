#ifndef _TRL_H
#define _TRL_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <time.h>
#include <math.h>
#include <vector>

#include <iostream>
#include <fstream>
#include "parset.h"

using namespace std;

class Dmat;
class Dvec;

//#define DEBUG_WEI_TRL

#define STD_DATA

/* The maximum number of string allowed for a information titile. */
#define STRING_LEN 132

#ifndef INT_MAX
#define INT_MAX       2147483647    /* maximum (signed) int value */
#endif

#ifdef STD_DATA
typedef int integer_;
typedef double doublereal_;
typedef float trl_real_;
typedef int logical_;

typedef int ftnlen_;
typedef int flag_;
typedef int ftnint_;

typedef struct
{	
	flag_ cierr;
	ftnint_ ciunit;
	flag_ ciend;
	char *cifmt;
	ftnint_ cirec;
} cilist_;

typedef struct {
	double r, i;
} trl_dcomplex_;
#define min(a,b) ( a < b ? a : b )
#define max(a,b) ( a > b ? a : b )
#define TRUE_ (1)
#define FALSE_ (0)
#endif

#define nint(x) (int)((x)+0.5)

/**
Prototype matrix-vector multiplication function.  Defined to be easier
to use with C/C++ code.

@arg nrow The number of rows locally on this processor.

@arg ncol The number of columns in vectors x and y.

@arg x Pointer to the elements of input vectors x.  It is assumed to
have [ldx * ncol] elements, with the ith column starting at element
ldx * i.

@arg ldx The leading dimension of x, the ith column of x starts at
element ldx * i.
Note that ldx must not be less than nrow.

@arg y Pointer to the elements of output vectors y.  It is assumed to
have [ldy * ncol] elements, with the ith column starting at element
ldy * i.

@arg ldy The leadying dimension of y.  The ith column of y starts at
position ldy * i.
Note that ldy must not be less than nrow.

@arg mvparam The extra parameter to be passed to the matrix-vector
multiplication through trl_info.  This parameter mvparam is used
exclusively in the matrix-vector multiplication function and nowhere
else in TRLan.
*/

/*
typedef void (*trl_matvec) (const int nrow, const int ncol,
							const double *x, const int ldx,
							double *y, const int ldy, void *mvparam);
*/

typedef void (*trl_matvec) (const int nrow, const int ncol,
							double *x, const int ldx,
							double *y, const int ldy, void *mvparam);

/**
The data structure to store the current information about
the eigenvalue problem and the progress of TRLAN.
*/
typedef struct strct_trl_info {
	int stat;			/* status  (error code) of TRLAN */

	/* specification of the eigenvalue problem */
	/** Which end of spectrum to compute.
	- lohi < 0 --> the smallest eigenvalues
	- lohi = 0 --> whichever converge first
	- lohi > 0 --> the largest eigenvalues
	*/
	int lohi;
	/** Number of eigenpairs wanted.  */
	int ned;
	/** Number of eigenpairs converged.  if the user has nec correct
	eigenvectors on input, then they are expected to be stored at the
	beginning of the eigenvector array.  */
	int nec;
	/* Convergence tolerance.  An eigenpair is declared converged if its
	residual norm is less than tol*||OP||.  */
	double tol;

	/* specification of resource allowed to use by TRLAN */
	/** The MPI communicator.  */
	int mpicom;
	/** The maximum basis size to be used.  */
	int maxlan;
	/** The actual basis size currently. This value may be smaller than
	maxlan.  It is set during restarting.  */
	int klan;
	/** The maximum number of MATVEC allowed. Note one MATVEC == one
	application of the operator on one vector.  */
	int maxmv;

	/** The restarting scheme to use.  */
	int restart;
	/** The number of eigenvalues locked.  */
	int locked;
	/** Option for handling initial guesses:
	- <= 0, user did not provide initial guess, use    
	[1,1,..,1].
	-  = 1, user has supplied initial guess, will only 
	use the first one.
	-  > 1, restart with previous check-point file.      */
	int guess;

	/* some information about the progress and resouce comsumption    */
	/** The number of MATVEC used by TRLAN.                    */
	int matvec;
	/** The number of restart of the Lanczos iterations.       */
	int nloop;
	/** The number of full orthogonalization invoked.          */
	int north;
	/** The number of times a random element is introduced. Random elements
	are introduced when an invariant subspace is found, but the number
	of converged eigen-pairs is less than desired.  */
	int nrand;
	/** Floating-point operations count (EXCLUDING MATVEC). */
	int flop;
	/** FLOPS used in re-orthogonalization.                */
	int flop_h;
	/** FLOPS used in restarting.                          */
	int flop_r;
	double rflp;
	double rflp_h;
	double rflp_r;

	/* variables to store timing results */
	/** system clock rate (SYSTEM_CLOCK)                */
	clock_t clk_rate;
	/** Maximum counter value                           */
	clock_t clk_max;
	/** Total time spent in TRLAN (in clock ticks)      */
	clock_t clk_tot;
	/** Time in applying the operator (MATVEC)          */
	clock_t clk_op;
	/** Time in re-orthogonalization                    */
	clock_t clk_orth;
	/** Time in restarting the Lanczos iterations       */
	clock_t clk_res;
	/** The sum of clk_tot and tick_t is the actual time */
	double tick_t;
	double tick_o;
	double tick_h;
	double tick_r;
	/** Time spent in reading input data file           */
	int clk_in;
	/** Number of real(8) words read                    */
	int wrds_in;
	/** Time spent in writing output data file          */
	int clk_out;
	/** Number of real(8) words written to file         */
	int wrds_out;

	/** Norm of the operator used.  This is an estimate based on the
	largest absolute value of a Rayleigh quotient.  */
	double anrm;

	/** The PE number of current processor (start with 0).  */
	int my_pe;
	/** number of PEs in the group                      */
	int npes;
	/** Local problem size                              */
	int nloc;
	/** Global problem size                             */
	int ntot;

	/** How much inforation to output during the execution of TRLAN.  By
	default, it only print information related to fatal errors.  If
	verbose > 0, more diagnostic messages are printed.  */
	int verbose;
	/** Variables needed to measure convergence factor (crat).  The
	convergence rate of the restarted Lanczos algorithm is measured by
	the reduction of residual norm per MATVEC.  The residual norm of the
	target is used.  */
	double crat;
	/** The Ritz value that might convege next.  */
	double trgt;
	/** The residual norm of the target.  */
	double tres;
	/** MATVEC used when target and tres were recorded  */
	int tmv;
	double avgm;

	/* Stores some convergence history for restart scheme 9           */
	double cfac;
	double ptres;
	double predicted_crate;
	double mgap;
	double mgamma, gamma0;
	double old_target;
	int target_id;
	int old_locked;
	int k1, k2, k;
	double rfact;

	/* Store "shift" */
	double ref;

	/* Fortran I/O unit number to be used for          */
	/* debug output. Used if verbose > 0.              */
	int log_io;
	FILE *log_fp;
	/** Base of the file names used by TRLAN to store debug info if verbose
	> 0, the filenames are computed by appending 'PE#' to this base.  */
	char log_file[128];

	/** check-pointing parameters.
	when cpflag_ is greater than 0, the basis vectors will be written
	out roughly 'cpflag_' times.  For simplicitly, each PE writes its
	own portion of the basis vectors to a file with cpfile followed by
	the processor number.  The file is written as unformatted fortran
	files with the following content:
	@pre
	nrow, kb(basis size)
	alpha(1:kb)
	beta(1:kb)
	1st basis vector, 2nd basis vector, ..., last basis vector
	the residual vector
	@endpre
	*/
	int cpflag_, cpio;
	FILE *cpt_fp;
	char cpfile[128], oldcpf[128];
	/**
	This parameter is only used by the matrix-vector multiplication
	function when the conditional macro TRL_FORTRAN_COMPATIBLE is not
	defined.  It is used nowhere else in TRLan.
	*/
	void* mvparam;
} trl_info;

namespace TRL{
	static int seed=0;

	void smg(vector< vector<double> > mat, vector<double> & egval, vector< vector<double> > & egvec);

	void sort_EGrst(int comp, int ned, int nrow, double * eval, double *evec, 
		vector<double> & egval, vector< vector<double> > & egvec, bool allned=false);
	void mt_op(const int nrow, const int ncol, double *xin, const int ldx,
		double *yout, const int ldy, void* mvparam);
	void mt_op2(const int nrow, const int ncol, const double *xin, const int ldx,
		double *yout, const int ldy, void* mvparam);
	void mt_op3(const int nrow, const int ncol, double *xin, const int ldx,
		double *yout, const int ldy, void* mvparam);
	void mt_op4(const int nrow, const int ncol, double *xin, const int ldx,
		double *yout, const int ldy, void* mvparam);
	//this mt_op2 not stable, why?
	//void mt_op2(const int nrow, const int ncol, const double *xin, const int ldx,
	//	double *yout, const int ldy, void* mvparam);
	void eg_nuTRan(trl_matvec op, vector< vector<double> > & mt,
		vector<double> & egval, vector< vector<double> > & egvec, 
		double pres ,const int lohi, int mev, int maxlan, int maxmv);
	void eg_nuTRan(trl_matvec op, vector< vector<double> > & mt,
		vector<double> & egval, vector< vector<double> > & egvec, 
		int ned, double pres ,const int lohi, int mev, int maxlan, int maxmv);
	void eg_nuTRan(trl_matvec op, void* mt, int nrow,
	vector<double> & egval, vector< vector<double> > & egvec, 
		int ned, double pres ,const int lohi, int mev, int maxlan, int maxmv);
	//void eg_nuTRan3(trl_matvec op, vector< vector<double> > & mtt,
	//	vector<double> & egval, vector< vector<double> > & egvec, 
	//	int ned, double pres ,const int lohi, int mev, int maxlan, int maxmv);
	//void eg_nuTRan2(trl_matvec op, double * mt, int nrow,
	//	vector<double> & egval, vector< vector<double> > & egvec, 
	//	int ned, double pres ,const int lohi, int mev, int maxlan, int maxmv);
	void indexx(vector<double> &arr, vector<int> &indx);

	inline void SWAP(int &a, int &b)
	{
		int temp=a;
		a=b;
		b=temp;
	}

	void trlan(trl_matvec op,
		trl_info * info, int nrow, int mev, double *eval,
		double *evec, int lde, int lwrk, double *wrk );

	void trl_check_ritz(trl_matvec op,
		trl_info * info, int nrow, int ncol, double *rvec,
		int ldrvec, double *alpha, int *check, double *beta,
		double *eval, double *wrk, int lwrk);


	void trl_init_info(trl_info *info, int nrow, int mxlan, int lohi,
		int ned, double tol, int restart, int maxmv,
		int mpicom, void *mvparam=0);


	void trl_set_restart(trl_info * info, double rfact);


	void trl_set_debug(trl_info * info, int msglvl, char *filename);

	void trl_set_checkpoint(trl_info * info, int cpflag_, char *file);

	void trl_set_iguess(trl_info * info, int nec, int iguess, int ncps,
		char *cpf );

	void trl_print_info(trl_info * info, int mvflop);

	void trl_terse_info(trl_info * info, FILE * iou);

	void
		trl_ritz_projection(trl_matvec op, trl_info *info, int mev,
		double *evec, int lde, double *eres, double *wrk, int lwrk,
		double *base);

	void
		trl_rayleigh_quotients(trl_matvec op, trl_info * info, int ncol, double *evec,
		int lde, double *eres, double *base);

/////////////////////////////////////////////////////////////////

	int indchar(char *a, char b);
	int close_file(FILE * fp, int err1, int err2);

	void trl_open_logfile(trl_info * info);

	void trl_reopen_logfile(trl_info * info);

	void trl_close_logfile(trl_info * info);

	void trl_open_cptfile(trl_info * info);

	void trl_close_cptfile(trl_info * info);

	void trl_print_int(trl_info * info, char *title, int size_array,
		int *array, int inc);

	void trl_print_real(trl_info * info, char *title, int size_array,
		double *array, int inc);

	void trl_print_progress(trl_info * info);

	void trl_check_orth(trl_info * info, int nrow, double *v1, int ld1,
		int j1, double *v2, int ld2, int j2, double *wrk,
		int lwrk);

	void
		trl_check_recurrence(trl_matvec op,
		trl_info * info, int nrow, int ncol, double *v1,
		int ld1, int m1, double *v2, int ld2, int m2,
		int kept, double *alpha, double *beta, double *wrk,
		int lwrk);

	int trl_write_checkpoint(char *filename, int nrow, double *alpha,
		double *beta, double *evec, int lde, int me,
		double *base, int ldb, int nb);

	int trl_read_checkpoint(char *filename, int nrow, double *evec, int lde,
		int mev, int *j1, double *base, int ldb, int nbas,
		int *j2, int nalpha, double *alpha, int nbeta,
		double *beta);

	void trl_pe_filename(int nlen, char *filename, char *base, int my_rank,
		int npe);

	void trl_time_stamp(FILE * iou);

/////////////////////////////////////////////////////////////////////
	void
		trlanczos(trl_matvec op, trl_info * info, int nrow, int mev, double *eval,
		double *evec, int lde, double *base, int ldb, int nbas,
		double *wrk, int lwrk);

	void trl_shuffle_eig(int nd, int mnd, double *lambda, double *res,
		trl_info * info, int *kept, int locked);


	void add_clock_ticks(trl_info * info, clock_t *time, double *rtime,
		clock_t clk1);

	void print_alpha_beta(trl_info * info, char *title, int i,
		double *alpha, double *beta);

	void print_all_alpha_beta(trl_info * info, char *title, int jnd,
		double *alfrot, double *betrot);

	void print_lambda_res(trl_info * info, int jnd, double *lambda,
		double *res);

	void trl_orth(int nrow, double *v1, int ld1, int m1, double *v2, int ld2,
		int m2, double *rr, int kept, double *alpha, double *beta,
		double *wrk, int lwrk, trl_info * info);

	void trl_initial_guess(int nrow, double *evec, int lde, int mev,
		double *base, int ldb, int nbas, double *alpha,
		double *beta, int *j1, int *j2, trl_info * info,
		double *wrk, int lwrk);

	void trl_tridiag(int nd, double *alpha, double *beta, double *rot,
		double *alfrot, double *betrot, double *wrk, int lwrk,
		int *ierr);

	void trl_sort_eig(int nd, int lohi, int nec, double ref, double *lambda,
		double *res);

	void trl_get_tvec(int nd, double *alpha, double *beta, int irot, int nrot,
		double *rot, int nlam, double *lambda, double *yy,
		int *iwrk, double *wrk, int lwrk, int *ierr);

	void trl_get_tvec_a(int nd, int kept, double *alpha, double *beta,
		int nlam, double *lambda, double *yy, double *wrk,
		int lwrk, int *iwrk, int *ierr);

	void trl_get_eval(int nd, int locked, double *alpha, double *beta,
		double *lambda, double *res, double *wrk, int lwrk,
		int *ierr);

	void trl_set_locking(int jnd, int nlam, double *lambda, double *res,
		double *yy, int anrm, int *locked);

	void trl_ritz_vectors(int nrow, int lck, int ny, double *yy, int ldy,
		double *vec1, int ld1, int m1, double *vec2,
		int ld2, int m2, double *wrk, int lwrk);

	int trl_cgs(trl_info * info, int nrow, double *v1, int ld1, int m1,
		double *v2, int ld2, int m2, double *rr, double *rnrm,
		double *alpha, int *north, double *wrk);

	void trl_convergence_test(int nd, double *lambda, double *res,
		trl_info * info, double *wrk);
	int trl_check_dgen(trl_info * info, int jnd, double *lambda,
		double *res);

///////////////////////////////////////////////////////////////////////

	static void
		trl_restart_fixed(int nd, int mnd, int tind, double *lambda,
		double *res, trl_info * info, int *kl, int *kr);


	static double trl_min_gap_ratio_(trl_info * info, int nd, int tind, double *res);
	inline double gap_ratio_(int i, int j, int tind, double *lambda);
	////
	static void
		trl_restart_small_res(int nd, int mnd, int tind, double *lambda,
		double *res, trl_info * info, int *kl, int *kr);
	static double maxval(int n, double *a);
	////
	static void
		trl_restart_max_gap_ratio(int nd, int tind, int kept, double *lambda,
		double *res, trl_info * info, int *kl, int *kr);
	static void
		trl_restart_max_progress(int nd, int tind, int kept, double *lambda,
		double *res, trl_info * info, int *kl, int *kr);

	static void
		trl_restart_max_reduction(int nd, int tind, int kept, double *lambda,
		double *res, trl_info * info, int *kl, int *kr);

	static void
		trl_restart_scan(int nd, double *res, trl_info * info, int kept,
		int *kl, int *kr);

	static void
		trl_restart_max_gap_cost_ratio(int n, int tind, trl_info * info,
		double *lambda, double *res, int *kl, int *kr);

	static void
		trl_restart_max_gap_cost_ratio_static(int n, int tind, trl_info * info,
		double *lambda, double *res,
		int *kl, int *kr);

	static void trl_restart_search_range_(int nd, double *lambda, double *res,
		trl_info * info, int ncl, int ncr,
		int *lohi, int tind, int *klm, int *krm);

///////////////////////////////////////////////////////////////////////////////

	void trl_dgemv(char *trans, int m, int n, double alpha, double *a, int lda,
		double *x, int incx, double beta, double *y, int incy);
	double trl_ddot(integer_ n, double *dx, integer_ incx,
		double *dy, integer_ incy);
	void trl_daxpy(int n, double da, double *dx, int incx, double *dy,
		int incy);
	void trl_dscal(int n, double da, double *dx, int incx);
	void trl_dgemm(char *transa, char *transb, int m, int n, int k,
		double alpha, double *a, int lda, double *b, int ldb,
		double beta, double *c, int ldc);
	void trl_dcopy(int n, double *dx, int incx, double *dy, int incy);
////////////////////////////////////////////////////////////////////////////////////


	int dgemv_(char *trans, integer_ *m, integer_ *n, doublereal_ *
		alpha, doublereal_ *a, integer_ *lda, doublereal_ *x, integer_ *incx, 
		doublereal_ *beta, doublereal_ *y, integer_ *incy);
	doublereal_ ddot_(integer_ *n, doublereal_ *dx, integer_ *incx, doublereal_ *dy, 
		integer_ *incy);
	int daxpy_(integer_ *n, doublereal_ *da, doublereal_ *dx, 
		integer_ *incx, doublereal_ *dy, integer_ *incy);
	int dscal_(integer_ *n, doublereal_ *da, doublereal_ *dx, 
		integer_ *incx);
	int dgemm_(char *transa, char *transb, integer_ *m, integer_ *
		n, integer_ *k, doublereal_ *alpha, doublereal_ *a, integer_ *lda, 
		doublereal_ *b, integer_ *ldb, doublereal_ *beta, doublereal_ *c__, 
		integer_ *ldc);
	int dcopy_(integer_ *n, doublereal_ *dx, integer_ *incx, 
		doublereal_ *dy, integer_ *incy);
	inline double d_sign(doublereal_ * a, doublereal_ * b);
	int dstqrb_(integer_ * n, doublereal_ * d__, doublereal_ * e,
		doublereal_ * z__, doublereal_ * work, integer_ * info);
	void dsort2(int N, double *ARRAY1, double *ARRAY2);
	void dsort2a(int N, double *ARRAY1, double *ARRAY2);
	void dsort2s(int N, double s, double *ARRAY1, double *ARRAY2);
	void dsort2su_(int N, double s, double *ARRAY1, double *ARRAY2);
	void dsort2sd(int N, double s, double *ARRAY1, double *ARRAY2);
	inline double d_sign(doublereal_ * a, doublereal_ * b);
	int dstqrb_(integer_ * n, doublereal_ * d__, doublereal_ * e,
		doublereal_ * z__, doublereal_ * work, integer_ * info);

/////////////////////////////////////////////////
	void trl_clear_counter(trl_info * info, int nrow, int mev, int lde);

	void trl_print_setup(trl_info * info, int lbas, int lmis, int lwrk);

	void trl_g_sum(int mpicom, int nelm, double *x, double *y);

	int trl_sync_flag_(int mpicom, int inflag_);

	void trl_g_dot_(int mpicom, int nrow, double *v1, int ld1, int m1,
		double *v2, int ld2, int m2, double *rr, double *wrk);


///////////////////////////////////////////
	double ran1(int &idum);
	double ran2(int &idum);
	void srand48_new();
	double drand48_new();
	void trl_smooth_rr(int n, double *rr);
	void write_checkpoint(trl_info * info, char *title, int nrow,
		double *alpha, double *beta, double *evec, int lde,
		double *base, int ldb, int j1n, int jnd, int j2n);
	void print_restart_state(trl_info * info, char *title, int nrow,
		int mev, double *alpha, double *beta,
		double *betrot, double *evec, double *yy,
		int kept, int locked, int *iwrk, double *wrk2,
		int i2, int jml);
	void print_final_state(trl_info * info, char *title, int nrow, int mev,
		double *eval, double *beta, double *evec,
		double *yy, int kept, int jml);
	void log_error_state(trl_info * info, int kept, int j1, int j2, int jnd,
		int nrow, int mev, double *eval, double *alpha,
		double *alfrot, double *beta, double *betrot,
		double *evec, double *base, double *qa, double *qb,
		double *rr, char *title, int *iwrk);
/////////////////////////////////////////////////////////////////////////////////////
	logical_ lsame_(char *ca, char *cb);
	//
	doublereal_ dlamch_(char *cmach);
	//
	int dlamc1_(integer_ *beta, integer_ *t, logical_ *rnd, logical_ 
		*ieee1);
	//
	int dlamc2_(integer_ *beta, integer_ *t, logical_ *rnd, 
		doublereal_ *eps, integer_ *emin, doublereal_ *rmin, integer_ *emax, 
		doublereal_ *rmax);
	//
	doublereal_ dlamc3_(doublereal_ *a, doublereal_ *b);
	//
	int dlamc4_(integer_ *emin, doublereal_ *start, integer_ *base);
	//
	int dlamc5_(integer_ *beta, integer_ *p, integer_ *emin, 
		logical_ *ieee, integer_ *emax, doublereal_ *rmax);
	//
	int dlascl_(char *type__, integer_ *kl, integer_ *ku, 
		doublereal_ *cfrom, doublereal_ *cto, integer_ *m, integer_ *n, 
		doublereal_ *a, integer_ *lda, integer_ *info);
	integer_ ilaenv_(integer_ *ispec, char *name__, char *opts, integer_ *n1, 
		integer_ *n2, integer_ *n3, integer_ *n4, ftnlen_ name_len, ftnlen_ 
		opts_len);
	int xerbla_(char *srname, integer_ *info);
	//
	int dsterf_(integer_ *n, doublereal_ *d__, doublereal_ *e, 
		integer_ *info);
	//
	doublereal_ dlansy_(char *norm, char *uplo, integer_ *n, doublereal_ *a, integer_ 
		*lda, doublereal_ *work);
	//
	int dtrsm_(char *side, char *uplo, char *transa, char *diag, 
		integer_ *m, integer_ *n, doublereal_ *alpha, doublereal_ *a, integer_ *
		lda, doublereal_ *b, integer_ *ldb);
	//
	int dsyrk_(char *uplo, char *trans, integer_ *n, integer_ *k, 
		doublereal_ *alpha, doublereal_ *a, integer_ *lda, doublereal_ *beta, 
		doublereal_ *c__, integer_ *ldc);
	//
	int dpotf2_(char *uplo, integer_ *n, doublereal_ *a, integer_ *
		lda, integer_ *info);
	//
	int dorgql_(integer_ *m, integer_ *n, integer_ *k, doublereal_ *
		a, integer_ *lda, doublereal_ *tau, doublereal_ *work, integer_ *lwork, 
		integer_ *info);
	//
	int dorgqr_(integer_ *m, integer_ *n, integer_ *k, doublereal_ *
		a, integer_ *lda, doublereal_ *tau, doublereal_ *work, integer_ *lwork, 
		integer_ *info);
	//
	integer_ ieeeck_(integer_ *ispec, trl_real_ *zero, trl_real_ *one);
	//
	integer_ iparmq_(integer_ *ispec, char *name__, char *opts, integer_ *n, integer_ 
		*ilo, integer_ *ihi, integer_ *lwork, ftnlen_ name_len, ftnlen_ 
		opts_len);
	//
	int dlae2_(doublereal_ *a, doublereal_ *b, doublereal_ *c__, 
		doublereal_ *rt1, doublereal_ *rt2);
	//
	doublereal_ dlapy2_(doublereal_ *x, doublereal_ *y);
	//
	doublereal_ dlanst_(char *norm, integer_ *n, doublereal_ *d__, doublereal_ *e);
	//
	int dlasrt_(char *id, integer_ *n, doublereal_ *d__, integer_ *
		info);
	//
	int dlassq_(integer_ *n, doublereal_ *x, integer_ *incx, 
		doublereal_ *scale, doublereal_ *sumsq);
	//
	int dorg2l_(integer_ *m, integer_ *n, integer_ *k, doublereal_ *
		a, integer_ *lda, doublereal_ *tau, doublereal_ *work, integer_ *info);
	//
	int dlarfb_(char *side, char *trans, char *direct, char *
		storev, integer_ *m, integer_ *n, integer_ *k, doublereal_ *v, integer_ *
		ldv, doublereal_ *t, integer_ *ldt, doublereal_ *c__, integer_ *ldc, 
		doublereal_ *work, integer_ *ldwork);
	int dlarft_(char *direct, char *storev, integer_ *n, integer_ *
		k, doublereal_ *v, integer_ *ldv, doublereal_ *tau, doublereal_ *t, 
		integer_ *ldt);
	//
	int dorg2r_(integer_ *m, integer_ *n, integer_ *k, doublereal_ *
		a, integer_ *lda, doublereal_ *tau, doublereal_ *work, integer_ *info);
	//
	int dlarf_(char *side, integer_ *m, integer_ *n, doublereal_ *v,
		integer_ *incv, doublereal_ *tau, doublereal_ *c__, integer_ *ldc, 
		doublereal_ *work);
	//
	int dtrmm_(char *side, char *uplo, char *transa, char *diag, 
		integer_ *m, integer_ *n, doublereal_ *alpha, doublereal_ *a, integer_ *
		lda, doublereal_ *b, integer_ *ldb);
	//
	int dtrmv_(char *uplo, char *trans, char *diag, integer_ *n, 
		doublereal_ *a, integer_ *lda, doublereal_ *x, integer_ *incx);
	//
	int dger_(integer_ *m, integer_ *n, doublereal_ *alpha, 
		doublereal_ *x, integer_ *incx, doublereal_ *y, integer_ *incy, 
		doublereal_ *a, integer_ *lda);
	/////////////////////////////////////
	///////////////////////////////
	///////////////////////////////////////
	//////////////////////////

	double pow_di(doublereal_ *ap, integer_ *bp);
	/*integer_ s_wsfe(cilist_ *a)*/

	void s_copy(register char *a, register char *b, ftnlen_ la, ftnlen_ lb);
	//
	integer_ s_cmp(char *a0, char *b0, ftnlen_ la, ftnlen_ lb);
	double d_sign(doublereal_ *a, doublereal_ *b);
	//
	integer_ i_nint(trl_real_ *x);
	double sqrt_(doublereal_ x);
	double log_(doublereal_ x);
	double cos_(doublereal_ x);


	////////////
	////////////////
	////////////////////////
	int dsteqr_(char *compz, integer_ *n, doublereal_ *d__, 
		doublereal_ *e, doublereal_ *z__, integer_ *ldz, doublereal_ *work, 
		integer_ *info);
	////
	int dlasr_(char *side, char *pivot, char *direct, integer_ *m,
		integer_ *n, doublereal_ *c__, doublereal_ *s, doublereal_ *a, integer_ *
		lda);
	//
	int dswap_(integer_ *n, doublereal_ *dx, integer_ *incx, 
		doublereal_ *dy, integer_ *incy);
	//
	int dlaev2_(doublereal_ *a, doublereal_ *b, doublereal_ *c__, 
		doublereal_ *rt1, doublereal_ *rt2, doublereal_ *cs1, doublereal_ *sn1);
	//

	int dlaset_(char *uplo, integer_ *m, integer_ *n, doublereal_ *
		alpha, doublereal_ *beta, doublereal_ *a, integer_ *lda);
	//

	int dlartg_(doublereal_ *f, doublereal_ *g, doublereal_ *cs, 
		doublereal_ *sn, doublereal_ *r__);
	//////////////////////////////////////////////////////////////////////////
	doublereal_ dnrm2_(integer_ *n, doublereal_ *x, integer_ *incx);

	//
	doublereal_ dasum_(integer_ *n, doublereal_ *dx, integer_ *incx);
	int dlagtf_(integer_ *n, doublereal_ *a, doublereal_ *lambda, 
		doublereal_ *b, doublereal_ *c__, doublereal_ *tol, doublereal_ *d__, 
		integer_ *in, integer_ *info);

	integer_ idamax_(integer_ *n, doublereal_ *dx, integer_ *incx);
	int dlagts_(integer_ *job, integer_ *n, doublereal_ *a, 
		doublereal_ *b, doublereal_ *c__, doublereal_ *d__, integer_ *in, 
		doublereal_ *y, doublereal_ *tol, integer_ *info);

	int dlarnv_(integer_ *idist, integer_ *iseed, integer_ *n, 
		doublereal_ *x);

	int dlaruv_(integer_ *iseed, integer_ *n, doublereal_ *x);


	//////////////////////////////////////////////////
	////////////////////////////////
	/////////////////////////////////////////

	int dsyev_(char *jobz, char *uplo, integer_ *n, doublereal_ *a,
		integer_ *lda, doublereal_ *w, doublereal_ *work, integer_ *lwork, 
		integer_ *info);


	int dpotrf_(char *uplo, integer_ *n, doublereal_ *a, integer_ *
		lda, integer_ *info);

	//
	int dtrtrs_(char *uplo, char *trans, char *diag, integer_ *n, 
		integer_ *nrhs, doublereal_ *a, integer_ *lda, doublereal_ *b, integer_ *
		ldb, integer_ *info);


	int dsytrd_(char *uplo, integer_ *n, doublereal_ *a, integer_ *
		lda, doublereal_ *d__, doublereal_ *e, doublereal_ *tau, doublereal_ *
		work, integer_ *lwork, integer_ *info);

	int dorgtr_(char *uplo, integer_ *n, doublereal_ *a, integer_ *
		lda, doublereal_ *tau, doublereal_ *work, integer_ *lwork, integer_ *info);

	int dstein_(integer_ *n, doublereal_ *d__, doublereal_ *e, 
		integer_ *m, doublereal_ *w, integer_ *iblock, integer_ *isplit, 
		doublereal_ *z__, integer_ *ldz, doublereal_ *work, integer_ *iwork, 
		integer_ *ifail, integer_ *info);



	//////////////////////////
	//////////////////////////////
	/////////////////////////////////////

	int dsytd2_(char *uplo, integer_ *n, doublereal_ *a, integer_ *
		lda, doublereal_ *d__, doublereal_ *e, doublereal_ *tau, integer_ *info);

	int dsyr2k_(char *uplo, char *trans, integer_ *n, integer_ *k, 
		doublereal_ *alpha, doublereal_ *a, integer_ *lda, doublereal_ *b, 
		integer_ *ldb, doublereal_ *beta, doublereal_ *c__, integer_ *ldc);



	/////////
	int dlatrd_(char *uplo, integer_ *n, integer_ *nb, doublereal_ *
		a, integer_ *lda, doublereal_ *e, doublereal_ *tau, doublereal_ *w, 
		integer_ *ldw);


	int dsyr2_(char *uplo, integer_ *n, doublereal_ *alpha, 
		doublereal_ *x, integer_ *incx, doublereal_ *y, integer_ *incy, 
		doublereal_ *a, integer_ *lda);

	//

	int dsymv_(char *uplo, integer_ *n, doublereal_ *alpha, 
		doublereal_ *a, integer_ *lda, doublereal_ *x, integer_ *incx, doublereal_ 
		*beta, doublereal_ *y, integer_ *incy);


	//

	int dlarfg_(integer_ *n, doublereal_ *alpha, doublereal_ *x, 
		integer_ *incx, doublereal_ *tau);
	//








































}


#endif