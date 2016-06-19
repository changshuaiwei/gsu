#include "trl.h"

int TRL::dgemv_(char *trans, integer_ *m, integer_ *n, doublereal_ *
							alpha, doublereal_ *a, integer_ *lda, doublereal_ *x, integer_ *incx, 
							doublereal_ *beta, doublereal_ *y, integer_ *incy)
{
	/* System generated locals */
	integer_ a_dim1, a_offset, i__1, i__2;
	/* Local variables */
	static integer_ info;
	static doublereal_ temp;
	static integer_ lenx, leny, i__, j;
	//extern logical_ lsame_(char *, char *);
	static integer_ ix, iy, jx, jy, kx, ky;
	//extern /* Subroutine */ int xerbla_(char *, integer_ *);
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]
	/*  Purpose   
	=======   
	DGEMV  performs one of the matrix-vector operations   
	y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,   
	where alpha and beta are scalars, x and y are vectors and A is an   
	m by n matrix.   
	Parameters   
	==========   
	TRANS  - CHARACTER*1.   
	On entry, TRANS specifies the operation to be performed as   
	follows:   
	TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.   
	TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.   
	TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.   
	Unchanged on exit.   
	M      - INTEGER.   
	On entry, M specifies the number of rows of the matrix A.   
	M must be at least zero.   
	Unchanged on exit.   
	N      - INTEGER.   
	On entry, N specifies the number of columns of the matrix A.   
	N must be at least zero.   
	Unchanged on exit.   
	ALPHA  - DOUBLE PRECISION.   
	On entry, ALPHA specifies the scalar alpha.   
	Unchanged on exit.   
	A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).   
	Before entry, the leading m by n part of the array A must   
	contain the matrix of coefficients.   
	Unchanged on exit.   
	LDA    - INTEGER.   
	On entry, LDA specifies the first dimension of A as declared   
	in the calling (sub) program. LDA must be at least   
	max( 1, m ).   
	Unchanged on exit.   
	X      - DOUBLE PRECISION array of DIMENSION at least   
	( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'   
	and at least   
	( 1 + ( m - 1 )*abs( INCX ) ) otherwise.   
	Before entry, the incremented array X must contain the   
	vector x.   
	Unchanged on exit.   
	INCX   - INTEGER.   
	On entry, INCX specifies the increment for the elements of   
	X. INCX must not be zero.   
	Unchanged on exit.   
	BETA   - DOUBLE PRECISION.   
	On entry, BETA specifies the scalar beta. When BETA is   
	supplied as zero then Y need not be set on input.   
	Unchanged on exit.   
	Y      - DOUBLE PRECISION array of DIMENSION at least   
	( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'   
	and at least   
	( 1 + ( n - 1 )*abs( INCY ) ) otherwise.   
	Before entry with BETA non-zero, the incremented array Y   
	must contain the vector y. On exit, Y is overwritten by the   
	updated vector y.   
	INCY   - INTEGER.   
	On entry, INCY specifies the increment for the elements of   
	Y. INCY must not be zero.   
	Unchanged on exit.   
	Level 2 Blas routine.   
	-- Written on 22-October-1986.   
	Jack Dongarra, Argonne National Lab.   
	Jeremy Du Croz, Nag Central Office.   
	Sven Hammarling, Nag Central Office.   
	Richard Hanson, Sandia National Labs.   
	Test the input parameters.   
	Parameter adjustments */
	a_dim1 = *lda;
	a_offset = 1 + a_dim1 * 1;
	a -= a_offset;
	--x;
	--y;
	/* Function Body */
	info = 0;
	if (! lsame_(trans, "N") && ! lsame_(trans, "T") && ! lsame_(trans, "C")
		) {
			info = 1;
	} else if (*m < 0) {
		info = 2;
	} else if (*n < 0) {
		info = 3;
	} else if (*lda < max(1,*m)) {
		info = 6;
	} else if (*incx == 0) {
		info = 8;
	} else if (*incy == 0) {
		info = 11;
	}
	if (info != 0) {
		xerbla_("DGEMV ", &info);
		return 0;
	}
	/*     Quick return if possible. */
	if (*m == 0 || *n == 0 || *alpha == 0. && *beta == 1.) {
		return 0;
	}
	/*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set   
	up the start points in  X  and  Y. */
	if (lsame_(trans, "N")) {
		lenx = *n;
		leny = *m;
	} else {
		lenx = *m;
		leny = *n;
	}
	if (*incx > 0) {
		kx = 1;
	} else {
		kx = 1 - (lenx - 1) * *incx;
	}
	if (*incy > 0) {
		ky = 1;
	} else {
		ky = 1 - (leny - 1) * *incy;
	}
	/*     Start the operations. In this version the elements of A are   
	accessed sequentially with one pass through A.   
	First form  y := beta*y. */
	if (*beta != 1.) {
		if (*incy == 1) {
			if (*beta == 0.) {
				i__1 = leny;
				for (i__ = 1; i__ <= i__1; ++i__) {
					y[i__] = 0.;
					/* L10: */
				}
			} else {
				i__1 = leny;
				for (i__ = 1; i__ <= i__1; ++i__) {
					y[i__] = *beta * y[i__];
					/* L20: */
				}
			}
		} else {
			iy = ky;
			if (*beta == 0.) {
				i__1 = leny;
				for (i__ = 1; i__ <= i__1; ++i__) {
					y[iy] = 0.;
					iy += *incy;
					/* L30: */
				}
			} else {
				i__1 = leny;
				for (i__ = 1; i__ <= i__1; ++i__) {
					y[iy] = *beta * y[iy];
					iy += *incy;
					/* L40: */
				}
			}
		}
	}
	if (*alpha == 0.) {
		return 0;
	}
	if (lsame_(trans, "N")) {
		/*        Form  y := alpha*A*x + y. */
		jx = kx;
		if (*incy == 1) {
			i__1 = *n;
			for (j = 1; j <= i__1; ++j) {
				if (x[jx] != 0.) {
					temp = *alpha * x[jx];
					i__2 = *m;
					for (i__ = 1; i__ <= i__2; ++i__) {
						y[i__] += temp * a_ref(i__, j);
						/* L50: */
					}
				}
				jx += *incx;
				/* L60: */
			}
		} else {
			i__1 = *n;
			for (j = 1; j <= i__1; ++j) {
				if (x[jx] != 0.) {
					temp = *alpha * x[jx];
					iy = ky;
					i__2 = *m;
					for (i__ = 1; i__ <= i__2; ++i__) {
						y[iy] += temp * a_ref(i__, j);
						iy += *incy;
						/* L70: */
					}
				}
				jx += *incx;
				/* L80: */
			}
		}
	} else {
		/*        Form  y := alpha*A'*x + y. */
		jy = ky;
		if (*incx == 1) {
			i__1 = *n;
			for (j = 1; j <= i__1; ++j) {
				temp = 0.;
				i__2 = *m;
				for (i__ = 1; i__ <= i__2; ++i__) {
					temp += a_ref(i__, j) * x[i__];
					/* L90: */
				}
				y[jy] += *alpha * temp;
				jy += *incy;
				/* L100: */
			}
		} else {
			i__1 = *n;
			for (j = 1; j <= i__1; ++j) {
				temp = 0.;
				ix = kx;
				i__2 = *m;
				for (i__ = 1; i__ <= i__2; ++i__) {
					temp += a_ref(i__, j) * x[ix];
					ix += *incx;
					/* L110: */
				}
				y[jy] += *alpha * temp;
				jy += *incy;
				/* L120: */
			}
		}
	}
	#undef a_ref
	return 0;
	/*     End of DGEMV . */
} /* dgemv_ */


doublereal_ TRL::ddot_(integer_ *n, doublereal_ *dx, integer_ *incx, doublereal_ *dy, 
				 integer_ *incy)
{
	/* System generated locals */
	integer_ i__1;
	doublereal_ ret_val;
	/* Local variables */
	static integer_ i__, m;
	static doublereal_ dtemp;
	static integer_ ix, iy, mp1;
	/*     forms the dot product of two vectors.   
	uses unrolled loops for increments equal to one.   
	jack dongarra, linpack, 3/11/78.   
	modified 12/3/93, array(1) declarations changed to array(*)   
	Parameter adjustments */
	--dy;
	--dx;
	/* Function Body */
	//printf( " ** in ddot **\n" );
	ret_val = 0.;
	dtemp = 0.;
	if (*n <= 0) {
		return ret_val;
	}
	if (*incx == 1 && *incy == 1) {
		goto L20;
	}
	/*        code for unequal increments or equal increments   
	not equal to 1 */
	ix = 1;
	iy = 1;
	if (*incx < 0) {
		ix = (-(*n) + 1) * *incx + 1;
	}
	if (*incy < 0) {
		iy = (-(*n) + 1) * *incy + 1;
	}
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		dtemp += dx[ix] * dy[iy];
		ix += *incx;
		iy += *incy;
		/* L10: */
	}
	ret_val = dtemp;
	return ret_val;
	/*        code for both increments equal to 1   
	clean-up loop */
L20:
	m = *n % 5;
	if (m == 0) {
		goto L40;
	}
	i__1 = m;
	for (i__ = 1; i__ <= i__1; ++i__) {
		dtemp += dx[i__] * dy[i__];
		/* L30: */
	}
	if (*n < 5) {
		goto L60;
	}
L40:
	mp1 = m + 1;
	i__1 = *n;

	//printf( "%d %d\n",mp1,i__1);
	for (i__ = mp1; i__ <= i__1; i__ += 5) {
		dtemp = dtemp + dx[i__] * dy[i__] + dx[i__ + 1] * dy[i__ + 1] + dx[
			i__ + 2] * dy[i__ + 2] + dx[i__ + 3] * dy[i__ + 3] + dx[i__ + 
				4] * dy[i__ + 4];
			//printf( "%e\n",dtemp );
			/* L50: */
	}
L60:
	ret_val = dtemp;
	return ret_val;
} /* ddot_ */




int TRL::daxpy_(integer_ *n, doublereal_ *da, doublereal_ *dx, 
							integer_ *incx, doublereal_ *dy, integer_ *incy)
{
	/* System generated locals */
	integer_ i__1;
	/* Local variables */
	static integer_ i__, m, ix, iy, mp1;
	/*     constant times a vector plus a vector.   
	uses unrolled loops for increments equal to one.   
	jack dongarra, linpack, 3/11/78.   
	modified 12/3/93, array(1) declarations changed to array(*)   
	Parameter adjustments */
	--dy;
	--dx;
	/* Function Body */
	if (*n <= 0) {
		return 0;
	}
	if (*da == 0.) {
		return 0;
	}
	if (*incx == 1 && *incy == 1) {
		goto L20;
	}
	/*        code for unequal increments or equal increments   
	not equal to 1 */
	ix = 1;
	iy = 1;
	if (*incx < 0) {
		ix = (-(*n) + 1) * *incx + 1;
	}
	if (*incy < 0) {
		iy = (-(*n) + 1) * *incy + 1;
	}
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		dy[iy] += *da * dx[ix];
		ix += *incx;
		iy += *incy;
		/* L10: */
	}
	return 0;
	/*        code for both increments equal to 1   
	clean-up loop */
L20:
	m = *n % 4;
	if (m == 0) {
		goto L40;
	}
	i__1 = m;
	for (i__ = 1; i__ <= i__1; ++i__) {
		dy[i__] += *da * dx[i__];
		/* L30: */
	}
	if (*n < 4) {
		return 0;
	}
L40:
	mp1 = m + 1;
	i__1 = *n;
	for (i__ = mp1; i__ <= i__1; i__ += 4) {
		dy[i__] += *da * dx[i__];
		dy[i__ + 1] += *da * dx[i__ + 1];
		dy[i__ + 2] += *da * dx[i__ + 2];
		dy[i__ + 3] += *da * dx[i__ + 3];
		/* L50: */
	}
	return 0;
} /* daxpy_ */



int TRL::dscal_(integer_ *n, doublereal_ *da, doublereal_ *dx, 
		   integer_ *incx)
{
	/* System generated locals */
	integer_ i__1, i__2;
	/* Local variables */
	static integer_ i__, m, nincx, mp1;
	/*     scales a vector by a constant.   
	uses unrolled loops for increment equal to one.   
	jack dongarra, linpack, 3/11/78.   
	modified 3/93 to return if incx .le. 0.   
	modified 12/3/93, array(1) declarations changed to array(*)   
	Parameter adjustments */
	--dx;
	/* Function Body */
	if (*n <= 0 || *incx <= 0) {
		return 0;
	}
	if (*incx == 1) {
		goto L20;
	}
	/*        code for increment not equal to 1 */
	nincx = *n * *incx;
	i__1 = nincx;
	i__2 = *incx;
	for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
		dx[i__] = *da * dx[i__];
		/* L10: */
	}
	return 0;
	/*        code for increment equal to 1   
	clean-up loop */
L20:
	m = *n % 5;
	if (m == 0) {
		goto L40;
	}
	i__2 = m;
	for (i__ = 1; i__ <= i__2; ++i__) {
		dx[i__] = *da * dx[i__];
		/* L30: */
	}
	if (*n < 5) {
		return 0;
	}
L40:
	mp1 = m + 1;
	i__2 = *n;
	for (i__ = mp1; i__ <= i__2; i__ += 5) {
		dx[i__] = *da * dx[i__];
		dx[i__ + 1] = *da * dx[i__ + 1];
		dx[i__ + 2] = *da * dx[i__ + 2];
		dx[i__ + 3] = *da * dx[i__ + 3];
		dx[i__ + 4] = *da * dx[i__ + 4];
		/* L50: */
	}
	return 0;
} /* dscal_ */





int TRL::dgemm_(char *transa, char *transb, integer_ *m, integer_ *
		   n, integer_ *k, doublereal_ *alpha, doublereal_ *a, integer_ *lda, 
		   doublereal_ *b, integer_ *ldb, doublereal_ *beta, doublereal_ *c__, 
		   integer_ *ldc)
{
	/* System generated locals */
	integer_ a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, 
		i__3;
	/* Local variables */
	static integer_ info;
	static logical_ nota, notb;
	static doublereal_ temp;
	static integer_ i__, j, l, ncola;
	//extern logical_ lsame_(char *, char *);
	static integer_ nrowa, nrowb;
	//extern /* Subroutine */ int xerbla_(char *, integer_ *);
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]
#define b_ref(a_1,a_2) b[(a_2)*b_dim1 + a_1]
#define c___ref(a_1,a_2) c__[(a_2)*c_dim1 + a_1]
	/*  Purpose   
	=======   
	DGEMM  performs one of the matrix-matrix operations   
	C := alpha*op( A )*op( B ) + beta*C,   
	where  op( X ) is one of   
	op( X ) = X   or   op( X ) = X',   
	alpha and beta are scalars, and A, B and C are matrices, with op( A )   
	an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.   
	Parameters   
	==========   
	TRANSA - CHARACTER*1.   
	On entry, TRANSA specifies the form of op( A ) to be used in   
	the matrix multiplication as follows:   
	TRANSA = 'N' or 'n',  op( A ) = A.   
	TRANSA = 'T' or 't',  op( A ) = A'.   
	TRANSA = 'C' or 'c',  op( A ) = A'.   
	Unchanged on exit.   
	TRANSB - CHARACTER*1.   
	On entry, TRANSB specifies the form of op( B ) to be used in   
	the matrix multiplication as follows:   
	TRANSB = 'N' or 'n',  op( B ) = B.   
	TRANSB = 'T' or 't',  op( B ) = B'.   
	TRANSB = 'C' or 'c',  op( B ) = B'.   
	Unchanged on exit.   
	M      - INTEGER.   
	On entry,  M  specifies  the number  of rows  of the  matrix   
	op( A )  and of the  matrix  C.  M  must  be at least  zero.   
	Unchanged on exit.   
	N      - INTEGER.   
	On entry,  N  specifies the number  of columns of the matrix   
	op( B ) and the number of columns of the matrix C. N must be   
	at least zero.   
	Unchanged on exit.   
	K      - INTEGER.   
	On entry,  K  specifies  the number of columns of the matrix   
	op( A ) and the number of rows of the matrix op( B ). K must   
	be at least  zero.   
	Unchanged on exit.   
	ALPHA  - DOUBLE PRECISION.   
	On entry, ALPHA specifies the scalar alpha.   
	Unchanged on exit.   
	A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is   
	k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.   
	Before entry with  TRANSA = 'N' or 'n',  the leading  m by k   
	part of the array  A  must contain the matrix  A,  otherwise   
	the leading  k by m  part of the array  A  must contain  the   
	matrix A.   
	Unchanged on exit.   
	LDA    - INTEGER.   
	On entry, LDA specifies the first dimension of A as declared   
	in the calling (sub) program. When  TRANSA = 'N' or 'n' then   
	LDA must be at least  max( 1, m ), otherwise  LDA must be at   
	least  max( 1, k ).   
	Unchanged on exit.   
	B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is   
	n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.   
	Before entry with  TRANSB = 'N' or 'n',  the leading  k by n   
	part of the array  B  must contain the matrix  B,  otherwise   
	the leading  n by k  part of the array  B  must contain  the   
	matrix B.   
	Unchanged on exit.   
	LDB    - INTEGER.   
	On entry, LDB specifies the first dimension of B as declared   
	in the calling (sub) program. When  TRANSB = 'N' or 'n' then   
	LDB must be at least  max( 1, k ), otherwise  LDB must be at   
	least  max( 1, n ).   
	Unchanged on exit.   
	BETA   - DOUBLE PRECISION.   
	On entry,  BETA  specifies the scalar  beta.  When  BETA  is   
	supplied as zero then C need not be set on input.   
	Unchanged on exit.   
	C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).   
	Before entry, the leading  m by n  part of the array  C must   
	contain the matrix  C,  except when  beta  is zero, in which   
	case C need not be set on entry.   
	On exit, the array  C  is overwritten by the  m by n  matrix   
	( alpha*op( A )*op( B ) + beta*C ).   
	LDC    - INTEGER.   
	On entry, LDC specifies the first dimension of C as declared   
	in  the  calling  (sub)  program.   LDC  must  be  at  least   
	max( 1, m ).   
	Unchanged on exit.   
	Level 3 Blas routine.   
	-- Written on 8-February-1989.   
	Jack Dongarra, Argonne National Laboratory.   
	Iain Duff, AERE Harwell.   
	Jeremy Du Croz, Numerical Algorithms Group Ltd.   
	Sven Hammarling, Numerical Algorithms Group Ltd.   
	Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not   
	transposed and set  NROWA, NCOLA and  NROWB  as the number of rows   
	and  columns of  A  and the  number of  rows  of  B  respectively.   
	Parameter adjustments */
	a_dim1 = *lda;
	a_offset = 1 + a_dim1 * 1;
	a -= a_offset;
	b_dim1 = *ldb;
	b_offset = 1 + b_dim1 * 1;
	b -= b_offset;
	c_dim1 = *ldc;
	c_offset = 1 + c_dim1 * 1;
	c__ -= c_offset;
	/* Function Body */
	nota = lsame_(transa, "N");
	notb = lsame_(transb, "N");
	if (nota) {
		nrowa = *m;
		ncola = *k;
	} else {
		nrowa = *k;
		ncola = *m;
	}
	if (notb) {
		nrowb = *k;
	} else {
		nrowb = *n;
	}
	/*     Test the input parameters. */
	info = 0;
	if (! nota && ! lsame_(transa, "C") && ! lsame_(
		transa, "T")) {
			info = 1;
	} else if (! notb && ! lsame_(transb, "C") && ! 
		lsame_(transb, "T")) {
			info = 2;
	} else if (*m < 0) {
		info = 3;
	} else if (*n < 0) {
		info = 4;
	} else if (*k < 0) {
		info = 5;
	} else if (*lda < max(1,nrowa)) {
		info = 8;
	} else if (*ldb < max(1,nrowb)) {
		info = 10;
	} else if (*ldc < max(1,*m)) {
		info = 13;
	}
	if (info != 0) {
		xerbla_("DGEMM ", &info);
		return 0;
	}
	/*     Quick return if possible. */
	if (*m == 0 || *n == 0 || (*alpha == 0. || *k == 0) && *beta == 1.) {
		return 0;
	}
	/*     And if  alpha.eq.zero. */
	if (*alpha == 0.) {
		if (*beta == 0.) {
			i__1 = *n;
			for (j = 1; j <= i__1; ++j) {
				i__2 = *m;
				for (i__ = 1; i__ <= i__2; ++i__) {
					c___ref(i__, j) = 0.;
					/* L10: */
				}
				/* L20: */
			}
		} else {
			i__1 = *n;
			for (j = 1; j <= i__1; ++j) {
				i__2 = *m;
				for (i__ = 1; i__ <= i__2; ++i__) {
					c___ref(i__, j) = *beta * c___ref(i__, j);
					/* L30: */
				}
				/* L40: */
			}
		}
		return 0;
	}
	/*     Start the operations. */
	if (notb) {
		if (nota) {
			/*           Form  C := alpha*A*B + beta*C. */
			i__1 = *n;
			for (j = 1; j <= i__1; ++j) {
				if (*beta == 0.) {
					i__2 = *m;
					for (i__ = 1; i__ <= i__2; ++i__) {
						c___ref(i__, j) = 0.;
						/* L50: */
					}
				} else if (*beta != 1.) {
					i__2 = *m;
					for (i__ = 1; i__ <= i__2; ++i__) {
						c___ref(i__, j) = *beta * c___ref(i__, j);
						/* L60: */
					}
				}
				i__2 = *k;
				for (l = 1; l <= i__2; ++l) {
					if (b_ref(l, j) != 0.) {
						temp = *alpha * b_ref(l, j);
						i__3 = *m;
						for (i__ = 1; i__ <= i__3; ++i__) {
							c___ref(i__, j) = c___ref(i__, j) + temp * a_ref(
								i__, l);
							/* L70: */
						}
					}
					/* L80: */
				}
				/* L90: */
			}
		} else {
			/*           Form  C := alpha*A'*B + beta*C */
			i__1 = *n;
			for (j = 1; j <= i__1; ++j) {
				i__2 = *m;
				for (i__ = 1; i__ <= i__2; ++i__) {
					temp = 0.;
					i__3 = *k;
					for (l = 1; l <= i__3; ++l) {
						temp += a_ref(l, i__) * b_ref(l, j);
						/* L100: */
					}
					if (*beta == 0.) {
						c___ref(i__, j) = *alpha * temp;
					} else {
						c___ref(i__, j) = *alpha * temp + *beta * c___ref(i__,
							j);
					}
					/* L110: */
				}
				/* L120: */
			}
		}
	} else {
		if (nota) {
			/*           Form  C := alpha*A*B' + beta*C */
			i__1 = *n;
			for (j = 1; j <= i__1; ++j) {
				if (*beta == 0.) {
					i__2 = *m;
					for (i__ = 1; i__ <= i__2; ++i__) {
						c___ref(i__, j) = 0.;
						/* L130: */
					}
				} else if (*beta != 1.) {
					i__2 = *m;
					for (i__ = 1; i__ <= i__2; ++i__) {
						c___ref(i__, j) = *beta * c___ref(i__, j);
						/* L140: */
					}
				}
				i__2 = *k;
				for (l = 1; l <= i__2; ++l) {
					if (b_ref(j, l) != 0.) {
						temp = *alpha * b_ref(j, l);
						i__3 = *m;
						for (i__ = 1; i__ <= i__3; ++i__) {
							c___ref(i__, j) = c___ref(i__, j) + temp * a_ref(
								i__, l);
							/* L150: */
						}
					}
					/* L160: */
				}
				/* L170: */
			}
		} else {
			/*           Form  C := alpha*A'*B' + beta*C */
			i__1 = *n;
			for (j = 1; j <= i__1; ++j) {
				i__2 = *m;
				for (i__ = 1; i__ <= i__2; ++i__) {
					temp = 0.;
					i__3 = *k;
					for (l = 1; l <= i__3; ++l) {
						temp += a_ref(l, i__) * b_ref(j, l);
						/* L180: */
					}
					if (*beta == 0.) {
						c___ref(i__, j) = *alpha * temp;
					} else {
						c___ref(i__, j) = *alpha * temp + *beta * c___ref(i__,
							j);
					}
					/* L190: */
				}
				/* L200: */
			}
		}
	}
#undef c___ref
#undef b_ref
#undef a_ref
	return 0;
	/*     End of DGEMM . */
} /* dgemm_ */








int TRL::dcopy_(integer_ *n, doublereal_ *dx, integer_ *incx, 
		   doublereal_ *dy, integer_ *incy)
{
	/* System generated locals */
	integer_ i__1;
	/* Local variables */
	static integer_ i__, m, ix, iy, mp1;
	/*     copies a vector, x, to a vector, y.   
	uses unrolled loops for increments equal to one.   
	jack dongarra, linpack, 3/11/78.   
	modified 12/3/93, array(1) declarations changed to array(*)   
	Parameter adjustments */
	--dy;
	--dx;
	/* Function Body */
	if (*n <= 0) {
		return 0;
	}
	if (*incx == 1 && *incy == 1) {
		goto L20;
	}
	/*        code for unequal increments or equal increments   
	not equal to 1 */
	ix = 1;
	iy = 1;
	if (*incx < 0) {
		ix = (-(*n) + 1) * *incx + 1;
	}
	if (*incy < 0) {
		iy = (-(*n) + 1) * *incy + 1;
	}
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		dy[iy] = dx[ix];
		ix += *incx;
		iy += *incy;
		/* L10: */
	}
	return 0;
	/*        code for both increments equal to 1   
	clean-up loop */
L20:
	m = *n % 7;
	if (m == 0) {
		goto L40;
	}
	i__1 = m;
	for (i__ = 1; i__ <= i__1; ++i__) {
		dy[i__] = dx[i__];
		/* L30: */
	}
	if (*n < 7) {
		return 0;
	}
L40:
	mp1 = m + 1;
	i__1 = *n;
	for (i__ = mp1; i__ <= i__1; i__ += 7) {
		dy[i__] = dx[i__];
		dy[i__ + 1] = dx[i__ + 1];
		dy[i__ + 2] = dx[i__ + 2];
		dy[i__ + 3] = dx[i__ + 3];
		dy[i__ + 4] = dx[i__ + 4];
		dy[i__ + 5] = dx[i__ + 5];
		dy[i__ + 6] = dx[i__ + 6];
		/* L50: */
	}
	return 0;
} /* dcopy_ */

inline double TRL::d_sign(doublereal_ * a, doublereal_ * b)
{
	double x;
	x = (*a >= 0 ? *a : -*a);
	return (*b >= 0 ? x : -x);
}

int TRL::dstqrb_(integer_ * n, doublereal_ * d__, doublereal_ * e,
			doublereal_ * z__, doublereal_ * work, integer_ * info)
{
	/*

	Purpose
	=======
	DSTQRB computes all eigenvalues and the last component of the eigenvectors of a
	symmetric tridiagonal matrix using the implicit QL or QR method.
	This is mainly a modification of the CLAPACK subroutine dsteqr.c

	Arguments
	=========

	N       (input) INTEGER
	The order of the matrix.  N >= 0.

	D       (input/output) DOUBLE PRECISION array, dimension (N)
	On entry, the diagonal elements of the tridiagonal matrix.
	On exit, if INFO = 0, the eigenvalues in ascending order.

	E       (input/output) DOUBLE PRECISION array, dimension (N-1)
	On entry, the (n-1) subdiagonal elements of the tridiagonal
	matrix.
	On exit, E has been destroyed.

	Z       (input/output) DOUBLE PRECISION array, dimension (LDZ, N)
	On entry, if  COMPZ = 'V', then Z contains the orthogonal
	matrix used in the reduction to tridiagonal form.
	On exit, if INFO = 0, then if  COMPZ = 'V', Z contains the
	orthonormal eigenvectors of the original symmetric matrix,
	and if COMPZ = 'I', Z contains the orthonormal eigenvectors
	of the symmetric tridiagonal matrix.
	If COMPZ = 'N', then Z is not referenced.

	WORK    (workspace) DOUBLE PRECISION array, dimension (max(1,2*N-2))
	If COMPZ = 'N', then WORK is not referenced.

	INFO    (output) INTEGER
	= 0:  successful exit
	< 0:  if INFO = -i, the i-th argument had an illegal value
	> 0:  the algorithm has failed to find all the eigenvalues in
	a total of 30*N iterations; if INFO = i, then i
	elements of E have not converged to zero; on exit, D
	and E contain the elements of a symmetric tridiagonal
	matrix which is orthogonally similar to the original
	matrix.

	=====================================================================
	*/
	/* Table of constant values */
	doublereal_ c_b10 = 1.;
	integer_ c__0 = 0;
	integer_ c__1 = 1;

	/* System generated locals */
	integer_ i__1, i__2;
	doublereal_ d__1, d__2;
	/* Builtin functions */
	//double sqrt_(doublereal_), d_sign(doublereal_ *, doublereal_ *);
	//extern /* Subroutine */ int dlae2_(doublereal_ *, doublereal_ *, doublereal_
	//	*, doublereal_ *, doublereal_ *);
	doublereal_ b, c__, f, g;
	integer_ i__, j, k, l, m;
	doublereal_ p, r__, s;
	//extern logical_ lsame_(char *, char *);
	//extern /* Subroutine */ int dlasr_(char *, char *, char *, integer_ *,
	//	integer_ *, doublereal_ *,
	//	doublereal_ *, doublereal_ *,
	//	integer_ *);
	doublereal_ anorm;
	//extern /* Subroutine */ int dswap_(integer_ *, doublereal_ *, integer_ *,
	//	doublereal_ *, integer_ *);
	integer_ l1;
	//extern /* Subroutine */ int dlaev2_(doublereal_ *, doublereal_ *,
	//	doublereal_ *, doublereal_ *,
	//	doublereal_ *, doublereal_ *,
	//	doublereal_ *);
	integer_ lendm1, lendp1;
	//extern doublereal_ dlapy2_(doublereal_ *, doublereal_ *);
	integer_ ii;
	//extern doublereal_ dlamch_(char *);
	integer_ mm, iscale;
	//extern /* Subroutine */ int dlascl_(char *, integer_ *, integer_ *,
	//	doublereal_ *, doublereal_ *,
	//	integer_ *, integer_ *, doublereal_ *,
	//	integer_ *, integer_ *),
	//	dlaset_(char *, integer_ *, integer_ *, doublereal_ *, doublereal_ *,
	//	doublereal_ *, integer_ *);
	doublereal_ safmin;
	//extern /* Subroutine */ int dlartg_(doublereal_ *, doublereal_ *,
	//	doublereal_ *, doublereal_ *,
	//	doublereal_ *);
	doublereal_ safmax;
	//extern /* Subroutine */ int xerbla_(char *, integer_ *);
	//extern doublereal_ dlanst_(char *, integer_ *, doublereal_ *,
	//	doublereal_ *);
	//extern /* Subroutine */ int dlasrt_(char *, integer_ *, doublereal_ *,
	//	integer_ *);
	/* Local variables */
	integer_ lend, jtot;
	integer_ lendsv;
	doublereal_ ssfmin;
	integer_ nmaxit;
	doublereal_ ssfmax;
	integer_ lm1, mm1, nm1;
	doublereal_ rt1, rt2, eps;
	integer_ lsv;
	doublereal_ tst, eps2;

	--d__;
	--e;
	--z__;
	/* z_dim1 = *ldz;             */
	/* z_offset = 1 + z_dim1 * 1; */
	/* z__ -= z_offset;           */
	--work;

	/* Function Body */
	*info = 0;
	/* Taken out for TRLan
	if (lsame_(compz, "N")) {
	icompz = 0;
	} else if (lsame_(compz, "V")) {
	icompz = 1;
	} else if (lsame_(compz, "I")) {
	icompz = 2;
	} else {
	icompz = -1;
	}
	if (icompz < 0) {
	*info = -1;
	} else if (*n < 0) {
	*info = -2;
	} else if (*ldz < 1 || icompz > 0 && *ldz < max(1,*n)) {
	*info = -6;
	}
	if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DSTEQR", &i__1);
	return 0;
	}
	*/
	/*  icompz = 2; */

	/*	Quick return if possible */

	if (*n == 0) {
		return 0;
	}

	if (*n == 1) {
		z__[1] = 1;
		return 0;
	}

	/*	Determine the unit roundoff and over/underflow thresholds. */

	eps = dlamch_("E");
	/*	Computing 2nd power */
	d__1 = eps;
	eps2 = d__1 * d__1;
	safmin = dlamch_("S");
	safmax = 1. / safmin;
	ssfmax = sqrt_(safmax) / 3.;
	ssfmin = sqrt_(safmin) / eps2;

	/*	Compute the eigenvalues and eigenvectors of the tridiagonal
	matrix. */
	/* Taken out for TRLan
	if (icompz == 2) {
	dlaset_("Full", n, n, &c_b9, &c_b10, &z__[z_offset], ldz);
	}
	*/
	for (j = 1; j < *n; j++) {
		z__[j] = 0.0;
	}
	z__[*n] = 1.0;
	nmaxit = *n * 30;
	jtot = 0;

	/*	Determine where the matrix splits and choose QL or QR iteration
	for each block, according to whether top or bottom diagonal
	element is smaller. */

	l1 = 1;
	nm1 = *n - 1;

L10:
	if (l1 > *n) {
		goto L160;
	}
	if (l1 > 1) {
		e[l1 - 1] = 0.;
	}
	if (l1 <= nm1) {
		i__1 = nm1;
		for (m = l1; m <= i__1; ++m) {
			tst = (d__1 = e[m], fabs(d__1));
			if (tst == 0.) {
				goto L30;
			}
			if (tst <=
				sqrt_((d__1 = d__[m], fabs(d__1))) * sqrt_((d__2 =
				d__[m + 1],
				fabs(d__2))) *
				eps) {
					e[m] = 0.;
					goto L30;
			}
			/* L20: */
		}
	}
	m = *n;

L30:
	l = l1;
	lsv = l;
	lend = m;
	lendsv = lend;
	l1 = m + 1;
	if (lend == l) {
		goto L10;
	}

	/*	Scale submatrix in rows and columns L to LEND */

	i__1 = lend - l + 1;
	anorm = dlanst_("I", &i__1, &d__[l], &e[l]);
	iscale = 0;
	if (anorm == 0.) {
		goto L10;
	}
	if (anorm > ssfmax) {
		iscale = 1;
		i__1 = lend - l + 1;
		dlascl_("G", &c__0, &c__0, &anorm, &ssfmax, &i__1, &c__1, &d__[l],
			n, info);
		i__1 = lend - l;
		dlascl_("G", &c__0, &c__0, &anorm, &ssfmax, &i__1, &c__1, &e[l], n,
			info);
	} else if (anorm < ssfmin) {
		iscale = 2;
		i__1 = lend - l + 1;
		dlascl_("G", &c__0, &c__0, &anorm, &ssfmin, &i__1, &c__1, &d__[l],
			n, info);
		i__1 = lend - l;
		dlascl_("G", &c__0, &c__0, &anorm, &ssfmin, &i__1, &c__1, &e[l], n,
			info);
	}

	/*	Choose between QL and QR iteration */

	if ((d__1 = d__[lend], fabs(d__1)) < (d__2 = d__[l], fabs(d__2))) {
		lend = lsv;
		l = lendsv;
	}

	if (lend > l) {

		/*	QL Iteration

		Look for small subdiagonal element. */

L40:

		if (l != lend) {
			lendm1 = lend - 1;
			i__1 = lendm1;
			for (m = l; m <= i__1; ++m) {
				/*		Computing 2nd power */
				d__2 = (d__1 = e[m], fabs(d__1));
				tst = d__2 * d__2;
				if (tst <=
					eps2 * (d__1 = d__[m], fabs(d__1)) * (d__2 =
					d__[m + 1],
					fabs(d__2)) +
					safmin) {
						goto L60;
				}
				/* L50: */
			}
		}

		m = lend;

L60:
		if (m < lend) {
			e[m] = 0.;
		}
		p = d__[l];
		if (m == l) {
			goto L80;
		}

		/*	If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
		to compute its eigensystem. */

		if (m == l + 1) {
			dlaev2_(&d__[l], &e[l], &d__[l + 1], &rt1, &rt2, &c__, &s);
			work[l] = c__;
			work[*n - 1 + l] = s;
			/* Taken out for TRLan
			dlasr_("R", "V", "B", n, &c__2, &work[l], &work[*n - 1 + l], &
			z___ref(1, l), ldz);
			*/
			tst = z__[l + 1];
			z__[l + 1] = c__ * tst - s * z__[l];
			z__[l] = s * tst + c__ * z__[l];
			d__[l] = rt1;
			d__[l + 1] = rt2;
			e[l] = 0.;
			l += 2;
			if (l <= lend) {
				goto L40;
			}
			goto L140;
		}

		if (jtot == nmaxit) {
			goto L140;
		}
		++jtot;

		/*	Form shift. */

		g = (d__[l + 1] - p) / (e[l] * 2.);
		r__ = dlapy2_(&g, &c_b10);
		g = d__[m] - p + e[l] / (g + d_sign(&r__, &g));

		s = 1.;
		c__ = 1.;
		p = 0.;

		/*	Inner loop */

		mm1 = m - 1;
		i__1 = l;
		for (i__ = mm1; i__ >= i__1; --i__) {
			f = s * e[i__];
			b = c__ * e[i__];
			dlartg_(&g, &f, &c__, &s, &r__);
			if (i__ != m - 1) {
				e[i__ + 1] = r__;
			}
			g = d__[i__ + 1] - p;
			r__ = (d__[i__] - g) * s + c__ * 2. * b;
			p = s * r__;
			d__[i__ + 1] = g + p;
			g = c__ * r__ - b;

			/*		If eigenvectors are desired, then save rotations. */

			work[i__] = c__;
			work[*n - 1 + i__] = -s;

			/* L70: */
		}

		/*	If eigenvectors are desired, then apply saved rotations. */

		mm = m - l + 1;
		/* Taken out for TRLan
		dlasr_("R", "V", "B", n, &mm, &work[l], &work[*n - 1 + l], &
		z___ref(1, l), ldz);
		*/
		dlasr_("R", "V", "B", &c__1, &mm, &work[l], &work[*n - 1 + l],
			&z__[l], &c__1);

		d__[l] -= p;
		e[l] = g;
		goto L40;

		/* 	Eigenvalue found. */

L80:
		d__[l] = p;

		++l;
		if (l <= lend) {
			goto L40;
		}
		goto L140;

	} else {

		/*	QR Iteration

		Look for small superdiagonal element. */

L90:
		if (l != lend) {
			lendp1 = lend + 1;
			i__1 = lendp1;
			for (m = l; m >= i__1; --m) {
				/*			Computing 2nd power */
				d__2 = (d__1 = e[m - 1], fabs(d__1));
				tst = d__2 * d__2;
				if (tst <=
					eps2 * (d__1 = d__[m], fabs(d__1)) * (d__2 =
					d__[m - 1],
					fabs(d__2)) +
					safmin) {
						goto L110;
				}
				/* L100: */
			}
		}

		m = lend;

L110:
		if (m > lend) {
			e[m - 1] = 0.;
		}
		p = d__[l];
		if (m == l) {
			goto L130;
		}

		/*	If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
		to compute its eigensystem. */

		if (m == l - 1) {
			dlaev2_(&d__[l - 1], &e[l - 1], &d__[l], &rt1, &rt2, &c__, &s);
			/* Taken out for TRLan
			work[m] = c__;
			work[*n - 1 + m] = s;
			dlasr_("R", "V", "F", n, &c__2, &work[m], &work[*n - 1 + m], &
			z___ref(1, l - 1), ldz);
			*/
			tst = z__[l];
			z__[l] = c__ * tst - s * z__[l - 1];
			z__[l - 1] = s * tst + c__ * z__[l - 1];

			d__[l - 1] = rt1;
			d__[l] = rt2;
			e[l - 1] = 0.;
			l += -2;
			if (l >= lend) {
				goto L90;
			}
			goto L140;
		}

		if (jtot == nmaxit) {
			goto L140;
		}
		++jtot;

		/*	Form shift. */

		g = (d__[l - 1] - p) / (e[l - 1] * 2.);
		r__ = dlapy2_(&g, &c_b10);
		g = d__[m] - p + e[l - 1] / (g + d_sign(&r__, &g));

		s = 1.;
		c__ = 1.;
		p = 0.;

		/*	Inner loop */

		lm1 = l - 1;
		i__1 = lm1;
		for (i__ = m; i__ <= i__1; ++i__) {
			f = s * e[i__];
			b = c__ * e[i__];
			dlartg_(&g, &f, &c__, &s, &r__);
			if (i__ != m) {
				e[i__ - 1] = r__;
			}
			g = d__[i__] - p;
			r__ = (d__[i__ + 1] - g) * s + c__ * 2. * b;
			p = s * r__;
			d__[i__] = g + p;
			g = c__ * r__ - b;

			/*		If eigenvectors are desired, then save rotations. */

			work[i__] = c__;
			work[*n - 1 + i__] = s;

			/* L120: */
		}

		/*	If eigenvectors are desired, then apply saved rotations. */

		mm = l - m + 1;
		/*
		dlasr_("R", "V", "F", n, &mm, &work[m], &work[*n - 1 + m], &
		z___ref(1, m), ldz);
		*/
		dlasr_("R", "V", "F", &c__1, &mm, &work[m], &work[*n - 1 + m],
			&z__[m], &c__1);

		d__[l] -= p;
		e[lm1] = g;
		goto L90;

		/*        Eigenvalue found. */

L130:
		d__[l] = p;

		--l;
		if (l >= lend) {
			goto L90;
		}
		goto L140;

	}

	/*     Undo scaling if necessary */

L140:
	if (iscale == 1) {
		i__1 = lendsv - lsv + 1;
		dlascl_("G", &c__0, &c__0, &ssfmax, &anorm, &i__1, &c__1,
			&d__[lsv], n, info);
		i__1 = lendsv - lsv;
		dlascl_("G", &c__0, &c__0, &ssfmax, &anorm, &i__1, &c__1, &e[lsv],
			n, info);
	} else if (iscale == 2) {
		i__1 = lendsv - lsv + 1;
		dlascl_("G", &c__0, &c__0, &ssfmin, &anorm, &i__1, &c__1,
			&d__[lsv], n, info);
		i__1 = lendsv - lsv;
		dlascl_("G", &c__0, &c__0, &ssfmin, &anorm, &i__1, &c__1, &e[lsv],
			n, info);
	}

	/*     Check for no convergence to an eigenvalue after a total
	of N*MAXIT iterations. */

	if (jtot < nmaxit) {
		goto L10;
	}
	i__1 = *n - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
		if (e[i__] != 0.) {
			++(*info);
		}
		/* L150: */
	}
	goto L190;

	/*     Order eigenvalues and eigenvectors. */

L160:

	/*        Use Selection Sort to minimize swaps of eigenvectors */

	i__1 = *n;
	for (ii = 2; ii <= i__1; ++ii) {
		i__ = ii - 1;
		k = i__;
		p = d__[i__];
		i__2 = *n;
		for (j = ii; j <= i__2; ++j) {
			if (d__[j] < p) {
				k = j;
				p = d__[j];
			}
			/* L170: */
		}
		if (k != i__) {
			d__[k] = d__[i__];
			d__[i__] = p;
			/* Taken out for TRLan
			dswap_(n, &z___ref(1, i__), &c__1, &z___ref(1, k), &c__1);
			*/
			p = z__[k];
			z__[k] = z__[i__];
			z__[i__] = p;
		}
		/* L180: */
	}

L190:
	return 0;
}	

void TRL::dsort2(int N, double *ARRAY1, double *ARRAY2)
{
	/*
	// Purpose
	// =======
	// Use quick sort algorithm to sort ARRAY1 and ARRAY2 in the increasing 
	// order of ARRAY1.
	//
	// Arguments
	// =========
	// N       (input) INTEGER
	//          On entry, specifies the size of the arrays.
	//
	// ARRAY1  (input/output) DOUBLE PRECISION ARRAY of LENGTH N
	//          On entry, contains the array to be sorted.
	//          On exit, contains the sorted array.
	//
	// ARRAY2  (input/output) DOUBLE PRECISION ARRAY of LENGTH N
	//          On entry, contains the array to be sorted.
	//          On exit, contains the sorted array.
	//
	// ..
	// .. local scalars ..
	*/
	int IGAP, I, J;
	double TEMP;
	/*
	// ..
	// .. executable statements ..
	*/
	IGAP = N / 2;
	while (IGAP > 0) {
		for (I = IGAP; I < N; I++) {
			J = I - IGAP;
			while (J >= 0) {
				if (ARRAY1[J] > ARRAY1[J + IGAP]) {
					TEMP = ARRAY1[J];
					ARRAY1[J] = ARRAY1[J + IGAP];
					ARRAY1[J + IGAP] = TEMP;
					TEMP = ARRAY2[J];
					ARRAY2[J] = ARRAY2[J + IGAP];
					ARRAY2[J + IGAP] = TEMP;
					J = J - IGAP;
				} else {
					break;
				}
			}
		}
		IGAP = IGAP / 2;
	}
	/*
	// .. end of dsort2_c_
	*/
}

void TRL::dsort2a(int N, double *ARRAY1, double *ARRAY2)
{
	/*
	// Purpose
	// =======
	// Use quick sort algorithm to sort ARRAY1 and ARRAY2 in the increasing 
	// order of abs(ARRAY1).
	//
	// Arguments
	// =========
	// N       (input) INTEGER
	//          On entry, specifies the size of the arrays.
	//
	// ARRAY1  (input/output) DOUBLE PRECISION ARRAY of LENGTH N
	//          On entry, contains the array to be sorted.
	//          On exit, contains the sorted array.
	//
	// ARRAY2  (input/output) DOUBLE PRECISION ARRAY of LENGTH N
	//          On entry, contains the array to be sorted.
	//          On exit, contains the sorted array.
	//
	// ..
	// .. local scalars ..
	*/
	int IGAP, I, J;
	double TEMP;
	/*
	// ..
	// .. executable statements ..
	*/
	IGAP = N / 2;
	while (IGAP > 0) {
		for (I = IGAP; I < N; I++) {
			J = I - IGAP;
			while (J >= 0) {
				if (fabs(ARRAY1[J]) > fabs(ARRAY1[J + IGAP])) {
					TEMP = ARRAY1[J];
					ARRAY1[J] = ARRAY1[J + IGAP];
					ARRAY1[J + IGAP] = TEMP;
					TEMP = ARRAY2[J];
					ARRAY2[J] = ARRAY2[J + IGAP];
					ARRAY2[J + IGAP] = TEMP;
					J = J - IGAP;
				} else {
					break;
				}
			}
		}
		IGAP = IGAP / 2;
	}
	/*
	// .. end of dsort2ac_
	*/
}

void TRL::dsort2s(int N, double s, double *ARRAY1, double *ARRAY2)
{
	/*
	// Purpose
	// =======
	// Use quick sort algorithm to sort ARRAY1 and ARRAY2 in the increasing 
	// order of abs(ARRAY1-s).
	//
	// Arguments
	// =========
	// N       (input) INTEGER
	//          On entry, specifies the size of the arrays.
	//
	// s       (input) DOUBLE PRECISION
	//          On entry, specifies the reference value s.
	//          
	// ARRAY1  (input/output) DOUBLE PRECISION ARRAY of LENGTH N
	//          On entry, contains the array to be sorted.
	//          On exit, contains the sorted array.
	//
	// ARRAY2  (input/output) DOUBLE PRECISION ARRAY of LENGTH N
	//          On entry, contains the array to be sorted.
	//          On exit, contains the sorted array.
	//
	// ..
	// .. local scalars ..
	*/
	int IGAP, I, J;
	double TEMP;
	/*
	// ..
	// .. executable statements ..
	*/
	IGAP = N / 2;
	while (IGAP > 0) {
		for (I = IGAP; I < N; I++) {
			J = I - IGAP;
			while (J >= 0) {
				if (fabs(ARRAY1[J] - s) > fabs(ARRAY1[J + IGAP] - s)) {
					TEMP = ARRAY1[J];
					ARRAY1[J] = ARRAY1[J + IGAP];
					ARRAY1[J + IGAP] = TEMP;
					TEMP = ARRAY2[J];
					ARRAY2[J] = ARRAY2[J + IGAP];
					ARRAY2[J + IGAP] = TEMP;
					J = J - IGAP;
				} else {
					break;
				}
			}
		}
		IGAP = IGAP / 2;
	}
	/*
	// .. end of dsort2ac_
	*/
}

void TRL::dsort2su_(int N, double s, double *ARRAY1, double *ARRAY2)
{
	/*
	// Purpose
	// =======
	// Use quick sort algorithm to sort ARRAY1 and ARRAY2 in the increasing order 
	// of ARRAY1-s if ARRAY1-s is non-negative. Negative ARRAY1-s are ordered after
	// those with non-negative ARRAY1-s and in the increasing order of ARRAY1.
	// 
	//
	// Arguments
	// =========
	// N       (input) INTEGER
	//          On entry, specifies the size of the arrays.
	//
	// ARRAY1  (input/output) DOUBLE PRECISION ARRAY of LENGTH N
	//          On entry, contains the array to be sorted.
	//          On exit, contains the sorted array.
	//
	// ARRAY2  (input/output) DOUBLE PRECISION ARRAY of LENGTH N
	//          On entry, contains the array to be sorted.
	//          On exit, contains the sorted array.
	//
	// ..
	// .. local scalars ..
	*/
	int IGAP, I, J;
	double TEMP, v1, v2, d1, d2, maxd;
	/*
	// ..
	// .. executable statements ..
	*/
	IGAP = N / 2;
	maxd = fabs(ARRAY1[0]);
	for (I = 1; I < N; I++) {
		if (maxd < fabs(ARRAY1[I])) {
			maxd = fabs(ARRAY1[I]);
		}
	}
	while (IGAP > 0) {
		for (I = IGAP; I < N; I++) {
			J = I - IGAP;
			while (J >= 0) {
				v1 = fabs(ARRAY1[J]);
				v2 = fabs(ARRAY1[J + IGAP]);
				d1 = v1 - s;
				d2 = v2 - s;
				if (d1 < 0.0) {
					d1 = maxd + v1;
				}
				if (d2 < 0.0) {
					d2 = maxd + v2;
				}
				if (d1 > d2) {
					TEMP = ARRAY1[J];
					ARRAY1[J] = ARRAY1[J + IGAP];
					ARRAY1[J + IGAP] = TEMP;
					TEMP = ARRAY2[J];
					ARRAY2[J] = ARRAY2[J + IGAP];
					ARRAY2[J + IGAP] = TEMP;
					J = J - IGAP;
				} else {
					break;
				}
			}
		}
		IGAP = IGAP / 2;
	}
	/*
	// .. end of dsort2ac_
	*/
}

void TRL::dsort2sd(int N, double s, double *ARRAY1, double *ARRAY2)
{
	/*
	// Purpose
	// =======
	// Use quick sort algorithm to sort ARRAY1 and ARRAY2 in the increasing order
	// of abs(ARRAY1-s) if ARRAY1-s is non-positive. Positive ARRAY1-s are ordered 
	// after those with non-positive ARRAY1-s and in the increasing order of -ARRAY1.
	//
	// Arguments
	// =========
	// N       (input) INTEGER
	//          On entry, specifies the size of the arrays.
	//
	// ARRAY1  (input/output) DOUBLE PRECISION ARRAY of LENGTH N
	//          On entry, contains the array to be sorted.
	//          On exit, contains the sorted array.
	//
	// ARRAY2  (input/output) DOUBLE PRECISION ARRAY of LENGTH N
	//          On entry, contains the array to be sorted.
	//          On exit, contains the sorted array.
	//
	// ..
	// .. local scalars ..
	*/
	int IGAP, I, J;
	double TEMP, v1, v2, d1, d2, maxd;
	/*
	// ..
	// .. executable statements ..
	*/
	IGAP = N / 2;
	maxd = fabs(ARRAY1[0]);
	for (I = 1; I < N; I++) {
		if (maxd < fabs(ARRAY1[I])) {
			maxd = fabs(ARRAY1[I]);
		}
	}
	maxd = maxd + 1.0;
	while (IGAP > 0) {
		for (I = IGAP; I < N; I++) {
			J = I - IGAP;
			while (J >= 0) {
				v1 = fabs(ARRAY1[J]);
				v2 = fabs(ARRAY1[J + IGAP]);
				d1 = s - v1;
				d2 = s - v2;
				if (d1 < 0.0) {
					d1 = s + (maxd - v1);
				}
				if (d2 < 0.0) {
					d2 = s + (maxd - v2);
				}
				if (d1 > d2) {
					TEMP = ARRAY1[J];
					ARRAY1[J] = ARRAY1[J + IGAP];
					ARRAY1[J + IGAP] = TEMP;
					TEMP = ARRAY2[J];
					ARRAY2[J] = ARRAY2[J + IGAP];
					ARRAY2[J + IGAP] = TEMP;
					J = J - IGAP;
				} else {
					break;
				}
			}
		}
		IGAP = IGAP / 2;
	}
	/*
	// .. end of dsort2ac_
	*/
}


//////////////////////////

double TRL::ran1(int &idum)
{
	const int IA=16807,IM=2147483647,IQ=127773,IR=2836,NTAB=32;
	const int NDIV=(1+(IM-1)/NTAB);
	const double EPS=3.0e-16,AM=1.0/IM,RNMX=(1.0-EPS);
	static int iy=0;
	static vector<int> iv(NTAB);
	int j,k;
	double temp;

	if (idum <= 0 || !iy) {
		if (-idum < 1) idum=1;
		else idum = -idum;
		for (j=NTAB+7;j>=0;j--) {
			k=idum/IQ;
			idum=IA*(idum-k*IQ)-IR*k;
			if (idum < 0) idum += IM;
			if (j < NTAB) iv[j] = idum;
		}
		iy=iv[0];
	}
	k=idum/IQ;
	idum=IA*(idum-k*IQ)-IR*k;
	if (idum < 0) idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = idum;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}

double TRL::ran2(int &idum)
{
	const int IM1=2147483563,IM2=2147483399;
	const int IA1=40014,IA2=40692,IQ1=53668,IQ2=52774;
	const int IR1=12211,IR2=3791,NTAB=32,IMM1=IM1-1;
	const int NDIV=1+IMM1/NTAB;
	const double EPS=3.0e-16,RNMX=1.0-EPS,AM=1.0/double(IM1);
	static int idum2=123456789,iy=0;
	static vector<int> iv(NTAB);
	int j,k;
	double temp;

	if (idum <= 0) {
		idum=(idum==0 ? 1 : -idum);
		idum2=idum;
		for (j=NTAB+7;j>=0;j--) {
			k=idum/IQ1;
			idum=IA1*(idum-k*IQ1)-k*IR1;
			if (idum < 0) idum += IM1;
			if (j < NTAB) iv[j] = idum;
		}
		iy=iv[0];
	}
	k=idum/IQ1;
	idum=IA1*(idum-k*IQ1)-k*IR1;
	if (idum < 0) idum += IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;
	if (idum2 < 0) idum2 += IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j] = idum;
	if (iy < 1) iy += IMM1;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}





void TRL::srand48_new()
{
	seed=time(NULL);
}

double TRL::drand48_new()
{
	return ran1(par::seed);
}

void TRL::indexx(vector<double> &arr, vector<int> &indx)
{
	const int M=7,NSTACK=50;
	int i,indxt,ir,j,k,jstack=-1,l=0;
	double a;
	vector<int> istack(NSTACK);

	int n=arr.size();
	ir=n-1;
	indx.resize(n);
	for (j=0;j<n;j++) indx[j]=j;
	for (;;) {
		if (ir-l < M) {
			for (j=l+1;j<=ir;j++) {
				indxt=indx[j];
				a=arr[indxt];
				for (i=j-1;i>=l;i--) {
					if (arr[indx[i]] <= a) break;
					indx[i+1]=indx[i];
				}
				indx[i+1]=indxt;
			}
			if (jstack < 0) break;
			ir=istack[jstack--];
			l=istack[jstack--];
		} else {
			k=(l+ir) >> 1;
			SWAP(indx[k],indx[l+1]);
			if (arr[indx[l]] > arr[indx[ir]]) {
				SWAP(indx[l],indx[ir]);
			}
			if (arr[indx[l+1]] > arr[indx[ir]]) {
				SWAP(indx[l+1],indx[ir]);
			}
			if (arr[indx[l]] > arr[indx[l+1]]) {
				SWAP(indx[l],indx[l+1]);
			}
			i=l+1;
			j=ir;
			indxt=indx[l+1];
			a=arr[indxt];
			for (;;) {
				do i++; while (arr[indx[i]] < a);
				do j--; while (arr[indx[j]] > a);
				if (j < i) break;
				SWAP(indx[i],indx[j]);
			}
			indx[l+1]=indx[j];
			indx[j]=indxt;
			jstack += 2;
			if (jstack >= NSTACK) cerr<<"NSTACK too small in indexx.";
			if (ir-i+1 >= j-l) {
				istack[jstack]=ir;
				istack[jstack-1]=i;
				ir=j-1;
			} else {
				istack[jstack]=j-1;
				istack[jstack-1]=l;
				l=i;
			}
		}
	}
}








