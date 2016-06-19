#include "trl.h"


//////////////////////////////

//dsyev_
//dpotrf_
//dtrtrs_
//dsytrd_
//dorgtr_
//dstein_


///////////////////////////

//lsame_
//dlamch_
//dlascl_
//ilaenv_
//xerbla_
//dsterf_
//dlansy_
//dtrsm_
//dsyrk_
//dpotf2_
//dorgql_
//dorgqr_
//ieeeck_
//iparmq_
//dlae2_
//dlapy2_
//dlanst_
//dlasrt_
//dlassq_
//dorg2l_
//dlarfb_
//dlarft_
//dorg2r_
//dlarf_
//dtrmm_
//dtrmv_
//dger_

//dsteqr_
//dlasr_
//dswap_
//dlaev2_
//dlapy2_
//dlaset_
//dlartg_

//dnrm2_
//dasum_
//dlagtf_
//idamax_
//dlagts_
//dlarnv_
//dlaruv_

//dsytd2_
//dsyr2k_
//dlatrd_

//dsyr2_
//dsymv_
//dlarfg_


//////////////////////////

//pow_di
//s_wsfe
//do_fio
//e_wsfe
//s_copy
//s_cmp
//sqrt
//d_sign
//log
//i_nint





























/////////////////////////////////////////////
///////////////////////////////////////////////////
/////////////////////////////////////////////////////////



logical_ TRL::lsame_(char *ca, char *cb)
{
	/*  -- LAPACK auxiliary routine (version 3.1) --   
	Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..   
	November 2006   


	Purpose   
	=======   

	LSAME returns .TRUE. if CA is the same letter as CB regardless of   
	case.   

	Arguments   
	=========   

	CA      (input) CHARACTER*1   

	CB      (input) CHARACTER*1   
	CA and CB specify the single characters to be compared.   

	=====================================================================   


	Test if the characters are equal */
	/* System generated locals */
	logical_ ret_val;
	/* Local variables */
	static integer_ inta, intb, zcode;


	ret_val = *(unsigned char *)ca == *(unsigned char *)cb;
	if (ret_val) {
		return ret_val;
	}

	/*     Now test for equivalence if both characters are alphabetic. */

	zcode = 'Z';

	/*     Use 'Z' rather than 'A' so that ASCII can be detected on Prime   
	machines, on which ICHAR returns a value with bit 8 set.   
	ICHAR('A') on Prime machines returns 193 which is the same as   
	ICHAR('A') on an EBCDIC machine. */

	inta = *(unsigned char *)ca;
	intb = *(unsigned char *)cb;

	if (zcode == 90 || zcode == 122) {

		/*        ASCII is assumed - ZCODE is the ASCII code of either lower or   
		upper case 'Z'. */

		if (inta >= 97 && inta <= 122) {
			inta += -32;
		}
		if (intb >= 97 && intb <= 122) {
			intb += -32;
		}

	} else if (zcode == 233 || zcode == 169) {

		/*        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or   
		upper case 'Z'. */

		if (inta >= 129 && inta <= 137 || inta >= 145 && inta <= 153 || inta 
			>= 162 && inta <= 169) {
				inta += 64;
		}
		if (intb >= 129 && intb <= 137 || intb >= 145 && intb <= 153 || intb 
			>= 162 && intb <= 169) {
				intb += 64;
		}

	} else if (zcode == 218 || zcode == 250) {

		/*        ASCII is assumed, on Prime machines - ZCODE is the ASCII code   
		plus 128 of either lower or upper case 'Z'. */

		if (inta >= 225 && inta <= 250) {
			inta += -32;
		}
		if (intb >= 225 && intb <= 250) {
			intb += -32;
		}
	}
	ret_val = inta == intb;

	/*     RETURN   

	End of LSAME */

	return ret_val;
} /* lsame_ */











//
doublereal_ TRL::dlamch_(char *cmach)
{
	static integer_ c__1 = 1;
	static doublereal_ c_b32 = 0.;
	/* Initialized data */

	static logical_ first = TRUE_;

	/* System generated locals */
	integer_ i__1;
	doublereal_ ret_val;

	/* Builtin functions */
	//double pow_di(doublereal_ *, integer_ *);

	/* Local variables */
	static doublereal_ t;
	static integer_ it;
	static doublereal_ rnd, eps, base;
	static integer_ beta;
	static doublereal_ emin, prec, emax;
	static integer_ imin, imax;
	static logical_ lrnd;
	static doublereal_ rmin, rmax, rmach;
	//extern logical_ lsame_(char *, char *);
	static doublereal_ small, sfmin;
	//extern /* Subroutine */ int dlamc2_(integer_ *, integer_ *, logical_ *, 
	//	doublereal_ *, integer_ *, doublereal_ *, integer_ *, doublereal_ *);


	/*  -- LAPACK auxiliary routine (version 3.1) --   
	Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..   
	November 2006   


	Purpose   
	=======   

	DLAMCH determines double precision machine parameters.   

	Arguments   
	=========   

	CMACH   (input) CHARACTER*1   
	Specifies the value to be returned by DLAMCH:   
	= 'E' or 'e',   DLAMCH := eps   
	= 'S' or 's ,   DLAMCH := sfmin   
	= 'B' or 'b',   DLAMCH := base   
	= 'P' or 'p',   DLAMCH := eps*base   
	= 'N' or 'n',   DLAMCH := t   
	= 'R' or 'r',   DLAMCH := rnd   
	= 'M' or 'm',   DLAMCH := emin   
	= 'U' or 'u',   DLAMCH := rmin   
	= 'L' or 'l',   DLAMCH := emax   
	= 'O' or 'o',   DLAMCH := rmax   

	where   

	eps   = relative machine precision   
	sfmin = safe minimum, such that 1/sfmin does not overflow   
	base  = base of the machine   
	prec  = eps*base   
	t     = number of (base) digits in the mantissa   
	rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise   
	emin  = minimum exponent before (gradual) underflow   
	rmin  = underflow threshold - base**(emin-1)   
	emax  = largest exponent before overflow   
	rmax  = overflow threshold  - (base**emax)*(1-eps)   

	===================================================================== */


	if (first) {
		dlamc2_(&beta, &it, &lrnd, &eps, &imin, &rmin, &imax, &rmax);
		base = (doublereal_) beta;
		t = (doublereal_) it;
		if (lrnd) {
			rnd = 1.;
			i__1 = 1 - it;
			eps = pow_di(&base, &i__1) / 2;
		} else {
			rnd = 0.;
			i__1 = 1 - it;
			eps = pow_di(&base, &i__1);
		}
		prec = eps * base;
		emin = (doublereal_) imin;
		emax = (doublereal_) imax;
		sfmin = rmin;
		small = 1. / rmax;
		if (small >= sfmin) {

			/*           Use SMALL plus a bit, to avoid the possibility of rounding   
			causing overflow when computing  1/sfmin. */

			sfmin = small * (eps + 1.);
		}
	}

	if (lsame_(cmach, "E")) {
		rmach = eps;
	} else if (lsame_(cmach, "S")) {
		rmach = sfmin;
	} else if (lsame_(cmach, "B")) {
		rmach = base;
	} else if (lsame_(cmach, "P")) {
		rmach = prec;
	} else if (lsame_(cmach, "N")) {
		rmach = t;
	} else if (lsame_(cmach, "R")) {
		rmach = rnd;
	} else if (lsame_(cmach, "M")) {
		rmach = emin;
	} else if (lsame_(cmach, "U")) {
		rmach = rmin;
	} else if (lsame_(cmach, "L")) {
		rmach = emax;
	} else if (lsame_(cmach, "O")) {
		rmach = rmax;
	}

	ret_val = rmach;
	first = FALSE_;
	return ret_val;

	/*     End of DLAMCH */

}

//
int TRL::dlamc1_(integer_ *beta, integer_ *t, logical_ *rnd, logical_ 
  *ieee1)
{
	static integer_ c__1 = 1;
	static doublereal_ c_b32 = 0.;
	/* Initialized data */

	static logical_ first = TRUE_;

	/* System generated locals */
	doublereal_ d__1, d__2;

	/* Local variables */
	static doublereal_ a, b, c__, f, t1, t2;
	static integer_ lt;
	static doublereal_ one, qtr;
	static logical_ lrnd;
	static integer_ lbeta;
	static doublereal_ savec;
	//extern doublereal_ dlamc3_(doublereal_ *, doublereal_ *);
	static logical_ lieee1;


	/*  -- LAPACK auxiliary routine (version 3.1) --   
	Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..   
	November 2006   


	Purpose   
	=======   

	DLAMC1 determines the machine parameters given by BETA, T, RND, and   
	IEEE1.   

	Arguments   
	=========   

	BETA    (output) INTEGER   
	The base of the machine.   

	T       (output) INTEGER   
	The number of ( BETA ) digits in the mantissa.   

	RND     (output) LOGICAL   
	Specifies whether proper rounding  ( RND = .TRUE. )  or   
	chopping  ( RND = .FALSE. )  occurs in addition. This may not   
	be a reliable guide to the way in which the machine performs   
	its arithmetic.   

	IEEE1   (output) LOGICAL   
	Specifies whether rounding appears to be done in the IEEE   
	'round to nearest' style.   

	Further Details   
	===============   

	The routine is based on the routine  ENVRON  by Malcolm and   
	incorporates suggestions by Gentleman and Marovich. See   

	Malcolm M. A. (1972) Algorithms to reveal properties of   
	floating-point arithmetic. Comms. of the ACM, 15, 949-951.   

	Gentleman W. M. and Marovich S. B. (1974) More on algorithms   
	that reveal properties of floating point arithmetic units.   
	Comms. of the ACM, 17, 276-277.   

	===================================================================== */


	if (first) {
		one = 1.;

		/*        LBETA,  LIEEE1,  LT and  LRND  are the  local values  of  BETA,   
		IEEE1, T and RND.   

		Throughout this routine  we use the function  DLAMC3  to ensure   
		that relevant values are  stored and not held in registers,  or   
		are not affected by optimizers.   

		Compute  a = 2.0**m  with the  smallest positive integer_ m such   
		that   

		fl( a + 1.0 ) = a. */

		a = 1.;
		c__ = 1.;

		/* +       WHILE( C.EQ.ONE )LOOP */
L10:
		if (c__ == one) {
			a *= 2;
			c__ = dlamc3_(&a, &one);
			d__1 = -a;
			c__ = dlamc3_(&c__, &d__1);
			goto L10;
		}
		/* +       END WHILE   

		Now compute  b = 2.0**m  with the smallest positive integer_ m   
		such that   

		fl( a + b ) .gt. a. */

		b = 1.;
		c__ = dlamc3_(&a, &b);

		/* +       WHILE( C.EQ.A )LOOP */
L20:
		if (c__ == a) {
			b *= 2;
			c__ = dlamc3_(&a, &b);
			goto L20;
		}
		/* +       END WHILE   

		Now compute the base.  a and c  are neighbouring floating point   
		numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and so   
		their difference is beta. Adding 0.25 to c is to ensure that it   
		is truncated to beta and not ( beta - 1 ). */

		qtr = one / 4;
		savec = c__;
		d__1 = -a;
		c__ = dlamc3_(&c__, &d__1);
		lbeta = (integer_) (c__ + qtr);

		/*        Now determine whether rounding or chopping occurs,  by adding a   
		bit  less  than  beta/2  and a  bit  more  than  beta/2  to  a. */

		b = (doublereal_) lbeta;
		d__1 = b / 2;
		d__2 = -b / 100;
		f = dlamc3_(&d__1, &d__2);
		c__ = dlamc3_(&f, &a);
		if (c__ == a) {
			lrnd = TRUE_;
		} else {
			lrnd = FALSE_;
		}
		d__1 = b / 2;
		d__2 = b / 100;
		f = dlamc3_(&d__1, &d__2);
		c__ = dlamc3_(&f, &a);
		if (lrnd && c__ == a) {
			lrnd = FALSE_;
		}

		/*        Try and decide whether rounding is done in the  IEEE  'round to   
		nearest' style. B/2 is half a unit in the last place of the two   
		numbers A and SAVEC. Furthermore, A is even, i.e. has last  bit   
		zero, and SAVEC is odd. Thus adding B/2 to A should not  change   
		A, but adding B/2 to SAVEC should change SAVEC. */

		d__1 = b / 2;
		t1 = dlamc3_(&d__1, &a);
		d__1 = b / 2;
		t2 = dlamc3_(&d__1, &savec);
		lieee1 = t1 == a && t2 > savec && lrnd;

		/*        Now find  the  mantissa, t.  It should  be the  integer_ part of   
		log to the base beta of a,  however it is safer to determine  t   
		by powering.  So we find t as the smallest positive integer_ for   
		which   

		fl( beta**t + 1.0 ) = 1.0. */

		lt = 0;
		a = 1.;
		c__ = 1.;

		/* +       WHILE( C.EQ.ONE )LOOP */
L30:
		if (c__ == one) {
			++lt;
			a *= lbeta;
			c__ = dlamc3_(&a, &one);
			d__1 = -a;
			c__ = dlamc3_(&c__, &d__1);
			goto L30;
		}
		/* +       END WHILE */

	}

	*beta = lbeta;
	*t = lt;
	*rnd = lrnd;
	*ieee1 = lieee1;
	first = FALSE_;
	return 0;

	/*     End of DLAMC1 */

}

//
int TRL::dlamc2_(integer_ *beta, integer_ *t, logical_ *rnd, 
  doublereal_ *eps, integer_ *emin, doublereal_ *rmin, integer_ *emax, 
  doublereal_ *rmax)
{
	static integer_ c__1 = 1;
	static doublereal_ c_b32 = 0.;
	/* Initialized data */

	static logical_ first = TRUE_;
	static logical_ iwarn = FALSE_;

	/* Format strings */
	static char fmt_9999[] = "(//\002 WARNING. The value EMIN may be incorre"
		"ct:-\002,\002  EMIN = \002,i8,/\002 If, after inspection, the va"
		"lue EMIN looks\002,\002 acceptable please comment out \002,/\002"
		" the IF block as marked within the code of routine\002,\002 DLAM"
		"C2,\002,/\002 otherwise supply EMIN explicitly.\002,/)";

	/* System generated locals */
	integer_ i__1;
	doublereal_ d__1, d__2, d__3, d__4, d__5;

	/* Builtin functions */
	//double pow_di(doublereal_ *, integer_ *);
	//integer_ s_wsfe(cilist_ *), do_fio(integer_ *, char *, ftnlen_), e_wsfe(void);

	/* Local variables */
	static doublereal_ a, b, c__;
	static integer_ i__, lt;
	static doublereal_ one, two;
	static logical_ ieee;
	static doublereal_ half;
	static logical_ lrnd;
	static doublereal_ leps, zero;
	static integer_ lbeta;
	static doublereal_ rbase;
	static integer_ lemin, lemax, gnmin;
	static doublereal_ small;
	static integer_ gpmin;
	static doublereal_ third, lrmin, lrmax, sixth;
	//extern /* Subroutine */ int dlamc1_(integer_ *, integer_ *, logical_ *, 
	//	logical_ *);
	//extern doublereal_ dlamc3_(doublereal_ *, doublereal_ *);
	static logical_ lieee1;
	//extern /* Subroutine */ int dlamc4_(integer_ *, doublereal_ *, integer_ *), 
	//	dlamc5_(integer_ *, integer_ *, integer_ *, logical_ *, integer_ *, 
	//	doublereal_ *);
	static integer_ ngnmin, ngpmin;

	/* Fortran I/O blocks */
	static cilist_ io___58 = { 0, 6, 0, fmt_9999, 0 };



	/*  -- LAPACK auxiliary routine (version 3.1) --   
	Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..   
	November 2006   


	Purpose   
	=======   

	DLAMC2 determines the machine parameters specified in its argument   
	list.   

	Arguments   
	=========   

	BETA    (output) INTEGER   
	The base of the machine.   

	T       (output) INTEGER   
	The number of ( BETA ) digits in the mantissa.   

	RND     (output) LOGICAL   
	Specifies whether proper rounding  ( RND = .TRUE. )  or   
	chopping  ( RND = .FALSE. )  occurs in addition. This may not   
	be a reliable guide to the way in which the machine performs   
	its arithmetic.   

	EPS     (output) DOUBLE PRECISION   
	The smallest positive number such that   

	fl( 1.0 - EPS ) .LT. 1.0,   

	where fl denotes the computed value.   

	EMIN    (output) INTEGER   
	The minimum exponent before (gradual) underflow occurs.   

	RMIN    (output) DOUBLE PRECISION   
	The smallest normalized number for the machine, given by   
	BASE**( EMIN - 1 ), where  BASE  is the floating point value   
	of BETA.   

	EMAX    (output) INTEGER   
	The maximum exponent before overflow occurs.   

	RMAX    (output) DOUBLE PRECISION   
	The largest positive number for the machine, given by   
	BASE**EMAX * ( 1 - EPS ), where  BASE  is the floating point   
	value of BETA.   

	Further Details   
	===============   

	The computation of  EPS  is based on a routine PARANOIA by   
	W. Kahan of the University of California at Berkeley.   

	===================================================================== */


	if (first) {
		zero = 0.;
		one = 1.;
		two = 2.;

		/*        LBETA, LT, LRND, LEPS, LEMIN and LRMIN  are the local values of   
		BETA, T, RND, EPS, EMIN and RMIN.   

		Throughout this routine  we use the function  DLAMC3  to ensure   
		that relevant values are stored  and not held in registers,  or   
		are not affected by optimizers.   

		DLAMC1 returns the parameters  LBETA, LT, LRND and LIEEE1. */

		dlamc1_(&lbeta, &lt, &lrnd, &lieee1);

		/*        Start to find EPS. */

		b = (doublereal_) lbeta;
		i__1 = -lt;
		a = pow_di(&b, &i__1);
		leps = a;

		/*        Try some tricks to see whether or not this is the correct  EPS. */

		b = two / 3;
		half = one / 2;
		d__1 = -half;
		sixth = dlamc3_(&b, &d__1);
		third = dlamc3_(&sixth, &sixth);
		d__1 = -half;
		b = dlamc3_(&third, &d__1);
		b = dlamc3_(&b, &sixth);
		b = abs(b);
		if (b < leps) {
			b = leps;
		}

		leps = 1.;

		/* +       WHILE( ( LEPS.GT.B ).AND.( B.GT.ZERO ) )LOOP */
L10:
		if (leps > b && b > zero) {
			leps = b;
			d__1 = half * leps;
			/* Computing 5th power */
			d__3 = two, d__4 = d__3, d__3 *= d__3;
			/* Computing 2nd power */
			d__5 = leps;
			d__2 = d__4 * (d__3 * d__3) * (d__5 * d__5);
			c__ = dlamc3_(&d__1, &d__2);
			d__1 = -c__;
			c__ = dlamc3_(&half, &d__1);
			b = dlamc3_(&half, &c__);
			d__1 = -b;
			c__ = dlamc3_(&half, &d__1);
			b = dlamc3_(&half, &c__);
			goto L10;
		}
		/* +       END WHILE */

		if (a < leps) {
			leps = a;
		}

		/*        Computation of EPS complete.   

		Now find  EMIN.  Let A = + or - 1, and + or - (1 + BASE**(-3)).   
		Keep dividing  A by BETA until (gradual) underflow occurs. This   
		is detected when we cannot recover the previous A. */

		rbase = one / lbeta;
		small = one;
		for (i__ = 1; i__ <= 3; ++i__) {
			d__1 = small * rbase;
			small = dlamc3_(&d__1, &zero);
			/* L20: */
		}
		a = dlamc3_(&one, &small);
		dlamc4_(&ngpmin, &one, &lbeta);
		d__1 = -one;
		dlamc4_(&ngnmin, &d__1, &lbeta);
		dlamc4_(&gpmin, &a, &lbeta);
		d__1 = -a;
		dlamc4_(&gnmin, &d__1, &lbeta);
		ieee = FALSE_;

		if (ngpmin == ngnmin && gpmin == gnmin) {
			if (ngpmin == gpmin) {
				lemin = ngpmin;
				/*            ( Non twos-complement machines, no gradual underflow;   
				e.g.,  VAX ) */
			} else if (gpmin - ngpmin == 3) {
				lemin = ngpmin - 1 + lt;
				ieee = TRUE_;
				/*            ( Non twos-complement machines, with gradual underflow;   
				e.g., IEEE standard followers ) */
			} else {
				lemin = min(ngpmin,gpmin);
				/*            ( A guess; no known machine ) */
				iwarn = TRUE_;
			}

		} else if (ngpmin == gpmin && ngnmin == gnmin) {
			if ((i__1 = ngpmin - ngnmin, abs(i__1)) == 1) {
				lemin = max(ngpmin,ngnmin);
				/*            ( Twos-complement machines, no gradual underflow;   
				e.g., CYBER 205 ) */
			} else {
				lemin = min(ngpmin,ngnmin);
				/*            ( A guess; no known machine ) */
				iwarn = TRUE_;
			}

		} else if ((i__1 = ngpmin - ngnmin, abs(i__1)) == 1 && gpmin == gnmin)
		{
			if (gpmin - min(ngpmin,ngnmin) == 3) {
				lemin = max(ngpmin,ngnmin) - 1 + lt;
				/*            ( Twos-complement machines with gradual underflow;   
				no known machine ) */
			} else {
				lemin = min(ngpmin,ngnmin);
				/*            ( A guess; no known machine ) */
				iwarn = TRUE_;
			}

		} else {
			/* Computing MIN */
			i__1 = min(ngpmin,ngnmin), i__1 = min(i__1,gpmin);
			lemin = min(i__1,gnmin);
			/*         ( A guess; no known machine ) */
			iwarn = TRUE_;
		}
		first = FALSE_;
		/* **   
		Comment out this if block if EMIN is ok */
		//if (iwarn) {
		//	first = TRUE_;
		//	s_wsfe(&io___58);
		//	do_fio(&c__1, (char *)&lemin, (ftnlen_)sizeof(integer_));
		//	e_wsfe();
		//}
		/* **   

		Assume IEEE arithmetic if we found denormalised  numbers above,   
		or if arithmetic seems to round in the  IEEE style,  determined   
		in routine DLAMC1. A true IEEE machine should have both  things   
		true; however, faulty machines may have one or the other. */

		ieee = ieee || lieee1;

		/*        Compute  RMIN by successive division by  BETA. We could compute   
		RMIN as BASE**( EMIN - 1 ),  but some machines underflow during   
		this computation. */

		lrmin = 1.;
		i__1 = 1 - lemin;
		for (i__ = 1; i__ <= i__1; ++i__) {
			d__1 = lrmin * rbase;
			lrmin = dlamc3_(&d__1, &zero);
			/* L30: */
		}

		/*        Finally, call DLAMC5 to compute EMAX and RMAX. */

		dlamc5_(&lbeta, &lt, &lemin, &ieee, &lemax, &lrmax);
	}

	*beta = lbeta;
	*t = lt;
	*rnd = lrnd;
	*eps = leps;
	*emin = lemin;
	*rmin = lrmin;
	*emax = lemax;
	*rmax = lrmax;

	return 0;


	/*     End of DLAMC2 */

} 

//
doublereal_ TRL::dlamc3_(doublereal_ *a, doublereal_ *b)
{
	static integer_ c__1 = 1;
	static doublereal_ c_b32 = 0.;
	/* System generated locals */
	doublereal_ ret_val;


	/*  -- LAPACK auxiliary routine (version 3.1) --   
	Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..   
	November 2006   


	Purpose   
	=======   

	DLAMC3  is intended to force  A  and  B  to be stored prior to doing   
	the addition of  A  and  B ,  for use in situations where optimizers   
	might hold one of these in a register.   

	Arguments   
	=========   

	A       (input) DOUBLE PRECISION   
	B       (input) DOUBLE PRECISION   
	The values A and B.   

	===================================================================== */


	ret_val = *a + *b;

	return ret_val;

	/*     End of DLAMC3 */

} 
//
int TRL::dlamc4_(integer_ *emin, doublereal_ *start, integer_ *base)
{
	static integer_ c__1 = 1;
	static doublereal_ c_b32 = 0.;
	/* System generated locals */
	integer_ i__1;
	doublereal_ d__1;

	/* Local variables */
	static doublereal_ a;
	static integer_ i__;
	static doublereal_ b1, b2, c1, c2, d1, d2, one, zero, rbase;
	//extern doublereal_ dlamc3_(doublereal_ *, doublereal_ *);


	/*  -- LAPACK auxiliary routine (version 3.1) --   
	Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..   
	November 2006   


	Purpose   
	=======   

	DLAMC4 is a service routine for DLAMC2.   

	Arguments   
	=========   

	EMIN    (output) INTEGER   
	The minimum exponent before (gradual) underflow, computed by   
	setting A = START and dividing by BASE until the previous A   
	can not be recovered.   

	START   (input) DOUBLE PRECISION   
	The starting point for determining EMIN.   

	BASE    (input) INTEGER   
	The base of the machine.   

	===================================================================== */


	a = *start;
	one = 1.;
	rbase = one / *base;
	zero = 0.;
	*emin = 1;
	d__1 = a * rbase;
	b1 = dlamc3_(&d__1, &zero);
	c1 = a;
	c2 = a;
	d1 = a;
	d2 = a;
	/* +    WHILE( ( C1.EQ.A ).AND.( C2.EQ.A ).AND.   
	$       ( D1.EQ.A ).AND.( D2.EQ.A )      )LOOP */
L10:
	if (c1 == a && c2 == a && d1 == a && d2 == a) {
		--(*emin);
		a = b1;
		d__1 = a / *base;
		b1 = dlamc3_(&d__1, &zero);
		d__1 = b1 * *base;
		c1 = dlamc3_(&d__1, &zero);
		d1 = zero;
		i__1 = *base;
		for (i__ = 1; i__ <= i__1; ++i__) {
			d1 += b1;
			/* L20: */
		}
		d__1 = a * rbase;
		b2 = dlamc3_(&d__1, &zero);
		d__1 = b2 / rbase;
		c2 = dlamc3_(&d__1, &zero);
		d2 = zero;
		i__1 = *base;
		for (i__ = 1; i__ <= i__1; ++i__) {
			d2 += b2;
			/* L30: */
		}
		goto L10;
	}
	/* +    END WHILE */

	return 0;

	/*     End of DLAMC4 */

}
//
int TRL::dlamc5_(integer_ *beta, integer_ *p, integer_ *emin, 
  logical_ *ieee, integer_ *emax, doublereal_ *rmax)
{
	static integer_ c__1 = 1;
	static doublereal_ c_b32 = 0.;
	/* System generated locals */
	integer_ i__1;
	doublereal_ d__1;

	/* Local variables */
	static integer_ i__;
	static doublereal_ y, z__;
	static integer_ try__, lexp;
	static doublereal_ oldy;
	static integer_ uexp, nbits;
	//extern doublereal_ dlamc3_(doublereal_ *, doublereal_ *);
	static doublereal_ recbas;
	static integer_ exbits, expsum;


	/*  -- LAPACK auxiliary routine (version 3.1) --   
	Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..   
	November 2006   


	Purpose   
	=======   

	DLAMC5 attempts to compute RMAX, the largest machine floating-point   
	number, without overflow.  It assumes that EMAX + abs(EMIN) sum   
	approximately to a power of 2.  It will fail on machines where this   
	assumption does not hold, for example, the Cyber 205 (EMIN = -28625,   
	EMAX = 28718).  It will also fail if the value supplied for EMIN is   
	too large (i.e. too close to zero), probably with overflow.   

	Arguments   
	=========   

	BETA    (input) INTEGER   
	The base of floating-point arithmetic.   

	P       (input) INTEGER   
	The number of base BETA digits in the mantissa of a   
	floating-point value.   

	EMIN    (input) INTEGER   
	The minimum exponent before (gradual) underflow.   

	IEEE    (input) LOGICAL   
	A logical_ flag_ specifying whether or not the arithmetic   
	system is thought to comply with the IEEE standard.   

	EMAX    (output) INTEGER   
	The largest exponent before overflow   

	RMAX    (output) DOUBLE PRECISION   
	The largest machine floating-point number.   

	=====================================================================   


	First compute LEXP and UEXP, two powers of 2 that bound   
	abs(EMIN). We then assume that EMAX + abs(EMIN) will sum   
	approximately to the bound that is closest to abs(EMIN).   
	(EMAX is the exponent of the required number RMAX). */

	lexp = 1;
	exbits = 1;
L10:
	try__ = lexp << 1;
	if (try__ <= -(*emin)) {
		lexp = try__;
		++exbits;
		goto L10;
	}
	if (lexp == -(*emin)) {
		uexp = lexp;
	} else {
		uexp = try__;
		++exbits;
	}

	/*     Now -LEXP is less than or equal to EMIN, and -UEXP is greater   
	than or equal to EMIN. EXBITS is the number of bits needed to   
	store the exponent. */

	if (uexp + *emin > -lexp - *emin) {
		expsum = lexp << 1;
	} else {
		expsum = uexp << 1;
	}

	/*     EXPSUM is the exponent range, approximately equal to   
	EMAX - EMIN + 1 . */

	*emax = expsum + *emin - 1;
	nbits = exbits + 1 + *p;

	/*     NBITS is the total number of bits needed to store a   
	floating-point number. */

	if (nbits % 2 == 1 && *beta == 2) {

		/*        Either there are an odd number of bits used to store a   
		floating-point number, which is unlikely, or some bits are   
		not used in the representation of numbers, which is possible,   
		(e.g. Cray machines) or the mantissa has an implicit bit,   
		(e.g. IEEE machines, Dec Vax machines), which is perhaps the   
		most likely. We have to assume the last alternative.   
		If this is true, then we need to reduce EMAX by one because   
		there must be some way of representing zero in an implicit-bit   
		system. On machines like Cray, we are reducing EMAX by one   
		unnecessarily. */

		--(*emax);
	}

	if (*ieee) {

		/*        Assume we are on an IEEE machine which reserves one exponent   
		for infinity and NaN. */

		--(*emax);
	}

	/*     Now create RMAX, the largest machine number, which should   
	be equal to (1.0 - BETA**(-P)) * BETA**EMAX .   

	First compute 1.0 - BETA**(-P), being careful that the   
	result is less than 1.0 . */

	recbas = 1. / *beta;
	z__ = *beta - 1.;
	y = 0.;
	i__1 = *p;
	for (i__ = 1; i__ <= i__1; ++i__) {
		z__ *= recbas;
		if (y < 1.) {
			oldy = y;
		}
		y = dlamc3_(&y, &z__);
		/* L20: */
	}
	if (y >= 1.) {
		y = oldy;
	}

	/*     Now multiply by BETA**EMAX to get RMAX. */

	i__1 = *emax;
	for (i__ = 1; i__ <= i__1; ++i__) {
		d__1 = y * *beta;
		y = dlamc3_(&d__1, &c_b32);
		/* L30: */
	}

	*rmax = y;
	return 0;

	/*     End of DLAMC5 */

} /* dlamc5_ */


//
int TRL::dlascl_(char *type__, integer_ *kl, integer_ *ku, 
			doublereal_ *cfrom, doublereal_ *cto, integer_ *m, integer_ *n, 
			doublereal_ *a, integer_ *lda, integer_ *info)
{
	/*  -- LAPACK auxiliary routine (version 3.1) --   
	Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..   
	November 2006   


	Purpose   
	=======   

	DLASCL multiplies the M by N trl_real_ matrix A by the trl_real_ scalar   
	CTO/CFROM.  This is done without over/underflow as long as the final   
	result CTO*A(I,J)/CFROM does not over/underflow. TYPE specifies that   
	A may be full, upper triangular, lower triangular, upper Hessenberg,   
	or banded.   

	Arguments   
	=========   

	TYPE    (input) CHARACTER*1   
	TYPE indices the storage type of the input matrix.   
	= 'G':  A is a full matrix.   
	= 'L':  A is a lower triangular matrix.   
	= 'U':  A is an upper triangular matrix.   
	= 'H':  A is an upper Hessenberg matrix.   
	= 'B':  A is a symmetric band matrix with lower bandwidth KL   
	and upper bandwidth KU and with the only the lower   
	half stored.   
	= 'Q':  A is a symmetric band matrix with lower bandwidth KL   
	and upper bandwidth KU and with the only the upper   
	half stored.   
	= 'Z':  A is a band matrix with lower bandwidth KL and upper   
	bandwidth KU.   

	KL      (input) INTEGER   
	The lower bandwidth of A.  Referenced only if TYPE = 'B',   
	'Q' or 'Z'.   

	KU      (input) INTEGER   
	The upper bandwidth of A.  Referenced only if TYPE = 'B',   
	'Q' or 'Z'.   

	CFROM   (input) DOUBLE PRECISION   
	CTO     (input) DOUBLE PRECISION   
	The matrix A is multiplied by CTO/CFROM. A(I,J) is computed   
	without over/underflow if the final result CTO*A(I,J)/CFROM   
	can be represented without over/underflow.  CFROM must be   
	nonzero.   

	M       (input) INTEGER   
	The number of rows of the matrix A.  M >= 0.   

	N       (input) INTEGER   
	The number of columns of the matrix A.  N >= 0.   

	A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
	The matrix to be multiplied by CTO/CFROM.  See TYPE for the   
	storage type.   

	LDA     (input) INTEGER   
	The leading dimension of the array A.  LDA >= max(1,M).   

	INFO    (output) INTEGER   
	0  - successful exit   
	<0 - if INFO = -i, the i-th argument had an illegal value.   

	=====================================================================   


	Test the input arguments   

	Parameter adjustments */
	/* System generated locals */
	integer_ a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
	/* Local variables */
	static integer_ i__, j, k1, k2, k3, k4;
	static doublereal_ mul, cto1;
	static logical_ done;
	static doublereal_ ctoc;
	//extern logical_ lsame_(char *, char *);
	static integer_ itype;
	static doublereal_ cfrom1;
	//extern doublereal_ dlamch_(char *);
	static doublereal_ cfromc;
	//extern /* Subroutine */ int xerbla_(char *, integer_ *);
	static doublereal_ bignum, smlnum;

	a_dim1 = *lda;
	a_offset = 1 + a_dim1;
	a -= a_offset;

	/* Function Body */
	*info = 0;

	if (lsame_(type__, "G")) {
		itype = 0;
	} else if (lsame_(type__, "L")) {
		itype = 1;
	} else if (lsame_(type__, "U")) {
		itype = 2;
	} else if (lsame_(type__, "H")) {
		itype = 3;
	} else if (lsame_(type__, "B")) {
		itype = 4;
	} else if (lsame_(type__, "Q")) {
		itype = 5;
	} else if (lsame_(type__, "Z")) {
		itype = 6;
	} else {
		itype = -1;
	}

	if (itype == -1) {
		*info = -1;
	} else if (*cfrom == 0.) {
		*info = -4;
	} else if (*m < 0) {
		*info = -6;
	} else if (*n < 0 || itype == 4 && *n != *m || itype == 5 && *n != *m) {
		*info = -7;
	} else if (itype <= 3 && *lda < max(1,*m)) {
		*info = -9;
	} else if (itype >= 4) {
		/* Computing MAX */
		i__1 = *m - 1;
		if (*kl < 0 || *kl > max(i__1,0)) {
			*info = -2;
		} else /* if(complicated condition) */ {
			/* Computing MAX */
			i__1 = *n - 1;
			if (*ku < 0 || *ku > max(i__1,0) || (itype == 4 || itype == 5) && 
				*kl != *ku) {
					*info = -3;
			} else if (itype == 4 && *lda < *kl + 1 || itype == 5 && *lda < *
				ku + 1 || itype == 6 && *lda < (*kl << 1) + *ku + 1) {
					*info = -9;
			}
		}
	}

	if (*info != 0) {
		i__1 = -(*info);
		xerbla_("DLASCL", &i__1);
		return 0;
	}

	/*     Quick return if possible */

	if (*n == 0 || *m == 0) {
		return 0;
	}

	/*     Get machine parameters */

	smlnum = dlamch_("S");
	bignum = 1. / smlnum;

	cfromc = *cfrom;
	ctoc = *cto;

L10:
	cfrom1 = cfromc * smlnum;
	cto1 = ctoc / bignum;
	if (abs(cfrom1) > abs(ctoc) && ctoc != 0.) {
		mul = smlnum;
		done = FALSE_;
		cfromc = cfrom1;
	} else if (abs(cto1) > abs(cfromc)) {
		mul = bignum;
		done = FALSE_;
		ctoc = cto1;
	} else {
		mul = ctoc / cfromc;
		done = TRUE_;
	}

	if (itype == 0) {

		/*        Full matrix */

		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
			i__2 = *m;
			for (i__ = 1; i__ <= i__2; ++i__) {
				a[i__ + j * a_dim1] *= mul;
				/* L20: */
			}
			/* L30: */
		}

	} else if (itype == 1) {

		/*        Lower triangular matrix */

		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
			i__2 = *m;
			for (i__ = j; i__ <= i__2; ++i__) {
				a[i__ + j * a_dim1] *= mul;
				/* L40: */
			}
			/* L50: */
		}

	} else if (itype == 2) {

		/*        Upper triangular matrix */

		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
			i__2 = min(j,*m);
			for (i__ = 1; i__ <= i__2; ++i__) {
				a[i__ + j * a_dim1] *= mul;
				/* L60: */
			}
			/* L70: */
		}

	} else if (itype == 3) {

		/*        Upper Hessenberg matrix */

		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
			/* Computing MIN */
			i__3 = j + 1;
			i__2 = min(i__3,*m);
			for (i__ = 1; i__ <= i__2; ++i__) {
				a[i__ + j * a_dim1] *= mul;
				/* L80: */
			}
			/* L90: */
		}

	} else if (itype == 4) {

		/*        Lower half of a symmetric band matrix */

		k3 = *kl + 1;
		k4 = *n + 1;
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
			/* Computing MIN */
			i__3 = k3, i__4 = k4 - j;
			i__2 = min(i__3,i__4);
			for (i__ = 1; i__ <= i__2; ++i__) {
				a[i__ + j * a_dim1] *= mul;
				/* L100: */
			}
			/* L110: */
		}

	} else if (itype == 5) {

		/*        Upper half of a symmetric band matrix */

		k1 = *ku + 2;
		k3 = *ku + 1;
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
			/* Computing MAX */
			i__2 = k1 - j;
			i__3 = k3;
			for (i__ = max(i__2,1); i__ <= i__3; ++i__) {
				a[i__ + j * a_dim1] *= mul;
				/* L120: */
			}
			/* L130: */
		}

	} else if (itype == 6) {

		/*        Band matrix */

		k1 = *kl + *ku + 2;
		k2 = *kl + 1;
		k3 = (*kl << 1) + *ku + 1;
		k4 = *kl + *ku + 1 + *m;
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
			/* Computing MAX */
			i__3 = k1 - j;
			/* Computing MIN */
			i__4 = k3, i__5 = k4 - j;
			i__2 = min(i__4,i__5);
			for (i__ = max(i__3,k2); i__ <= i__2; ++i__) {
				a[i__ + j * a_dim1] *= mul;
				/* L140: */
			}
			/* L150: */
		}

	}

	if (! done) {
		goto L10;
	}

	return 0;

	/*     End of DLASCL */

}




integer_ TRL::ilaenv_(integer_ *ispec, char *name__, char *opts, integer_ *n1, 
				integer_ *n2, integer_ *n3, integer_ *n4, ftnlen_ name_len, ftnlen_ 
				opts_len)
{
	static integer_ c__0 = 0;
	static trl_real_ c_b163 = 0.f;
	static trl_real_ c_b164 = 1.f;
	static integer_ c__1 = 1;
	/* System generated locals */
	integer_ ret_val;

	/* Builtin functions   
	Subroutine */ //void s_copy(char *, char *, ftnlen_, ftnlen_);
	//integer_ s_cmp(char *, char *, ftnlen_, ftnlen_);

	/* Local variables */
	static integer_ i__;
	static char c1[1], c2[2], c3[3], c4[2];
	static integer_ ic, nb, iz, nx;
	static logical_ cname;
	static integer_ nbmin;
	static logical_ sname;
	//extern integer_ ieeeck_(integer_ *, trl_real_ *, trl_real_ *);
	static char subnam[6];
	//extern integer_ iparmq_(integer_ *, char *, char *, integer_ *, integer_ *, integer_ *, integer_ *, ftnlen_, ftnlen_);


	/*  -- LAPACK auxiliary routine (version 3.1.1) --   
	Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..   
	January 2007   


	Purpose   
	=======   

	ILAENV is called from the LAPACK routines to choose problem-dependent   
	parameters for the local environment.  See ISPEC for a description of   
	the parameters.   

	ILAENV returns an INTEGER   
	if ILAENV >= 0: ILAENV returns the value of the parameter specified by ISPEC   
	if ILAENV < 0:  if ILAENV = -k, the k-th argument had an illegal value.   

	This version provides a set of parameters which should give good,   
	but not optimal, performance on many of the currently available   
	computers.  Users are encouraged to modify this subroutine to set   
	the tuning parameters for their particular machine using the option   
	and problem size information in the arguments.   

	This routine will not function correctly if it is converted to all   
	lower case.  Converting it to all upper case is allowed.   

	Arguments   
	=========   

	ISPEC   (input) INTEGER   
	Specifies the parameter to be returned as the value of   
	ILAENV.   
	= 1: the optimal blocksize; if this value is 1, an unblocked   
	algorithm will give the best performance.   
	= 2: the minimum block size for which the block routine   
	should be used; if the usable block size is less than   
	this value, an unblocked routine should be used.   
	= 3: the crossover point (in a block routine, for N less   
	than this value, an unblocked routine should be used)   
	= 4: the number of shifts, used in the nonsymmetric   
	eigenvalue routines (DEPRECATED)   
	= 5: the minimum column dimension for blocking to be used;   
	rectangular blocks must have dimension at least k by m,   
	where k is given by ILAENV(2,...) and m by ILAENV(5,...)   
	= 6: the crossover point for the SVD (when reducing an m by n   
	matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds   
	this value, a QR factorization is used first to reduce   
	the matrix to a triangular form.)   
	= 7: the number of processors   
	= 8: the crossover point for the multishift QR method   
	for nonsymmetric eigenvalue problems (DEPRECATED)   
	= 9: maximum size of the subproblems at the bottom of the   
	computation tree in the divide-and-conquer algorithm   
	(used by xGELSD and xGESDD)   
	=10: ieee NaN arithmetic can be trusted not to trap   
	=11: infinity arithmetic can be trusted not to trap   
	12 <= ISPEC <= 16:   
	xHSEQR or one of its subroutines,   
	see IPARMQ for detailed explanation   

	NAME    (input) CHARACTER*(*)   
	The name of the calling subroutine, in either upper case or   
	lower case.   

	OPTS    (input) CHARACTER*(*)   
	The character options to the subroutine NAME, concatenated   
	into a single character string.  For example, UPLO = 'U',   
	TRANS = 'T', and DIAG = 'N' for a triangular routine would   
	be specified as OPTS = 'UTN'.   

	N1      (input) INTEGER   
	N2      (input) INTEGER   
	N3      (input) INTEGER   
	N4      (input) INTEGER   
	Problem dimensions for the subroutine NAME; these may not all   
	be required.   

	Further Details   
	===============   

	The following conventions have been used when calling ILAENV from the   
	LAPACK routines:   
	1)  OPTS is a concatenation of all of the character options to   
	subroutine NAME, in the same order that they appear in the   
	argument list for NAME, even if they are not used in determining   
	the value of the parameter specified by ISPEC.   
	2)  The problem dimensions N1, N2, N3, N4 are specified in the order   
	that they appear in the argument list for NAME.  N1 is used   
	first, N2 second, and so on, and unused problem dimensions are   
	passed a value of -1.   
	3)  The parameter value returned by ILAENV is checked for validity in   
	the calling subroutine.  For example, ILAENV is used to retrieve   
	the optimal blocksize for STRTRI as follows:   

	NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )   
	IF( NB.LE.1 ) NB = MAX( 1, N )   

	===================================================================== */


	switch (*ispec) {
	case 1:  goto L10;
	case 2:  goto L10;
	case 3:  goto L10;
	case 4:  goto L80;
	case 5:  goto L90;
	case 6:  goto L100;
	case 7:  goto L110;
	case 8:  goto L120;
	case 9:  goto L130;
	case 10:  goto L140;
	case 11:  goto L150;
	case 12:  goto L160;
	case 13:  goto L160;
	case 14:  goto L160;
	case 15:  goto L160;
	case 16:  goto L160;
	}

	/*     Invalid value for ISPEC */

	ret_val = -1;
	return ret_val;

L10:

	/*     Convert NAME to upper case if the first character is lower case. */

	ret_val = 1;
	s_copy(subnam, name__, (ftnlen_)6, name_len);
	ic = *(unsigned char *)subnam;
	iz = 'Z';
	if (iz == 90 || iz == 122) {

		/*        ASCII character set */

		if (ic >= 97 && ic <= 122) {
			*(unsigned char *)subnam = (char) (ic - 32);
			for (i__ = 2; i__ <= 6; ++i__) {
				ic = *(unsigned char *)&subnam[i__ - 1];
				if (ic >= 97 && ic <= 122) {
					*(unsigned char *)&subnam[i__ - 1] = (char) (ic - 32);
				}
				/* L20: */
			}
		}

	} else if (iz == 233 || iz == 169) {

		/*        EBCDIC character set */

		if (ic >= 129 && ic <= 137 || ic >= 145 && ic <= 153 || ic >= 162 && 
			ic <= 169) {
				*(unsigned char *)subnam = (char) (ic + 64);
				for (i__ = 2; i__ <= 6; ++i__) {
					ic = *(unsigned char *)&subnam[i__ - 1];
					if (ic >= 129 && ic <= 137 || ic >= 145 && ic <= 153 || ic >= 
						162 && ic <= 169) {
							*(unsigned char *)&subnam[i__ - 1] = (char) (ic + 64);
					}
					/* L30: */
				}
		}

	} else if (iz == 218 || iz == 250) {

		/*        Prime machines:  ASCII+128 */

		if (ic >= 225 && ic <= 250) {
			*(unsigned char *)subnam = (char) (ic - 32);
			for (i__ = 2; i__ <= 6; ++i__) {
				ic = *(unsigned char *)&subnam[i__ - 1];
				if (ic >= 225 && ic <= 250) {
					*(unsigned char *)&subnam[i__ - 1] = (char) (ic - 32);
				}
				/* L40: */
			}
		}
	}

	*(unsigned char *)c1 = *(unsigned char *)subnam;
	sname = *(unsigned char *)c1 == 'S' || *(unsigned char *)c1 == 'D';
	cname = *(unsigned char *)c1 == 'C' || *(unsigned char *)c1 == 'Z';
	if (! (cname || sname)) {
		return ret_val;
	}
	s_copy(c2, subnam + 1, (ftnlen_)2, (ftnlen_)2);
	s_copy(c3, subnam + 3, (ftnlen_)3, (ftnlen_)3);
	s_copy(c4, c3 + 1, (ftnlen_)2, (ftnlen_)2);

	switch (*ispec) {
	case 1:  goto L50;
	case 2:  goto L60;
	case 3:  goto L70;
	}

L50:

	/*     ISPEC = 1:  block size   

	In these examples, separate code is provided for setting NB for   
	real and complex.  We assume that NB will take the same value in   
	single or double precision. */

	nb = 1;

	if (s_cmp(c2, "GE", (ftnlen_)2, (ftnlen_)2) == 0) {
		if (s_cmp(c3, "TRF", (ftnlen_)3, (ftnlen_)3) == 0) {
			if (sname) {
				nb = 64;
			} else {
				nb = 64;
			}
		} else if (s_cmp(c3, "QRF", (ftnlen_)3, (ftnlen_)3) == 0 || s_cmp(c3, 
			"RQF", (ftnlen_)3, (ftnlen_)3) == 0 || s_cmp(c3, "LQF", (ftnlen_)
			3, (ftnlen_)3) == 0 || s_cmp(c3, "QLF", (ftnlen_)3, (ftnlen_)3) 
			== 0) {
				if (sname) {
					nb = 32;
				} else {
					nb = 32;
				}
		} else if (s_cmp(c3, "HRD", (ftnlen_)3, (ftnlen_)3) == 0) {
			if (sname) {
				nb = 32;
			} else {
				nb = 32;
			}
		} else if (s_cmp(c3, "BRD", (ftnlen_)3, (ftnlen_)3) == 0) {
			if (sname) {
				nb = 32;
			} else {
				nb = 32;
			}
		} else if (s_cmp(c3, "TRI", (ftnlen_)3, (ftnlen_)3) == 0) {
			if (sname) {
				nb = 64;
			} else {
				nb = 64;
			}
		}
	} else if (s_cmp(c2, "PO", (ftnlen_)2, (ftnlen_)2) == 0) {
		if (s_cmp(c3, "TRF", (ftnlen_)3, (ftnlen_)3) == 0) {
			if (sname) {
				nb = 64;
			} else {
				nb = 64;
			}
		}
	} else if (s_cmp(c2, "SY", (ftnlen_)2, (ftnlen_)2) == 0) {
		if (s_cmp(c3, "TRF", (ftnlen_)3, (ftnlen_)3) == 0) {
			if (sname) {
				nb = 64;
			} else {
				nb = 64;
			}
		} else if (sname && s_cmp(c3, "TRD", (ftnlen_)3, (ftnlen_)3) == 0) {
			nb = 32;
		} else if (sname && s_cmp(c3, "GST", (ftnlen_)3, (ftnlen_)3) == 0) {
			nb = 64;
		}
	} else if (cname && s_cmp(c2, "HE", (ftnlen_)2, (ftnlen_)2) == 0) {
		if (s_cmp(c3, "TRF", (ftnlen_)3, (ftnlen_)3) == 0) {
			nb = 64;
		} else if (s_cmp(c3, "TRD", (ftnlen_)3, (ftnlen_)3) == 0) {
			nb = 32;
		} else if (s_cmp(c3, "GST", (ftnlen_)3, (ftnlen_)3) == 0) {
			nb = 64;
		}
	} else if (sname && s_cmp(c2, "OR", (ftnlen_)2, (ftnlen_)2) == 0) {
		if (*(unsigned char *)c3 == 'G') {
			if (s_cmp(c4, "QR", (ftnlen_)2, (ftnlen_)2) == 0 || s_cmp(c4, "RQ", 
				(ftnlen_)2, (ftnlen_)2) == 0 || s_cmp(c4, "LQ", (ftnlen_)2, (
				ftnlen_)2) == 0 || s_cmp(c4, "QL", (ftnlen_)2, (ftnlen_)2) ==
				0 || s_cmp(c4, "HR", (ftnlen_)2, (ftnlen_)2) == 0 || s_cmp(
				c4, "TR", (ftnlen_)2, (ftnlen_)2) == 0 || s_cmp(c4, "BR", (
				ftnlen_)2, (ftnlen_)2) == 0) {
					nb = 32;
			}
		} else if (*(unsigned char *)c3 == 'M') {
			if (s_cmp(c4, "QR", (ftnlen_)2, (ftnlen_)2) == 0 || s_cmp(c4, "RQ", 
				(ftnlen_)2, (ftnlen_)2) == 0 || s_cmp(c4, "LQ", (ftnlen_)2, (
				ftnlen_)2) == 0 || s_cmp(c4, "QL", (ftnlen_)2, (ftnlen_)2) ==
				0 || s_cmp(c4, "HR", (ftnlen_)2, (ftnlen_)2) == 0 || s_cmp(
				c4, "TR", (ftnlen_)2, (ftnlen_)2) == 0 || s_cmp(c4, "BR", (
				ftnlen_)2, (ftnlen_)2) == 0) {
					nb = 32;
			}
		}
	} else if (cname && s_cmp(c2, "UN", (ftnlen_)2, (ftnlen_)2) == 0) {
		if (*(unsigned char *)c3 == 'G') {
			if (s_cmp(c4, "QR", (ftnlen_)2, (ftnlen_)2) == 0 || s_cmp(c4, "RQ", 
				(ftnlen_)2, (ftnlen_)2) == 0 || s_cmp(c4, "LQ", (ftnlen_)2, (
				ftnlen_)2) == 0 || s_cmp(c4, "QL", (ftnlen_)2, (ftnlen_)2) ==
				0 || s_cmp(c4, "HR", (ftnlen_)2, (ftnlen_)2) == 0 || s_cmp(
				c4, "TR", (ftnlen_)2, (ftnlen_)2) == 0 || s_cmp(c4, "BR", (
				ftnlen_)2, (ftnlen_)2) == 0) {
					nb = 32;
			}
		} else if (*(unsigned char *)c3 == 'M') {
			if (s_cmp(c4, "QR", (ftnlen_)2, (ftnlen_)2) == 0 || s_cmp(c4, "RQ", 
				(ftnlen_)2, (ftnlen_)2) == 0 || s_cmp(c4, "LQ", (ftnlen_)2, (
				ftnlen_)2) == 0 || s_cmp(c4, "QL", (ftnlen_)2, (ftnlen_)2) ==
				0 || s_cmp(c4, "HR", (ftnlen_)2, (ftnlen_)2) == 0 || s_cmp(
				c4, "TR", (ftnlen_)2, (ftnlen_)2) == 0 || s_cmp(c4, "BR", (
				ftnlen_)2, (ftnlen_)2) == 0) {
					nb = 32;
			}
		}
	} else if (s_cmp(c2, "GB", (ftnlen_)2, (ftnlen_)2) == 0) {
		if (s_cmp(c3, "TRF", (ftnlen_)3, (ftnlen_)3) == 0) {
			if (sname) {
				if (*n4 <= 64) {
					nb = 1;
				} else {
					nb = 32;
				}
			} else {
				if (*n4 <= 64) {
					nb = 1;
				} else {
					nb = 32;
				}
			}
		}
	} else if (s_cmp(c2, "PB", (ftnlen_)2, (ftnlen_)2) == 0) {
		if (s_cmp(c3, "TRF", (ftnlen_)3, (ftnlen_)3) == 0) {
			if (sname) {
				if (*n2 <= 64) {
					nb = 1;
				} else {
					nb = 32;
				}
			} else {
				if (*n2 <= 64) {
					nb = 1;
				} else {
					nb = 32;
				}
			}
		}
	} else if (s_cmp(c2, "TR", (ftnlen_)2, (ftnlen_)2) == 0) {
		if (s_cmp(c3, "TRI", (ftnlen_)3, (ftnlen_)3) == 0) {
			if (sname) {
				nb = 64;
			} else {
				nb = 64;
			}
		}
	} else if (s_cmp(c2, "LA", (ftnlen_)2, (ftnlen_)2) == 0) {
		if (s_cmp(c3, "UUM", (ftnlen_)3, (ftnlen_)3) == 0) {
			if (sname) {
				nb = 64;
			} else {
				nb = 64;
			}
		}
	} else if (sname && s_cmp(c2, "ST", (ftnlen_)2, (ftnlen_)2) == 0) {
		if (s_cmp(c3, "EBZ", (ftnlen_)3, (ftnlen_)3) == 0) {
			nb = 1;
		}
	}
	ret_val = nb;
	return ret_val;

L60:

	/*     ISPEC = 2:  minimum block size */

	nbmin = 2;
	if (s_cmp(c2, "GE", (ftnlen_)2, (ftnlen_)2) == 0) {
		if (s_cmp(c3, "QRF", (ftnlen_)3, (ftnlen_)3) == 0 || s_cmp(c3, "RQF", (
			ftnlen_)3, (ftnlen_)3) == 0 || s_cmp(c3, "LQF", (ftnlen_)3, (
			ftnlen_)3) == 0 || s_cmp(c3, "QLF", (ftnlen_)3, (ftnlen_)3) == 0)
		{
			if (sname) {
				nbmin = 2;
			} else {
				nbmin = 2;
			}
		} else if (s_cmp(c3, "HRD", (ftnlen_)3, (ftnlen_)3) == 0) {
			if (sname) {
				nbmin = 2;
			} else {
				nbmin = 2;
			}
		} else if (s_cmp(c3, "BRD", (ftnlen_)3, (ftnlen_)3) == 0) {
			if (sname) {
				nbmin = 2;
			} else {
				nbmin = 2;
			}
		} else if (s_cmp(c3, "TRI", (ftnlen_)3, (ftnlen_)3) == 0) {
			if (sname) {
				nbmin = 2;
			} else {
				nbmin = 2;
			}
		}
	} else if (s_cmp(c2, "SY", (ftnlen_)2, (ftnlen_)2) == 0) {
		if (s_cmp(c3, "TRF", (ftnlen_)3, (ftnlen_)3) == 0) {
			if (sname) {
				nbmin = 8;
			} else {
				nbmin = 8;
			}
		} else if (sname && s_cmp(c3, "TRD", (ftnlen_)3, (ftnlen_)3) == 0) {
			nbmin = 2;
		}
	} else if (cname && s_cmp(c2, "HE", (ftnlen_)2, (ftnlen_)2) == 0) {
		if (s_cmp(c3, "TRD", (ftnlen_)3, (ftnlen_)3) == 0) {
			nbmin = 2;
		}
	} else if (sname && s_cmp(c2, "OR", (ftnlen_)2, (ftnlen_)2) == 0) {
		if (*(unsigned char *)c3 == 'G') {
			if (s_cmp(c4, "QR", (ftnlen_)2, (ftnlen_)2) == 0 || s_cmp(c4, "RQ", 
				(ftnlen_)2, (ftnlen_)2) == 0 || s_cmp(c4, "LQ", (ftnlen_)2, (
				ftnlen_)2) == 0 || s_cmp(c4, "QL", (ftnlen_)2, (ftnlen_)2) ==
				0 || s_cmp(c4, "HR", (ftnlen_)2, (ftnlen_)2) == 0 || s_cmp(
				c4, "TR", (ftnlen_)2, (ftnlen_)2) == 0 || s_cmp(c4, "BR", (
				ftnlen_)2, (ftnlen_)2) == 0) {
					nbmin = 2;
			}
		} else if (*(unsigned char *)c3 == 'M') {
			if (s_cmp(c4, "QR", (ftnlen_)2, (ftnlen_)2) == 0 || s_cmp(c4, "RQ", 
				(ftnlen_)2, (ftnlen_)2) == 0 || s_cmp(c4, "LQ", (ftnlen_)2, (
				ftnlen_)2) == 0 || s_cmp(c4, "QL", (ftnlen_)2, (ftnlen_)2) ==
				0 || s_cmp(c4, "HR", (ftnlen_)2, (ftnlen_)2) == 0 || s_cmp(
				c4, "TR", (ftnlen_)2, (ftnlen_)2) == 0 || s_cmp(c4, "BR", (
				ftnlen_)2, (ftnlen_)2) == 0) {
					nbmin = 2;
			}
		}
	} else if (cname && s_cmp(c2, "UN", (ftnlen_)2, (ftnlen_)2) == 0) {
		if (*(unsigned char *)c3 == 'G') {
			if (s_cmp(c4, "QR", (ftnlen_)2, (ftnlen_)2) == 0 || s_cmp(c4, "RQ", 
				(ftnlen_)2, (ftnlen_)2) == 0 || s_cmp(c4, "LQ", (ftnlen_)2, (
				ftnlen_)2) == 0 || s_cmp(c4, "QL", (ftnlen_)2, (ftnlen_)2) ==
				0 || s_cmp(c4, "HR", (ftnlen_)2, (ftnlen_)2) == 0 || s_cmp(
				c4, "TR", (ftnlen_)2, (ftnlen_)2) == 0 || s_cmp(c4, "BR", (
				ftnlen_)2, (ftnlen_)2) == 0) {
					nbmin = 2;
			}
		} else if (*(unsigned char *)c3 == 'M') {
			if (s_cmp(c4, "QR", (ftnlen_)2, (ftnlen_)2) == 0 || s_cmp(c4, "RQ", 
				(ftnlen_)2, (ftnlen_)2) == 0 || s_cmp(c4, "LQ", (ftnlen_)2, (
				ftnlen_)2) == 0 || s_cmp(c4, "QL", (ftnlen_)2, (ftnlen_)2) ==
				0 || s_cmp(c4, "HR", (ftnlen_)2, (ftnlen_)2) == 0 || s_cmp(
				c4, "TR", (ftnlen_)2, (ftnlen_)2) == 0 || s_cmp(c4, "BR", (
				ftnlen_)2, (ftnlen_)2) == 0) {
					nbmin = 2;
			}
		}
	}
	ret_val = nbmin;
	return ret_val;

L70:

	/*     ISPEC = 3:  crossover point */

	nx = 0;
	if (s_cmp(c2, "GE", (ftnlen_)2, (ftnlen_)2) == 0) {
		if (s_cmp(c3, "QRF", (ftnlen_)3, (ftnlen_)3) == 0 || s_cmp(c3, "RQF", (
			ftnlen_)3, (ftnlen_)3) == 0 || s_cmp(c3, "LQF", (ftnlen_)3, (
			ftnlen_)3) == 0 || s_cmp(c3, "QLF", (ftnlen_)3, (ftnlen_)3) == 0)
		{
			if (sname) {
				nx = 128;
			} else {
				nx = 128;
			}
		} else if (s_cmp(c3, "HRD", (ftnlen_)3, (ftnlen_)3) == 0) {
			if (sname) {
				nx = 128;
			} else {
				nx = 128;
			}
		} else if (s_cmp(c3, "BRD", (ftnlen_)3, (ftnlen_)3) == 0) {
			if (sname) {
				nx = 128;
			} else {
				nx = 128;
			}
		}
	} else if (s_cmp(c2, "SY", (ftnlen_)2, (ftnlen_)2) == 0) {
		if (sname && s_cmp(c3, "TRD", (ftnlen_)3, (ftnlen_)3) == 0) {
			nx = 32;
		}
	} else if (cname && s_cmp(c2, "HE", (ftnlen_)2, (ftnlen_)2) == 0) {
		if (s_cmp(c3, "TRD", (ftnlen_)3, (ftnlen_)3) == 0) {
			nx = 32;
		}
	} else if (sname && s_cmp(c2, "OR", (ftnlen_)2, (ftnlen_)2) == 0) {
		if (*(unsigned char *)c3 == 'G') {
			if (s_cmp(c4, "QR", (ftnlen_)2, (ftnlen_)2) == 0 || s_cmp(c4, "RQ", 
				(ftnlen_)2, (ftnlen_)2) == 0 || s_cmp(c4, "LQ", (ftnlen_)2, (
				ftnlen_)2) == 0 || s_cmp(c4, "QL", (ftnlen_)2, (ftnlen_)2) ==
				0 || s_cmp(c4, "HR", (ftnlen_)2, (ftnlen_)2) == 0 || s_cmp(
				c4, "TR", (ftnlen_)2, (ftnlen_)2) == 0 || s_cmp(c4, "BR", (
				ftnlen_)2, (ftnlen_)2) == 0) {
					nx = 128;
			}
		}
	} else if (cname && s_cmp(c2, "UN", (ftnlen_)2, (ftnlen_)2) == 0) {
		if (*(unsigned char *)c3 == 'G') {
			if (s_cmp(c4, "QR", (ftnlen_)2, (ftnlen_)2) == 0 || s_cmp(c4, "RQ", 
				(ftnlen_)2, (ftnlen_)2) == 0 || s_cmp(c4, "LQ", (ftnlen_)2, (
				ftnlen_)2) == 0 || s_cmp(c4, "QL", (ftnlen_)2, (ftnlen_)2) ==
				0 || s_cmp(c4, "HR", (ftnlen_)2, (ftnlen_)2) == 0 || s_cmp(
				c4, "TR", (ftnlen_)2, (ftnlen_)2) == 0 || s_cmp(c4, "BR", (
				ftnlen_)2, (ftnlen_)2) == 0) {
					nx = 128;
			}
		}
	}
	ret_val = nx;
	return ret_val;

L80:

	/*     ISPEC = 4:  number of shifts (used by xHSEQR) */

	ret_val = 6;
	return ret_val;

L90:

	/*     ISPEC = 5:  minimum column dimension (not used) */

	ret_val = 2;
	return ret_val;

L100:

	/*     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD) */

	ret_val = (integer_) ((trl_real_) min(*n1,*n2) * 1.6f);
	return ret_val;

L110:

	/*     ISPEC = 7:  number of processors (not used) */

	ret_val = 1;
	return ret_val;

L120:

	/*     ISPEC = 8:  crossover point for multishift (used by xHSEQR) */

	ret_val = 50;
	return ret_val;

L130:

	/*     ISPEC = 9:  maximum size of the subproblems at the bottom of the   
	computation tree in the divide-and-conquer algorithm   
	(used by xGELSD and xGESDD) */

	ret_val = 25;
	return ret_val;

L140:

	/*     ISPEC = 10: ieee NaN arithmetic can be trusted not to trap   

	ILAENV = 0 */
	ret_val = 1;
	if (ret_val == 1) {
		ret_val = ieeeck_(&c__0, &c_b163, &c_b164);
	}
	return ret_val;

L150:

	/*     ISPEC = 11: infinity arithmetic can be trusted not to trap   

	ILAENV = 0 */
	ret_val = 1;
	if (ret_val == 1) {
		ret_val = ieeeck_(&c__1, &c_b163, &c_b164);
	}
	return ret_val;

L160:

	/*     12 <= ISPEC <= 16: xHSEQR or one of its subroutines. */

	ret_val = iparmq_(ispec, name__, opts, n1, n2, n3, n4, name_len, opts_len)
		;
	return ret_val;

	/*     End of ILAENV */

} /* ilaenv_ */

int TRL::xerbla_(char *srname, integer_ *info)
{
	/*  -- LAPACK auxiliary routine (version 3.1) --   
	Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
	Courant Institute, Argonne National Lab, and Rice University   
	November 2006


	Purpose   
	=======   

	XERBLA  is an error handler for the LAPACK routines.   
	It is called by an LAPACK routine if an input parameter has an   
	invalid value.  A message is printed and execution stops.   

	Installers may consider modifying the STOP statement in order to   
	call system-specific exception-handling facilities.   

	Arguments   
	=========   

	SRNAME  (input) CHARACTER*6   
	The name of the routine which called XERBLA.   

	INFO    (input) INTEGER   
	The position of the invalid parameter in the parameter list   
	of the calling routine.   

	===================================================================== 
	*/

	printf("** On entry to %6s, parameter number %2i had an illegal value\n",
		srname, *info);

	/*     End of XERBLA */

	return 0;
}

//
int TRL::dsterf_(integer_ *n, doublereal_ *d__, doublereal_ *e, 
			integer_ *info)
{
	/*  -- LAPACK routine (version 3.1) --   
	Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..   
	November 2006   


	Purpose   
	=======   

	DSTERF computes all eigenvalues of a symmetric tridiagonal matrix   
	using the Pal-Walker-Kahan variant of the QL or QR algorithm.   

	Arguments   
	=========   

	N       (input) INTEGER   
	The order of the matrix.  N >= 0.   

	D       (input/output) DOUBLE PRECISION array, dimension (N)   
	On entry, the n diagonal elements of the tridiagonal matrix.   
	On exit, if INFO = 0, the eigenvalues in ascending order.   

	E       (input/output) DOUBLE PRECISION array, dimension (N-1)   
	On entry, the (n-1) subdiagonal elements of the tridiagonal   
	matrix.   
	On exit, E has been destroyed.   

	INFO    (output) INTEGER   
	= 0:  successful exit   
	< 0:  if INFO = -i, the i-th argument had an illegal value   
	> 0:  the algorithm failed to find all of the eigenvalues in   
	a total of 30*N iterations; if INFO = i, then i   
	elements of E have not converged to zero.   

	=====================================================================   


	Test the input parameters.   

	Parameter adjustments */
	/* Table of constant values */
	static integer_ c__0 = 0;
	static integer_ c__1 = 1;
	static doublereal_ c_b32 = 1.;

	/* System generated locals */
	integer_ i__1;
	doublereal_ d__1, d__2, d__3;
	/* Builtin functions */
	//double sqrt_(doublereal_), d_sign(doublereal_ *, doublereal_ *);
	/* Local variables */
	static doublereal_ c__;
	static integer_ i__, l, m;
	static doublereal_ p, r__, s;
	static integer_ l1;
	static doublereal_ bb, rt1, rt2, eps, rte;
	static integer_ lsv;
	static doublereal_ eps2, oldc;
	static integer_ lend, jtot;
	//extern /* Subroutine */ int dlae2_(doublereal_ *, doublereal_ *, doublereal_ 
	//	*, doublereal_ *, doublereal_ *);
	static doublereal_ gamma, alpha, sigma, anorm;
	//extern doublereal_ dlapy2_(doublereal_ *, doublereal_ *), dlamch_(char *);
	static integer_ iscale;
	//extern /* Subroutine */ int dlascl_(char *, integer_ *, integer_ *, 
	//	doublereal_ *, doublereal_ *, integer_ *, integer_ *, doublereal_ *, 
	//	integer_ *, integer_ *);
	static doublereal_ oldgam, safmin;
	//extern /* Subroutine */ int xerbla_(char *, integer_ *);
	static doublereal_ safmax;
	//extern doublereal_ dlanst_(char *, integer_ *, doublereal_ *, doublereal_ *);
	//extern /* Subroutine */ int dlasrt_(char *, integer_ *, doublereal_ *, 
	//	integer_ *);
	static integer_ lendsv;
	static doublereal_ ssfmin;
	static integer_ nmaxit;
	static doublereal_ ssfmax;


	--e;
	--d__;

	/* Function Body */
	*info = 0;

	/*     Quick return if possible */

	if (*n < 0) {
		*info = -1;
		i__1 = -(*info);
		xerbla_("DSTERF", &i__1);
		return 0;
	}
	if (*n <= 1) {
		return 0;
	}

	/*     Determine the unit roundoff for this environment. */

	eps = dlamch_("E");
	/* Computing 2nd power */
	d__1 = eps;
	eps2 = d__1 * d__1;
	safmin = dlamch_("S");
	safmax = 1. / safmin;
	ssfmax = sqrt_(safmax) / 3.;
	ssfmin = sqrt_(safmin) / eps2;

	/*     Compute the eigenvalues of the tridiagonal matrix. */

	nmaxit = *n * 30;
	sigma = 0.;
	jtot = 0;

	/*     Determine where the matrix splits and choose QL or QR iteration   
	for each block, according to whether top or bottom diagonal   
	element is smaller. */

	l1 = 1;

L10:
	if (l1 > *n) {
		goto L170;
	}
	if (l1 > 1) {
		e[l1 - 1] = 0.;
	}
	i__1 = *n - 1;
	for (m = l1; m <= i__1; ++m) {
		if ((d__3 = e[m], abs(d__3)) <= sqrt_((d__1 = d__[m], abs(d__1))) * 
			sqrt_((d__2 = d__[m + 1], abs(d__2))) * eps) {
				e[m] = 0.;
				goto L30;
		}
		/* L20: */
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

	/*     Scale submatrix in rows and columns L to LEND */

	i__1 = lend - l + 1;
	anorm = dlanst_("I", &i__1, &d__[l], &e[l]);
	iscale = 0;
	if (anorm > ssfmax) {
		iscale = 1;
		i__1 = lend - l + 1;
		dlascl_("G", &c__0, &c__0, &anorm, &ssfmax, &i__1, &c__1, &d__[l], n, 
			info);
		i__1 = lend - l;
		dlascl_("G", &c__0, &c__0, &anorm, &ssfmax, &i__1, &c__1, &e[l], n, 
			info);
	} else if (anorm < ssfmin) {
		iscale = 2;
		i__1 = lend - l + 1;
		dlascl_("G", &c__0, &c__0, &anorm, &ssfmin, &i__1, &c__1, &d__[l], n, 
			info);
		i__1 = lend - l;
		dlascl_("G", &c__0, &c__0, &anorm, &ssfmin, &i__1, &c__1, &e[l], n, 
			info);
	}

	i__1 = lend - 1;
	for (i__ = l; i__ <= i__1; ++i__) {
		/* Computing 2nd power */
		d__1 = e[i__];
		e[i__] = d__1 * d__1;
		/* L40: */
	}

	/*     Choose between QL and QR iteration */

	if ((d__1 = d__[lend], abs(d__1)) < (d__2 = d__[l], abs(d__2))) {
		lend = lsv;
		l = lendsv;
	}

	if (lend >= l) {

		/*        QL Iteration   

		Look for small subdiagonal element. */

L50:
		if (l != lend) {
			i__1 = lend - 1;
			for (m = l; m <= i__1; ++m) {
				if ((d__2 = e[m], abs(d__2)) <= eps2 * (d__1 = d__[m] * d__[m 
					+ 1], abs(d__1))) {
						goto L70;
				}
				/* L60: */
			}
		}
		m = lend;

L70:
		if (m < lend) {
			e[m] = 0.;
		}
		p = d__[l];
		if (m == l) {
			goto L90;
		}

		/*        If remaining matrix is 2 by 2, use DLAE2 to compute its   
		eigenvalues. */

		if (m == l + 1) {
			rte = sqrt_(e[l]);
			dlae2_(&d__[l], &rte, &d__[l + 1], &rt1, &rt2);
			d__[l] = rt1;
			d__[l + 1] = rt2;
			e[l] = 0.;
			l += 2;
			if (l <= lend) {
				goto L50;
			}
			goto L150;
		}

		if (jtot == nmaxit) {
			goto L150;
		}
		++jtot;

		/*        Form shift. */

		rte = sqrt_(e[l]);
		sigma = (d__[l + 1] - p) / (rte * 2.);
		r__ = dlapy2_(&sigma, &c_b32);
		sigma = p - rte / (sigma + d_sign(&r__, &sigma));

		c__ = 1.;
		s = 0.;
		gamma = d__[m] - sigma;
		p = gamma * gamma;

		/*        Inner loop */

		i__1 = l;
		for (i__ = m - 1; i__ >= i__1; --i__) {
			bb = e[i__];
			r__ = p + bb;
			if (i__ != m - 1) {
				e[i__ + 1] = s * r__;
			}
			oldc = c__;
			c__ = p / r__;
			s = bb / r__;
			oldgam = gamma;
			alpha = d__[i__];
			gamma = c__ * (alpha - sigma) - s * oldgam;
			d__[i__ + 1] = oldgam + (alpha - gamma);
			if (c__ != 0.) {
				p = gamma * gamma / c__;
			} else {
				p = oldc * bb;
			}
			/* L80: */
		}

		e[l] = s * p;
		d__[l] = sigma + gamma;
		goto L50;

		/*        Eigenvalue found. */

L90:
		d__[l] = p;

		++l;
		if (l <= lend) {
			goto L50;
		}
		goto L150;

	} else {

		/*        QR Iteration   

		Look for small superdiagonal element. */

L100:
		i__1 = lend + 1;
		for (m = l; m >= i__1; --m) {
			if ((d__2 = e[m - 1], abs(d__2)) <= eps2 * (d__1 = d__[m] * d__[m 
				- 1], abs(d__1))) {
					goto L120;
			}
			/* L110: */
		}
		m = lend;

L120:
		if (m > lend) {
			e[m - 1] = 0.;
		}
		p = d__[l];
		if (m == l) {
			goto L140;
		}

		/*        If remaining matrix is 2 by 2, use DLAE2 to compute its   
		eigenvalues. */

		if (m == l - 1) {
			rte = sqrt_(e[l - 1]);
			dlae2_(&d__[l], &rte, &d__[l - 1], &rt1, &rt2);
			d__[l] = rt1;
			d__[l - 1] = rt2;
			e[l - 1] = 0.;
			l += -2;
			if (l >= lend) {
				goto L100;
			}
			goto L150;
		}

		if (jtot == nmaxit) {
			goto L150;
		}
		++jtot;

		/*        Form shift. */

		rte = sqrt_(e[l - 1]);
		sigma = (d__[l - 1] - p) / (rte * 2.);
		r__ = dlapy2_(&sigma, &c_b32);
		sigma = p - rte / (sigma + d_sign(&r__, &sigma));

		c__ = 1.;
		s = 0.;
		gamma = d__[m] - sigma;
		p = gamma * gamma;

		/*        Inner loop */

		i__1 = l - 1;
		for (i__ = m; i__ <= i__1; ++i__) {
			bb = e[i__];
			r__ = p + bb;
			if (i__ != m) {
				e[i__ - 1] = s * r__;
			}
			oldc = c__;
			c__ = p / r__;
			s = bb / r__;
			oldgam = gamma;
			alpha = d__[i__ + 1];
			gamma = c__ * (alpha - sigma) - s * oldgam;
			d__[i__] = oldgam + (alpha - gamma);
			if (c__ != 0.) {
				p = gamma * gamma / c__;
			} else {
				p = oldc * bb;
			}
			/* L130: */
		}

		e[l - 1] = s * p;
		d__[l] = sigma + gamma;
		goto L100;

		/*        Eigenvalue found. */

L140:
		d__[l] = p;

		--l;
		if (l >= lend) {
			goto L100;
		}
		goto L150;

	}

	/*     Undo scaling if necessary */

L150:
	if (iscale == 1) {
		i__1 = lendsv - lsv + 1;
		dlascl_("G", &c__0, &c__0, &ssfmax, &anorm, &i__1, &c__1, &d__[lsv], 
			n, info);
	}
	if (iscale == 2) {
		i__1 = lendsv - lsv + 1;
		dlascl_("G", &c__0, &c__0, &ssfmin, &anorm, &i__1, &c__1, &d__[lsv], 
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
		/* L160: */
	}
	goto L180;

	/*     Sort eigenvalues in increasing order. */

L170:
	dlasrt_("I", n, &d__[1], info);

L180:
	return 0;

	/*     End of DSTERF */

} /* dsterf_ */

//
doublereal_ TRL::dlansy_(char *norm, char *uplo, integer_ *n, doublereal_ *a, integer_ 
				   *lda, doublereal_ *work)
{
	/*  -- LAPACK auxiliary routine (version 3.1) --   
	Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..   
	November 2006   


	Purpose   
	=======   

	DLANSY  returns the value of the one norm,  or the Frobenius norm, or   
	the  infinity norm,  or the  element of  largest absolute value  of a   
	real symmetric matrix A.   

	Description   
	===========   

	DLANSY returns the value   

	DLANSY = ( max(abs(A(i,j))), NORM = 'M' or 'm'   
	(   
	( norm1(A),         NORM = '1', 'O' or 'o'   
	(   
	( normI(A),         NORM = 'I' or 'i'   
	(   
	( normF(A),         NORM = 'F', 'f', 'E' or 'e'   

	where  norm1  denotes the  one norm of a matrix (maximum column sum),   
	normI  denotes the  infinity norm  of a matrix  (maximum row sum) and   
	normF  denotes the  Frobenius norm of a matrix (square root of sum of   
	squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.   

	Arguments   
	=========   

	NORM    (input) CHARACTER*1   
	Specifies the value to be returned in DLANSY as described   
	above.   

	UPLO    (input) CHARACTER*1   
	Specifies whether the upper or lower triangular part of the   
	symmetric matrix A is to be referenced.   
	= 'U':  Upper triangular part of A is referenced   
	= 'L':  Lower triangular part of A is referenced   

	N       (input) INTEGER   
	The order of the matrix A.  N >= 0.  When N = 0, DLANSY is   
	set to zero.   

	A       (input) DOUBLE PRECISION array, dimension (LDA,N)   
	The symmetric matrix A.  If UPLO = 'U', the leading n by n   
	upper triangular part of A contains the upper triangular part   
	of the matrix A, and the strictly lower triangular part of A   
	is not referenced.  If UPLO = 'L', the leading n by n lower   
	triangular part of A contains the lower triangular part of   
	the matrix A, and the strictly upper triangular part of A is   
	not referenced.   

	LDA     (input) INTEGER   
	The leading dimension of the array A.  LDA >= max(N,1).   

	WORK    (workspace) DOUBLE PRECISION array, dimension (MAX(1,LWORK)),   
	where LWORK >= N when NORM = 'I' or '1' or 'O'; otherwise,   
	WORK is not referenced.   

	=====================================================================   


	Parameter adjustments */
	/* Table of constant values */
	static integer_ c__1 = 1;

	/* System generated locals */
	integer_ a_dim1, a_offset, i__1, i__2;
	doublereal_ ret_val, d__1, d__2, d__3;
	/* Builtin functions */
	//double sqrt_(doublereal_);
	/* Local variables */
	static integer_ i__, j;
	static doublereal_ sum, absa, scale;
	//extern logical_ lsame_(char *, char *);
	static doublereal_ value;
	//extern /* Subroutine */ int dlassq_(integer_ *, doublereal_ *, integer_ *, 
	//	doublereal_ *, doublereal_ *);


	a_dim1 = *lda;
	a_offset = 1 + a_dim1;
	a -= a_offset;
	--work;

	/* Function Body */
	if (*n == 0) {
		value = 0.;
	} else if (lsame_(norm, "M")) {

		/*        Find max(abs(A(i,j))). */

		value = 0.;
		if (lsame_(uplo, "U")) {
			i__1 = *n;
			for (j = 1; j <= i__1; ++j) {
				i__2 = j;
				for (i__ = 1; i__ <= i__2; ++i__) {
					/* Computing MAX */
					d__2 = value, d__3 = (d__1 = a[i__ + j * a_dim1], abs(
						d__1));
					value = max(d__2,d__3);
					/* L10: */
				}
				/* L20: */
			}
		} else {
			i__1 = *n;
			for (j = 1; j <= i__1; ++j) {
				i__2 = *n;
				for (i__ = j; i__ <= i__2; ++i__) {
					/* Computing MAX */
					d__2 = value, d__3 = (d__1 = a[i__ + j * a_dim1], abs(
						d__1));
					value = max(d__2,d__3);
					/* L30: */
				}
				/* L40: */
			}
		}
	} else if (lsame_(norm, "I") || lsame_(norm, "O") || *(unsigned char *)norm == '1') {

		/*        Find normI(A) ( = norm1(A), since A is symmetric). */

		value = 0.;
		if (lsame_(uplo, "U")) {
			i__1 = *n;
			for (j = 1; j <= i__1; ++j) {
				sum = 0.;
				i__2 = j - 1;
				for (i__ = 1; i__ <= i__2; ++i__) {
					absa = (d__1 = a[i__ + j * a_dim1], abs(d__1));
					sum += absa;
					work[i__] += absa;
					/* L50: */
				}
				work[j] = sum + (d__1 = a[j + j * a_dim1], abs(d__1));
				/* L60: */
			}
			i__1 = *n;
			for (i__ = 1; i__ <= i__1; ++i__) {
				/* Computing MAX */
				d__1 = value, d__2 = work[i__];
				value = max(d__1,d__2);
				/* L70: */
			}
		} else {
			i__1 = *n;
			for (i__ = 1; i__ <= i__1; ++i__) {
				work[i__] = 0.;
				/* L80: */
			}
			i__1 = *n;
			for (j = 1; j <= i__1; ++j) {
				sum = work[j] + (d__1 = a[j + j * a_dim1], abs(d__1));
				i__2 = *n;
				for (i__ = j + 1; i__ <= i__2; ++i__) {
					absa = (d__1 = a[i__ + j * a_dim1], abs(d__1));
					sum += absa;
					work[i__] += absa;
					/* L90: */
				}
				value = max(value,sum);
				/* L100: */
			}
		}
	} else if (lsame_(norm, "F") || lsame_(norm, "E")) {

		/*        Find normF(A). */

		scale = 0.;
		sum = 1.;
		if (lsame_(uplo, "U")) {
			i__1 = *n;
			for (j = 2; j <= i__1; ++j) {
				i__2 = j - 1;
				dlassq_(&i__2, &a[j * a_dim1 + 1], &c__1, &scale, &sum);
				/* L110: */
			}
		} else {
			i__1 = *n - 1;
			for (j = 1; j <= i__1; ++j) {
				i__2 = *n - j;
				dlassq_(&i__2, &a[j + 1 + j * a_dim1], &c__1, &scale, &sum);
				/* L120: */
			}
		}
		sum *= 2;
		i__1 = *lda + 1;
		dlassq_(n, &a[a_offset], &i__1, &scale, &sum);
		value = scale * sqrt_(sum);
	}

	ret_val = value;
	return ret_val;

	/*     End of DLANSY */

} /* dlansy_ */


//
int TRL::dtrsm_(char *side, char *uplo, char *transa, char *diag, 
		   integer_ *m, integer_ *n, doublereal_ *alpha, doublereal_ *a, integer_ *
		   lda, doublereal_ *b, integer_ *ldb)
{
	/* System generated locals */
	integer_ a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;
	/* Local variables */
	static integer_ i__, j, k, info;
	static doublereal_ temp;
	static logical_ lside;
	//extern logical_ lsame_(char *, char *);
	static integer_ nrowa;
	static logical_ upper;
	//extern /* Subroutine */ int xerbla_(char *, integer_ *);
	static logical_ nounit;
	/*  Purpose   
	=======   
	DTRSM  solves one of the matrix equations   
	op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,   
	where alpha is a scalar, X and B are m by n matrices, A is a unit, or   
	non-unit,  upper or lower triangular matrix  and  op( A )  is one  of   
	op( A ) = A   or   op( A ) = A'.   
	The matrix X is overwritten on B.   
	Arguments   
	==========   
	SIDE   - CHARACTER*1.   
	On entry, SIDE specifies whether op( A ) appears on the left   
	or right of X as follows:   
	SIDE = 'L' or 'l'   op( A )*X = alpha*B.   
	SIDE = 'R' or 'r'   X*op( A ) = alpha*B.   
	Unchanged on exit.   
	UPLO   - CHARACTER*1.   
	On entry, UPLO specifies whether the matrix A is an upper or   
	lower triangular matrix as follows:   
	UPLO = 'U' or 'u'   A is an upper triangular matrix.   
	UPLO = 'L' or 'l'   A is a lower triangular matrix.   
	Unchanged on exit.   
	TRANSA - CHARACTER*1.   
	On entry, TRANSA specifies the form of op( A ) to be used in   
	the matrix multiplication as follows:   
	TRANSA = 'N' or 'n'   op( A ) = A.   
	TRANSA = 'T' or 't'   op( A ) = A'.   
	TRANSA = 'C' or 'c'   op( A ) = A'.   
	Unchanged on exit.   
	DIAG   - CHARACTER*1.   
	On entry, DIAG specifies whether or not A is unit triangular   
	as follows:   
	DIAG = 'U' or 'u'   A is assumed to be unit triangular.   
	DIAG = 'N' or 'n'   A is not assumed to be unit   
	triangular.   
	Unchanged on exit.   
	M      - INTEGER.   
	On entry, M specifies the number of rows of B. M must be at   
	least zero.   
	Unchanged on exit.   
	N      - INTEGER.   
	On entry, N specifies the number of columns of B.  N must be   
	at least zero.   
	Unchanged on exit.   
	ALPHA  - DOUBLE PRECISION.   
	On entry,  ALPHA specifies the scalar  alpha. When  alpha is   
	zero then  A is not referenced and  B need not be set before   
	entry.   
	Unchanged on exit.   
	A      - DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m   
	when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.   
	Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k   
	upper triangular part of the array  A must contain the upper   
	triangular matrix  and the strictly lower triangular part of   
	A is not referenced.   
	Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k   
	lower triangular part of the array  A must contain the lower   
	triangular matrix  and the strictly upper triangular part of   
	A is not referenced.   
	Note that when  DIAG = 'U' or 'u',  the diagonal elements of   
	A  are not referenced either,  but are assumed to be  unity.   
	Unchanged on exit.   
	LDA    - INTEGER.   
	On entry, LDA specifies the first dimension of A as declared   
	in the calling (sub) program.  When  SIDE = 'L' or 'l'  then   
	LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'   
	then LDA must be at least max( 1, n ).   
	Unchanged on exit.   
	B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ).   
	Before entry,  the leading  m by n part of the array  B must   
	contain  the  right-hand  side  matrix  B,  and  on exit  is   
	overwritten by the solution matrix  X.   
	LDB    - INTEGER.   
	On entry, LDB specifies the first dimension of B as declared   
	in  the  calling  (sub)  program.   LDB  must  be  at  least   
	max( 1, m ).   
	Unchanged on exit.   
	Level 3 Blas routine.   
	-- Written on 8-February-1989.   
	Jack Dongarra, Argonne National Laboratory.   
	Iain Duff, AERE Harwell.   
	Jeremy Du Croz, Numerical Algorithms Group Ltd.   
	Sven Hammarling, Numerical Algorithms Group Ltd.   
	Test the input parameters.   
	Parameter adjustments */
	a_dim1 = *lda;
	a_offset = 1 + a_dim1;
	a -= a_offset;
	b_dim1 = *ldb;
	b_offset = 1 + b_dim1;
	b -= b_offset;
	/* Function Body */
	lside = lsame_(side, "L");
	if (lside) {
		nrowa = *m;
	} else {
		nrowa = *n;
	}
	nounit = lsame_(diag, "N");
	upper = lsame_(uplo, "U");
	info = 0;
	if (! lside && ! lsame_(side, "R")) {
		info = 1;
	} else if (! upper && ! lsame_(uplo, "L")) {
		info = 2;
	} else if (! lsame_(transa, "N") && ! lsame_(transa, 
		"T") && ! lsame_(transa, "C")) {
			info = 3;
	} else if (! lsame_(diag, "U") && ! lsame_(diag, 
		"N")) {
			info = 4;
	} else if (*m < 0) {
		info = 5;
	} else if (*n < 0) {
		info = 6;
	} else if (*lda < max(1,nrowa)) {
		info = 9;
	} else if (*ldb < max(1,*m)) {
		info = 11;
	}
	if (info != 0) {
		xerbla_("DTRSM ", &info);
		return 0;
	}
	/*     Quick return if possible. */
	if (*n == 0) {
		return 0;
	}
	/*     And when  alpha.eq.zero. */
	if (*alpha == 0.) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
			i__2 = *m;
			for (i__ = 1; i__ <= i__2; ++i__) {
				b[i__ + j * b_dim1] = 0.;
				/* L10: */
			}
			/* L20: */
		}
		return 0;
	}
	/*     Start the operations. */
	if (lside) {
		if (lsame_(transa, "N")) {
			/*           Form  B := alpha*inv( A )*B. */
			if (upper) {
				i__1 = *n;
				for (j = 1; j <= i__1; ++j) {
					if (*alpha != 1.) {
						i__2 = *m;
						for (i__ = 1; i__ <= i__2; ++i__) {
							b[i__ + j * b_dim1] = *alpha * b[i__ + j * b_dim1]
							;
							/* L30: */
						}
					}
					for (k = *m; k >= 1; --k) {
						if (b[k + j * b_dim1] != 0.) {
							if (nounit) {
								b[k + j * b_dim1] /= a[k + k * a_dim1];
							}
							i__2 = k - 1;
							for (i__ = 1; i__ <= i__2; ++i__) {
								b[i__ + j * b_dim1] -= b[k + j * b_dim1] * a[
									i__ + k * a_dim1];
									/* L40: */
							}
						}
						/* L50: */
					}
					/* L60: */
				}
			} else {
				i__1 = *n;
				for (j = 1; j <= i__1; ++j) {
					if (*alpha != 1.) {
						i__2 = *m;
						for (i__ = 1; i__ <= i__2; ++i__) {
							b[i__ + j * b_dim1] = *alpha * b[i__ + j * b_dim1]
							;
							/* L70: */
						}
					}
					i__2 = *m;
					for (k = 1; k <= i__2; ++k) {
						if (b[k + j * b_dim1] != 0.) {
							if (nounit) {
								b[k + j * b_dim1] /= a[k + k * a_dim1];
							}
							i__3 = *m;
							for (i__ = k + 1; i__ <= i__3; ++i__) {
								b[i__ + j * b_dim1] -= b[k + j * b_dim1] * a[
									i__ + k * a_dim1];
									/* L80: */
							}
						}
						/* L90: */
					}
					/* L100: */
				}
			}
		} else {
			/*           Form  B := alpha*inv( A' )*B. */
			if (upper) {
				i__1 = *n;
				for (j = 1; j <= i__1; ++j) {
					i__2 = *m;
					for (i__ = 1; i__ <= i__2; ++i__) {
						temp = *alpha * b[i__ + j * b_dim1];
						i__3 = i__ - 1;
						for (k = 1; k <= i__3; ++k) {
							temp -= a[k + i__ * a_dim1] * b[k + j * b_dim1];
							/* L110: */
						}
						if (nounit) {
							temp /= a[i__ + i__ * a_dim1];
						}
						b[i__ + j * b_dim1] = temp;
						/* L120: */
					}
					/* L130: */
				}
			} else {
				i__1 = *n;
				for (j = 1; j <= i__1; ++j) {
					for (i__ = *m; i__ >= 1; --i__) {
						temp = *alpha * b[i__ + j * b_dim1];
						i__2 = *m;
						for (k = i__ + 1; k <= i__2; ++k) {
							temp -= a[k + i__ * a_dim1] * b[k + j * b_dim1];
							/* L140: */
						}
						if (nounit) {
							temp /= a[i__ + i__ * a_dim1];
						}
						b[i__ + j * b_dim1] = temp;
						/* L150: */
					}
					/* L160: */
				}
			}
		}
	} else {
		if (lsame_(transa, "N")) {
			/*           Form  B := alpha*B*inv( A ). */
			if (upper) {
				i__1 = *n;
				for (j = 1; j <= i__1; ++j) {
					if (*alpha != 1.) {
						i__2 = *m;
						for (i__ = 1; i__ <= i__2; ++i__) {
							b[i__ + j * b_dim1] = *alpha * b[i__ + j * b_dim1]
							;
							/* L170: */
						}
					}
					i__2 = j - 1;
					for (k = 1; k <= i__2; ++k) {
						if (a[k + j * a_dim1] != 0.) {
							i__3 = *m;
							for (i__ = 1; i__ <= i__3; ++i__) {
								b[i__ + j * b_dim1] -= a[k + j * a_dim1] * b[
									i__ + k * b_dim1];
									/* L180: */
							}
						}
						/* L190: */
					}
					if (nounit) {
						temp = 1. / a[j + j * a_dim1];
						i__2 = *m;
						for (i__ = 1; i__ <= i__2; ++i__) {
							b[i__ + j * b_dim1] = temp * b[i__ + j * b_dim1];
							/* L200: */
						}
					}
					/* L210: */
				}
			} else {
				for (j = *n; j >= 1; --j) {
					if (*alpha != 1.) {
						i__1 = *m;
						for (i__ = 1; i__ <= i__1; ++i__) {
							b[i__ + j * b_dim1] = *alpha * b[i__ + j * b_dim1]
							;
							/* L220: */
						}
					}
					i__1 = *n;
					for (k = j + 1; k <= i__1; ++k) {
						if (a[k + j * a_dim1] != 0.) {
							i__2 = *m;
							for (i__ = 1; i__ <= i__2; ++i__) {
								b[i__ + j * b_dim1] -= a[k + j * a_dim1] * b[
									i__ + k * b_dim1];
									/* L230: */
							}
						}
						/* L240: */
					}
					if (nounit) {
						temp = 1. / a[j + j * a_dim1];
						i__1 = *m;
						for (i__ = 1; i__ <= i__1; ++i__) {
							b[i__ + j * b_dim1] = temp * b[i__ + j * b_dim1];
							/* L250: */
						}
					}
					/* L260: */
				}
			}
		} else {
			/*           Form  B := alpha*B*inv( A' ). */
			if (upper) {
				for (k = *n; k >= 1; --k) {
					if (nounit) {
						temp = 1. / a[k + k * a_dim1];
						i__1 = *m;
						for (i__ = 1; i__ <= i__1; ++i__) {
							b[i__ + k * b_dim1] = temp * b[i__ + k * b_dim1];
							/* L270: */
						}
					}
					i__1 = k - 1;
					for (j = 1; j <= i__1; ++j) {
						if (a[j + k * a_dim1] != 0.) {
							temp = a[j + k * a_dim1];
							i__2 = *m;
							for (i__ = 1; i__ <= i__2; ++i__) {
								b[i__ + j * b_dim1] -= temp * b[i__ + k * 
									b_dim1];
								/* L280: */
							}
						}
						/* L290: */
					}
					if (*alpha != 1.) {
						i__1 = *m;
						for (i__ = 1; i__ <= i__1; ++i__) {
							b[i__ + k * b_dim1] = *alpha * b[i__ + k * b_dim1]
							;
							/* L300: */
						}
					}
					/* L310: */
				}
			} else {
				i__1 = *n;
				for (k = 1; k <= i__1; ++k) {
					if (nounit) {
						temp = 1. / a[k + k * a_dim1];
						i__2 = *m;
						for (i__ = 1; i__ <= i__2; ++i__) {
							b[i__ + k * b_dim1] = temp * b[i__ + k * b_dim1];
							/* L320: */
						}
					}
					i__2 = *n;
					for (j = k + 1; j <= i__2; ++j) {
						if (a[j + k * a_dim1] != 0.) {
							temp = a[j + k * a_dim1];
							i__3 = *m;
							for (i__ = 1; i__ <= i__3; ++i__) {
								b[i__ + j * b_dim1] -= temp * b[i__ + k * 
									b_dim1];
								/* L330: */
							}
						}
						/* L340: */
					}
					if (*alpha != 1.) {
						i__2 = *m;
						for (i__ = 1; i__ <= i__2; ++i__) {
							b[i__ + k * b_dim1] = *alpha * b[i__ + k * b_dim1]
							;
							/* L350: */
						}
					}
					/* L360: */
				}
			}
		}
	}
	return 0;
	/*     End of DTRSM . */
}


//
int TRL::dsyrk_(char *uplo, char *trans, integer_ *n, integer_ *k, 
		   doublereal_ *alpha, doublereal_ *a, integer_ *lda, doublereal_ *beta, 
		   doublereal_ *c__, integer_ *ldc)
{
	/* System generated locals */
	integer_ a_dim1, a_offset, c_dim1, c_offset, i__1, i__2, i__3;
	/* Local variables */
	static integer_ i__, j, l, info;
	static doublereal_ temp;
	//extern logical_ lsame_(char *, char *);
	static integer_ nrowa;
	static logical_ upper;
	//extern /* Subroutine */ int xerbla_(char *, integer_ *);
	/*  Purpose   
	=======   
	DSYRK  performs one of the symmetric rank k operations   
	C := alpha*A*A' + beta*C,   
	or   
	C := alpha*A'*A + beta*C,   
	where  alpha and beta  are scalars, C is an  n by n  symmetric matrix   
	and  A  is an  n by k  matrix in the first case and a  k by n  matrix   
	in the second case.   
	Arguments   
	==========   
	UPLO   - CHARACTER*1.   
	On  entry,   UPLO  specifies  whether  the  upper  or  lower   
	triangular  part  of the  array  C  is to be  referenced  as   
	follows:   
	UPLO = 'U' or 'u'   Only the  upper triangular part of  C   
	is to be referenced.   
	UPLO = 'L' or 'l'   Only the  lower triangular part of  C   
	is to be referenced.   
	Unchanged on exit.   
	TRANS  - CHARACTER*1.   
	On entry,  TRANS  specifies the operation to be performed as   
	follows:   
	TRANS = 'N' or 'n'   C := alpha*A*A' + beta*C.   
	TRANS = 'T' or 't'   C := alpha*A'*A + beta*C.   
	TRANS = 'C' or 'c'   C := alpha*A'*A + beta*C.   
	Unchanged on exit.   
	N      - INTEGER.   
	On entry,  N specifies the order of the matrix C.  N must be   
	at least zero.   
	Unchanged on exit.   
	K      - INTEGER.   
	On entry with  TRANS = 'N' or 'n',  K  specifies  the number   
	of  columns   of  the   matrix   A,   and  on   entry   with   
	TRANS = 'T' or 't' or 'C' or 'c',  K  specifies  the  number   
	of rows of the matrix  A.  K must be at least zero.   
	Unchanged on exit.   
	ALPHA  - DOUBLE PRECISION.   
	On entry, ALPHA specifies the scalar alpha.   
	Unchanged on exit.   
	A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is   
	k  when  TRANS = 'N' or 'n',  and is  n  otherwise.   
	Before entry with  TRANS = 'N' or 'n',  the  leading  n by k   
	part of the array  A  must contain the matrix  A,  otherwise   
	the leading  k by n  part of the array  A  must contain  the   
	matrix A.   
	Unchanged on exit.   
	LDA    - INTEGER.   
	On entry, LDA specifies the first dimension of A as declared   
	in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'   
	then  LDA must be at least  max( 1, n ), otherwise  LDA must   
	be at least  max( 1, k ).   
	Unchanged on exit.   
	BETA   - DOUBLE PRECISION.   
	On entry, BETA specifies the scalar beta.   
	Unchanged on exit.   
	C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).   
	Before entry  with  UPLO = 'U' or 'u',  the leading  n by n   
	upper triangular part of the array C must contain the upper   
	triangular part  of the  symmetric matrix  and the strictly   
	lower triangular part of C is not referenced.  On exit, the   
	upper triangular part of the array  C is overwritten by the   
	upper triangular part of the updated matrix.   
	Before entry  with  UPLO = 'L' or 'l',  the leading  n by n   
	lower triangular part of the array C must contain the lower   
	triangular part  of the  symmetric matrix  and the strictly   
	upper triangular part of C is not referenced.  On exit, the   
	lower triangular part of the array  C is overwritten by the   
	lower triangular part of the updated matrix.   
	LDC    - INTEGER.   
	On entry, LDC specifies the first dimension of C as declared   
	in  the  calling  (sub)  program.   LDC  must  be  at  least   
	max( 1, n ).   
	Unchanged on exit.   
	Level 3 Blas routine.   
	-- Written on 8-February-1989.   
	Jack Dongarra, Argonne National Laboratory.   
	Iain Duff, AERE Harwell.   
	Jeremy Du Croz, Numerical Algorithms Group Ltd.   
	Sven Hammarling, Numerical Algorithms Group Ltd.   
	Test the input parameters.   
	Parameter adjustments */
	a_dim1 = *lda;
	a_offset = 1 + a_dim1;
	a -= a_offset;
	c_dim1 = *ldc;
	c_offset = 1 + c_dim1;
	c__ -= c_offset;
	/* Function Body */
	if (lsame_(trans, "N")) {
		nrowa = *n;
	} else {
		nrowa = *k;
	}
	upper = lsame_(uplo, "U");
	info = 0;
	if (! upper && ! lsame_(uplo, "L")) {
		info = 1;
	} else if (! lsame_(trans, "N") && ! lsame_(trans, 
		"T") && ! lsame_(trans, "C")) {
			info = 2;
	} else if (*n < 0) {
		info = 3;
	} else if (*k < 0) {
		info = 4;
	} else if (*lda < max(1,nrowa)) {
		info = 7;
	} else if (*ldc < max(1,*n)) {
		info = 10;
	}
	if (info != 0) {
		xerbla_("DSYRK ", &info);
		return 0;
	}
	/*     Quick return if possible. */
	if (*n == 0 || (*alpha == 0. || *k == 0) && *beta == 1.) {
		return 0;
	}
	/*     And when  alpha.eq.zero. */
	if (*alpha == 0.) {
		if (upper) {
			if (*beta == 0.) {
				i__1 = *n;
				for (j = 1; j <= i__1; ++j) {
					i__2 = j;
					for (i__ = 1; i__ <= i__2; ++i__) {
						c__[i__ + j * c_dim1] = 0.;
						/* L10: */
					}
					/* L20: */
				}
			} else {
				i__1 = *n;
				for (j = 1; j <= i__1; ++j) {
					i__2 = j;
					for (i__ = 1; i__ <= i__2; ++i__) {
						c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
						/* L30: */
					}
					/* L40: */
				}
			}
		} else {
			if (*beta == 0.) {
				i__1 = *n;
				for (j = 1; j <= i__1; ++j) {
					i__2 = *n;
					for (i__ = j; i__ <= i__2; ++i__) {
						c__[i__ + j * c_dim1] = 0.;
						/* L50: */
					}
					/* L60: */
				}
			} else {
				i__1 = *n;
				for (j = 1; j <= i__1; ++j) {
					i__2 = *n;
					for (i__ = j; i__ <= i__2; ++i__) {
						c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
						/* L70: */
					}
					/* L80: */
				}
			}
		}
		return 0;
	}
	/*     Start the operations. */
	if (lsame_(trans, "N")) {
		/*        Form  C := alpha*A*A' + beta*C. */
		if (upper) {
			i__1 = *n;
			for (j = 1; j <= i__1; ++j) {
				if (*beta == 0.) {
					i__2 = j;
					for (i__ = 1; i__ <= i__2; ++i__) {
						c__[i__ + j * c_dim1] = 0.;
						/* L90: */
					}
				} else if (*beta != 1.) {
					i__2 = j;
					for (i__ = 1; i__ <= i__2; ++i__) {
						c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
						/* L100: */
					}
				}
				i__2 = *k;
				for (l = 1; l <= i__2; ++l) {
					if (a[j + l * a_dim1] != 0.) {
						temp = *alpha * a[j + l * a_dim1];
						i__3 = j;
						for (i__ = 1; i__ <= i__3; ++i__) {
							c__[i__ + j * c_dim1] += temp * a[i__ + l * 
								a_dim1];
							/* L110: */
						}
					}
					/* L120: */
				}
				/* L130: */
			}
		} else {
			i__1 = *n;
			for (j = 1; j <= i__1; ++j) {
				if (*beta == 0.) {
					i__2 = *n;
					for (i__ = j; i__ <= i__2; ++i__) {
						c__[i__ + j * c_dim1] = 0.;
						/* L140: */
					}
				} else if (*beta != 1.) {
					i__2 = *n;
					for (i__ = j; i__ <= i__2; ++i__) {
						c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
						/* L150: */
					}
				}
				i__2 = *k;
				for (l = 1; l <= i__2; ++l) {
					if (a[j + l * a_dim1] != 0.) {
						temp = *alpha * a[j + l * a_dim1];
						i__3 = *n;
						for (i__ = j; i__ <= i__3; ++i__) {
							c__[i__ + j * c_dim1] += temp * a[i__ + l * 
								a_dim1];
							/* L160: */
						}
					}
					/* L170: */
				}
				/* L180: */
			}
		}
	} else {
		/*        Form  C := alpha*A'*A + beta*C. */
		if (upper) {
			i__1 = *n;
			for (j = 1; j <= i__1; ++j) {
				i__2 = j;
				for (i__ = 1; i__ <= i__2; ++i__) {
					temp = 0.;
					i__3 = *k;
					for (l = 1; l <= i__3; ++l) {
						temp += a[l + i__ * a_dim1] * a[l + j * a_dim1];
						/* L190: */
					}
					if (*beta == 0.) {
						c__[i__ + j * c_dim1] = *alpha * temp;
					} else {
						c__[i__ + j * c_dim1] = *alpha * temp + *beta * c__[
							i__ + j * c_dim1];
					}
					/* L200: */
				}
				/* L210: */
			}
		} else {
			i__1 = *n;
			for (j = 1; j <= i__1; ++j) {
				i__2 = *n;
				for (i__ = j; i__ <= i__2; ++i__) {
					temp = 0.;
					i__3 = *k;
					for (l = 1; l <= i__3; ++l) {
						temp += a[l + i__ * a_dim1] * a[l + j * a_dim1];
						/* L220: */
					}
					if (*beta == 0.) {
						c__[i__ + j * c_dim1] = *alpha * temp;
					} else {
						c__[i__ + j * c_dim1] = *alpha * temp + *beta * c__[
							i__ + j * c_dim1];
					}
					/* L230: */
				}
				/* L240: */
			}
		}
	}
	return 0;
	/*     End of DSYRK . */
} /* dsyrk_ */


//
int TRL::dpotf2_(char *uplo, integer_ *n, doublereal_ *a, integer_ *
			lda, integer_ *info)
{
	/*  -- LAPACK routine (version 3.1) --   
	Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..   
	November 2006   


	Purpose   
	=======   

	DPOTF2 computes the Cholesky factorization of a trl_real_ symmetric   
	positive definite matrix A.   

	The factorization has the form   
	A = U' * U ,  if UPLO = 'U', or   
	A = L  * L',  if UPLO = 'L',   
	where U is an upper triangular matrix and L is lower triangular.   

	This is the unblocked version of the algorithm, calling Level 2 BLAS.   

	Arguments   
	=========   

	UPLO    (input) CHARACTER*1   
	Specifies whether the upper or lower triangular part of the   
	symmetric matrix A is stored.   
	= 'U':  Upper triangular   
	= 'L':  Lower triangular   

	N       (input) INTEGER   
	The order of the matrix A.  N >= 0.   

	A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
	On entry, the symmetric matrix A.  If UPLO = 'U', the leading   
	n by n upper triangular part of A contains the upper   
	triangular part of the matrix A, and the strictly lower   
	triangular part of A is not referenced.  If UPLO = 'L', the   
	leading n by n lower triangular part of A contains the lower   
	triangular part of the matrix A, and the strictly upper   
	triangular part of A is not referenced.   

	On exit, if INFO = 0, the factor U or L from the Cholesky   
	factorization A = U'*U  or A = L*L'.   

	LDA     (input) INTEGER   
	The leading dimension of the array A.  LDA >= max(1,N).   

	INFO    (output) INTEGER   
	= 0: successful exit   
	< 0: if INFO = -k, the k-th argument had an illegal value   
	> 0: if INFO = k, the leading minor of order k is not   
	positive definite, and the factorization could not be   
	completed.   

	=====================================================================   


	Test the input parameters.   

	Parameter adjustments */
	/* Table of constant values */
	static integer_ c__1 = 1;
	static doublereal_ c_b10 = -1.;
	static doublereal_ c_b12 = 1.;

	/* System generated locals */
	integer_ a_dim1, a_offset, i__1, i__2, i__3;
	doublereal_ d__1;
	/* Builtin functions */
	//double sqrt_(doublereal_);
	/* Local variables */
	static integer_ j;
	static doublereal_ ajj;
	//extern doublereal_ ddot_(integer_ *, doublereal_ *, integer_ *, doublereal_ *, 
	//	integer_ *);
	//extern /* Subroutine */ int dscal_(integer_ *, doublereal_ *, doublereal_ *, 
	//	integer_ *);
	//extern logical_ lsame_(char *, char *);
	//extern /* Subroutine */ int dgemv_(char *, integer_ *, integer_ *, 
	//	doublereal_ *, doublereal_ *, integer_ *, doublereal_ *, integer_ *, 
	//	doublereal_ *, doublereal_ *, integer_ *);
	static logical_ upper;
	//extern /* Subroutine */ int xerbla_(char *, integer_ *);


	a_dim1 = *lda;
	a_offset = 1 + a_dim1;
	a -= a_offset;

	/* Function Body */
	*info = 0;
	upper = lsame_(uplo, "U");
	if (! upper && ! lsame_(uplo, "L")) {
		*info = -1;
	} else if (*n < 0) {
		*info = -2;
	} else if (*lda < max(1,*n)) {
		*info = -4;
	}
	if (*info != 0) {
		i__1 = -(*info);
		xerbla_("DPOTF2", &i__1);
		return 0;
	}

	/*     Quick return if possible */

	if (*n == 0) {
		return 0;
	}

	if (upper) {

		/*        Compute the Cholesky factorization A = U'*U. */

		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {

			/*           Compute U(J,J) and test for non-positive-definiteness. */

			i__2 = j - 1;
			ajj = a[j + j * a_dim1] - ddot_(&i__2, &a[j * a_dim1 + 1], &c__1, 
				&a[j * a_dim1 + 1], &c__1);
			if (ajj <= 0.) {
				a[j + j * a_dim1] = ajj;
				goto L30;
			}
			ajj = sqrt_(ajj);
			a[j + j * a_dim1] = ajj;

			/*           Compute elements J+1:N of row J. */

			if (j < *n) {
				i__2 = j - 1;
				i__3 = *n - j;
				dgemv_("Transpose", &i__2, &i__3, &c_b10, &a[(j + 1) * a_dim1 
					+ 1], lda, &a[j * a_dim1 + 1], &c__1, &c_b12, &a[j + (
					j + 1) * a_dim1], lda);
				i__2 = *n - j;
				d__1 = 1. / ajj;
				dscal_(&i__2, &d__1, &a[j + (j + 1) * a_dim1], lda);
			}
			/* L10: */
		}
	} else {

		/*        Compute the Cholesky factorization A = L*L'. */

		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {

			/*           Compute L(J,J) and test for non-positive-definiteness. */

			i__2 = j - 1;
			ajj = a[j + j * a_dim1] - ddot_(&i__2, &a[j + a_dim1], lda, &a[j 
				+ a_dim1], lda);
			if (ajj <= 0.) {
				a[j + j * a_dim1] = ajj;
				goto L30;
			}
			ajj = sqrt_(ajj);
			a[j + j * a_dim1] = ajj;

			/*           Compute elements J+1:N of column J. */

			if (j < *n) {
				i__2 = *n - j;
				i__3 = j - 1;
				dgemv_("No transpose", &i__2, &i__3, &c_b10, &a[j + 1 + 
					a_dim1], lda, &a[j + a_dim1], lda, &c_b12, &a[j + 1 + 
					j * a_dim1], &c__1);
				i__2 = *n - j;
				d__1 = 1. / ajj;
				dscal_(&i__2, &d__1, &a[j + 1 + j * a_dim1], &c__1);
			}
			/* L20: */
		}
	}
	goto L40;

L30:
	*info = j;

L40:
	return 0;

	/*     End of DPOTF2 */

} /* dpotf2_ */

//
int TRL::dorgql_(integer_ *m, integer_ *n, integer_ *k, doublereal_ *
			a, integer_ *lda, doublereal_ *tau, doublereal_ *work, integer_ *lwork, 
			integer_ *info)
{
	/*  -- LAPACK routine (version 3.1) --   
	Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..   
	November 2006   


	Purpose   
	=======   

	DORGQL generates an M-by-N trl_real_ matrix Q with orthonormal columns,   
	which is defined as the last N columns of a product of K elementary   
	reflectors of order M   

	Q  =  H(k) . . . H(2) H(1)   

	as returned by DGEQLF.   

	Arguments   
	=========   

	M       (input) INTEGER   
	The number of rows of the matrix Q. M >= 0.   

	N       (input) INTEGER   
	The number of columns of the matrix Q. M >= N >= 0.   

	K       (input) INTEGER   
	The number of elementary reflectors whose product defines the   
	matrix Q. N >= K >= 0.   

	A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
	On entry, the (n-k+i)-th column must contain the vector which   
	defines the elementary reflector H(i), for i = 1,2,...,k, as   
	returned by DGEQLF in the last k columns of its array   
	argument A.   
	On exit, the M-by-N matrix Q.   

	LDA     (input) INTEGER   
	The first dimension of the array A. LDA >= max(1,M).   

	TAU     (input) DOUBLE PRECISION array, dimension (K)   
	TAU(i) must contain the scalar factor of the elementary   
	reflector H(i), as returned by DGEQLF.   

	WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))   
	On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

	LWORK   (input) INTEGER   
	The dimension of the array WORK. LWORK >= max(1,N).   
	For optimum performance LWORK >= N*NB, where NB is the   
	optimal blocksize.   

	If LWORK = -1, then a workspace query is assumed; the routine   
	only calculates the optimal size of the WORK array, returns   
	this value as the first entry of the WORK array, and no error   
	message related to LWORK is issued by XERBLA.   

	INFO    (output) INTEGER   
	= 0:  successful exit   
	< 0:  if INFO = -i, the i-th argument has an illegal value   

	=====================================================================   


	Test the input arguments   

	Parameter adjustments */
	/* Table of constant values */
	static integer_ c__1 = 1;
	static integer_ c_n1 = -1;
	static integer_ c__3 = 3;
	static integer_ c__2 = 2;

	/* System generated locals */
	integer_ a_dim1, a_offset, i__1, i__2, i__3, i__4;
	/* Local variables */
	static integer_ i__, j, l, ib, nb, kk, nx, iws, nbmin, iinfo;
	//extern /* Subroutine */ int dorg2l_(integer_ *, integer_ *, integer_ *, 
	//	doublereal_ *, integer_ *, doublereal_ *, doublereal_ *, integer_ *), 
	//	dlarfb_(char *, char *, char *, char *, integer_ *, integer_ *, 
	//	integer_ *, doublereal_ *, integer_ *, doublereal_ *, integer_ *, 
	//	doublereal_ *, integer_ *, doublereal_ *, integer_ *), dlarft_(char *, char *, integer_ *, integer_ *, 
	//	doublereal_ *, integer_ *, doublereal_ *, doublereal_ *, integer_ *), xerbla_(char *, integer_ *);
	//extern integer_ ilaenv_(integer_ *, char *, char *, integer_ *, integer_ *, 
	//	integer_ *, integer_ *, ftnlen_, ftnlen_);
	static integer_ ldwork, lwkopt;
	static logical_ lquery;


	a_dim1 = *lda;
	a_offset = 1 + a_dim1;
	a -= a_offset;
	--tau;
	--work;

	/* Function Body */
	*info = 0;
	lquery = *lwork == -1;
	if (*m < 0) {
		*info = -1;
	} else if (*n < 0 || *n > *m) {
		*info = -2;
	} else if (*k < 0 || *k > *n) {
		*info = -3;
	} else if (*lda < max(1,*m)) {
		*info = -5;
	}

	if (*info == 0) {
		if (*n == 0) {
			lwkopt = 1;
		} else {
			nb = ilaenv_(&c__1, "DORGQL", " ", m, n, k, &c_n1, (ftnlen_)6, (
				ftnlen_)1);
			lwkopt = *n * nb;
		}
		work[1] = (doublereal_) lwkopt;

		if (*lwork < max(1,*n) && ! lquery) {
			*info = -8;
		}
	}

	if (*info != 0) {
		i__1 = -(*info);
		xerbla_("DORGQL", &i__1);
		return 0;
	} else if (lquery) {
		return 0;
	}

	/*     Quick return if possible */

	if (*n <= 0) {
		return 0;
	}

	nbmin = 2;
	nx = 0;
	iws = *n;
	if (nb > 1 && nb < *k) {

		/*        Determine when to cross over from blocked to unblocked code.   

		Computing MAX */
		i__1 = 0, i__2 = ilaenv_(&c__3, "DORGQL", " ", m, n, k, &c_n1, (
			ftnlen_)6, (ftnlen_)1);
		nx = max(i__1,i__2);
		if (nx < *k) {

			/*           Determine if workspace is large enough for blocked code. */

			ldwork = *n;
			iws = ldwork * nb;
			if (*lwork < iws) {

				/*              Not enough workspace to use optimal NB:  reduce NB and   
				determine the minimum value of NB. */

				nb = *lwork / ldwork;
				/* Computing MAX */
				i__1 = 2, i__2 = ilaenv_(&c__2, "DORGQL", " ", m, n, k, &c_n1,
					(ftnlen_)6, (ftnlen_)1);
				nbmin = max(i__1,i__2);
			}
		}
	}

	if (nb >= nbmin && nb < *k && nx < *k) {

		/*        Use blocked code after the first block.   
		The last kk columns are handled by the block method.   

		Computing MIN */
		i__1 = *k, i__2 = (*k - nx + nb - 1) / nb * nb;
		kk = min(i__1,i__2);

		/*        Set A(m-kk+1:m,1:n-kk) to zero. */

		i__1 = *n - kk;
		for (j = 1; j <= i__1; ++j) {
			i__2 = *m;
			for (i__ = *m - kk + 1; i__ <= i__2; ++i__) {
				a[i__ + j * a_dim1] = 0.;
				/* L10: */
			}
			/* L20: */
		}
	} else {
		kk = 0;
	}

	/*     Use unblocked code for the first or only block. */

	i__1 = *m - kk;
	i__2 = *n - kk;
	i__3 = *k - kk;
	dorg2l_(&i__1, &i__2, &i__3, &a[a_offset], lda, &tau[1], &work[1], &iinfo)
		;

	if (kk > 0) {

		/*        Use blocked code */

		i__1 = *k;
		i__2 = nb;
		for (i__ = *k - kk + 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += 
			i__2) {
				/* Computing MIN */
				i__3 = nb, i__4 = *k - i__ + 1;
				ib = min(i__3,i__4);
				if (*n - *k + i__ > 1) {

					/*              Form the triangular factor of the block reflector   
					H = H(i+ib-1) . . . H(i+1) H(i) */

					i__3 = *m - *k + i__ + ib - 1;
					dlarft_("Backward", "Columnwise", &i__3, &ib, &a[(*n - *k + 
						i__) * a_dim1 + 1], lda, &tau[i__], &work[1], &ldwork);

					/*              Apply H to A(1:m-k+i+ib-1,1:n-k+i-1) from the left */

					i__3 = *m - *k + i__ + ib - 1;
					i__4 = *n - *k + i__ - 1;
					dlarfb_("Left", "No transpose", "Backward", "Columnwise", &
						i__3, &i__4, &ib, &a[(*n - *k + i__) * a_dim1 + 1], 
						lda, &work[1], &ldwork, &a[a_offset], lda, &work[ib + 
						1], &ldwork);
				}

				/*           Apply H to rows 1:m-k+i+ib-1 of current block */

				i__3 = *m - *k + i__ + ib - 1;
				dorg2l_(&i__3, &ib, &ib, &a[(*n - *k + i__) * a_dim1 + 1], lda, &
					tau[i__], &work[1], &iinfo);

				/*           Set rows m-k+i+ib:m of current block to zero */

				i__3 = *n - *k + i__ + ib - 1;
				for (j = *n - *k + i__; j <= i__3; ++j) {
					i__4 = *m;
					for (l = *m - *k + i__ + ib; l <= i__4; ++l) {
						a[l + j * a_dim1] = 0.;
						/* L30: */
					}
					/* L40: */
				}
				/* L50: */
		}
	}

	work[1] = (doublereal_) iws;
	return 0;

	/*     End of DORGQL */

} /* dorgql_ */

//
int TRL::dorgqr_(integer_ *m, integer_ *n, integer_ *k, doublereal_ *
			a, integer_ *lda, doublereal_ *tau, doublereal_ *work, integer_ *lwork, 
			integer_ *info)
{
	/*  -- LAPACK routine (version 3.1) --   
	Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..   
	November 2006   


	Purpose   
	=======   

	DORGQR generates an M-by-N trl_real_ matrix Q with orthonormal columns,   
	which is defined as the first N columns of a product of K elementary   
	reflectors of order M   

	Q  =  H(1) H(2) . . . H(k)   

	as returned by DGEQRF.   

	Arguments   
	=========   

	M       (input) INTEGER   
	The number of rows of the matrix Q. M >= 0.   

	N       (input) INTEGER   
	The number of columns of the matrix Q. M >= N >= 0.   

	K       (input) INTEGER   
	The number of elementary reflectors whose product defines the   
	matrix Q. N >= K >= 0.   

	A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
	On entry, the i-th column must contain the vector which   
	defines the elementary reflector H(i), for i = 1,2,...,k, as   
	returned by DGEQRF in the first k columns of its array   
	argument A.   
	On exit, the M-by-N matrix Q.   

	LDA     (input) INTEGER   
	The first dimension of the array A. LDA >= max(1,M).   

	TAU     (input) DOUBLE PRECISION array, dimension (K)   
	TAU(i) must contain the scalar factor of the elementary   
	reflector H(i), as returned by DGEQRF.   

	WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))   
	On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

	LWORK   (input) INTEGER   
	The dimension of the array WORK. LWORK >= max(1,N).   
	For optimum performance LWORK >= N*NB, where NB is the   
	optimal blocksize.   

	If LWORK = -1, then a workspace query is assumed; the routine   
	only calculates the optimal size of the WORK array, returns   
	this value as the first entry of the WORK array, and no error   
	message related to LWORK is issued by XERBLA.   

	INFO    (output) INTEGER   
	= 0:  successful exit   
	< 0:  if INFO = -i, the i-th argument has an illegal value   

	=====================================================================   


	Test the input arguments   

	Parameter adjustments */
	/* Table of constant values */
	static integer_ c__1 = 1;
	static integer_ c_n1 = -1;
	static integer_ c__3 = 3;
	static integer_ c__2 = 2;

	/* System generated locals */
	integer_ a_dim1, a_offset, i__1, i__2, i__3;
	/* Local variables */
	static integer_ i__, j, l, ib, nb, ki, kk, nx, iws, nbmin, iinfo;
	//extern /* Subroutine */ int dorg2r_(integer_ *, integer_ *, integer_ *, 
	//	doublereal_ *, integer_ *, doublereal_ *, doublereal_ *, integer_ *), 
	//	dlarfb_(char *, char *, char *, char *, integer_ *, integer_ *, 
	//	integer_ *, doublereal_ *, integer_ *, doublereal_ *, integer_ *, 
	//	doublereal_ *, integer_ *, doublereal_ *, integer_ *), dlarft_(char *, char *, integer_ *, integer_ *, 
	//	doublereal_ *, integer_ *, doublereal_ *, doublereal_ *, integer_ *), xerbla_(char *, integer_ *);
	//extern integer_ ilaenv_(integer_ *, char *, char *, integer_ *, integer_ *, 
	//	integer_ *, integer_ *, ftnlen_, ftnlen_);
	static integer_ ldwork, lwkopt;
	static logical_ lquery;


	a_dim1 = *lda;
	a_offset = 1 + a_dim1;
	a -= a_offset;
	--tau;
	--work;

	/* Function Body */
	*info = 0;
	nb = ilaenv_(&c__1, "DORGQR", " ", m, n, k, &c_n1, (ftnlen_)6, (ftnlen_)1);
	lwkopt = max(1,*n) * nb;
	work[1] = (doublereal_) lwkopt;
	lquery = *lwork == -1;
	if (*m < 0) {
		*info = -1;
	} else if (*n < 0 || *n > *m) {
		*info = -2;
	} else if (*k < 0 || *k > *n) {
		*info = -3;
	} else if (*lda < max(1,*m)) {
		*info = -5;
	} else if (*lwork < max(1,*n) && ! lquery) {
		*info = -8;
	}
	if (*info != 0) {
		i__1 = -(*info);
		xerbla_("DORGQR", &i__1);
		return 0;
	} else if (lquery) {
		return 0;
	}

	/*     Quick return if possible */

	if (*n <= 0) {
		work[1] = 1.;
		return 0;
	}

	nbmin = 2;
	nx = 0;
	iws = *n;
	if (nb > 1 && nb < *k) {

		/*        Determine when to cross over from blocked to unblocked code.   

		Computing MAX */
		i__1 = 0, i__2 = ilaenv_(&c__3, "DORGQR", " ", m, n, k, &c_n1, (
			ftnlen_)6, (ftnlen_)1);
		nx = max(i__1,i__2);
		if (nx < *k) {

			/*           Determine if workspace is large enough for blocked code. */

			ldwork = *n;
			iws = ldwork * nb;
			if (*lwork < iws) {

				/*              Not enough workspace to use optimal NB:  reduce NB and   
				determine the minimum value of NB. */

				nb = *lwork / ldwork;
				/* Computing MAX */
				i__1 = 2, i__2 = ilaenv_(&c__2, "DORGQR", " ", m, n, k, &c_n1,
					(ftnlen_)6, (ftnlen_)1);
				nbmin = max(i__1,i__2);
			}
		}
	}

	if (nb >= nbmin && nb < *k && nx < *k) {

		/*        Use blocked code after the last block.   
		The first kk columns are handled by the block method. */

		ki = (*k - nx - 1) / nb * nb;
		/* Computing MIN */
		i__1 = *k, i__2 = ki + nb;
		kk = min(i__1,i__2);

		/*        Set A(1:kk,kk+1:n) to zero. */

		i__1 = *n;
		for (j = kk + 1; j <= i__1; ++j) {
			i__2 = kk;
			for (i__ = 1; i__ <= i__2; ++i__) {
				a[i__ + j * a_dim1] = 0.;
				/* L10: */
			}
			/* L20: */
		}
	} else {
		kk = 0;
	}

	/*     Use unblocked code for the last or only block. */

	if (kk < *n) {
		i__1 = *m - kk;
		i__2 = *n - kk;
		i__3 = *k - kk;
		dorg2r_(&i__1, &i__2, &i__3, &a[kk + 1 + (kk + 1) * a_dim1], lda, &
			tau[kk + 1], &work[1], &iinfo);
	}

	if (kk > 0) {

		/*        Use blocked code */

		i__1 = -nb;
		for (i__ = ki + 1; i__1 < 0 ? i__ >= 1 : i__ <= 1; i__ += i__1) {
			/* Computing MIN */
			i__2 = nb, i__3 = *k - i__ + 1;
			ib = min(i__2,i__3);
			if (i__ + ib <= *n) {

				/*              Form the triangular factor of the block reflector   
				H = H(i) H(i+1) . . . H(i+ib-1) */

				i__2 = *m - i__ + 1;
				dlarft_("Forward", "Columnwise", &i__2, &ib, &a[i__ + i__ * 
					a_dim1], lda, &tau[i__], &work[1], &ldwork);

				/*              Apply H to A(i:m,i+ib:n) from the left */

				i__2 = *m - i__ + 1;
				i__3 = *n - i__ - ib + 1;
				dlarfb_("Left", "No transpose", "Forward", "Columnwise", &
					i__2, &i__3, &ib, &a[i__ + i__ * a_dim1], lda, &work[
						1], &ldwork, &a[i__ + (i__ + ib) * a_dim1], lda, &
							work[ib + 1], &ldwork);
			}

			/*           Apply H to rows i:m of current block */

			i__2 = *m - i__ + 1;
			dorg2r_(&i__2, &ib, &ib, &a[i__ + i__ * a_dim1], lda, &tau[i__], &
				work[1], &iinfo);

			/*           Set rows 1:i-1 of current block to zero */

			i__2 = i__ + ib - 1;
			for (j = i__; j <= i__2; ++j) {
				i__3 = i__ - 1;
				for (l = 1; l <= i__3; ++l) {
					a[l + j * a_dim1] = 0.;
					/* L30: */
				}
				/* L40: */
			}
			/* L50: */
		}
	}

	work[1] = (doublereal_) iws;
	return 0;

	/*     End of DORGQR */

} /* dorgqr_ */

//
integer_ TRL::ieeeck_(integer_ *ispec, trl_real_ *zero, trl_real_ *one)
{
	/*  -- LAPACK auxiliary routine (version 3.1) --   
	Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..   
	November 2006   


	Purpose   
	=======   

	IEEECK is called from the ILAENV to verify that Infinity and   
	possibly NaN arithmetic is safe (i.e. will not trap).   

	Arguments   
	=========   

	ISPEC   (input) INTEGER   
	Specifies whether to test just for inifinity arithmetic   
	or whether to test for infinity and NaN arithmetic.   
	= 0: Verify infinity arithmetic only.   
	= 1: Verify infinity and NaN arithmetic.   

	ZERO    (input) REAL   
	Must contain the value 0.0   
	This is passed to prevent the compiler from optimizing   
	away this code.   

	ONE     (input) REAL   
	Must contain the value 1.0   
	This is passed to prevent the compiler from optimizing   
	away this code.   

	RETURN VALUE:  INTEGER   
	= 0:  Arithmetic failed to produce the correct answers   
	= 1:  Arithmetic produced the correct answers */
	/* System generated locals */
	integer_ ret_val;
	/* Local variables */
	static trl_real_ nan1, nan2, nan3, nan4, nan5, nan6, neginf, posinf, negzro, 
		newzro;


	ret_val = 1;

	posinf = *one / *zero;
	if (posinf <= *one) {
		ret_val = 0;
		return ret_val;
	}

	neginf = -(*one) / *zero;
	if (neginf >= *zero) {
		ret_val = 0;
		return ret_val;
	}

	negzro = *one / (neginf + *one);
	if (negzro != *zero) {
		ret_val = 0;
		return ret_val;
	}

	neginf = *one / negzro;
	if (neginf >= *zero) {
		ret_val = 0;
		return ret_val;
	}

	newzro = negzro + *zero;
	if (newzro != *zero) {
		ret_val = 0;
		return ret_val;
	}

	posinf = *one / newzro;
	if (posinf <= *one) {
		ret_val = 0;
		return ret_val;
	}

	neginf *= posinf;
	if (neginf >= *zero) {
		ret_val = 0;
		return ret_val;
	}

	posinf *= posinf;
	if (posinf <= *one) {
		ret_val = 0;
		return ret_val;
	}




	/*     Return if we were only asked to check infinity arithmetic */

	if (*ispec == 0) {
		return ret_val;
	}

	nan1 = posinf + neginf;

	nan2 = posinf / neginf;

	nan3 = posinf / posinf;

	nan4 = posinf * *zero;

	nan5 = neginf * negzro;

	nan6 = nan5 * 0.f;

	if (nan1 == nan1) {
		ret_val = 0;
		return ret_val;
	}

	if (nan2 == nan2) {
		ret_val = 0;
		return ret_val;
	}

	if (nan3 == nan3) {
		ret_val = 0;
		return ret_val;
	}

	if (nan4 == nan4) {
		ret_val = 0;
		return ret_val;
	}

	if (nan5 == nan5) {
		ret_val = 0;
		return ret_val;
	}

	if (nan6 == nan6) {
		ret_val = 0;
		return ret_val;
	}

	return ret_val;
} /* ieeeck_ */

//
integer_ TRL::iparmq_(integer_ *ispec, char *name__, char *opts, integer_ *n, integer_ 
					 *ilo, integer_ *ihi, integer_ *lwork, ftnlen_ name_len, ftnlen_ 
					 opts_len)
{
	/*  -- LAPACK auxiliary routine (version 3.1) --   
	Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..   
	November 2006   


	Purpose   
	=======   

	This program sets problem and machine dependent parameters   
	useful for xHSEQR and its subroutines. It is called whenever   
	ILAENV is called with 12 <= ISPEC <= 16   

	Arguments   
	=========   

	ISPEC  (input) integer_ scalar   
	ISPEC specifies which tunable parameter IPARMQ should   
	return.   

	ISPEC=12: (INMIN)  Matrices of order nmin or less   
	are sent directly to xLAHQR, the implicit   
	double shift QR algorithm.  NMIN must be   
	at least 11.   

	ISPEC=13: (INWIN)  Size of the deflation window.   
	This is best set greater than or equal to   
	the number of simultaneous shifts NS.   
	Larger matrices benefit from larger deflation   
	windows.   

	ISPEC=14: (INIBL) Determines when to stop nibbling and   
	invest in an (expensive) multi-shift QR sweep.   
	If the aggressive early deflation subroutine   
	finds LD converged eigenvalues from an order   
	NW deflation window and LD.GT.(NW*NIBBLE)/100,   
	then the next QR sweep is skipped and early   
	deflation is applied immediately to the   
	remaining active diagonal block.  Setting   
	IPARMQ(ISPEC=14) = 0 causes TTQRE to skip a   
	multi-shift QR sweep whenever early deflation   
	finds a converged eigenvalue.  Setting   
	IPARMQ(ISPEC=14) greater than or equal to 100   
	prevents TTQRE from skipping a multi-shift   
	QR sweep.   

	ISPEC=15: (NSHFTS) The number of simultaneous shifts in   
	a multi-shift QR iteration.   

	ISPEC=16: (IACC22) IPARMQ is set to 0, 1 or 2 with the   
	following meanings.   
	0:  During the multi-shift QR sweep,   
	xLAQR5 does not accumulate reflections and   
	does not use matrix-matrix multiply to   
	update the far-from-diagonal matrix   
	entries.   
	1:  During the multi-shift QR sweep,   
	xLAQR5 and/or xLAQRaccumulates reflections and uses   
	matrix-matrix multiply to update the   
	far-from-diagonal matrix entries.   
	2:  During the multi-shift QR sweep.   
	xLAQR5 accumulates reflections and takes   
	advantage of 2-by-2 block structure during   
	matrix-matrix multiplies.   
	(If xTRMM is slower than xGEMM, then   
	IPARMQ(ISPEC=16)=1 may be more efficient than   
	IPARMQ(ISPEC=16)=2 despite the greater level of   
	arithmetic work implied by the latter choice.)   

	NAME    (input) character string   
	Name of the calling subroutine   

	OPTS    (input) character string   
	This is a concatenation of the string arguments to   
	TTQRE.   

	N       (input) integer_ scalar   
	N is the order of the Hessenberg matrix H.   

	ILO     (input) INTEGER   
	IHI     (input) INTEGER   
	It is assumed that H is already upper triangular   
	in rows and columns 1:ILO-1 and IHI+1:N.   

	LWORK   (input) integer_ scalar   
	The amount of workspace available.   

	Further Details   
	===============   

	Little is known about how best to choose these parameters.   
	It is possible to use different values of the parameters   
	for each of CHSEQR, DHSEQR, SHSEQR and ZHSEQR.   

	It is probably best to choose different parameters for   
	different matrices and different parameters at different   
	times during the iteration, but this has not been   
	implemented --- yet.   


	The best choices of most of the parameters depend   
	in an ill-understood way on the relative execution   
	rate of xLAQR3 and xLAQR5 and on the nature of each   
	particular eigenvalue problem.  Experiment may be the   
	only practical way to determine which choices are most   
	effective.   

	Following is a list of default values supplied by IPARMQ.   
	These defaults may be adjusted in order to attain better   
	performance in any particular computational environment.   

	IPARMQ(ISPEC=12) The xLAHQR vs xLAQR0 crossover point.   
	Default: 75. (Must be at least 11.)   

	IPARMQ(ISPEC=13) Recommended deflation window size.   
	This depends on ILO, IHI and NS, the   
	number of simultaneous shifts returned   
	by IPARMQ(ISPEC=15).  The default for   
	(IHI-ILO+1).LE.500 is NS.  The default   
	for (IHI-ILO+1).GT.500 is 3*NS/2.   

	IPARMQ(ISPEC=14) Nibble crossover point.  Default: 14.   

	IPARMQ(ISPEC=15) Number of simultaneous shifts, NS.   
	a multi-shift QR iteration.   

	If IHI-ILO+1 is ...   

	greater than      ...but less    ... the   
	or equal to ...      than        default is   

	0               30       NS =   2+   
	30               60       NS =   4+   
	60              150       NS =  10   
	150              590       NS =  **   
	590             3000       NS =  64   
	3000             6000       NS = 128   
	6000             infinity   NS = 256   

	(+)  By default matrices of this order are   
	passed to the implicit double shift routine   
	xLAHQR.  See IPARMQ(ISPEC=12) above.   These   
	values of NS are used only in case of a rare   
	xLAHQR failure.   

	(**) The asterisks (**) indicate an ad-hoc   
	function increasing from 10 to 64.   

	IPARMQ(ISPEC=16) Select structured matrix multiply.   
	(See ISPEC=16 above for details.)   
	Default: 3.   

	================================================================ */
	/* System generated locals */
	integer_ ret_val, i__1, i__2;
	trl_real_ r__1;
	/* Builtin functions */
	//double log_(doublereal_);
	//integer_ i_nint(trl_real_ *);
	/* Local variables */
	static integer_ nh, ns;

	if (*ispec == 15 || *ispec == 13 || *ispec == 16) {

		/*        ==== Set the number simultaneous shifts ==== */

		nh = *ihi - *ilo + 1;
		ns = 2;
		if (nh >= 30) {
			ns = 4;
		}
		if (nh >= 60) {
			ns = 10;
		}
		if (nh >= 150) {
			/* Computing MAX */
			r__1 = log_((trl_real_) nh) / log_(2.f);
			i__1 = 10, i__2 = nh / i_nint(&r__1);
			ns = max(i__1,i__2);
		}
		if (nh >= 590) {
			ns = 64;
		}
		if (nh >= 3000) {
			ns = 128;
		}
		if (nh >= 6000) {
			ns = 256;
		}
		/* Computing MAX */
		i__1 = 2, i__2 = ns - ns % 2;
		ns = max(i__1,i__2);
	}

	if (*ispec == 12) {


		/*        ===== Matrices of order smaller than NMIN get sent   
		.     to xLAHQR, the classic double shift algorithm.   
		.     This must be at least 11. ==== */

		ret_val = 75;

	} else if (*ispec == 14) {

		/*        ==== INIBL: skip a multi-shift qr iteration and   
		.    whenever aggressive early deflation finds   
		.    at least (NIBBLE*(window size)/100) deflations. ==== */

		ret_val = 14;

	} else if (*ispec == 15) {

		/*        ==== NSHFTS: The number of simultaneous shifts ===== */

		ret_val = ns;

	} else if (*ispec == 13) {

		/*        ==== NW: deflation window size.  ==== */

		if (nh <= 500) {
			ret_val = ns;
		} else {
			ret_val = ns * 3 / 2;
		}

	} else if (*ispec == 16) {

		/*        ==== IACC22: Whether to accumulate reflections   
		.     before updating the far-from-diagonal elements   
		.     and whether to use 2-by-2 block structure while   
		.     doing it.  A small amount of work could be saved   
		.     by making this choice dependent also upon the   
		.     NH=IHI-ILO+1. */

		ret_val = 0;
		if (ns >= 14) {
			ret_val = 1;
		}
		if (ns >= 14) {
			ret_val = 2;
		}

	} else {
		/*        ===== invalid value of ispec ===== */
		ret_val = -1;

	}

	/*     ==== End of IPARMQ ==== */

	return ret_val;
} /* iparmq_ */

//
int TRL::dlae2_(doublereal_ *a, doublereal_ *b, doublereal_ *c__, 
		   doublereal_ *rt1, doublereal_ *rt2)
{
	/*  -- LAPACK auxiliary routine (version 3.1) --   
	Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..   
	November 2006   


	Purpose   
	=======   

	DLAE2  computes the eigenvalues of a 2-by-2 symmetric matrix   
	[  A   B  ]   
	[  B   C  ].   
	On return, RT1 is the eigenvalue of larger absolute value, and RT2   
	is the eigenvalue of smaller absolute value.   

	Arguments   
	=========   

	A       (input) DOUBLE PRECISION   
	The (1,1) element of the 2-by-2 matrix.   

	B       (input) DOUBLE PRECISION   
	The (1,2) and (2,1) elements of the 2-by-2 matrix.   

	C       (input) DOUBLE PRECISION   
	The (2,2) element of the 2-by-2 matrix.   

	RT1     (output) DOUBLE PRECISION   
	The eigenvalue of larger absolute value.   

	RT2     (output) DOUBLE PRECISION   
	The eigenvalue of smaller absolute value.   

	Further Details   
	===============   

	RT1 is accurate to a few ulps barring over/underflow.   

	RT2 may be inaccurate if there is massive cancellation in the   
	determinant A*C-B*B; higher precision or correctly rounded or   
	correctly truncated arithmetic would be needed to compute RT2   
	accurately in all cases.   

	Overflow is possible only if RT1 is within a factor of 5 of overflow.   
	Underflow is harmless if the input data is 0 or exceeds   
	underflow_threshold / macheps.   

	=====================================================================   


	Compute the eigenvalues */
	/* System generated locals */
	doublereal_ d__1;
	/* Builtin functions */
	//double sqrt_(doublereal_);
	/* Local variables */
	static doublereal_ ab, df, tb, sm, rt, adf, acmn, acmx;


	sm = *a + *c__;
	df = *a - *c__;
	adf = abs(df);
	tb = *b + *b;
	ab = abs(tb);
	if (abs(*a) > abs(*c__)) {
		acmx = *a;
		acmn = *c__;
	} else {
		acmx = *c__;
		acmn = *a;
	}
	if (adf > ab) {
		/* Computing 2nd power */
		d__1 = ab / adf;
		rt = adf * sqrt_(d__1 * d__1 + 1.);
	} else if (adf < ab) {
		/* Computing 2nd power */
		d__1 = adf / ab;
		rt = ab * sqrt_(d__1 * d__1 + 1.);
	} else {

		/*        Includes case AB=ADF=0 */

		rt = ab * sqrt_(2.);
	}
	if (sm < 0.) {
		*rt1 = (sm - rt) * .5;

		/*        Order of execution important.   
		To get fully accurate smaller eigenvalue,   
		next line needs to be executed in higher precision. */

		*rt2 = acmx / *rt1 * acmn - *b / *rt1 * *b;
	} else if (sm > 0.) {
		*rt1 = (sm + rt) * .5;

		/*        Order of execution important.   
		To get fully accurate smaller eigenvalue,   
		next line needs to be executed in higher precision. */

		*rt2 = acmx / *rt1 * acmn - *b / *rt1 * *b;
	} else {

		/*        Includes case RT1 = RT2 = 0 */

		*rt1 = rt * .5;
		*rt2 = rt * -.5;
	}
	return 0;

	/*     End of DLAE2 */

} /* dlae2_ */




//
doublereal_ TRL::dlapy2_(doublereal_ *x, doublereal_ *y)
{
	/*  -- LAPACK auxiliary routine (version 3.1) --   
	Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..   
	November 2006   


	Purpose   
	=======   

	DLAPY2 returns sqrt_(x**2+y**2), taking care not to cause unnecessary   
	overflow.   

	Arguments   
	=========   

	X       (input) DOUBLE PRECISION   
	Y       (input) DOUBLE PRECISION   
	X and Y specify the values x and y.   

	===================================================================== */
	/* System generated locals */
	doublereal_ ret_val, d__1;
	/* Builtin functions */
	//double sqrt_(doublereal_);
	/* Local variables */
	static doublereal_ w, z__, xabs, yabs;



	xabs = abs(*x);
	yabs = abs(*y);
	w = max(xabs,yabs);
	z__ = min(xabs,yabs);
	if (z__ == 0.) {
		ret_val = w;
	} else {
		/* Computing 2nd power */
		d__1 = z__ / w;
		ret_val = w * sqrt_(d__1 * d__1 + 1.);
	}
	return ret_val;

	/*     End of DLAPY2 */

} /* dlapy2_ */

//
doublereal_ TRL::dlanst_(char *norm, integer_ *n, doublereal_ *d__, doublereal_ *e)
{
	/*  -- LAPACK auxiliary routine (version 3.1) --   
	Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..   
	November 2006   


	Purpose   
	=======   

	DLANST  returns the value of the one norm,  or the Frobenius norm, or   
	the  infinity norm,  or the  element of  largest absolute value  of a   
	real symmetric tridiagonal matrix A.   

	Description   
	===========   

	DLANST returns the value   

	DLANST = ( max(abs(A(i,j))), NORM = 'M' or 'm'   
	(   
	( norm1(A),         NORM = '1', 'O' or 'o'   
	(   
	( normI(A),         NORM = 'I' or 'i'   
	(   
	( normF(A),         NORM = 'F', 'f', 'E' or 'e'   

	where  norm1  denotes the  one norm of a matrix (maximum column sum),   
	normI  denotes the  infinity norm  of a matrix  (maximum row sum) and   
	normF  denotes the  Frobenius norm of a matrix (square root of sum of   
	squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.   

	Arguments   
	=========   

	NORM    (input) CHARACTER*1   
	Specifies the value to be returned in DLANST as described   
	above.   

	N       (input) INTEGER   
	The order of the matrix A.  N >= 0.  When N = 0, DLANST is   
	set to zero.   

	D       (input) DOUBLE PRECISION array, dimension (N)   
	The diagonal elements of A.   

	E       (input) DOUBLE PRECISION array, dimension (N-1)   
	The (n-1) sub-diagonal or super-diagonal elements of A.   

	=====================================================================   


	Parameter adjustments */
	/* Table of constant values */
	static integer_ c__1 = 1;

	/* System generated locals */
	integer_ i__1;
	doublereal_ ret_val, d__1, d__2, d__3, d__4, d__5;
	/* Builtin functions */
	//double sqrt_(doublereal_);
	/* Local variables */
	static integer_ i__;
	static doublereal_ sum, scale;
	//extern logical_ lsame_(char *, char *);
	static doublereal_ anorm;
	//extern /* Subroutine */ int dlassq_(integer_ *, doublereal_ *, integer_ *, 
	//	doublereal_ *, doublereal_ *);


	--e;
	--d__;

	/* Function Body */
	if (*n <= 0) {
		anorm = 0.;
	} else if (lsame_(norm, "M")) {

		/*        Find max(abs(A(i,j))). */

		anorm = (d__1 = d__[*n], abs(d__1));
		i__1 = *n - 1;
		for (i__ = 1; i__ <= i__1; ++i__) {
			/* Computing MAX */
			d__2 = anorm, d__3 = (d__1 = d__[i__], abs(d__1));
			anorm = max(d__2,d__3);
			/* Computing MAX */
			d__2 = anorm, d__3 = (d__1 = e[i__], abs(d__1));
			anorm = max(d__2,d__3);
			/* L10: */
		}
	} else if (lsame_(norm, "O") || *(unsigned char *)
		norm == '1' || lsame_(norm, "I")) {

			/*        Find norm1(A). */

			if (*n == 1) {
				anorm = abs(d__[1]);
			} else {
				/* Computing MAX */
				d__3 = abs(d__[1]) + abs(e[1]), d__4 = (d__1 = e[*n - 1], abs(
					d__1)) + (d__2 = d__[*n], abs(d__2));
				anorm = max(d__3,d__4);
				i__1 = *n - 1;
				for (i__ = 2; i__ <= i__1; ++i__) {
					/* Computing MAX */
					d__4 = anorm, d__5 = (d__1 = d__[i__], abs(d__1)) + (d__2 = e[
						i__], abs(d__2)) + (d__3 = e[i__ - 1], abs(d__3));
						anorm = max(d__4,d__5);
						/* L20: */
				}
			}
	} else if (lsame_(norm, "F") || lsame_(norm, "E")) {

		/*        Find normF(A). */

		scale = 0.;
		sum = 1.;
		if (*n > 1) {
			i__1 = *n - 1;
			dlassq_(&i__1, &e[1], &c__1, &scale, &sum);
			sum *= 2;
		}
		dlassq_(n, &d__[1], &c__1, &scale, &sum);
		anorm = scale * sqrt_(sum);
	}

	ret_val = anorm;
	return ret_val;

	/*     End of DLANST */

} /* dlanst_ */

//
int TRL::dlasrt_(char *id, integer_ *n, doublereal_ *d__, integer_ *
			info)
{
	/*  -- LAPACK routine (version 3.1) --   
	Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..   
	November 2006   


	Purpose   
	=======   

	Sort the numbers in D in increasing order (if ID = 'I') or   
	in decreasing order (if ID = 'D' ).   

	Use Quick Sort, reverting to Insertion sort on arrays of   
	size <= 20. Dimension of STACK limits N to about 2**32.   

	Arguments   
	=========   

	ID      (input) CHARACTER*1   
	= 'I': sort D in increasing order;   
	= 'D': sort D in decreasing order.   

	N       (input) INTEGER   
	The length of the array D.   

	D       (input/output) DOUBLE PRECISION array, dimension (N)   
	On entry, the array to be sorted.   
	On exit, D has been sorted into increasing order   
	(D(1) <= ... <= D(N) ) or into decreasing order   
	(D(1) >= ... >= D(N) ), depending on ID.   

	INFO    (output) INTEGER   
	= 0:  successful exit   
	< 0:  if INFO = -i, the i-th argument had an illegal value   

	=====================================================================   


	Test the input paramters.   

	Parameter adjustments */
	/* System generated locals */
	integer_ i__1, i__2;
	/* Local variables */
	static integer_ i__, j;
	static doublereal_ d1, d2, d3;
	static integer_ dir;
	static doublereal_ tmp;
	static integer_ endd;
	//extern logical_ lsame_(char *, char *);
	static integer_ stack[64]	/* was [2][32] */;
	static doublereal_ dmnmx;
	static integer_ start;
	//extern /* Subroutine */ int xerbla_(char *, integer_ *);
	static integer_ stkpnt;

	--d__;

	/* Function Body */
	*info = 0;
	dir = -1;
	if (lsame_(id, "D")) {
		dir = 0;
	} else if (lsame_(id, "I")) {
		dir = 1;
	}
	if (dir == -1) {
		*info = -1;
	} else if (*n < 0) {
		*info = -2;
	}
	if (*info != 0) {
		i__1 = -(*info);
		xerbla_("DLASRT", &i__1);
		return 0;
	}

	/*     Quick return if possible */

	if (*n <= 1) {
		return 0;
	}

	stkpnt = 1;
	stack[0] = 1;
	stack[1] = *n;
L10:
	start = stack[(stkpnt << 1) - 2];
	endd = stack[(stkpnt << 1) - 1];
	--stkpnt;
	if (endd - start <= 20 && endd - start > 0) {

		/*        Do Insertion sort on D( START:ENDD ) */

		if (dir == 0) {

			/*           Sort into decreasing order */

			i__1 = endd;
			for (i__ = start + 1; i__ <= i__1; ++i__) {
				i__2 = start + 1;
				for (j = i__; j >= i__2; --j) {
					if (d__[j] > d__[j - 1]) {
						dmnmx = d__[j];
						d__[j] = d__[j - 1];
						d__[j - 1] = dmnmx;
					} else {
						goto L30;
					}
					/* L20: */
				}
L30:
				;
			}

		} else {

			/*           Sort into increasing order */

			i__1 = endd;
			for (i__ = start + 1; i__ <= i__1; ++i__) {
				i__2 = start + 1;
				for (j = i__; j >= i__2; --j) {
					if (d__[j] < d__[j - 1]) {
						dmnmx = d__[j];
						d__[j] = d__[j - 1];
						d__[j - 1] = dmnmx;
					} else {
						goto L50;
					}
					/* L40: */
				}
L50:
				;
			}

		}

	} else if (endd - start > 20) {

		/*        Partition D( START:ENDD ) and stack parts, largest one first   

		Choose partition entry as median of 3 */

		d1 = d__[start];
		d2 = d__[endd];
		i__ = (start + endd) / 2;
		d3 = d__[i__];
		if (d1 < d2) {
			if (d3 < d1) {
				dmnmx = d1;
			} else if (d3 < d2) {
				dmnmx = d3;
			} else {
				dmnmx = d2;
			}
		} else {
			if (d3 < d2) {
				dmnmx = d2;
			} else if (d3 < d1) {
				dmnmx = d3;
			} else {
				dmnmx = d1;
			}
		}

		if (dir == 0) {

			/*           Sort into decreasing order */

			i__ = start - 1;
			j = endd + 1;
L60:
L70:
			--j;
			if (d__[j] < dmnmx) {
				goto L70;
			}
L80:
			++i__;
			if (d__[i__] > dmnmx) {
				goto L80;
			}
			if (i__ < j) {
				tmp = d__[i__];
				d__[i__] = d__[j];
				d__[j] = tmp;
				goto L60;
			}
			if (j - start > endd - j - 1) {
				++stkpnt;
				stack[(stkpnt << 1) - 2] = start;
				stack[(stkpnt << 1) - 1] = j;
				++stkpnt;
				stack[(stkpnt << 1) - 2] = j + 1;
				stack[(stkpnt << 1) - 1] = endd;
			} else {
				++stkpnt;
				stack[(stkpnt << 1) - 2] = j + 1;
				stack[(stkpnt << 1) - 1] = endd;
				++stkpnt;
				stack[(stkpnt << 1) - 2] = start;
				stack[(stkpnt << 1) - 1] = j;
			}
		} else {

			/*           Sort into increasing order */

			i__ = start - 1;
			j = endd + 1;
L90:
L100:
			--j;
			if (d__[j] > dmnmx) {
				goto L100;
			}
L110:
			++i__;
			if (d__[i__] < dmnmx) {
				goto L110;
			}
			if (i__ < j) {
				tmp = d__[i__];
				d__[i__] = d__[j];
				d__[j] = tmp;
				goto L90;
			}
			if (j - start > endd - j - 1) {
				++stkpnt;
				stack[(stkpnt << 1) - 2] = start;
				stack[(stkpnt << 1) - 1] = j;
				++stkpnt;
				stack[(stkpnt << 1) - 2] = j + 1;
				stack[(stkpnt << 1) - 1] = endd;
			} else {
				++stkpnt;
				stack[(stkpnt << 1) - 2] = j + 1;
				stack[(stkpnt << 1) - 1] = endd;
				++stkpnt;
				stack[(stkpnt << 1) - 2] = start;
				stack[(stkpnt << 1) - 1] = j;
			}
		}
	}
	if (stkpnt > 0) {
		goto L10;
	}
	return 0;

	/*     End of DLASRT */

} /* dlasrt_ */

//
int TRL::dlassq_(integer_ *n, doublereal_ *x, integer_ *incx, 
			doublereal_ *scale, doublereal_ *sumsq)
{
	/*  -- LAPACK auxiliary routine (version 3.1) --   
	Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..   
	November 2006   


	Purpose   
	=======   

	DLASSQ  returns the values  scl  and  smsq  such that   

	( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,   

	where  x( i ) = X( 1 + ( i - 1 )*INCX ). The value of  sumsq  is   
	assumed to be non-negative and  scl  returns the value   

	scl = max( scale, abs( x( i ) ) ).   

	scale and sumsq must be supplied in SCALE and SUMSQ and   
	scl and smsq are overwritten on SCALE and SUMSQ respectively.   

	The routine makes only one pass through the vector x.   

	Arguments   
	=========   

	N       (input) INTEGER   
	The number of elements to be used from the vector X.   

	X       (input) DOUBLE PRECISION array, dimension (N)   
	The vector for which a scaled sum of squares is computed.   
	x( i )  = X( 1 + ( i - 1 )*INCX ), 1 <= i <= n.   

	INCX    (input) INTEGER   
	The increment between successive values of the vector X.   
	INCX > 0.   

	SCALE   (input/output) DOUBLE PRECISION   
	On entry, the value  scale  in the equation above.   
	On exit, SCALE is overwritten with  scl , the scaling factor   
	for the sum of squares.   

	SUMSQ   (input/output) DOUBLE PRECISION   
	On entry, the value  sumsq  in the equation above.   
	On exit, SUMSQ is overwritten with  smsq , the basic sum of   
	squares from which  scl  has been factored out.   

	=====================================================================   


	Parameter adjustments */
	/* System generated locals */
	integer_ i__1, i__2;
	doublereal_ d__1;
	/* Local variables */
	static integer_ ix;
	static doublereal_ absxi;

	--x;

	/* Function Body */
	if (*n > 0) {
		i__1 = (*n - 1) * *incx + 1;
		i__2 = *incx;
		for (ix = 1; i__2 < 0 ? ix >= i__1 : ix <= i__1; ix += i__2) {
			if (x[ix] != 0.) {
				absxi = (d__1 = x[ix], abs(d__1));
				if (*scale < absxi) {
					/* Computing 2nd power */
					d__1 = *scale / absxi;
					*sumsq = *sumsq * (d__1 * d__1) + 1;
					*scale = absxi;
				} else {
					/* Computing 2nd power */
					d__1 = absxi / *scale;
					*sumsq += d__1 * d__1;
				}
			}
			/* L10: */
		}
	}
	return 0;

	/*     End of DLASSQ */

} /* dlassq_ */

//
int TRL::dorg2l_(integer_ *m, integer_ *n, integer_ *k, doublereal_ *
			a, integer_ *lda, doublereal_ *tau, doublereal_ *work, integer_ *info)
{
	/*  -- LAPACK routine (version 3.1) --   
	Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..   
	November 2006   


	Purpose   
	=======   

	DORG2L generates an m by n trl_real_ matrix Q with orthonormal columns,   
	which is defined as the last n columns of a product of k elementary   
	reflectors of order m   

	Q  =  H(k) . . . H(2) H(1)   

	as returned by DGEQLF.   

	Arguments   
	=========   

	M       (input) INTEGER   
	The number of rows of the matrix Q. M >= 0.   

	N       (input) INTEGER   
	The number of columns of the matrix Q. M >= N >= 0.   

	K       (input) INTEGER   
	The number of elementary reflectors whose product defines the   
	matrix Q. N >= K >= 0.   

	A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
	On entry, the (n-k+i)-th column must contain the vector which   
	defines the elementary reflector H(i), for i = 1,2,...,k, as   
	returned by DGEQLF in the last k columns of its array   
	argument A.   
	On exit, the m by n matrix Q.   

	LDA     (input) INTEGER   
	The first dimension of the array A. LDA >= max(1,M).   

	TAU     (input) DOUBLE PRECISION array, dimension (K)   
	TAU(i) must contain the scalar factor of the elementary   
	reflector H(i), as returned by DGEQLF.   

	WORK    (workspace) DOUBLE PRECISION array, dimension (N)   

	INFO    (output) INTEGER   
	= 0: successful exit   
	< 0: if INFO = -i, the i-th argument has an illegal value   

	=====================================================================   


	Test the input arguments   

	Parameter adjustments */
	/* Table of constant values */
	static integer_ c__1 = 1;

	/* System generated locals */
	integer_ a_dim1, a_offset, i__1, i__2, i__3;
	doublereal_ d__1;
	/* Local variables */
	static integer_ i__, j, l, ii;
	//extern /* Subroutine */ int dscal_(integer_ *, doublereal_ *, doublereal_ *, 
	//	integer_ *), dlarf_(char *, integer_ *, integer_ *, doublereal_ *, 
	//	integer_ *, doublereal_ *, doublereal_ *, integer_ *, doublereal_ *), xerbla_(char *, integer_ *);


	a_dim1 = *lda;
	a_offset = 1 + a_dim1;
	a -= a_offset;
	--tau;
	--work;

	/* Function Body */
	*info = 0;
	if (*m < 0) {
		*info = -1;
	} else if (*n < 0 || *n > *m) {
		*info = -2;
	} else if (*k < 0 || *k > *n) {
		*info = -3;
	} else if (*lda < max(1,*m)) {
		*info = -5;
	}
	if (*info != 0) {
		i__1 = -(*info);
		xerbla_("DORG2L", &i__1);
		return 0;
	}

	/*     Quick return if possible */

	if (*n <= 0) {
		return 0;
	}

	/*     Initialise columns 1:n-k to columns of the unit matrix */

	i__1 = *n - *k;
	for (j = 1; j <= i__1; ++j) {
		i__2 = *m;
		for (l = 1; l <= i__2; ++l) {
			a[l + j * a_dim1] = 0.;
			/* L10: */
		}
		a[*m - *n + j + j * a_dim1] = 1.;
		/* L20: */
	}

	i__1 = *k;
	for (i__ = 1; i__ <= i__1; ++i__) {
		ii = *n - *k + i__;

		/*        Apply H(i) to A(1:m-k+i,1:n-k+i) from the left */

		a[*m - *n + ii + ii * a_dim1] = 1.;
		i__2 = *m - *n + ii;
		i__3 = ii - 1;
		dlarf_("Left", &i__2, &i__3, &a[ii * a_dim1 + 1], &c__1, &tau[i__], &
			a[a_offset], lda, &work[1]);
		i__2 = *m - *n + ii - 1;
		d__1 = -tau[i__];
		dscal_(&i__2, &d__1, &a[ii * a_dim1 + 1], &c__1);
		a[*m - *n + ii + ii * a_dim1] = 1. - tau[i__];

		/*        Set A(m-k+i+1:m,n-k+i) to zero */

		i__2 = *m;
		for (l = *m - *n + ii + 1; l <= i__2; ++l) {
			a[l + ii * a_dim1] = 0.;
			/* L30: */
		}
		/* L40: */
	}
	return 0;

	/*     End of DORG2L */

} /* dorg2l_ */


//
int TRL::dlarfb_(char *side, char *trans, char *direct, char *
			storev, integer_ *m, integer_ *n, integer_ *k, doublereal_ *v, integer_ *
			ldv, doublereal_ *t, integer_ *ldt, doublereal_ *c__, integer_ *ldc, 
			doublereal_ *work, integer_ *ldwork)
{
	/*  -- LAPACK auxiliary routine (version 3.1) --   
	Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..   
	November 2006   


	Purpose   
	=======   

	DLARFB applies a trl_real_ block reflector H or its transpose H' to a   
	real m by n matrix C, from either the left or the right.   

	Arguments   
	=========   

	SIDE    (input) CHARACTER*1   
	= 'L': apply H or H' from the Left   
	= 'R': apply H or H' from the Right   

	TRANS   (input) CHARACTER*1   
	= 'N': apply H (No transpose)   
	= 'T': apply H' (Transpose)   

	DIRECT  (input) CHARACTER*1   
	Indicates how H is formed from a product of elementary   
	reflectors   
	= 'F': H = H(1) H(2) . . . H(k) (Forward)   
	= 'B': H = H(k) . . . H(2) H(1) (Backward)   

	STOREV  (input) CHARACTER*1   
	Indicates how the vectors which define the elementary   
	reflectors are stored:   
	= 'C': Columnwise   
	= 'R': Rowwise   

	M       (input) INTEGER   
	The number of rows of the matrix C.   

	N       (input) INTEGER   
	The number of columns of the matrix C.   

	K       (input) INTEGER   
	The order of the matrix T (= the number of elementary   
	reflectors whose product defines the block reflector).   

	V       (input) DOUBLE PRECISION array, dimension   
	(LDV,K) if STOREV = 'C'   
	(LDV,M) if STOREV = 'R' and SIDE = 'L'   
	(LDV,N) if STOREV = 'R' and SIDE = 'R'   
	The matrix V. See further details.   

	LDV     (input) INTEGER   
	The leading dimension of the array V.   
	If STOREV = 'C' and SIDE = 'L', LDV >= max(1,M);   
	if STOREV = 'C' and SIDE = 'R', LDV >= max(1,N);   
	if STOREV = 'R', LDV >= K.   

	T       (input) DOUBLE PRECISION array, dimension (LDT,K)   
	The triangular k by k matrix T in the representation of the   
	block reflector.   

	LDT     (input) INTEGER   
	The leading dimension of the array T. LDT >= K.   

	C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)   
	On entry, the m by n matrix C.   
	On exit, C is overwritten by H*C or H'*C or C*H or C*H'.   

	LDC     (input) INTEGER   
	The leading dimension of the array C. LDA >= max(1,M).   

	WORK    (workspace) DOUBLE PRECISION array, dimension (LDWORK,K)   

	LDWORK  (input) INTEGER   
	The leading dimension of the array WORK.   
	If SIDE = 'L', LDWORK >= max(1,N);   
	if SIDE = 'R', LDWORK >= max(1,M).   

	=====================================================================   


	Quick return if possible   

	Parameter adjustments */
	/* Table of constant values */
	static integer_ c__1 = 1;
	static doublereal_ c_b14 = 1.;
	static doublereal_ c_b25 = -1.;

	/* System generated locals */
	integer_ c_dim1, c_offset, t_dim1, t_offset, v_dim1, v_offset, work_dim1, 
		work_offset, i__1, i__2;
	/* Local variables */
	static integer_ i__, j;
	//extern /* Subroutine */ int dgemm_(char *, char *, integer_ *, integer_ *, 
	//	integer_ *, doublereal_ *, doublereal_ *, integer_ *, doublereal_ *, 
	//	integer_ *, doublereal_ *, doublereal_ *, integer_ *);
	//extern logical_ lsame_(char *, char *);
	//extern /* Subroutine */ int dcopy_(integer_ *, doublereal_ *, integer_ *, 
	//	doublereal_ *, integer_ *), dtrmm_(char *, char *, char *, char *, 
	//	integer_ *, integer_ *, doublereal_ *, doublereal_ *, integer_ *, 
	//	doublereal_ *, integer_ *);
	static char transt[1];


	v_dim1 = *ldv;
	v_offset = 1 + v_dim1;
	v -= v_offset;
	t_dim1 = *ldt;
	t_offset = 1 + t_dim1;
	t -= t_offset;
	c_dim1 = *ldc;
	c_offset = 1 + c_dim1;
	c__ -= c_offset;
	work_dim1 = *ldwork;
	work_offset = 1 + work_dim1;
	work -= work_offset;

	/* Function Body */
	if (*m <= 0 || *n <= 0) {
		return 0;
	}

	if (lsame_(trans, "N")) {
		*(unsigned char *)transt = 'T';
	} else {
		*(unsigned char *)transt = 'N';
	}

	if (lsame_(storev, "C")) {

		if (lsame_(direct, "F")) {

			/*           Let  V =  ( V1 )    (first K rows)   
			( V2 )   
			where  V1  is unit lower triangular. */

			if (lsame_(side, "L")) {

				/*              Form  H * C  or  H' * C  where  C = ( C1 )   
				( C2 )   

				W := C' * V  =  (C1'*V1 + C2'*V2)  (stored in WORK)   

				W := C1' */

				i__1 = *k;
				for (j = 1; j <= i__1; ++j) {
					dcopy_(n, &c__[j + c_dim1], ldc, &work[j * work_dim1 + 1], 
						&c__1);
					/* L10: */
				}

				/*              W := W * V1 */

				dtrmm_("Right", "Lower", "No transpose", "Unit", n, k, &c_b14, 
					&v[v_offset], ldv, &work[work_offset], ldwork);
				if (*m > *k) {

					/*                 W := W + C2'*V2 */

					i__1 = *m - *k;
					dgemm_("Transpose", "No transpose", n, k, &i__1, &c_b14, &
						c__[*k + 1 + c_dim1], ldc, &v[*k + 1 + v_dim1], 
						ldv, &c_b14, &work[work_offset], ldwork);
				}

				/*              W := W * T'  or  W * T */

				dtrmm_("Right", "Upper", transt, "Non-unit", n, k, &c_b14, &t[
					t_offset], ldt, &work[work_offset], ldwork);

					/*              C := C - V * W' */

					if (*m > *k) {

						/*                 C2 := C2 - V2 * W' */

						i__1 = *m - *k;
						dgemm_("No transpose", "Transpose", &i__1, n, k, &c_b25, &
							v[*k + 1 + v_dim1], ldv, &work[work_offset], 
							ldwork, &c_b14, &c__[*k + 1 + c_dim1], ldc);
					}

					/*              W := W * V1' */

					dtrmm_("Right", "Lower", "Transpose", "Unit", n, k, &c_b14, &
						v[v_offset], ldv, &work[work_offset], ldwork);

					/*              C1 := C1 - W' */

					i__1 = *k;
					for (j = 1; j <= i__1; ++j) {
						i__2 = *n;
						for (i__ = 1; i__ <= i__2; ++i__) {
							c__[j + i__ * c_dim1] -= work[i__ + j * work_dim1];
							/* L20: */
						}
						/* L30: */
					}

			} else if (lsame_(side, "R")) {

				/*              Form  C * H  or  C * H'  where  C = ( C1  C2 )   

				W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)   

				W := C1 */

				i__1 = *k;
				for (j = 1; j <= i__1; ++j) {
					dcopy_(m, &c__[j * c_dim1 + 1], &c__1, &work[j * 
						work_dim1 + 1], &c__1);
					/* L40: */
				}

				/*              W := W * V1 */

				dtrmm_("Right", "Lower", "No transpose", "Unit", m, k, &c_b14, 
					&v[v_offset], ldv, &work[work_offset], ldwork);
				if (*n > *k) {

					/*                 W := W + C2 * V2 */

					i__1 = *n - *k;
					dgemm_("No transpose", "No transpose", m, k, &i__1, &
						c_b14, &c__[(*k + 1) * c_dim1 + 1], ldc, &v[*k + 
						1 + v_dim1], ldv, &c_b14, &work[work_offset], 
						ldwork);
				}

				/*              W := W * T  or  W * T' */

				dtrmm_("Right", "Upper", trans, "Non-unit", m, k, &c_b14, &t[
					t_offset], ldt, &work[work_offset], ldwork);

					/*              C := C - W * V' */

					if (*n > *k) {

						/*                 C2 := C2 - W * V2' */

						i__1 = *n - *k;
						dgemm_("No transpose", "Transpose", m, &i__1, k, &c_b25, &
							work[work_offset], ldwork, &v[*k + 1 + v_dim1], 
							ldv, &c_b14, &c__[(*k + 1) * c_dim1 + 1], ldc);
					}

					/*              W := W * V1' */

					dtrmm_("Right", "Lower", "Transpose", "Unit", m, k, &c_b14, &
						v[v_offset], ldv, &work[work_offset], ldwork);

					/*              C1 := C1 - W */

					i__1 = *k;
					for (j = 1; j <= i__1; ++j) {
						i__2 = *m;
						for (i__ = 1; i__ <= i__2; ++i__) {
							c__[i__ + j * c_dim1] -= work[i__ + j * work_dim1];
							/* L50: */
						}
						/* L60: */
					}
			}

		} else {

			/*           Let  V =  ( V1 )   
			( V2 )    (last K rows)   
			where  V2  is unit upper triangular. */

			if (lsame_(side, "L")) {

				/*              Form  H * C  or  H' * C  where  C = ( C1 )   
				( C2 )   

				W := C' * V  =  (C1'*V1 + C2'*V2)  (stored in WORK)   

				W := C2' */

				i__1 = *k;
				for (j = 1; j <= i__1; ++j) {
					dcopy_(n, &c__[*m - *k + j + c_dim1], ldc, &work[j * 
						work_dim1 + 1], &c__1);
					/* L70: */
				}

				/*              W := W * V2 */

				dtrmm_("Right", "Upper", "No transpose", "Unit", n, k, &c_b14, 
					&v[*m - *k + 1 + v_dim1], ldv, &work[work_offset], 
					ldwork);
				if (*m > *k) {

					/*                 W := W + C1'*V1 */

					i__1 = *m - *k;
					dgemm_("Transpose", "No transpose", n, k, &i__1, &c_b14, &
						c__[c_offset], ldc, &v[v_offset], ldv, &c_b14, &
						work[work_offset], ldwork);
				}

				/*              W := W * T'  or  W * T */

				dtrmm_("Right", "Lower", transt, "Non-unit", n, k, &c_b14, &t[
					t_offset], ldt, &work[work_offset], ldwork);

					/*              C := C - V * W' */

					if (*m > *k) {

						/*                 C1 := C1 - V1 * W' */

						i__1 = *m - *k;
						dgemm_("No transpose", "Transpose", &i__1, n, k, &c_b25, &
							v[v_offset], ldv, &work[work_offset], ldwork, &
							c_b14, &c__[c_offset], ldc)
							;
					}

					/*              W := W * V2' */

					dtrmm_("Right", "Upper", "Transpose", "Unit", n, k, &c_b14, &
						v[*m - *k + 1 + v_dim1], ldv, &work[work_offset], 
						ldwork);

					/*              C2 := C2 - W' */

					i__1 = *k;
					for (j = 1; j <= i__1; ++j) {
						i__2 = *n;
						for (i__ = 1; i__ <= i__2; ++i__) {
							c__[*m - *k + j + i__ * c_dim1] -= work[i__ + j * 
								work_dim1];
							/* L80: */
						}
						/* L90: */
					}

			} else if (lsame_(side, "R")) {

				/*              Form  C * H  or  C * H'  where  C = ( C1  C2 )   

				W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)   

				W := C2 */

				i__1 = *k;
				for (j = 1; j <= i__1; ++j) {
					dcopy_(m, &c__[(*n - *k + j) * c_dim1 + 1], &c__1, &work[
						j * work_dim1 + 1], &c__1);
						/* L100: */
				}

				/*              W := W * V2 */

				dtrmm_("Right", "Upper", "No transpose", "Unit", m, k, &c_b14, 
					&v[*n - *k + 1 + v_dim1], ldv, &work[work_offset], 
					ldwork);
				if (*n > *k) {

					/*                 W := W + C1 * V1 */

					i__1 = *n - *k;
					dgemm_("No transpose", "No transpose", m, k, &i__1, &
						c_b14, &c__[c_offset], ldc, &v[v_offset], ldv, &
						c_b14, &work[work_offset], ldwork);
				}

				/*              W := W * T  or  W * T' */

				dtrmm_("Right", "Lower", trans, "Non-unit", m, k, &c_b14, &t[
					t_offset], ldt, &work[work_offset], ldwork);

					/*              C := C - W * V' */

					if (*n > *k) {

						/*                 C1 := C1 - W * V1' */

						i__1 = *n - *k;
						dgemm_("No transpose", "Transpose", m, &i__1, k, &c_b25, &
							work[work_offset], ldwork, &v[v_offset], ldv, &
							c_b14, &c__[c_offset], ldc)
							;
					}

					/*              W := W * V2' */

					dtrmm_("Right", "Upper", "Transpose", "Unit", m, k, &c_b14, &
						v[*n - *k + 1 + v_dim1], ldv, &work[work_offset], 
						ldwork);

					/*              C2 := C2 - W */

					i__1 = *k;
					for (j = 1; j <= i__1; ++j) {
						i__2 = *m;
						for (i__ = 1; i__ <= i__2; ++i__) {
							c__[i__ + (*n - *k + j) * c_dim1] -= work[i__ + j * 
								work_dim1];
							/* L110: */
						}
						/* L120: */
					}
			}
		}

	} else if (lsame_(storev, "R")) {

		if (lsame_(direct, "F")) {

			/*           Let  V =  ( V1  V2 )    (V1: first K columns)   
			where  V1  is unit upper triangular. */

			if (lsame_(side, "L")) {

				/*              Form  H * C  or  H' * C  where  C = ( C1 )   
				( C2 )   

				W := C' * V'  =  (C1'*V1' + C2'*V2') (stored in WORK)   

				W := C1' */

				i__1 = *k;
				for (j = 1; j <= i__1; ++j) {
					dcopy_(n, &c__[j + c_dim1], ldc, &work[j * work_dim1 + 1], 
						&c__1);
					/* L130: */
				}

				/*              W := W * V1' */

				dtrmm_("Right", "Upper", "Transpose", "Unit", n, k, &c_b14, &
					v[v_offset], ldv, &work[work_offset], ldwork);
				if (*m > *k) {

					/*                 W := W + C2'*V2' */

					i__1 = *m - *k;
					dgemm_("Transpose", "Transpose", n, k, &i__1, &c_b14, &
						c__[*k + 1 + c_dim1], ldc, &v[(*k + 1) * v_dim1 + 
						1], ldv, &c_b14, &work[work_offset], ldwork);
				}

				/*              W := W * T'  or  W * T */

				dtrmm_("Right", "Upper", transt, "Non-unit", n, k, &c_b14, &t[
					t_offset], ldt, &work[work_offset], ldwork);

					/*              C := C - V' * W' */

					if (*m > *k) {

						/*                 C2 := C2 - V2' * W' */

						i__1 = *m - *k;
						dgemm_("Transpose", "Transpose", &i__1, n, k, &c_b25, &v[(
							*k + 1) * v_dim1 + 1], ldv, &work[work_offset], 
							ldwork, &c_b14, &c__[*k + 1 + c_dim1], ldc);
					}

					/*              W := W * V1 */

					dtrmm_("Right", "Upper", "No transpose", "Unit", n, k, &c_b14, 
						&v[v_offset], ldv, &work[work_offset], ldwork);

					/*              C1 := C1 - W' */

					i__1 = *k;
					for (j = 1; j <= i__1; ++j) {
						i__2 = *n;
						for (i__ = 1; i__ <= i__2; ++i__) {
							c__[j + i__ * c_dim1] -= work[i__ + j * work_dim1];
							/* L140: */
						}
						/* L150: */
					}

			} else if (lsame_(side, "R")) {

				/*              Form  C * H  or  C * H'  where  C = ( C1  C2 )   

				W := C * V'  =  (C1*V1' + C2*V2')  (stored in WORK)   

				W := C1 */

				i__1 = *k;
				for (j = 1; j <= i__1; ++j) {
					dcopy_(m, &c__[j * c_dim1 + 1], &c__1, &work[j * 
						work_dim1 + 1], &c__1);
					/* L160: */
				}

				/*              W := W * V1' */

				dtrmm_("Right", "Upper", "Transpose", "Unit", m, k, &c_b14, &
					v[v_offset], ldv, &work[work_offset], ldwork);
				if (*n > *k) {

					/*                 W := W + C2 * V2' */

					i__1 = *n - *k;
					dgemm_("No transpose", "Transpose", m, k, &i__1, &c_b14, &
						c__[(*k + 1) * c_dim1 + 1], ldc, &v[(*k + 1) * 
						v_dim1 + 1], ldv, &c_b14, &work[work_offset], 
						ldwork);
				}

				/*              W := W * T  or  W * T' */

				dtrmm_("Right", "Upper", trans, "Non-unit", m, k, &c_b14, &t[
					t_offset], ldt, &work[work_offset], ldwork);

					/*              C := C - W * V */

					if (*n > *k) {

						/*                 C2 := C2 - W * V2 */

						i__1 = *n - *k;
						dgemm_("No transpose", "No transpose", m, &i__1, k, &
							c_b25, &work[work_offset], ldwork, &v[(*k + 1) * 
							v_dim1 + 1], ldv, &c_b14, &c__[(*k + 1) * c_dim1 
							+ 1], ldc);
					}

					/*              W := W * V1 */

					dtrmm_("Right", "Upper", "No transpose", "Unit", m, k, &c_b14, 
						&v[v_offset], ldv, &work[work_offset], ldwork);

					/*              C1 := C1 - W */

					i__1 = *k;
					for (j = 1; j <= i__1; ++j) {
						i__2 = *m;
						for (i__ = 1; i__ <= i__2; ++i__) {
							c__[i__ + j * c_dim1] -= work[i__ + j * work_dim1];
							/* L170: */
						}
						/* L180: */
					}

			}

		} else {

			/*           Let  V =  ( V1  V2 )    (V2: last K columns)   
			where  V2  is unit lower triangular. */

			if (lsame_(side, "L")) {

				/*              Form  H * C  or  H' * C  where  C = ( C1 )   
				( C2 )   

				W := C' * V'  =  (C1'*V1' + C2'*V2') (stored in WORK)   

				W := C2' */

				i__1 = *k;
				for (j = 1; j <= i__1; ++j) {
					dcopy_(n, &c__[*m - *k + j + c_dim1], ldc, &work[j * 
						work_dim1 + 1], &c__1);
					/* L190: */
				}

				/*              W := W * V2' */

				dtrmm_("Right", "Lower", "Transpose", "Unit", n, k, &c_b14, &
					v[(*m - *k + 1) * v_dim1 + 1], ldv, &work[work_offset]
				, ldwork);
				if (*m > *k) {

					/*                 W := W + C1'*V1' */

					i__1 = *m - *k;
					dgemm_("Transpose", "Transpose", n, k, &i__1, &c_b14, &
						c__[c_offset], ldc, &v[v_offset], ldv, &c_b14, &
						work[work_offset], ldwork);
				}

				/*              W := W * T'  or  W * T */

				dtrmm_("Right", "Lower", transt, "Non-unit", n, k, &c_b14, &t[
					t_offset], ldt, &work[work_offset], ldwork);

					/*              C := C - V' * W' */

					if (*m > *k) {

						/*                 C1 := C1 - V1' * W' */

						i__1 = *m - *k;
						dgemm_("Transpose", "Transpose", &i__1, n, k, &c_b25, &v[
							v_offset], ldv, &work[work_offset], ldwork, &
								c_b14, &c__[c_offset], ldc);
					}

					/*              W := W * V2 */

					dtrmm_("Right", "Lower", "No transpose", "Unit", n, k, &c_b14, 
						&v[(*m - *k + 1) * v_dim1 + 1], ldv, &work[
							work_offset], ldwork);

							/*              C2 := C2 - W' */

							i__1 = *k;
							for (j = 1; j <= i__1; ++j) {
								i__2 = *n;
								for (i__ = 1; i__ <= i__2; ++i__) {
									c__[*m - *k + j + i__ * c_dim1] -= work[i__ + j * 
										work_dim1];
									/* L200: */
								}
								/* L210: */
							}

			} else if (lsame_(side, "R")) {

				/*              Form  C * H  or  C * H'  where  C = ( C1  C2 )   

				W := C * V'  =  (C1*V1' + C2*V2')  (stored in WORK)   

				W := C2 */

				i__1 = *k;
				for (j = 1; j <= i__1; ++j) {
					dcopy_(m, &c__[(*n - *k + j) * c_dim1 + 1], &c__1, &work[
						j * work_dim1 + 1], &c__1);
						/* L220: */
				}

				/*              W := W * V2' */

				dtrmm_("Right", "Lower", "Transpose", "Unit", m, k, &c_b14, &
					v[(*n - *k + 1) * v_dim1 + 1], ldv, &work[work_offset]
				, ldwork);
				if (*n > *k) {

					/*                 W := W + C1 * V1' */

					i__1 = *n - *k;
					dgemm_("No transpose", "Transpose", m, k, &i__1, &c_b14, &
						c__[c_offset], ldc, &v[v_offset], ldv, &c_b14, &
						work[work_offset], ldwork);
				}

				/*              W := W * T  or  W * T' */

				dtrmm_("Right", "Lower", trans, "Non-unit", m, k, &c_b14, &t[
					t_offset], ldt, &work[work_offset], ldwork);

					/*              C := C - W * V */

					if (*n > *k) {

						/*                 C1 := C1 - W * V1 */

						i__1 = *n - *k;
						dgemm_("No transpose", "No transpose", m, &i__1, k, &
							c_b25, &work[work_offset], ldwork, &v[v_offset], 
							ldv, &c_b14, &c__[c_offset], ldc);
					}

					/*              W := W * V2 */

					dtrmm_("Right", "Lower", "No transpose", "Unit", m, k, &c_b14, 
						&v[(*n - *k + 1) * v_dim1 + 1], ldv, &work[
							work_offset], ldwork);

							/*              C1 := C1 - W */

							i__1 = *k;
							for (j = 1; j <= i__1; ++j) {
								i__2 = *m;
								for (i__ = 1; i__ <= i__2; ++i__) {
									c__[i__ + (*n - *k + j) * c_dim1] -= work[i__ + j * 
										work_dim1];
									/* L230: */
								}
								/* L240: */
							}

			}

		}
	}

	return 0;

	/*     End of DLARFB */

} /* dlarfb_ */

int TRL::dlarft_(char *direct, char *storev, integer_ *n, integer_ *
			k, doublereal_ *v, integer_ *ldv, doublereal_ *tau, doublereal_ *t, 
			integer_ *ldt)
{
	/*  -- LAPACK auxiliary routine (version 3.1) --   
	Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..   
	November 2006   


	Purpose   
	=======   

	DLARFT forms the triangular factor T of a trl_real_ block reflector H   
	of order n, which is defined as a product of k elementary reflectors.   

	If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;   

	If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.   

	If STOREV = 'C', the vector which defines the elementary reflector   
	H(i) is stored in the i-th column of the array V, and   

	H  =  I - V * T * V'   

	If STOREV = 'R', the vector which defines the elementary reflector   
	H(i) is stored in the i-th row of the array V, and   

	H  =  I - V' * T * V   

	Arguments   
	=========   

	DIRECT  (input) CHARACTER*1   
	Specifies the order in which the elementary reflectors are   
	multiplied to form the block reflector:   
	= 'F': H = H(1) H(2) . . . H(k) (Forward)   
	= 'B': H = H(k) . . . H(2) H(1) (Backward)   

	STOREV  (input) CHARACTER*1   
	Specifies how the vectors which define the elementary   
	reflectors are stored (see also Further Details):   
	= 'C': columnwise   
	= 'R': rowwise   

	N       (input) INTEGER   
	The order of the block reflector H. N >= 0.   

	K       (input) INTEGER   
	The order of the triangular factor T (= the number of   
	elementary reflectors). K >= 1.   

	V       (input/output) DOUBLE PRECISION array, dimension   
	(LDV,K) if STOREV = 'C'   
	(LDV,N) if STOREV = 'R'   
	The matrix V. See further details.   

	LDV     (input) INTEGER   
	The leading dimension of the array V.   
	If STOREV = 'C', LDV >= max(1,N); if STOREV = 'R', LDV >= K.   

	TAU     (input) DOUBLE PRECISION array, dimension (K)   
	TAU(i) must contain the scalar factor of the elementary   
	reflector H(i).   

	T       (output) DOUBLE PRECISION array, dimension (LDT,K)   
	The k by k triangular factor T of the block reflector.   
	If DIRECT = 'F', T is upper triangular; if DIRECT = 'B', T is   
	lower triangular. The rest of the array is not used.   

	LDT     (input) INTEGER   
	The leading dimension of the array T. LDT >= K.   

	Further Details   
	===============   

	The shape of the matrix V and the storage of the vectors which define   
	the H(i) is best illustrated by the following example with n = 5 and   
	k = 3. The elements equal to 1 are not stored; the corresponding   
	array elements are modified but restored on exit. The rest of the   
	array is not used.   

	DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R':   

	V = (  1       )                 V = (  1 v1 v1 v1 v1 )   
	( v1  1    )                     (     1 v2 v2 v2 )   
	( v1 v2  1 )                     (        1 v3 v3 )   
	( v1 v2 v3 )   
	( v1 v2 v3 )   

	DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R':   

	V = ( v1 v2 v3 )                 V = ( v1 v1  1       )   
	( v1 v2 v3 )                     ( v2 v2 v2  1    )   
	(  1 v2 v3 )                     ( v3 v3 v3 v3  1 )   
	(     1 v3 )   
	(        1 )   

	=====================================================================   


	Quick return if possible   

	Parameter adjustments */
	/* Table of constant values */
	static integer_ c__1 = 1;
	static doublereal_ c_b8 = 0.;

	/* System generated locals */
	integer_ t_dim1, t_offset, v_dim1, v_offset, i__1, i__2, i__3;
	doublereal_ d__1;
	/* Local variables */
	static integer_ i__, j;
	static doublereal_ vii;
	//extern logical_ lsame_(char *, char *);
	//extern /* Subroutine */ int dgemv_(char *, integer_ *, integer_ *, 
	//	doublereal_ *, doublereal_ *, integer_ *, doublereal_ *, integer_ *, 
	//	doublereal_ *, doublereal_ *, integer_ *), dtrmv_(char *, 
	//	char *, char *, integer_ *, doublereal_ *, integer_ *, doublereal_ *, 
	//	integer_ *);


	v_dim1 = *ldv;
	v_offset = 1 + v_dim1;
	v -= v_offset;
	--tau;
	t_dim1 = *ldt;
	t_offset = 1 + t_dim1;
	t -= t_offset;

	/* Function Body */
	if (*n == 0) {
		return 0;
	}

	if (lsame_(direct, "F")) {
		i__1 = *k;
		for (i__ = 1; i__ <= i__1; ++i__) {
			if (tau[i__] == 0.) {

				/*              H(i)  =  I */

				i__2 = i__;
				for (j = 1; j <= i__2; ++j) {
					t[j + i__ * t_dim1] = 0.;
					/* L10: */
				}
			} else {

				/*              general case */

				vii = v[i__ + i__ * v_dim1];
				v[i__ + i__ * v_dim1] = 1.;
				if (lsame_(storev, "C")) {

					/*                 T(1:i-1,i) := - tau(i) * V(i:n,1:i-1)' * V(i:n,i) */

					i__2 = *n - i__ + 1;
					i__3 = i__ - 1;
					d__1 = -tau[i__];
					dgemv_("Transpose", &i__2, &i__3, &d__1, &v[i__ + v_dim1], 
						ldv, &v[i__ + i__ * v_dim1], &c__1, &c_b8, &t[
							i__ * t_dim1 + 1], &c__1);
				} else {

					/*                 T(1:i-1,i) := - tau(i) * V(1:i-1,i:n) * V(i,i:n)' */

					i__2 = i__ - 1;
					i__3 = *n - i__ + 1;
					d__1 = -tau[i__];
					dgemv_("No transpose", &i__2, &i__3, &d__1, &v[i__ * 
						v_dim1 + 1], ldv, &v[i__ + i__ * v_dim1], ldv, &
						c_b8, &t[i__ * t_dim1 + 1], &c__1);
				}
				v[i__ + i__ * v_dim1] = vii;

				/*              T(1:i-1,i) := T(1:i-1,1:i-1) * T(1:i-1,i) */

				i__2 = i__ - 1;
				dtrmv_("Upper", "No transpose", "Non-unit", &i__2, &t[
					t_offset], ldt, &t[i__ * t_dim1 + 1], &c__1);
					t[i__ + i__ * t_dim1] = tau[i__];
			}
			/* L20: */
		}
	} else {
		for (i__ = *k; i__ >= 1; --i__) {
			if (tau[i__] == 0.) {

				/*              H(i)  =  I */

				i__1 = *k;
				for (j = i__; j <= i__1; ++j) {
					t[j + i__ * t_dim1] = 0.;
					/* L30: */
				}
			} else {

				/*              general case */

				if (i__ < *k) {
					if (lsame_(storev, "C")) {
						vii = v[*n - *k + i__ + i__ * v_dim1];
						v[*n - *k + i__ + i__ * v_dim1] = 1.;

						/*                    T(i+1:k,i) :=   
						- tau(i) * V(1:n-k+i,i+1:k)' * V(1:n-k+i,i) */

						i__1 = *n - *k + i__;
						i__2 = *k - i__;
						d__1 = -tau[i__];
						dgemv_("Transpose", &i__1, &i__2, &d__1, &v[(i__ + 1) 
							* v_dim1 + 1], ldv, &v[i__ * v_dim1 + 1], &
							c__1, &c_b8, &t[i__ + 1 + i__ * t_dim1], &
							c__1);
						v[*n - *k + i__ + i__ * v_dim1] = vii;
					} else {
						vii = v[i__ + (*n - *k + i__) * v_dim1];
						v[i__ + (*n - *k + i__) * v_dim1] = 1.;

						/*                    T(i+1:k,i) :=   
						- tau(i) * V(i+1:k,1:n-k+i) * V(i,1:n-k+i)' */

						i__1 = *k - i__;
						i__2 = *n - *k + i__;
						d__1 = -tau[i__];
						dgemv_("No transpose", &i__1, &i__2, &d__1, &v[i__ + 
							1 + v_dim1], ldv, &v[i__ + v_dim1], ldv, &
							c_b8, &t[i__ + 1 + i__ * t_dim1], &c__1);
						v[i__ + (*n - *k + i__) * v_dim1] = vii;
					}

					/*                 T(i+1:k,i) := T(i+1:k,i+1:k) * T(i+1:k,i) */

					i__1 = *k - i__;
					dtrmv_("Lower", "No transpose", "Non-unit", &i__1, &t[i__ 
						+ 1 + (i__ + 1) * t_dim1], ldt, &t[i__ + 1 + i__ *
						t_dim1], &c__1)
						;
				}
				t[i__ + i__ * t_dim1] = tau[i__];
			}
			/* L40: */
		}
	}
	return 0;

	/*     End of DLARFT */

} /* dlarft_ */

//
int TRL::dorg2r_(integer_ *m, integer_ *n, integer_ *k, doublereal_ *
			a, integer_ *lda, doublereal_ *tau, doublereal_ *work, integer_ *info)
{
	/*  -- LAPACK routine (version 3.1) --   
	Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..   
	November 2006   


	Purpose   
	=======   

	DORG2R generates an m by n trl_real_ matrix Q with orthonormal columns,   
	which is defined as the first n columns of a product of k elementary   
	reflectors of order m   

	Q  =  H(1) H(2) . . . H(k)   

	as returned by DGEQRF.   

	Arguments   
	=========   

	M       (input) INTEGER   
	The number of rows of the matrix Q. M >= 0.   

	N       (input) INTEGER   
	The number of columns of the matrix Q. M >= N >= 0.   

	K       (input) INTEGER   
	The number of elementary reflectors whose product defines the   
	matrix Q. N >= K >= 0.   

	A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
	On entry, the i-th column must contain the vector which   
	defines the elementary reflector H(i), for i = 1,2,...,k, as   
	returned by DGEQRF in the first k columns of its array   
	argument A.   
	On exit, the m-by-n matrix Q.   

	LDA     (input) INTEGER   
	The first dimension of the array A. LDA >= max(1,M).   

	TAU     (input) DOUBLE PRECISION array, dimension (K)   
	TAU(i) must contain the scalar factor of the elementary   
	reflector H(i), as returned by DGEQRF.   

	WORK    (workspace) DOUBLE PRECISION array, dimension (N)   

	INFO    (output) INTEGER   
	= 0: successful exit   
	< 0: if INFO = -i, the i-th argument has an illegal value   

	=====================================================================   


	Test the input arguments   

	Parameter adjustments */
	/* Table of constant values */
	static integer_ c__1 = 1;

	/* System generated locals */
	integer_ a_dim1, a_offset, i__1, i__2;
	doublereal_ d__1;
	/* Local variables */
	static integer_ i__, j, l;
	//extern /* Subroutine */ int dscal_(integer_ *, doublereal_ *, doublereal_ *, 
	//	integer_ *), dlarf_(char *, integer_ *, integer_ *, doublereal_ *, 
	//	integer_ *, doublereal_ *, doublereal_ *, integer_ *, doublereal_ *), xerbla_(char *, integer_ *);


	a_dim1 = *lda;
	a_offset = 1 + a_dim1;
	a -= a_offset;
	--tau;
	--work;

	/* Function Body */
	*info = 0;
	if (*m < 0) {
		*info = -1;
	} else if (*n < 0 || *n > *m) {
		*info = -2;
	} else if (*k < 0 || *k > *n) {
		*info = -3;
	} else if (*lda < max(1,*m)) {
		*info = -5;
	}
	if (*info != 0) {
		i__1 = -(*info);
		xerbla_("DORG2R", &i__1);
		return 0;
	}

	/*     Quick return if possible */

	if (*n <= 0) {
		return 0;
	}

	/*     Initialise columns k+1:n to columns of the unit matrix */

	i__1 = *n;
	for (j = *k + 1; j <= i__1; ++j) {
		i__2 = *m;
		for (l = 1; l <= i__2; ++l) {
			a[l + j * a_dim1] = 0.;
			/* L10: */
		}
		a[j + j * a_dim1] = 1.;
		/* L20: */
	}

	for (i__ = *k; i__ >= 1; --i__) {

		/*        Apply H(i) to A(i:m,i:n) from the left */

		if (i__ < *n) {
			a[i__ + i__ * a_dim1] = 1.;
			i__1 = *m - i__ + 1;
			i__2 = *n - i__;
			dlarf_("Left", &i__1, &i__2, &a[i__ + i__ * a_dim1], &c__1, &tau[
				i__], &a[i__ + (i__ + 1) * a_dim1], lda, &work[1]);
		}
		if (i__ < *m) {
			i__1 = *m - i__;
			d__1 = -tau[i__];
			dscal_(&i__1, &d__1, &a[i__ + 1 + i__ * a_dim1], &c__1);
		}
		a[i__ + i__ * a_dim1] = 1. - tau[i__];

		/*        Set A(1:i-1,i) to zero */

		i__1 = i__ - 1;
		for (l = 1; l <= i__1; ++l) {
			a[l + i__ * a_dim1] = 0.;
			/* L30: */
		}
		/* L40: */
	}
	return 0;

	/*     End of DORG2R */

} /* dorg2r_ */

//
int TRL::dlarf_(char *side, integer_ *m, integer_ *n, doublereal_ *v,
		   integer_ *incv, doublereal_ *tau, doublereal_ *c__, integer_ *ldc, 
		   doublereal_ *work)
{
	/*  -- LAPACK auxiliary routine (version 3.1) --   
	Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..   
	November 2006   


	Purpose   
	=======   

	DLARF applies a trl_real_ elementary reflector H to a trl_real_ m by n matrix   
	C, from either the left or the right. H is represented in the form   

	H = I - tau * v * v'   

	where tau is a trl_real_ scalar and v is a trl_real_ vector.   

	If tau = 0, then H is taken to be the unit matrix.   

	Arguments   
	=========   

	SIDE    (input) CHARACTER*1   
	= 'L': form  H * C   
	= 'R': form  C * H   

	M       (input) INTEGER   
	The number of rows of the matrix C.   

	N       (input) INTEGER   
	The number of columns of the matrix C.   

	V       (input) DOUBLE PRECISION array, dimension   
	(1 + (M-1)*abs(INCV)) if SIDE = 'L'   
	or (1 + (N-1)*abs(INCV)) if SIDE = 'R'   
	The vector v in the representation of H. V is not used if   
	TAU = 0.   

	INCV    (input) INTEGER   
	The increment between elements of v. INCV <> 0.   

	TAU     (input) DOUBLE PRECISION   
	The value tau in the representation of H.   

	C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)   
	On entry, the m by n matrix C.   
	On exit, C is overwritten by the matrix H * C if SIDE = 'L',   
	or C * H if SIDE = 'R'.   

	LDC     (input) INTEGER   
	The leading dimension of the array C. LDC >= max(1,M).   

	WORK    (workspace) DOUBLE PRECISION array, dimension   
	(N) if SIDE = 'L'   
	or (M) if SIDE = 'R'   

	=====================================================================   


	Parameter adjustments */
	/* Table of constant values */
	static doublereal_ c_b4 = 1.;
	static doublereal_ c_b5 = 0.;
	static integer_ c__1 = 1;

	/* System generated locals */
	integer_ c_dim1, c_offset;
	doublereal_ d__1;
	/* Local variables */
	//extern /* Subroutine */ int dger_(integer_ *, integer_ *, doublereal_ *, 
	//	doublereal_ *, integer_ *, doublereal_ *, integer_ *, doublereal_ *, 
	//	integer_ *);
	//extern logical_ lsame_(char *, char *);
	//extern /* Subroutine */ int dgemv_(char *, integer_ *, integer_ *, 
	//	doublereal_ *, doublereal_ *, integer_ *, doublereal_ *, integer_ *, 
	//	doublereal_ *, doublereal_ *, integer_ *);


	--v;
	c_dim1 = *ldc;
	c_offset = 1 + c_dim1;
	c__ -= c_offset;
	--work;

	/* Function Body */
	if (lsame_(side, "L")) {

		/*        Form  H * C */

		if (*tau != 0.) {

			/*           w := C' * v */

			dgemv_("Transpose", m, n, &c_b4, &c__[c_offset], ldc, &v[1], incv, 
				&c_b5, &work[1], &c__1);

			/*           C := C - v * w' */

			d__1 = -(*tau);
			dger_(m, n, &d__1, &v[1], incv, &work[1], &c__1, &c__[c_offset], 
				ldc);
		}
	} else {

		/*        Form  C * H */

		if (*tau != 0.) {

			/*           w := C * v */

			dgemv_("No transpose", m, n, &c_b4, &c__[c_offset], ldc, &v[1], 
				incv, &c_b5, &work[1], &c__1);

			/*           C := C - w * v' */

			d__1 = -(*tau);
			dger_(m, n, &d__1, &work[1], &c__1, &v[1], incv, &c__[c_offset], 
				ldc);
		}
	}
	return 0;

	/*     End of DLARF */

} /* dlarf_ */

//
int TRL::dtrmm_(char *side, char *uplo, char *transa, char *diag, 
		   integer_ *m, integer_ *n, doublereal_ *alpha, doublereal_ *a, integer_ *
		   lda, doublereal_ *b, integer_ *ldb)
{
	/* System generated locals */
	integer_ a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;
	/* Local variables */
	static integer_ i__, j, k, info;
	static doublereal_ temp;
	static logical_ lside;
	//extern logical_ lsame_(char *, char *);
	static integer_ nrowa;
	static logical_ upper;
	//extern /* Subroutine */ int xerbla_(char *, integer_ *);
	static logical_ nounit;
	/*  Purpose   
	=======   
	DTRMM  performs one of the matrix-matrix operations   
	B := alpha*op( A )*B,   or   B := alpha*B*op( A ),   
	where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or   
	non-unit,  upper or lower triangular matrix  and  op( A )  is one  of   
	op( A ) = A   or   op( A ) = A'.   
	Arguments   
	==========   
	SIDE   - CHARACTER*1.   
	On entry,  SIDE specifies whether  op( A ) multiplies B from   
	the left or right as follows:   
	SIDE = 'L' or 'l'   B := alpha*op( A )*B.   
	SIDE = 'R' or 'r'   B := alpha*B*op( A ).   
	Unchanged on exit.   
	UPLO   - CHARACTER*1.   
	On entry, UPLO specifies whether the matrix A is an upper or   
	lower triangular matrix as follows:   
	UPLO = 'U' or 'u'   A is an upper triangular matrix.   
	UPLO = 'L' or 'l'   A is a lower triangular matrix.   
	Unchanged on exit.   
	TRANSA - CHARACTER*1.   
	On entry, TRANSA specifies the form of op( A ) to be used in   
	the matrix multiplication as follows:   
	TRANSA = 'N' or 'n'   op( A ) = A.   
	TRANSA = 'T' or 't'   op( A ) = A'.   
	TRANSA = 'C' or 'c'   op( A ) = A'.   
	Unchanged on exit.   
	DIAG   - CHARACTER*1.   
	On entry, DIAG specifies whether or not A is unit triangular   
	as follows:   
	DIAG = 'U' or 'u'   A is assumed to be unit triangular.   
	DIAG = 'N' or 'n'   A is not assumed to be unit   
	triangular.   
	Unchanged on exit.   
	M      - INTEGER.   
	On entry, M specifies the number of rows of B. M must be at   
	least zero.   
	Unchanged on exit.   
	N      - INTEGER.   
	On entry, N specifies the number of columns of B.  N must be   
	at least zero.   
	Unchanged on exit.   
	ALPHA  - DOUBLE PRECISION.   
	On entry,  ALPHA specifies the scalar  alpha. When  alpha is   
	zero then  A is not referenced and  B need not be set before   
	entry.   
	Unchanged on exit.   
	A      - DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m   
	when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.   
	Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k   
	upper triangular part of the array  A must contain the upper   
	triangular matrix  and the strictly lower triangular part of   
	A is not referenced.   
	Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k   
	lower triangular part of the array  A must contain the lower   
	triangular matrix  and the strictly upper triangular part of   
	A is not referenced.   
	Note that when  DIAG = 'U' or 'u',  the diagonal elements of   
	A  are not referenced either,  but are assumed to be  unity.   
	Unchanged on exit.   
	LDA    - INTEGER.   
	On entry, LDA specifies the first dimension of A as declared   
	in the calling (sub) program.  When  SIDE = 'L' or 'l'  then   
	LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'   
	then LDA must be at least max( 1, n ).   
	Unchanged on exit.   
	B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ).   
	Before entry,  the leading  m by n part of the array  B must   
	contain the matrix  B,  and  on exit  is overwritten  by the   
	transformed matrix.   
	LDB    - INTEGER.   
	On entry, LDB specifies the first dimension of B as declared   
	in  the  calling  (sub)  program.   LDB  must  be  at  least   
	max( 1, m ).   
	Unchanged on exit.   
	Level 3 Blas routine.   
	-- Written on 8-February-1989.   
	Jack Dongarra, Argonne National Laboratory.   
	Iain Duff, AERE Harwell.   
	Jeremy Du Croz, Numerical Algorithms Group Ltd.   
	Sven Hammarling, Numerical Algorithms Group Ltd.   
	Test the input parameters.   
	Parameter adjustments */
	a_dim1 = *lda;
	a_offset = 1 + a_dim1;
	a -= a_offset;
	b_dim1 = *ldb;
	b_offset = 1 + b_dim1;
	b -= b_offset;
	/* Function Body */
	lside = lsame_(side, "L");
	if (lside) {
		nrowa = *m;
	} else {
		nrowa = *n;
	}
	nounit = lsame_(diag, "N");
	upper = lsame_(uplo, "U");
	info = 0;
	if (! lside && ! lsame_(side, "R")) {
		info = 1;
	} else if (! upper && ! lsame_(uplo, "L")) {
		info = 2;
	} else if (! lsame_(transa, "N") && ! lsame_(transa, 
		"T") && ! lsame_(transa, "C")) {
			info = 3;
	} else if (! lsame_(diag, "U") && ! lsame_(diag, 
		"N")) {
			info = 4;
	} else if (*m < 0) {
		info = 5;
	} else if (*n < 0) {
		info = 6;
	} else if (*lda < max(1,nrowa)) {
		info = 9;
	} else if (*ldb < max(1,*m)) {
		info = 11;
	}
	if (info != 0) {
		xerbla_("DTRMM ", &info);
		return 0;
	}
	/*     Quick return if possible. */
	if (*n == 0) {
		return 0;
	}
	/*     And when  alpha.eq.zero. */
	if (*alpha == 0.) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
			i__2 = *m;
			for (i__ = 1; i__ <= i__2; ++i__) {
				b[i__ + j * b_dim1] = 0.;
				/* L10: */
			}
			/* L20: */
		}
		return 0;
	}
	/*     Start the operations. */
	if (lside) {
		if (lsame_(transa, "N")) {
			/*           Form  B := alpha*A*B. */
			if (upper) {
				i__1 = *n;
				for (j = 1; j <= i__1; ++j) {
					i__2 = *m;
					for (k = 1; k <= i__2; ++k) {
						if (b[k + j * b_dim1] != 0.) {
							temp = *alpha * b[k + j * b_dim1];
							i__3 = k - 1;
							for (i__ = 1; i__ <= i__3; ++i__) {
								b[i__ + j * b_dim1] += temp * a[i__ + k * 
									a_dim1];
								/* L30: */
							}
							if (nounit) {
								temp *= a[k + k * a_dim1];
							}
							b[k + j * b_dim1] = temp;
						}
						/* L40: */
					}
					/* L50: */
				}
			} else {
				i__1 = *n;
				for (j = 1; j <= i__1; ++j) {
					for (k = *m; k >= 1; --k) {
						if (b[k + j * b_dim1] != 0.) {
							temp = *alpha * b[k + j * b_dim1];
							b[k + j * b_dim1] = temp;
							if (nounit) {
								b[k + j * b_dim1] *= a[k + k * a_dim1];
							}
							i__2 = *m;
							for (i__ = k + 1; i__ <= i__2; ++i__) {
								b[i__ + j * b_dim1] += temp * a[i__ + k * 
									a_dim1];
								/* L60: */
							}
						}
						/* L70: */
					}
					/* L80: */
				}
			}
		} else {
			/*           Form  B := alpha*A'*B. */
			if (upper) {
				i__1 = *n;
				for (j = 1; j <= i__1; ++j) {
					for (i__ = *m; i__ >= 1; --i__) {
						temp = b[i__ + j * b_dim1];
						if (nounit) {
							temp *= a[i__ + i__ * a_dim1];
						}
						i__2 = i__ - 1;
						for (k = 1; k <= i__2; ++k) {
							temp += a[k + i__ * a_dim1] * b[k + j * b_dim1];
							/* L90: */
						}
						b[i__ + j * b_dim1] = *alpha * temp;
						/* L100: */
					}
					/* L110: */
				}
			} else {
				i__1 = *n;
				for (j = 1; j <= i__1; ++j) {
					i__2 = *m;
					for (i__ = 1; i__ <= i__2; ++i__) {
						temp = b[i__ + j * b_dim1];
						if (nounit) {
							temp *= a[i__ + i__ * a_dim1];
						}
						i__3 = *m;
						for (k = i__ + 1; k <= i__3; ++k) {
							temp += a[k + i__ * a_dim1] * b[k + j * b_dim1];
							/* L120: */
						}
						b[i__ + j * b_dim1] = *alpha * temp;
						/* L130: */
					}
					/* L140: */
				}
			}
		}
	} else {
		if (lsame_(transa, "N")) {
			/*           Form  B := alpha*B*A. */
			if (upper) {
				for (j = *n; j >= 1; --j) {
					temp = *alpha;
					if (nounit) {
						temp *= a[j + j * a_dim1];
					}
					i__1 = *m;
					for (i__ = 1; i__ <= i__1; ++i__) {
						b[i__ + j * b_dim1] = temp * b[i__ + j * b_dim1];
						/* L150: */
					}
					i__1 = j - 1;
					for (k = 1; k <= i__1; ++k) {
						if (a[k + j * a_dim1] != 0.) {
							temp = *alpha * a[k + j * a_dim1];
							i__2 = *m;
							for (i__ = 1; i__ <= i__2; ++i__) {
								b[i__ + j * b_dim1] += temp * b[i__ + k * 
									b_dim1];
								/* L160: */
							}
						}
						/* L170: */
					}
					/* L180: */
				}
			} else {
				i__1 = *n;
				for (j = 1; j <= i__1; ++j) {
					temp = *alpha;
					if (nounit) {
						temp *= a[j + j * a_dim1];
					}
					i__2 = *m;
					for (i__ = 1; i__ <= i__2; ++i__) {
						b[i__ + j * b_dim1] = temp * b[i__ + j * b_dim1];
						/* L190: */
					}
					i__2 = *n;
					for (k = j + 1; k <= i__2; ++k) {
						if (a[k + j * a_dim1] != 0.) {
							temp = *alpha * a[k + j * a_dim1];
							i__3 = *m;
							for (i__ = 1; i__ <= i__3; ++i__) {
								b[i__ + j * b_dim1] += temp * b[i__ + k * 
									b_dim1];
								/* L200: */
							}
						}
						/* L210: */
					}
					/* L220: */
				}
			}
		} else {
			/*           Form  B := alpha*B*A'. */
			if (upper) {
				i__1 = *n;
				for (k = 1; k <= i__1; ++k) {
					i__2 = k - 1;
					for (j = 1; j <= i__2; ++j) {
						if (a[j + k * a_dim1] != 0.) {
							temp = *alpha * a[j + k * a_dim1];
							i__3 = *m;
							for (i__ = 1; i__ <= i__3; ++i__) {
								b[i__ + j * b_dim1] += temp * b[i__ + k * 
									b_dim1];
								/* L230: */
							}
						}
						/* L240: */
					}
					temp = *alpha;
					if (nounit) {
						temp *= a[k + k * a_dim1];
					}
					if (temp != 1.) {
						i__2 = *m;
						for (i__ = 1; i__ <= i__2; ++i__) {
							b[i__ + k * b_dim1] = temp * b[i__ + k * b_dim1];
							/* L250: */
						}
					}
					/* L260: */
				}
			} else {
				for (k = *n; k >= 1; --k) {
					i__1 = *n;
					for (j = k + 1; j <= i__1; ++j) {
						if (a[j + k * a_dim1] != 0.) {
							temp = *alpha * a[j + k * a_dim1];
							i__2 = *m;
							for (i__ = 1; i__ <= i__2; ++i__) {
								b[i__ + j * b_dim1] += temp * b[i__ + k * 
									b_dim1];
								/* L270: */
							}
						}
						/* L280: */
					}
					temp = *alpha;
					if (nounit) {
						temp *= a[k + k * a_dim1];
					}
					if (temp != 1.) {
						i__1 = *m;
						for (i__ = 1; i__ <= i__1; ++i__) {
							b[i__ + k * b_dim1] = temp * b[i__ + k * b_dim1];
							/* L290: */
						}
					}
					/* L300: */
				}
			}
		}
	}
	return 0;
	/*     End of DTRMM . */
} /* dtrmm_ */



//
int TRL::dtrmv_(char *uplo, char *trans, char *diag, integer_ *n, 
		   doublereal_ *a, integer_ *lda, doublereal_ *x, integer_ *incx)
{
	/* System generated locals */
	integer_ a_dim1, a_offset, i__1, i__2;
	/* Local variables */
	static integer_ i__, j, ix, jx, kx, info;
	static doublereal_ temp;
	//extern logical_ lsame_(char *, char *);
	//extern /* Subroutine */ int xerbla_(char *, integer_ *);
	static logical_ nounit;
	/*  Purpose   
	=======   
	DTRMV  performs one of the matrix-vector operations   
	x := A*x,   or   x := A'*x,   
	where x is an n element vector and  A is an n by n unit, or non-unit,   
	upper or lower triangular matrix.   
	Arguments   
	==========   
	UPLO   - CHARACTER*1.   
	On entry, UPLO specifies whether the matrix is an upper or   
	lower triangular matrix as follows:   
	UPLO = 'U' or 'u'   A is an upper triangular matrix.   
	UPLO = 'L' or 'l'   A is a lower triangular matrix.   
	Unchanged on exit.   
	TRANS  - CHARACTER*1.   
	On entry, TRANS specifies the operation to be performed as   
	follows:   
	TRANS = 'N' or 'n'   x := A*x.   
	TRANS = 'T' or 't'   x := A'*x.   
	TRANS = 'C' or 'c'   x := A'*x.   
	Unchanged on exit.   
	DIAG   - CHARACTER*1.   
	On entry, DIAG specifies whether or not A is unit   
	triangular as follows:   
	DIAG = 'U' or 'u'   A is assumed to be unit triangular.   
	DIAG = 'N' or 'n'   A is not assumed to be unit   
	triangular.   
	Unchanged on exit.   
	N      - INTEGER.   
	On entry, N specifies the order of the matrix A.   
	N must be at least zero.   
	Unchanged on exit.   
	A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).   
	Before entry with  UPLO = 'U' or 'u', the leading n by n   
	upper triangular part of the array A must contain the upper   
	triangular matrix and the strictly lower triangular part of   
	A is not referenced.   
	Before entry with UPLO = 'L' or 'l', the leading n by n   
	lower triangular part of the array A must contain the lower   
	triangular matrix and the strictly upper triangular part of   
	A is not referenced.   
	Note that when  DIAG = 'U' or 'u', the diagonal elements of   
	A are not referenced either, but are assumed to be unity.   
	Unchanged on exit.   
	LDA    - INTEGER.   
	On entry, LDA specifies the first dimension of A as declared   
	in the calling (sub) program. LDA must be at least   
	max( 1, n ).   
	Unchanged on exit.   
	X      - DOUBLE PRECISION array of dimension at least   
	( 1 + ( n - 1 )*abs( INCX ) ).   
	Before entry, the incremented array X must contain the n   
	element vector x. On exit, X is overwritten with the   
	tranformed vector x.   
	INCX   - INTEGER.   
	On entry, INCX specifies the increment for the elements of   
	X. INCX must not be zero.   
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
	a_offset = 1 + a_dim1;
	a -= a_offset;
	--x;
	/* Function Body */
	info = 0;
	if (! lsame_(uplo, "U") && ! lsame_(uplo, "L")) {
		info = 1;
	} else if (! lsame_(trans, "N") && ! lsame_(trans, 
		"T") && ! lsame_(trans, "C")) {
			info = 2;
	} else if (! lsame_(diag, "U") && ! lsame_(diag, 
		"N")) {
			info = 3;
	} else if (*n < 0) {
		info = 4;
	} else if (*lda < max(1,*n)) {
		info = 6;
	} else if (*incx == 0) {
		info = 8;
	}
	if (info != 0) {
		xerbla_("DTRMV ", &info);
		return 0;
	}
	/*     Quick return if possible. */
	if (*n == 0) {
		return 0;
	}
	nounit = lsame_(diag, "N");
	/*     Set up the start point in X if the increment is not unity. This   
	will be  ( N - 1 )*INCX  too small for descending loops. */
	if (*incx <= 0) {
		kx = 1 - (*n - 1) * *incx;
	} else if (*incx != 1) {
		kx = 1;
	}
	/*     Start the operations. In this version the elements of A are   
	accessed sequentially with one pass through A. */
	if (lsame_(trans, "N")) {
		/*        Form  x := A*x. */
		if (lsame_(uplo, "U")) {
			if (*incx == 1) {
				i__1 = *n;
				for (j = 1; j <= i__1; ++j) {
					if (x[j] != 0.) {
						temp = x[j];
						i__2 = j - 1;
						for (i__ = 1; i__ <= i__2; ++i__) {
							x[i__] += temp * a[i__ + j * a_dim1];
							/* L10: */
						}
						if (nounit) {
							x[j] *= a[j + j * a_dim1];
						}
					}
					/* L20: */
				}
			} else {
				jx = kx;
				i__1 = *n;
				for (j = 1; j <= i__1; ++j) {
					if (x[jx] != 0.) {
						temp = x[jx];
						ix = kx;
						i__2 = j - 1;
						for (i__ = 1; i__ <= i__2; ++i__) {
							x[ix] += temp * a[i__ + j * a_dim1];
							ix += *incx;
							/* L30: */
						}
						if (nounit) {
							x[jx] *= a[j + j * a_dim1];
						}
					}
					jx += *incx;
					/* L40: */
				}
			}
		} else {
			if (*incx == 1) {
				for (j = *n; j >= 1; --j) {
					if (x[j] != 0.) {
						temp = x[j];
						i__1 = j + 1;
						for (i__ = *n; i__ >= i__1; --i__) {
							x[i__] += temp * a[i__ + j * a_dim1];
							/* L50: */
						}
						if (nounit) {
							x[j] *= a[j + j * a_dim1];
						}
					}
					/* L60: */
				}
			} else {
				kx += (*n - 1) * *incx;
				jx = kx;
				for (j = *n; j >= 1; --j) {
					if (x[jx] != 0.) {
						temp = x[jx];
						ix = kx;
						i__1 = j + 1;
						for (i__ = *n; i__ >= i__1; --i__) {
							x[ix] += temp * a[i__ + j * a_dim1];
							ix -= *incx;
							/* L70: */
						}
						if (nounit) {
							x[jx] *= a[j + j * a_dim1];
						}
					}
					jx -= *incx;
					/* L80: */
				}
			}
		}
	} else {
		/*        Form  x := A'*x. */
		if (lsame_(uplo, "U")) {
			if (*incx == 1) {
				for (j = *n; j >= 1; --j) {
					temp = x[j];
					if (nounit) {
						temp *= a[j + j * a_dim1];
					}
					for (i__ = j - 1; i__ >= 1; --i__) {
						temp += a[i__ + j * a_dim1] * x[i__];
						/* L90: */
					}
					x[j] = temp;
					/* L100: */
				}
			} else {
				jx = kx + (*n - 1) * *incx;
				for (j = *n; j >= 1; --j) {
					temp = x[jx];
					ix = jx;
					if (nounit) {
						temp *= a[j + j * a_dim1];
					}
					for (i__ = j - 1; i__ >= 1; --i__) {
						ix -= *incx;
						temp += a[i__ + j * a_dim1] * x[ix];
						/* L110: */
					}
					x[jx] = temp;
					jx -= *incx;
					/* L120: */
				}
			}
		} else {
			if (*incx == 1) {
				i__1 = *n;
				for (j = 1; j <= i__1; ++j) {
					temp = x[j];
					if (nounit) {
						temp *= a[j + j * a_dim1];
					}
					i__2 = *n;
					for (i__ = j + 1; i__ <= i__2; ++i__) {
						temp += a[i__ + j * a_dim1] * x[i__];
						/* L130: */
					}
					x[j] = temp;
					/* L140: */
				}
			} else {
				jx = kx;
				i__1 = *n;
				for (j = 1; j <= i__1; ++j) {
					temp = x[jx];
					ix = jx;
					if (nounit) {
						temp *= a[j + j * a_dim1];
					}
					i__2 = *n;
					for (i__ = j + 1; i__ <= i__2; ++i__) {
						ix += *incx;
						temp += a[i__ + j * a_dim1] * x[ix];
						/* L150: */
					}
					x[jx] = temp;
					jx += *incx;
					/* L160: */
				}
			}
		}
	}
	return 0;
	/*     End of DTRMV . */
} /* dtrmv_ */


//
int TRL::dger_(integer_ *m, integer_ *n, doublereal_ *alpha, 
		  doublereal_ *x, integer_ *incx, doublereal_ *y, integer_ *incy, 
		  doublereal_ *a, integer_ *lda)
{
	/* System generated locals */
	integer_ a_dim1, a_offset, i__1, i__2;
	/* Local variables */
	static integer_ i__, j, ix, jy, kx, info;
	static doublereal_ temp;
	//extern /* Subroutine */ int xerbla_(char *, integer_ *);
	/*  Purpose   
	=======   
	DGER   performs the rank 1 operation   
	A := alpha*x*y' + A,   
	where alpha is a scalar, x is an m element vector, y is an n element   
	vector and A is an m by n matrix.   
	Arguments   
	==========   
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
	X      - DOUBLE PRECISION array of dimension at least   
	( 1 + ( m - 1 )*abs( INCX ) ).   
	Before entry, the incremented array X must contain the m   
	element vector x.   
	Unchanged on exit.   
	INCX   - INTEGER.   
	On entry, INCX specifies the increment for the elements of   
	X. INCX must not be zero.   
	Unchanged on exit.   
	Y      - DOUBLE PRECISION array of dimension at least   
	( 1 + ( n - 1 )*abs( INCY ) ).   
	Before entry, the incremented array Y must contain the n   
	element vector y.   
	Unchanged on exit.   
	INCY   - INTEGER.   
	On entry, INCY specifies the increment for the elements of   
	Y. INCY must not be zero.   
	Unchanged on exit.   
	A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).   
	Before entry, the leading m by n part of the array A must   
	contain the matrix of coefficients. On exit, A is   
	overwritten by the updated matrix.   
	LDA    - INTEGER.   
	On entry, LDA specifies the first dimension of A as declared   
	in the calling (sub) program. LDA must be at least   
	max( 1, m ).   
	Unchanged on exit.   
	Level 2 Blas routine.   
	-- Written on 22-October-1986.   
	Jack Dongarra, Argonne National Lab.   
	Jeremy Du Croz, Nag Central Office.   
	Sven Hammarling, Nag Central Office.   
	Richard Hanson, Sandia National Labs.   
	Test the input parameters.   
	Parameter adjustments */
	--x;
	--y;
	a_dim1 = *lda;
	a_offset = 1 + a_dim1;
	a -= a_offset;
	/* Function Body */
	info = 0;
	if (*m < 0) {
		info = 1;
	} else if (*n < 0) {
		info = 2;
	} else if (*incx == 0) {
		info = 5;
	} else if (*incy == 0) {
		info = 7;
	} else if (*lda < max(1,*m)) {
		info = 9;
	}
	if (info != 0) {
		xerbla_("DGER  ", &info);
		return 0;
	}
	/*     Quick return if possible. */
	if (*m == 0 || *n == 0 || *alpha == 0.) {
		return 0;
	}
	/*     Start the operations. In this version the elements of A are   
	accessed sequentially with one pass through A. */
	if (*incy > 0) {
		jy = 1;
	} else {
		jy = 1 - (*n - 1) * *incy;
	}
	if (*incx == 1) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
			if (y[jy] != 0.) {
				temp = *alpha * y[jy];
				i__2 = *m;
				for (i__ = 1; i__ <= i__2; ++i__) {
					a[i__ + j * a_dim1] += x[i__] * temp;
					/* L10: */
				}
			}
			jy += *incy;
			/* L20: */
		}
	} else {
		if (*incx > 0) {
			kx = 1;
		} else {
			kx = 1 - (*m - 1) * *incx;
		}
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
			if (y[jy] != 0.) {
				temp = *alpha * y[jy];
				ix = kx;
				i__2 = *m;
				for (i__ = 1; i__ <= i__2; ++i__) {
					a[i__ + j * a_dim1] += x[ix] * temp;
					ix += *incx;
					/* L30: */
				}
			}
			jy += *incy;
			/* L40: */
		}
	}
	return 0;
	/*     End of DGER  . */
} /* dger_ */



/////////////////////////////////////
///////////////////////////////
///////////////////////////////////////
//////////////////////////

double TRL::pow_di(doublereal_ *ap, integer_ *bp)
{
	double pow, x;
	integer_ n;
	unsigned long u;

	pow = 1;
	x = *ap;
	n = *bp;

	if(n != 0)
	{
		if(n < 0)
		{
			n = -n;
			x = 1/x;
		}
		for(u = n; ; )
		{
			if(u & 01)
				pow *= x;
			if(u >>= 1)
				x *= x;
			else
				break;
		}
	}
	return(pow);
}

/*integer_ s_wsfe(cilist_ *a)*/

void TRL::s_copy(register char *a, register char *b, ftnlen_ la, ftnlen_ lb)
{
	register char *aend, *bend;

	aend = a + la;

	if(la <= lb)
#ifndef NO_OVERWRITE
		if (a <= b || a >= b + la)
#endif
			while(a < aend)
				*a++ = *b++;
#ifndef NO_OVERWRITE
		else
			for(b += la; a < aend; )
				*--aend = *--b;
#endif

	else {
		bend = b + lb;
#ifndef NO_OVERWRITE
		if (a <= b || a >= bend)
#endif
			while(b < bend)
				*a++ = *b++;
#ifndef NO_OVERWRITE
		else {
			a += lb;
			while(b < bend)
				*--a = *--bend;
			a += lb;
		}
#endif
		while(a < aend)
			*a++ = ' ';
	}
}
//
integer_ TRL::s_cmp(char *a0, char *b0, ftnlen_ la, ftnlen_ lb)
{
	register unsigned char *a, *aend, *b, *bend;
	a = (unsigned char *)a0;
	b = (unsigned char *)b0;
	aend = a + la;
	bend = b + lb;

	if(la <= lb)
	{
		while(a < aend)
			if(*a != *b)
				return( *a - *b );
			else
			{ ++a; ++b; }

			while(b < bend)
				if(*b != ' ')
					return( ' ' - *b );
				else	++b;
	}

	else
	{
		while(b < bend)
			if(*a == *b)
			{ ++a; ++b; }
			else
				return( *a - *b );
		while(a < aend)
			if(*a != ' ')
				return(*a - ' ');
			else	++a;
	}
	return(0);
}

double TRL::d_sign(doublereal_ *a, doublereal_ *b)
{
	double x;
	x = (*a >= 0 ? *a : - *a);
	return( *b >= 0 ? x : -x);
}
//

integer_ TRL::i_nint(trl_real_ *x)
{
	return (integer_)(*x >= 0 ? floor(*x + .5) : -floor(.5 - *x));
}

double TRL::sqrt_(doublereal_ x)
{
	return sqrt(double(x));
}

double TRL::log_(doublereal_ x)
{
	return log(double(x));
}
double TRL::cos_(doublereal_ x)
{
	return cos(double(x));
}



//////////////////////////
////////////////
////////////////////////
int TRL::dsteqr_(char *compz, integer_ *n, doublereal_ *d__, 
			doublereal_ *e, doublereal_ *z__, integer_ *ldz, doublereal_ *work, 
			integer_ *info)
{
	/*  -- LAPACK routine (version 3.1) --   
	Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..   
	November 2006   


	Purpose   
	=======   

	DSTEQR computes all eigenvalues and, optionally, eigenvectors of a   
	symmetric tridiagonal matrix using the implicit QL or QR method.   
	The eigenvectors of a full or band symmetric matrix can also be found   
	if DSYTRD or DSPTRD or DSBTRD has been used to reduce this matrix to   
	tridiagonal form.   

	Arguments   
	=========   

	COMPZ   (input) CHARACTER*1   
	= 'N':  Compute eigenvalues only.   
	= 'V':  Compute eigenvalues and eigenvectors of the original   
	symmetric matrix.  On entry, Z must contain the   
	orthogonal matrix used to reduce the original matrix   
	to tridiagonal form.   
	= 'I':  Compute eigenvalues and eigenvectors of the   
	tridiagonal matrix.  Z is initialized to the identity   
	matrix.   

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

	LDZ     (input) INTEGER   
	The leading dimension of the array Z.  LDZ >= 1, and if   
	eigenvectors are desired, then  LDZ >= max(1,N).   

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


	Test the input parameters.   

	Parameter adjustments */
	/* Table of constant values */
	static doublereal_ c_b9 = 0.;
	static doublereal_ c_b10 = 1.;
	static integer_ c__0 = 0;
	static integer_ c__1 = 1;
	static integer_ c__2 = 2;

	/* System generated locals */
	integer_ z_dim1, z_offset, i__1, i__2;
	doublereal_ d__1, d__2;
	/* Builtin functions */
	//double sqrt_(doublereal_), d_sign(doublereal_ *, doublereal_ *);
	/* Local variables */
	static doublereal_ b, c__, f, g;
	static integer_ i__, j, k, l, m;
	static doublereal_ p, r__, s;
	static integer_ l1, ii, mm, lm1, mm1, nm1;
	static doublereal_ rt1, rt2, eps;
	static integer_ lsv;
	static doublereal_ tst, eps2;
	static integer_ lend, jtot;
	//extern /* Subroutine */ int dlae2_(doublereal_ *, doublereal_ *, doublereal_ 
	//	*, doublereal_ *, doublereal_ *);
	//extern logical_ lsame_(char *, char *);
	//extern /* Subroutine */ int dlasr_(char *, char *, char *, integer_ *, 
	//	integer_ *, doublereal_ *, doublereal_ *, doublereal_ *, integer_ *);
	static doublereal_ anorm;
	//extern /* Subroutine */ int dswap_(integer_ *, doublereal_ *, integer_ *, 
	//	doublereal_ *, integer_ *), dlaev2_(doublereal_ *, doublereal_ *, 
	//	doublereal_ *, doublereal_ *, doublereal_ *, doublereal_ *, 
	//	doublereal_ *);
	static integer_ lendm1, lendp1;
	//extern doublereal_ dlapy2_(doublereal_ *, doublereal_ *), dlamch_(char *);
	static integer_ iscale;
	//extern /* Subroutine */ int dlascl_(char *, integer_ *, integer_ *, 
	//	doublereal_ *, doublereal_ *, integer_ *, integer_ *, doublereal_ *, 
	//	integer_ *, integer_ *), dlaset_(char *, integer_ *, integer_ 
	//	*, doublereal_ *, doublereal_ *, doublereal_ *, integer_ *);
	static doublereal_ safmin;
	//extern /* Subroutine */ int dlartg_(doublereal_ *, doublereal_ *, 
	//	doublereal_ *, doublereal_ *, doublereal_ *);
	static doublereal_ safmax;
	//extern /* Subroutine */ int xerbla_(char *, integer_ *);
	//extern doublereal_ dlanst_(char *, integer_ *, doublereal_ *, doublereal_ *);
	//extern /* Subroutine */ int dlasrt_(char *, integer_ *, doublereal_ *, 
	//	integer_ *);
	static integer_ lendsv;
	static doublereal_ ssfmin;
	static integer_ nmaxit, icompz;
	static doublereal_ ssfmax;


	--d__;
	--e;
	z_dim1 = *ldz;
	z_offset = 1 + z_dim1;
	z__ -= z_offset;
	--work;

	/* Function Body */
	*info = 0;

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

	/*     Quick return if possible */

	if (*n == 0) {
		return 0;
	}

	if (*n == 1) {
		if (icompz == 2) {
			z__[z_dim1 + 1] = 1.;
		}
		return 0;
	}

	/*     Determine the unit roundoff and over/underflow thresholds. */

	eps = dlamch_("E");
	/* Computing 2nd power */
	d__1 = eps;
	eps2 = d__1 * d__1;
	safmin = dlamch_("S");
	safmax = 1. / safmin;
	ssfmax = sqrt_(safmax) / 3.;
	ssfmin = sqrt_(safmin) / eps2;

	/*     Compute the eigenvalues and eigenvectors of the tridiagonal   
	matrix. */

	if (icompz == 2) {
		dlaset_("Full", n, n, &c_b9, &c_b10, &z__[z_offset], ldz);
	}

	nmaxit = *n * 30;
	jtot = 0;

	/*     Determine where the matrix splits and choose QL or QR iteration   
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
			tst = (d__1 = e[m], abs(d__1));
			if (tst == 0.) {
				goto L30;
			}
			if (tst <= sqrt_((d__1 = d__[m], abs(d__1))) * sqrt_((d__2 = d__[m 
				+ 1], abs(d__2))) * eps) {
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

	/*     Scale submatrix in rows and columns L to LEND */

	i__1 = lend - l + 1;
	anorm = dlanst_("I", &i__1, &d__[l], &e[l]);
	iscale = 0;
	if (anorm == 0.) {
		goto L10;
	}
	if (anorm > ssfmax) {
		iscale = 1;
		i__1 = lend - l + 1;
		dlascl_("G", &c__0, &c__0, &anorm, &ssfmax, &i__1, &c__1, &d__[l], n, 
			info);
		i__1 = lend - l;
		dlascl_("G", &c__0, &c__0, &anorm, &ssfmax, &i__1, &c__1, &e[l], n, 
			info);
	} else if (anorm < ssfmin) {
		iscale = 2;
		i__1 = lend - l + 1;
		dlascl_("G", &c__0, &c__0, &anorm, &ssfmin, &i__1, &c__1, &d__[l], n, 
			info);
		i__1 = lend - l;
		dlascl_("G", &c__0, &c__0, &anorm, &ssfmin, &i__1, &c__1, &e[l], n, 
			info);
	}

	/*     Choose between QL and QR iteration */

	if ((d__1 = d__[lend], abs(d__1)) < (d__2 = d__[l], abs(d__2))) {
		lend = lsv;
		l = lendsv;
	}

	if (lend > l) {

		/*        QL Iteration   

		Look for small subdiagonal element. */

L40:
		if (l != lend) {
			lendm1 = lend - 1;
			i__1 = lendm1;
			for (m = l; m <= i__1; ++m) {
				/* Computing 2nd power */
				d__2 = (d__1 = e[m], abs(d__1));
				tst = d__2 * d__2;
				if (tst <= eps2 * (d__1 = d__[m], abs(d__1)) * (d__2 = d__[m 
					+ 1], abs(d__2)) + safmin) {
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

		/*        If remaining matrix is 2-by-2, use DLAE2 or SLAEV2   
		to compute its eigensystem. */

		if (m == l + 1) {
			if (icompz > 0) {
				dlaev2_(&d__[l], &e[l], &d__[l + 1], &rt1, &rt2, &c__, &s);
				work[l] = c__;
				work[*n - 1 + l] = s;
				dlasr_("R", "V", "B", n, &c__2, &work[l], &work[*n - 1 + l], &
					z__[l * z_dim1 + 1], ldz);
			} else {
				dlae2_(&d__[l], &e[l], &d__[l + 1], &rt1, &rt2);
			}
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

		/*        Form shift. */

		g = (d__[l + 1] - p) / (e[l] * 2.);
		r__ = dlapy2_(&g, &c_b10);
		g = d__[m] - p + e[l] / (g + d_sign(&r__, &g));

		s = 1.;
		c__ = 1.;
		p = 0.;

		/*        Inner loop */

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

			/*           If eigenvectors are desired, then save rotations. */

			if (icompz > 0) {
				work[i__] = c__;
				work[*n - 1 + i__] = -s;
			}

			/* L70: */
		}

		/*        If eigenvectors are desired, then apply saved rotations. */

		if (icompz > 0) {
			mm = m - l + 1;
			dlasr_("R", "V", "B", n, &mm, &work[l], &work[*n - 1 + l], &z__[l 
				* z_dim1 + 1], ldz);
		}

		d__[l] -= p;
		e[l] = g;
		goto L40;

		/*        Eigenvalue found. */

L80:
		d__[l] = p;

		++l;
		if (l <= lend) {
			goto L40;
		}
		goto L140;

	} else {

		/*        QR Iteration   

		Look for small superdiagonal element. */

L90:
		if (l != lend) {
			lendp1 = lend + 1;
			i__1 = lendp1;
			for (m = l; m >= i__1; --m) {
				/* Computing 2nd power */
				d__2 = (d__1 = e[m - 1], abs(d__1));
				tst = d__2 * d__2;
				if (tst <= eps2 * (d__1 = d__[m], abs(d__1)) * (d__2 = d__[m 
					- 1], abs(d__2)) + safmin) {
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

		/*        If remaining matrix is 2-by-2, use DLAE2 or SLAEV2   
		to compute its eigensystem. */

		if (m == l - 1) {
			if (icompz > 0) {
				dlaev2_(&d__[l - 1], &e[l - 1], &d__[l], &rt1, &rt2, &c__, &s)
					;
				work[m] = c__;
				work[*n - 1 + m] = s;
				dlasr_("R", "V", "F", n, &c__2, &work[m], &work[*n - 1 + m], &
					z__[(l - 1) * z_dim1 + 1], ldz);
			} else {
				dlae2_(&d__[l - 1], &e[l - 1], &d__[l], &rt1, &rt2);
			}
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

		/*        Form shift. */

		g = (d__[l - 1] - p) / (e[l - 1] * 2.);
		r__ = dlapy2_(&g, &c_b10);
		g = d__[m] - p + e[l - 1] / (g + d_sign(&r__, &g));

		s = 1.;
		c__ = 1.;
		p = 0.;

		/*        Inner loop */

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

			/*           If eigenvectors are desired, then save rotations. */

			if (icompz > 0) {
				work[i__] = c__;
				work[*n - 1 + i__] = s;
			}

			/* L120: */
		}

		/*        If eigenvectors are desired, then apply saved rotations. */

		if (icompz > 0) {
			mm = l - m + 1;
			dlasr_("R", "V", "F", n, &mm, &work[m], &work[*n - 1 + m], &z__[m 
				* z_dim1 + 1], ldz);
		}

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
		dlascl_("G", &c__0, &c__0, &ssfmax, &anorm, &i__1, &c__1, &d__[lsv], 
			n, info);
		i__1 = lendsv - lsv;
		dlascl_("G", &c__0, &c__0, &ssfmax, &anorm, &i__1, &c__1, &e[lsv], n, 
			info);
	} else if (iscale == 2) {
		i__1 = lendsv - lsv + 1;
		dlascl_("G", &c__0, &c__0, &ssfmin, &anorm, &i__1, &c__1, &d__[lsv], 
			n, info);
		i__1 = lendsv - lsv;
		dlascl_("G", &c__0, &c__0, &ssfmin, &anorm, &i__1, &c__1, &e[lsv], n, 
			info);
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
	if (icompz == 0) {

		/*        Use Quick Sort */

		dlasrt_("I", n, &d__[1], info);

	} else {

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
				dswap_(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[k * z_dim1 + 1], 
					&c__1);
			}
			/* L180: */
		}
	}

L190:
	return 0;

	/*     End of DSTEQR */

} /* dsteqr_ */

////
int TRL::dlasr_(char *side, char *pivot, char *direct, integer_ *m,
		   integer_ *n, doublereal_ *c__, doublereal_ *s, doublereal_ *a, integer_ *
		   lda)
{
	/*  -- LAPACK auxiliary routine (version 3.1) --   
	Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..   
	November 2006   


	Purpose   
	=======   

	DLASR applies a sequence of plane rotations to a trl_real_ matrix A,   
	from either the left or the right.   

	When SIDE = 'L', the transformation takes the form   

	A := P*A   

	and when SIDE = 'R', the transformation takes the form   

	A := A*P**T   

	where P is an orthogonal matrix consisting of a sequence of z plane   
	rotations, with z = M when SIDE = 'L' and z = N when SIDE = 'R',   
	and P**T is the transpose of P.   

	When DIRECT = 'F' (Forward sequence), then   

	P = P(z-1) * ... * P(2) * P(1)   

	and when DIRECT = 'B' (Backward sequence), then   

	P = P(1) * P(2) * ... * P(z-1)   

	where P(k) is a plane rotation matrix defined by the 2-by-2 rotation   

	R(k) = (  c(k)  s(k) )   
	= ( -s(k)  c(k) ).   

	When PIVOT = 'V' (Variable pivot), the rotation is performed   
	for the plane (k,k+1), i.e., P(k) has the form   

	P(k) = (  1                                            )   
	(       ...                                     )   
	(              1                                )   
	(                   c(k)  s(k)                  )   
	(                  -s(k)  c(k)                  )   
	(                                1              )   
	(                                     ...       )   
	(                                            1  )   

	where R(k) appears as a rank-2 modification to the identity matrix in   
	rows and columns k and k+1.   

	When PIVOT = 'T' (Top pivot), the rotation is performed for the   
	plane (1,k+1), so P(k) has the form   

	P(k) = (  c(k)                    s(k)                 )   
	(         1                                     )   
	(              ...                              )   
	(                     1                         )   
	( -s(k)                    c(k)                 )   
	(                                 1             )   
	(                                      ...      )   
	(                                             1 )   

	where R(k) appears in rows and columns 1 and k+1.   

	Similarly, when PIVOT = 'B' (Bottom pivot), the rotation is   
	performed for the plane (k,z), giving P(k) the form   

	P(k) = ( 1                                             )   
	(      ...                                      )   
	(             1                                 )   
	(                  c(k)                    s(k) )   
	(                         1                     )   
	(                              ...              )   
	(                                     1         )   
	(                 -s(k)                    c(k) )   

	where R(k) appears in rows and columns k and z.  The rotations are   
	performed without ever forming P(k) explicitly.   

	Arguments   
	=========   

	SIDE    (input) CHARACTER*1   
	Specifies whether the plane rotation matrix P is applied to   
	A on the left or the right.   
	= 'L':  Left, compute A := P*A   
	= 'R':  Right, compute A:= A*P**T   

	PIVOT   (input) CHARACTER*1   
	Specifies the plane for which P(k) is a plane rotation   
	matrix.   
	= 'V':  Variable pivot, the plane (k,k+1)   
	= 'T':  Top pivot, the plane (1,k+1)   
	= 'B':  Bottom pivot, the plane (k,z)   

	DIRECT  (input) CHARACTER*1   
	Specifies whether P is a forward or backward sequence of   
	plane rotations.   
	= 'F':  Forward, P = P(z-1)*...*P(2)*P(1)   
	= 'B':  Backward, P = P(1)*P(2)*...*P(z-1)   

	M       (input) INTEGER   
	The number of rows of the matrix A.  If m <= 1, an immediate   
	return is effected.   

	N       (input) INTEGER   
	The number of columns of the matrix A.  If n <= 1, an   
	immediate return is effected.   

	C       (input) DOUBLE PRECISION array, dimension   
	(M-1) if SIDE = 'L'   
	(N-1) if SIDE = 'R'   
	The cosines c(k) of the plane rotations.   

	S       (input) DOUBLE PRECISION array, dimension   
	(M-1) if SIDE = 'L'   
	(N-1) if SIDE = 'R'   
	The sines s(k) of the plane rotations.  The 2-by-2 plane   
	rotation part of the matrix P(k), R(k), has the form   
	R(k) = (  c(k)  s(k) )   
	( -s(k)  c(k) ).   

	A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
	The M-by-N matrix A.  On exit, A is overwritten by P*A if   
	SIDE = 'R' or by A*P**T if SIDE = 'L'.   

	LDA     (input) INTEGER   
	The leading dimension of the array A.  LDA >= max(1,M).   

	=====================================================================   


	Test the input parameters   

	Parameter adjustments */
	/* System generated locals */
	integer_ a_dim1, a_offset, i__1, i__2;
	/* Local variables */
	static integer_ i__, j, info;
	static doublereal_ temp;
	//extern logical_ lsame_(char *, char *);
	static doublereal_ ctemp, stemp;
	//extern /* Subroutine */ int xerbla_(char *, integer_ *);

	--c__;
	--s;
	a_dim1 = *lda;
	a_offset = 1 + a_dim1;
	a -= a_offset;

	/* Function Body */
	info = 0;
	if (! (lsame_(side, "L") || lsame_(side, "R"))) {
		info = 1;
	} else if (! (lsame_(pivot, "V") || lsame_(pivot, 
		"T") || lsame_(pivot, "B"))) {
			info = 2;
	} else if (! (lsame_(direct, "F") || lsame_(direct, 
		"B"))) {
			info = 3;
	} else if (*m < 0) {
		info = 4;
	} else if (*n < 0) {
		info = 5;
	} else if (*lda < max(1,*m)) {
		info = 9;
	}
	if (info != 0) {
		xerbla_("DLASR ", &info);
		return 0;
	}

	/*     Quick return if possible */

	if (*m == 0 || *n == 0) {
		return 0;
	}
	if (lsame_(side, "L")) {

		/*        Form  P * A */

		if (lsame_(pivot, "V")) {
			if (lsame_(direct, "F")) {
				i__1 = *m - 1;
				for (j = 1; j <= i__1; ++j) {
					ctemp = c__[j];
					stemp = s[j];
					if (ctemp != 1. || stemp != 0.) {
						i__2 = *n;
						for (i__ = 1; i__ <= i__2; ++i__) {
							temp = a[j + 1 + i__ * a_dim1];
							a[j + 1 + i__ * a_dim1] = ctemp * temp - stemp * 
								a[j + i__ * a_dim1];
							a[j + i__ * a_dim1] = stemp * temp + ctemp * a[j 
								+ i__ * a_dim1];
							/* L10: */
						}
					}
					/* L20: */
				}
			} else if (lsame_(direct, "B")) {
				for (j = *m - 1; j >= 1; --j) {
					ctemp = c__[j];
					stemp = s[j];
					if (ctemp != 1. || stemp != 0.) {
						i__1 = *n;
						for (i__ = 1; i__ <= i__1; ++i__) {
							temp = a[j + 1 + i__ * a_dim1];
							a[j + 1 + i__ * a_dim1] = ctemp * temp - stemp * 
								a[j + i__ * a_dim1];
							a[j + i__ * a_dim1] = stemp * temp + ctemp * a[j 
								+ i__ * a_dim1];
							/* L30: */
						}
					}
					/* L40: */
				}
			}
		} else if (lsame_(pivot, "T")) {
			if (lsame_(direct, "F")) {
				i__1 = *m;
				for (j = 2; j <= i__1; ++j) {
					ctemp = c__[j - 1];
					stemp = s[j - 1];
					if (ctemp != 1. || stemp != 0.) {
						i__2 = *n;
						for (i__ = 1; i__ <= i__2; ++i__) {
							temp = a[j + i__ * a_dim1];
							a[j + i__ * a_dim1] = ctemp * temp - stemp * a[
								i__ * a_dim1 + 1];
								a[i__ * a_dim1 + 1] = stemp * temp + ctemp * a[
									i__ * a_dim1 + 1];
									/* L50: */
						}
					}
					/* L60: */
				}
			} else if (lsame_(direct, "B")) {
				for (j = *m; j >= 2; --j) {
					ctemp = c__[j - 1];
					stemp = s[j - 1];
					if (ctemp != 1. || stemp != 0.) {
						i__1 = *n;
						for (i__ = 1; i__ <= i__1; ++i__) {
							temp = a[j + i__ * a_dim1];
							a[j + i__ * a_dim1] = ctemp * temp - stemp * a[
								i__ * a_dim1 + 1];
								a[i__ * a_dim1 + 1] = stemp * temp + ctemp * a[
									i__ * a_dim1 + 1];
									/* L70: */
						}
					}
					/* L80: */
				}
			}
		} else if (lsame_(pivot, "B")) {
			if (lsame_(direct, "F")) {
				i__1 = *m - 1;
				for (j = 1; j <= i__1; ++j) {
					ctemp = c__[j];
					stemp = s[j];
					if (ctemp != 1. || stemp != 0.) {
						i__2 = *n;
						for (i__ = 1; i__ <= i__2; ++i__) {
							temp = a[j + i__ * a_dim1];
							a[j + i__ * a_dim1] = stemp * a[*m + i__ * a_dim1]
							+ ctemp * temp;
							a[*m + i__ * a_dim1] = ctemp * a[*m + i__ * 
								a_dim1] - stemp * temp;
							/* L90: */
						}
					}
					/* L100: */
				}
			} else if (lsame_(direct, "B")) {
				for (j = *m - 1; j >= 1; --j) {
					ctemp = c__[j];
					stemp = s[j];
					if (ctemp != 1. || stemp != 0.) {
						i__1 = *n;
						for (i__ = 1; i__ <= i__1; ++i__) {
							temp = a[j + i__ * a_dim1];
							a[j + i__ * a_dim1] = stemp * a[*m + i__ * a_dim1]
							+ ctemp * temp;
							a[*m + i__ * a_dim1] = ctemp * a[*m + i__ * 
								a_dim1] - stemp * temp;
							/* L110: */
						}
					}
					/* L120: */
				}
			}
		}
	} else if (lsame_(side, "R")) {

		/*        Form A * P' */

		if (lsame_(pivot, "V")) {
			if (lsame_(direct, "F")) {
				i__1 = *n - 1;
				for (j = 1; j <= i__1; ++j) {
					ctemp = c__[j];
					stemp = s[j];
					if (ctemp != 1. || stemp != 0.) {
						i__2 = *m;
						for (i__ = 1; i__ <= i__2; ++i__) {
							temp = a[i__ + (j + 1) * a_dim1];
							a[i__ + (j + 1) * a_dim1] = ctemp * temp - stemp *
								a[i__ + j * a_dim1];
							a[i__ + j * a_dim1] = stemp * temp + ctemp * a[
								i__ + j * a_dim1];
								/* L130: */
						}
					}
					/* L140: */
				}
			} else if (lsame_(direct, "B")) {
				for (j = *n - 1; j >= 1; --j) {
					ctemp = c__[j];
					stemp = s[j];
					if (ctemp != 1. || stemp != 0.) {
						i__1 = *m;
						for (i__ = 1; i__ <= i__1; ++i__) {
							temp = a[i__ + (j + 1) * a_dim1];
							a[i__ + (j + 1) * a_dim1] = ctemp * temp - stemp *
								a[i__ + j * a_dim1];
							a[i__ + j * a_dim1] = stemp * temp + ctemp * a[
								i__ + j * a_dim1];
								/* L150: */
						}
					}
					/* L160: */
				}
			}
		} else if (lsame_(pivot, "T")) {
			if (lsame_(direct, "F")) {
				i__1 = *n;
				for (j = 2; j <= i__1; ++j) {
					ctemp = c__[j - 1];
					stemp = s[j - 1];
					if (ctemp != 1. || stemp != 0.) {
						i__2 = *m;
						for (i__ = 1; i__ <= i__2; ++i__) {
							temp = a[i__ + j * a_dim1];
							a[i__ + j * a_dim1] = ctemp * temp - stemp * a[
								i__ + a_dim1];
								a[i__ + a_dim1] = stemp * temp + ctemp * a[i__ + 
									a_dim1];
								/* L170: */
						}
					}
					/* L180: */
				}
			} else if (lsame_(direct, "B")) {
				for (j = *n; j >= 2; --j) {
					ctemp = c__[j - 1];
					stemp = s[j - 1];
					if (ctemp != 1. || stemp != 0.) {
						i__1 = *m;
						for (i__ = 1; i__ <= i__1; ++i__) {
							temp = a[i__ + j * a_dim1];
							a[i__ + j * a_dim1] = ctemp * temp - stemp * a[
								i__ + a_dim1];
								a[i__ + a_dim1] = stemp * temp + ctemp * a[i__ + 
									a_dim1];
								/* L190: */
						}
					}
					/* L200: */
				}
			}
		} else if (lsame_(pivot, "B")) {
			if (lsame_(direct, "F")) {
				i__1 = *n - 1;
				for (j = 1; j <= i__1; ++j) {
					ctemp = c__[j];
					stemp = s[j];
					if (ctemp != 1. || stemp != 0.) {
						i__2 = *m;
						for (i__ = 1; i__ <= i__2; ++i__) {
							temp = a[i__ + j * a_dim1];
							a[i__ + j * a_dim1] = stemp * a[i__ + *n * a_dim1]
							+ ctemp * temp;
							a[i__ + *n * a_dim1] = ctemp * a[i__ + *n * 
								a_dim1] - stemp * temp;
							/* L210: */
						}
					}
					/* L220: */
				}
			} else if (lsame_(direct, "B")) {
				for (j = *n - 1; j >= 1; --j) {
					ctemp = c__[j];
					stemp = s[j];
					if (ctemp != 1. || stemp != 0.) {
						i__1 = *m;
						for (i__ = 1; i__ <= i__1; ++i__) {
							temp = a[i__ + j * a_dim1];
							a[i__ + j * a_dim1] = stemp * a[i__ + *n * a_dim1]
							+ ctemp * temp;
							a[i__ + *n * a_dim1] = ctemp * a[i__ + *n * 
								a_dim1] - stemp * temp;
							/* L230: */
						}
					}
					/* L240: */
				}
			}
		}
	}

	return 0;

	/*     End of DLASR */

} /* dlasr_ */

//
int TRL::dswap_(integer_ *n, doublereal_ *dx, integer_ *incx, 
		   doublereal_ *dy, integer_ *incy)
{
	/* System generated locals */
	integer_ i__1;
	/* Local variables */
	static integer_ i__, m, ix, iy, mp1;
	static doublereal_ dtemp;
	/*  Purpose   
	=======   
	interchanges two vectors.   
	uses unrolled loops for increments equal one.   
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
	/*       code for unequal increments or equal increments not equal   
	to 1 */
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
		dtemp = dx[ix];
		dx[ix] = dy[iy];
		dy[iy] = dtemp;
		ix += *incx;
		iy += *incy;
		/* L10: */
	}
	return 0;
	/*       code for both increments equal to 1   
	clean-up loop */
L20:
	m = *n % 3;
	if (m == 0) {
		goto L40;
	}
	i__1 = m;
	for (i__ = 1; i__ <= i__1; ++i__) {
		dtemp = dx[i__];
		dx[i__] = dy[i__];
		dy[i__] = dtemp;
		/* L30: */
	}
	if (*n < 3) {
		return 0;
	}
L40:
	mp1 = m + 1;
	i__1 = *n;
	for (i__ = mp1; i__ <= i__1; i__ += 3) {
		dtemp = dx[i__];
		dx[i__] = dy[i__];
		dy[i__] = dtemp;
		dtemp = dx[i__ + 1];
		dx[i__ + 1] = dy[i__ + 1];
		dy[i__ + 1] = dtemp;
		dtemp = dx[i__ + 2];
		dx[i__ + 2] = dy[i__ + 2];
		dy[i__ + 2] = dtemp;
		/* L50: */
	}
	return 0;
} /* dswap_ */


//
int TRL::dlaev2_(doublereal_ *a, doublereal_ *b, doublereal_ *c__, 
			doublereal_ *rt1, doublereal_ *rt2, doublereal_ *cs1, doublereal_ *sn1)
{
	/*  -- LAPACK auxiliary routine (version 3.1) --   
	Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..   
	November 2006   


	Purpose   
	=======   

	DLAEV2 computes the eigendecomposition of a 2-by-2 symmetric matrix   
	[  A   B  ]   
	[  B   C  ].   
	On return, RT1 is the eigenvalue of larger absolute value, RT2 is the   
	eigenvalue of smaller absolute value, and (CS1,SN1) is the unit right   
	eigenvector for RT1, giving the decomposition   

	[ CS1  SN1 ] [  A   B  ] [ CS1 -SN1 ]  =  [ RT1  0  ]   
	[-SN1  CS1 ] [  B   C  ] [ SN1  CS1 ]     [  0  RT2 ].   

	Arguments   
	=========   

	A       (input) DOUBLE PRECISION   
	The (1,1) element of the 2-by-2 matrix.   

	B       (input) DOUBLE PRECISION   
	The (1,2) element and the conjugate of the (2,1) element of   
	the 2-by-2 matrix.   

	C       (input) DOUBLE PRECISION   
	The (2,2) element of the 2-by-2 matrix.   

	RT1     (output) DOUBLE PRECISION   
	The eigenvalue of larger absolute value.   

	RT2     (output) DOUBLE PRECISION   
	The eigenvalue of smaller absolute value.   

	CS1     (output) DOUBLE PRECISION   
	SN1     (output) DOUBLE PRECISION   
	The vector (CS1, SN1) is a unit right eigenvector for RT1.   

	Further Details   
	===============   

	RT1 is accurate to a few ulps barring over/underflow.   

	RT2 may be inaccurate if there is massive cancellation in the   
	determinant A*C-B*B; higher precision or correctly rounded or   
	correctly truncated arithmetic would be needed to compute RT2   
	accurately in all cases.   

	CS1 and SN1 are accurate to a few ulps barring over/underflow.   

	Overflow is possible only if RT1 is within a factor of 5 of overflow.   
	Underflow is harmless if the input data is 0 or exceeds   
	underflow_threshold / macheps.   

	=====================================================================   


	Compute the eigenvalues */
	/* System generated locals */
	doublereal_ d__1;
	/* Builtin functions */
	//double sqrt_(doublereal_);
	/* Local variables */
	static doublereal_ ab, df, cs, ct, tb, sm, tn, rt, adf, acs;
	static integer_ sgn1, sgn2;
	static doublereal_ acmn, acmx;


	sm = *a + *c__;
	df = *a - *c__;
	adf = abs(df);
	tb = *b + *b;
	ab = abs(tb);
	if (abs(*a) > abs(*c__)) {
		acmx = *a;
		acmn = *c__;
	} else {
		acmx = *c__;
		acmn = *a;
	}
	if (adf > ab) {
		/* Computing 2nd power */
		d__1 = ab / adf;
		rt = adf * sqrt_(d__1 * d__1 + 1.);
	} else if (adf < ab) {
		/* Computing 2nd power */
		d__1 = adf / ab;
		rt = ab * sqrt_(d__1 * d__1 + 1.);
	} else {

		/*        Includes case AB=ADF=0 */

		rt = ab * sqrt_(2.);
	}
	if (sm < 0.) {
		*rt1 = (sm - rt) * .5;
		sgn1 = -1;

		/*        Order of execution important.   
		To get fully accurate smaller eigenvalue,   
		next line needs to be executed in higher precision. */

		*rt2 = acmx / *rt1 * acmn - *b / *rt1 * *b;
	} else if (sm > 0.) {
		*rt1 = (sm + rt) * .5;
		sgn1 = 1;

		/*        Order of execution important.   
		To get fully accurate smaller eigenvalue,   
		next line needs to be executed in higher precision. */

		*rt2 = acmx / *rt1 * acmn - *b / *rt1 * *b;
	} else {

		/*        Includes case RT1 = RT2 = 0 */

		*rt1 = rt * .5;
		*rt2 = rt * -.5;
		sgn1 = 1;
	}

	/*     Compute the eigenvector */

	if (df >= 0.) {
		cs = df + rt;
		sgn2 = 1;
	} else {
		cs = df - rt;
		sgn2 = -1;
	}
	acs = abs(cs);
	if (acs > ab) {
		ct = -tb / cs;
		*sn1 = 1. / sqrt_(ct * ct + 1.);
		*cs1 = ct * *sn1;
	} else {
		if (ab == 0.) {
			*cs1 = 1.;
			*sn1 = 0.;
		} else {
			tn = -cs / tb;
			*cs1 = 1. / sqrt_(tn * tn + 1.);
			*sn1 = tn * *cs1;
		}
	}
	if (sgn1 == sgn2) {
		tn = *cs1;
		*cs1 = -(*sn1);
		*sn1 = tn;
	}
	return 0;

	/*     End of DLAEV2 */

} /* dlaev2_ */


//

int TRL::dlaset_(char *uplo, integer_ *m, integer_ *n, doublereal_ *
			alpha, doublereal_ *beta, doublereal_ *a, integer_ *lda)
{
	/*  -- LAPACK auxiliary routine (version 3.1) --   
	Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..   
	November 2006   


	Purpose   
	=======   

	DLASET initializes an m-by-n matrix A to BETA on the diagonal and   
	ALPHA on the offdiagonals.   

	Arguments   
	=========   

	UPLO    (input) CHARACTER*1   
	Specifies the part of the matrix A to be set.   
	= 'U':      Upper triangular part is set; the strictly lower   
	triangular part of A is not changed.   
	= 'L':      Lower triangular part is set; the strictly upper   
	triangular part of A is not changed.   
	Otherwise:  All of the matrix A is set.   

	M       (input) INTEGER   
	The number of rows of the matrix A.  M >= 0.   

	N       (input) INTEGER   
	The number of columns of the matrix A.  N >= 0.   

	ALPHA   (input) DOUBLE PRECISION   
	The constant to which the offdiagonal elements are to be set.   

	BETA    (input) DOUBLE PRECISION   
	The constant to which the diagonal elements are to be set.   

	A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
	On exit, the leading m-by-n submatrix of A is set as follows:   

	if UPLO = 'U', A(i,j) = ALPHA, 1<=i<=j-1, 1<=j<=n,   
	if UPLO = 'L', A(i,j) = ALPHA, j+1<=i<=m, 1<=j<=n,   
	otherwise,     A(i,j) = ALPHA, 1<=i<=m, 1<=j<=n, i.ne.j,   

	and, for all UPLO, A(i,i) = BETA, 1<=i<=min(m,n).   

	LDA     (input) INTEGER   
	The leading dimension of the array A.  LDA >= max(1,M).   

	=====================================================================   


	Parameter adjustments */
	/* System generated locals */
	integer_ a_dim1, a_offset, i__1, i__2, i__3;
	/* Local variables */
	static integer_ i__, j;
	//extern logical_ lsame_(char *, char *);

	a_dim1 = *lda;
	a_offset = 1 + a_dim1;
	a -= a_offset;

	/* Function Body */
	if (lsame_(uplo, "U")) {

		/*        Set the strictly upper triangular or trapezoidal part of the   
		array to ALPHA. */

		i__1 = *n;
		for (j = 2; j <= i__1; ++j) {
			/* Computing MIN */
			i__3 = j - 1;
			i__2 = min(i__3,*m);
			for (i__ = 1; i__ <= i__2; ++i__) {
				a[i__ + j * a_dim1] = *alpha;
				/* L10: */
			}
			/* L20: */
		}

	} else if (lsame_(uplo, "L")) {

		/*        Set the strictly lower triangular or trapezoidal part of the   
		array to ALPHA. */

		i__1 = min(*m,*n);
		for (j = 1; j <= i__1; ++j) {
			i__2 = *m;
			for (i__ = j + 1; i__ <= i__2; ++i__) {
				a[i__ + j * a_dim1] = *alpha;
				/* L30: */
			}
			/* L40: */
		}

	} else {

		/*        Set the leading m-by-n submatrix to ALPHA. */

		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
			i__2 = *m;
			for (i__ = 1; i__ <= i__2; ++i__) {
				a[i__ + j * a_dim1] = *alpha;
				/* L50: */
			}
			/* L60: */
		}
	}

	/*     Set the first min(M,N) diagonal elements to BETA. */

	i__1 = min(*m,*n);
	for (i__ = 1; i__ <= i__1; ++i__) {
		a[i__ + i__ * a_dim1] = *beta;
		/* L70: */
	}

	return 0;

	/*     End of DLASET */

} /* dlaset_ */

//

int TRL::dlartg_(doublereal_ *f, doublereal_ *g, doublereal_ *cs, 
			doublereal_ *sn, doublereal_ *r__)
{
	/*  -- LAPACK auxiliary routine (version 3.1) --   
	Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..   
	November 2006   


	Purpose   
	=======   

	DLARTG generate a plane rotation so that   

	[  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1.   
	[ -SN  CS  ]     [ G ]     [ 0 ]   

	This is a slower, more accurate version of the BLAS1 routine DROTG,   
	with the following other differences:   
	F and G are unchanged on return.   
	If G=0, then CS=1 and SN=0.   
	If F=0 and (G .ne. 0), then CS=0 and SN=1 without doing any   
	floating point operations (saves work in DBDSQR when   
	there are zeros on the diagonal).   

	If F exceeds G in magnitude, CS will be positive.   

	Arguments   
	=========   

	F       (input) DOUBLE PRECISION   
	The first component of vector to be rotated.   

	G       (input) DOUBLE PRECISION   
	The second component of vector to be rotated.   

	CS      (output) DOUBLE PRECISION   
	The cosine of the rotation.   

	SN      (output) DOUBLE PRECISION   
	The sine of the rotation.   

	R       (output) DOUBLE PRECISION   
	The nonzero component of the rotated vector.   

	This version has a few statements commented out for thread safety   
	(machine parameters are computed on each entry). 10 feb 03, SJH.   

	=====================================================================   

	LOGICAL            FIRST   
	SAVE               FIRST, SAFMX2, SAFMIN, SAFMN2   
	DATA               FIRST / .TRUE. /   

	IF( FIRST ) THEN */
	/* System generated locals */
	integer_ i__1;
	doublereal_ d__1, d__2;
	/* Builtin functions */
	//double log_(doublereal_), pow_di(doublereal_ *, integer_ *), sqrt_(doublereal_);
	/* Local variables */
	static integer_ i__;
	static doublereal_ f1, g1, eps, scale;
	static integer_ count;
	static doublereal_ safmn2, safmx2;
	//extern doublereal_ dlamch_(char *);
	static doublereal_ safmin;

	safmin = dlamch_("S");
	eps = dlamch_("E");
	d__1 = dlamch_("B");
	i__1 = (integer_) (log_(safmin / eps) / log_(dlamch_("B")) / 2.);
	safmn2 = pow_di(&d__1, &i__1);
	safmx2 = 1. / safmn2;
	/*        FIRST = .FALSE.   
	END IF */
	if (*g == 0.) {
		*cs = 1.;
		*sn = 0.;
		*r__ = *f;
	} else if (*f == 0.) {
		*cs = 0.;
		*sn = 1.;
		*r__ = *g;
	} else {
		f1 = *f;
		g1 = *g;
		/* Computing MAX */
		d__1 = abs(f1), d__2 = abs(g1);
		scale = max(d__1,d__2);
		if (scale >= safmx2) {
			count = 0;
L10:
			++count;
			f1 *= safmn2;
			g1 *= safmn2;
			/* Computing MAX */
			d__1 = abs(f1), d__2 = abs(g1);
			scale = max(d__1,d__2);
			if (scale >= safmx2) {
				goto L10;
			}
			/* Computing 2nd power */
			d__1 = f1;
			/* Computing 2nd power */
			d__2 = g1;
			*r__ = sqrt_(d__1 * d__1 + d__2 * d__2);
			*cs = f1 / *r__;
			*sn = g1 / *r__;
			i__1 = count;
			for (i__ = 1; i__ <= i__1; ++i__) {
				*r__ *= safmx2;
				/* L20: */
			}
		} else if (scale <= safmn2) {
			count = 0;
L30:
			++count;
			f1 *= safmx2;
			g1 *= safmx2;
			/* Computing MAX */
			d__1 = abs(f1), d__2 = abs(g1);
			scale = max(d__1,d__2);
			if (scale <= safmn2) {
				goto L30;
			}
			/* Computing 2nd power */
			d__1 = f1;
			/* Computing 2nd power */
			d__2 = g1;
			*r__ = sqrt_(d__1 * d__1 + d__2 * d__2);
			*cs = f1 / *r__;
			*sn = g1 / *r__;
			i__1 = count;
			for (i__ = 1; i__ <= i__1; ++i__) {
				*r__ *= safmn2;
				/* L40: */
			}
		} else {
			/* Computing 2nd power */
			d__1 = f1;
			/* Computing 2nd power */
			d__2 = g1;
			*r__ = sqrt_(d__1 * d__1 + d__2 * d__2);
			*cs = f1 / *r__;
			*sn = g1 / *r__;
		}
		if (abs(*f) > abs(*g) && *cs < 0.) {
			*cs = -(*cs);
			*sn = -(*sn);
			*r__ = -(*r__);
		}
	}
	return 0;

	/*     End of DLARTG */

} /* dlartg_ */


//////////////////////////////////////////////////////////////////////////
doublereal_ TRL::dnrm2_(integer_ *n, doublereal_ *x, integer_ *incx)
{
	/*        The following loop is equivalent to this call to the LAPACK   
	auxiliary routine:   
	CALL DLASSQ( N, X, INCX, SCALE, SSQ ) */
	/* System generated locals */
	integer_ i__1, i__2;
	doublereal_ ret_val, d__1;
	/* Builtin functions */
	//double sqrt_(doublereal_);
	/* Local variables */
	static integer_ ix;
	static doublereal_ ssq, norm, scale, absxi;
	/*  Purpose   
	=======   
	DNRM2 returns the euclidean norm of a vector via the function   
	name, so that   
	DNRM2 := sqrt_( x'*x )   
	-- This version written on 25-October-1982.   
	Modified on 14-October-1993 to inline the call to DLASSQ.   
	Sven Hammarling, Nag Ltd.   
	Parameter adjustments */
	--x;
	/* Function Body */
	if (*n < 1 || *incx < 1) {
		norm = 0.;
	} else if (*n == 1) {
		norm = abs(x[1]);
	} else {
		scale = 0.;
		ssq = 1.;


		i__1 = (*n - 1) * *incx + 1;
		i__2 = *incx;
		for (ix = 1; i__2 < 0 ? ix >= i__1 : ix <= i__1; ix += i__2) {
			if (x[ix] != 0.) {
				absxi = (d__1 = x[ix], abs(d__1));
				if (scale < absxi) {
					/* Computing 2nd power */
					d__1 = scale / absxi;
					ssq = ssq * (d__1 * d__1) + 1.;
					scale = absxi;
				} else {
					/* Computing 2nd power */
					d__1 = absxi / scale;
					ssq += d__1 * d__1;
				}
			}
			/* L10: */
		}
		norm = scale * sqrt_(ssq);
	}

	ret_val = norm;
	return ret_val;

	/*     End of DNRM2. */

} /* dnrm2_ */


//
doublereal_ TRL::dasum_(integer_ *n, doublereal_ *dx, integer_ *incx)
{
	/* System generated locals */
	integer_ i__1, i__2;
	doublereal_ ret_val, d__1, d__2, d__3, d__4, d__5, d__6;
	/* Local variables */
	static integer_ i__, m, mp1;
	static doublereal_ dtemp;
	static integer_ nincx;
	/*  Purpose   
	=======   
	takes the sum of the absolute values.   
	jack dongarra, linpack, 3/11/78.   
	modified 3/93 to return if incx .le. 0.   
	modified 12/3/93, array(1) declarations changed to array(*)   
	Parameter adjustments */
	--dx;
	/* Function Body */
	ret_val = 0.;
	dtemp = 0.;
	if (*n <= 0 || *incx <= 0) {
		return ret_val;
	}
	if (*incx == 1) {
		goto L20;
	}
	/*        code for increment not equal to 1 */
	nincx = *n * *incx;
	i__1 = nincx;
	i__2 = *incx;
	for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
		dtemp += (d__1 = dx[i__], abs(d__1));
		/* L10: */
	}
	ret_val = dtemp;
	return ret_val;
	/*        code for increment equal to 1   
	clean-up loop */
L20:
	m = *n % 6;
	if (m == 0) {
		goto L40;
	}
	i__2 = m;
	for (i__ = 1; i__ <= i__2; ++i__) {
		dtemp += (d__1 = dx[i__], abs(d__1));
		/* L30: */
	}
	if (*n < 6) {
		goto L60;
	}
L40:
	mp1 = m + 1;
	i__2 = *n;
	for (i__ = mp1; i__ <= i__2; i__ += 6) {
		dtemp = dtemp + (d__1 = dx[i__], abs(d__1)) + (d__2 = dx[i__ + 1], 
			abs(d__2)) + (d__3 = dx[i__ + 2], abs(d__3)) + (d__4 = dx[i__ 
			+ 3], abs(d__4)) + (d__5 = dx[i__ + 4], abs(d__5)) + (d__6 = 
			dx[i__ + 5], abs(d__6));
		/* L50: */
	}
L60:
	ret_val = dtemp;
	return ret_val;
} /* dasum_ */
//

//

int TRL::dlagtf_(integer_ *n, doublereal_ *a, doublereal_ *lambda, 
			doublereal_ *b, doublereal_ *c__, doublereal_ *tol, doublereal_ *d__, 
			integer_ *in, integer_ *info)
{
	/*  -- LAPACK routine (version 3.1) --   
	Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..   
	November 2006   


	Purpose   
	=======   

	DLAGTF factorizes the matrix (T - lambda*I), where T is an n by n   
	tridiagonal matrix and lambda is a scalar, as   

	T - lambda*I = PLU,   

	where P is a permutation matrix, L is a unit lower tridiagonal matrix   
	with at most one non-zero sub-diagonal elements per column and U is   
	an upper triangular matrix with at most two non-zero super-diagonal   
	elements per column.   

	The factorization is obtained by Gaussian elimination with partial   
	pivoting and implicit row scaling.   

	The parameter LAMBDA is included in the routine so that DLAGTF may   
	be used, in conjunction with DLAGTS, to obtain eigenvectors of T by   
	inverse iteration.   

	Arguments   
	=========   

	N       (input) INTEGER   
	The order of the matrix T.   

	A       (input/output) DOUBLE PRECISION array, dimension (N)   
	On entry, A must contain the diagonal elements of T.   

	On exit, A is overwritten by the n diagonal elements of the   
	upper triangular matrix U of the factorization of T.   

	LAMBDA  (input) DOUBLE PRECISION   
	On entry, the scalar lambda.   

	B       (input/output) DOUBLE PRECISION array, dimension (N-1)   
	On entry, B must contain the (n-1) super-diagonal elements of   
	T.   

	On exit, B is overwritten by the (n-1) super-diagonal   
	elements of the matrix U of the factorization of T.   

	C       (input/output) DOUBLE PRECISION array, dimension (N-1)   
	On entry, C must contain the (n-1) sub-diagonal elements of   
	T.   

	On exit, C is overwritten by the (n-1) sub-diagonal elements   
	of the matrix L of the factorization of T.   

	TOL     (input) DOUBLE PRECISION   
	On entry, a relative tolerance used to indicate whether or   
	not the matrix (T - lambda*I) is nearly singular. TOL should   
	normally be chose as approximately the largest relative error   
	in the elements of T. For example, if the elements of T are   
	correct to about 4 significant figures, then TOL should be   
	set to about 5*10**(-4). If TOL is supplied as less than eps,   
	where eps is the relative machine precision, then the value   
	eps is used in place of TOL.   

	D       (output) DOUBLE PRECISION array, dimension (N-2)   
	On exit, D is overwritten by the (n-2) second super-diagonal   
	elements of the matrix U of the factorization of T.   

	IN      (output) INTEGER array, dimension (N)   
	On exit, IN contains details of the permutation matrix P. If   
	an interchange occurred at the kth step of the elimination,   
	then IN(k) = 1, otherwise IN(k) = 0. The element IN(n)   
	returns the smallest positive integer_ j such that   

	abs( u(j,j) ).le. norm( (T - lambda*I)(j) )*TOL,   

	where norm( A(j) ) denotes the sum of the absolute values of   
	the jth row of the matrix A. If no such j exists then IN(n)   
	is returned as zero. If IN(n) is returned as positive, then a   
	diagonal element of U is small, indicating that   
	(T - lambda*I) is singular or nearly singular,   

	INFO    (output) INTEGER   
	= 0   : successful exit   
	.lt. 0: if INFO = -k, the kth argument had an illegal value   

	=====================================================================   


	Parameter adjustments */
	/* System generated locals */
	integer_ i__1;
	doublereal_ d__1, d__2;
	/* Local variables */
	static integer_ k;
	static doublereal_ tl, eps, piv1, piv2, temp, mult, scale1, scale2;
	//extern doublereal_ dlamch_(char *);
	//extern /* Subroutine */ int xerbla_(char *, integer_ *);

	--in;
	--d__;
	--c__;
	--b;
	--a;

	/* Function Body */
	*info = 0;
	if (*n < 0) {
		*info = -1;
		i__1 = -(*info);
		xerbla_("DLAGTF", &i__1);
		return 0;
	}

	if (*n == 0) {
		return 0;
	}

	a[1] -= *lambda;
	in[*n] = 0;
	if (*n == 1) {
		if (a[1] == 0.) {
			in[1] = 1;
		}
		return 0;
	}

	eps = dlamch_("Epsilon");

	tl = max(*tol,eps);
	scale1 = abs(a[1]) + abs(b[1]);
	i__1 = *n - 1;
	for (k = 1; k <= i__1; ++k) {
		a[k + 1] -= *lambda;
		scale2 = (d__1 = c__[k], abs(d__1)) + (d__2 = a[k + 1], abs(d__2));
		if (k < *n - 1) {
			scale2 += (d__1 = b[k + 1], abs(d__1));
		}
		if (a[k] == 0.) {
			piv1 = 0.;
		} else {
			piv1 = (d__1 = a[k], abs(d__1)) / scale1;
		}
		if (c__[k] == 0.) {
			in[k] = 0;
			piv2 = 0.;
			scale1 = scale2;
			if (k < *n - 1) {
				d__[k] = 0.;
			}
		} else {
			piv2 = (d__1 = c__[k], abs(d__1)) / scale2;
			if (piv2 <= piv1) {
				in[k] = 0;
				scale1 = scale2;
				c__[k] /= a[k];
				a[k + 1] -= c__[k] * b[k];
				if (k < *n - 1) {
					d__[k] = 0.;
				}
			} else {
				in[k] = 1;
				mult = a[k] / c__[k];
				a[k] = c__[k];
				temp = a[k + 1];
				a[k + 1] = b[k] - mult * temp;
				if (k < *n - 1) {
					d__[k] = b[k + 1];
					b[k + 1] = -mult * d__[k];
				}
				b[k] = temp;
				c__[k] = mult;
			}
		}
		if (max(piv1,piv2) <= tl && in[*n] == 0) {
			in[*n] = k;
		}
		/* L10: */
	}
	if ((d__1 = a[*n], abs(d__1)) <= scale1 * tl && in[*n] == 0) {
		in[*n] = *n;
	}

	return 0;

	/*     End of DLAGTF */

} /* dlagtf_ */


integer_ TRL::idamax_(integer_ *n, doublereal_ *dx, integer_ *incx)
{
	/* System generated locals */
	integer_ ret_val, i__1;
	doublereal_ d__1;
	/* Local variables */
	static integer_ i__, ix;
	static doublereal_ dmax__;
	/*  Purpose   
	=======   
	finds the index of element having max. absolute value.   
	jack dongarra, linpack, 3/11/78.   
	modified 3/93 to return if incx .le. 0.   
	modified 12/3/93, array(1) declarations changed to array(*)   
	Parameter adjustments */
	--dx;
	/* Function Body */
	ret_val = 0;
	if (*n < 1 || *incx <= 0) {
		return ret_val;
	}
	ret_val = 1;
	if (*n == 1) {
		return ret_val;
	}
	if (*incx == 1) {
		goto L20;
	}
	/*        code for increment not equal to 1 */
	ix = 1;
	dmax__ = abs(dx[1]);
	ix += *incx;
	i__1 = *n;
	for (i__ = 2; i__ <= i__1; ++i__) {
		if ((d__1 = dx[ix], abs(d__1)) <= dmax__) {
			goto L5;
		}
		ret_val = i__;
		dmax__ = (d__1 = dx[ix], abs(d__1));
L5:
		ix += *incx;
		/* L10: */
	}
	return ret_val;
	/*        code for increment equal to 1 */
L20:
	dmax__ = abs(dx[1]);
	i__1 = *n;
	for (i__ = 2; i__ <= i__1; ++i__) {
		if ((d__1 = dx[i__], abs(d__1)) <= dmax__) {
			goto L30;
		}
		ret_val = i__;
		dmax__ = (d__1 = dx[i__], abs(d__1));
L30:
		;
	}
	return ret_val;
} /* idamax_ */


int TRL::dlagts_(integer_ *job, integer_ *n, doublereal_ *a, 
			doublereal_ *b, doublereal_ *c__, doublereal_ *d__, integer_ *in, 
			doublereal_ *y, doublereal_ *tol, integer_ *info)
{
	/*  -- LAPACK auxiliary routine (version 3.1) --   
	Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..   
	November 2006   


	Purpose   
	=======   

	DLAGTS may be used to solve one of the systems of equations   

	(T - lambda*I)*x = y   or   (T - lambda*I)'*x = y,   

	where T is an n by n tridiagonal matrix, for x, following the   
	factorization of (T - lambda*I) as   

	(T - lambda*I) = P*L*U ,   

	by routine DLAGTF. The choice of equation to be solved is   
	controlled by the argument JOB, and in each case there is an option   
	to perturb zero or very small diagonal elements of U, this option   
	being intended for use in applications such as inverse iteration.   

	Arguments   
	=========   

	JOB     (input) INTEGER   
	Specifies the job to be performed by DLAGTS as follows:   
	=  1: The equations  (T - lambda*I)x = y  are to be solved,   
	but diagonal elements of U are not to be perturbed.   
	= -1: The equations  (T - lambda*I)x = y  are to be solved   
	and, if overflow would otherwise occur, the diagonal   
	elements of U are to be perturbed. See argument TOL   
	below.   
	=  2: The equations  (T - lambda*I)'x = y  are to be solved,   
	but diagonal elements of U are not to be perturbed.   
	= -2: The equations  (T - lambda*I)'x = y  are to be solved   
	and, if overflow would otherwise occur, the diagonal   
	elements of U are to be perturbed. See argument TOL   
	below.   

	N       (input) INTEGER   
	The order of the matrix T.   

	A       (input) DOUBLE PRECISION array, dimension (N)   
	On entry, A must contain the diagonal elements of U as   
	returned from DLAGTF.   

	B       (input) DOUBLE PRECISION array, dimension (N-1)   
	On entry, B must contain the first super-diagonal elements of   
	U as returned from DLAGTF.   

	C       (input) DOUBLE PRECISION array, dimension (N-1)   
	On entry, C must contain the sub-diagonal elements of L as   
	returned from DLAGTF.   

	D       (input) DOUBLE PRECISION array, dimension (N-2)   
	On entry, D must contain the second super-diagonal elements   
	of U as returned from DLAGTF.   

	IN      (input) INTEGER array, dimension (N)   
	On entry, IN must contain details of the matrix P as returned   
	from DLAGTF.   

	Y       (input/output) DOUBLE PRECISION array, dimension (N)   
	On entry, the right hand side vector y.   
	On exit, Y is overwritten by the solution vector x.   

	TOL     (input/output) DOUBLE PRECISION   
	On entry, with  JOB .lt. 0, TOL should be the minimum   
	perturbation to be made to very small diagonal elements of U.   
	TOL should normally be chosen as about eps*norm(U), where eps   
	is the relative machine precision, but if TOL is supplied as   
	non-positive, then it is reset to eps*max( abs( u(i,j) ) ).   
	If  JOB .gt. 0  then TOL is not referenced.   

	On exit, TOL is changed as described above, only if TOL is   
	non-positive on entry. Otherwise TOL is unchanged.   

	INFO    (output) INTEGER   
	= 0   : successful exit   
	.lt. 0: if INFO = -i, the i-th argument had an illegal value   
	.gt. 0: overflow would occur when computing the INFO(th)   
	element of the solution vector x. This can only occur   
	when JOB is supplied as positive and either means   
	that a diagonal element of U is very small, or that   
	the elements of the right-hand side vector y are very   
	large.   

	=====================================================================   


	Parameter adjustments */
	/* System generated locals */
	integer_ i__1;
	doublereal_ d__1, d__2, d__3, d__4, d__5;
	/* Builtin functions */
	//double d_sign(doublereal_ *, doublereal_ *);
	/* Local variables */
	static integer_ k;
	static doublereal_ ak, eps, temp, pert, absak, sfmin;
	//extern doublereal_ dlamch_(char *);
	//extern /* Subroutine */ int xerbla_(char *, integer_ *);
	static doublereal_ bignum;

	--y;
	--in;
	--d__;
	--c__;
	--b;
	--a;

	/* Function Body */
	*info = 0;
	if (abs(*job) > 2 || *job == 0) {
		*info = -1;
	} else if (*n < 0) {
		*info = -2;
	}
	if (*info != 0) {
		i__1 = -(*info);
		xerbla_("DLAGTS", &i__1);
		return 0;
	}

	if (*n == 0) {
		return 0;
	}

	eps = dlamch_("Epsilon");
	sfmin = dlamch_("Safe minimum");
	bignum = 1. / sfmin;

	if (*job < 0) {
		if (*tol <= 0.) {
			*tol = abs(a[1]);
			if (*n > 1) {
				/* Computing MAX */
				d__1 = *tol, d__2 = abs(a[2]), d__1 = max(d__1,d__2), d__2 = 
					abs(b[1]);
				*tol = max(d__1,d__2);
			}
			i__1 = *n;
			for (k = 3; k <= i__1; ++k) {
				/* Computing MAX */
				d__4 = *tol, d__5 = (d__1 = a[k], abs(d__1)), d__4 = max(d__4,
					d__5), d__5 = (d__2 = b[k - 1], abs(d__2)), d__4 = 
					max(d__4,d__5), d__5 = (d__3 = d__[k - 2], abs(d__3));
				*tol = max(d__4,d__5);
				/* L10: */
			}
			*tol *= eps;
			if (*tol == 0.) {
				*tol = eps;
			}
		}
	}

	if (abs(*job) == 1) {
		i__1 = *n;
		for (k = 2; k <= i__1; ++k) {
			if (in[k - 1] == 0) {
				y[k] -= c__[k - 1] * y[k - 1];
			} else {
				temp = y[k - 1];
				y[k - 1] = y[k];
				y[k] = temp - c__[k - 1] * y[k];
			}
			/* L20: */
		}
		if (*job == 1) {
			for (k = *n; k >= 1; --k) {
				if (k <= *n - 2) {
					temp = y[k] - b[k] * y[k + 1] - d__[k] * y[k + 2];
				} else if (k == *n - 1) {
					temp = y[k] - b[k] * y[k + 1];
				} else {
					temp = y[k];
				}
				ak = a[k];
				absak = abs(ak);
				if (absak < 1.) {
					if (absak < sfmin) {
						if (absak == 0. || abs(temp) * sfmin > absak) {
							*info = k;
							return 0;
						} else {
							temp *= bignum;
							ak *= bignum;
						}
					} else if (abs(temp) > absak * bignum) {
						*info = k;
						return 0;
					}
				}
				y[k] = temp / ak;
				/* L30: */
			}
		} else {
			for (k = *n; k >= 1; --k) {
				if (k <= *n - 2) {
					temp = y[k] - b[k] * y[k + 1] - d__[k] * y[k + 2];
				} else if (k == *n - 1) {
					temp = y[k] - b[k] * y[k + 1];
				} else {
					temp = y[k];
				}
				ak = a[k];
				pert = d_sign(tol, &ak);
L40:
				absak = abs(ak);
				if (absak < 1.) {
					if (absak < sfmin) {
						if (absak == 0. || abs(temp) * sfmin > absak) {
							ak += pert;
							pert *= 2;
							goto L40;
						} else {
							temp *= bignum;
							ak *= bignum;
						}
					} else if (abs(temp) > absak * bignum) {
						ak += pert;
						pert *= 2;
						goto L40;
					}
				}
				y[k] = temp / ak;
				/* L50: */
			}
		}
	} else {

		/*        Come to here if  JOB = 2 or -2 */

		if (*job == 2) {
			i__1 = *n;
			for (k = 1; k <= i__1; ++k) {
				if (k >= 3) {
					temp = y[k] - b[k - 1] * y[k - 1] - d__[k - 2] * y[k - 2];
				} else if (k == 2) {
					temp = y[k] - b[k - 1] * y[k - 1];
				} else {
					temp = y[k];
				}
				ak = a[k];
				absak = abs(ak);
				if (absak < 1.) {
					if (absak < sfmin) {
						if (absak == 0. || abs(temp) * sfmin > absak) {
							*info = k;
							return 0;
						} else {
							temp *= bignum;
							ak *= bignum;
						}
					} else if (abs(temp) > absak * bignum) {
						*info = k;
						return 0;
					}
				}
				y[k] = temp / ak;
				/* L60: */
			}
		} else {
			i__1 = *n;
			for (k = 1; k <= i__1; ++k) {
				if (k >= 3) {
					temp = y[k] - b[k - 1] * y[k - 1] - d__[k - 2] * y[k - 2];
				} else if (k == 2) {
					temp = y[k] - b[k - 1] * y[k - 1];
				} else {
					temp = y[k];
				}
				ak = a[k];
				pert = d_sign(tol, &ak);
L70:
				absak = abs(ak);
				if (absak < 1.) {
					if (absak < sfmin) {
						if (absak == 0. || abs(temp) * sfmin > absak) {
							ak += pert;
							pert *= 2;
							goto L70;
						} else {
							temp *= bignum;
							ak *= bignum;
						}
					} else if (abs(temp) > absak * bignum) {
						ak += pert;
						pert *= 2;
						goto L70;
					}
				}
				y[k] = temp / ak;
				/* L80: */
			}
		}

		for (k = *n; k >= 2; --k) {
			if (in[k - 1] == 0) {
				y[k - 1] -= c__[k - 1] * y[k];
			} else {
				temp = y[k - 1];
				y[k - 1] = y[k];
				y[k] = temp - c__[k - 1] * y[k];
			}
			/* L90: */
		}
	}

	/*     End of DLAGTS */

	return 0;
} /* dlagts_ */


int TRL::dlarnv_(integer_ *idist, integer_ *iseed, integer_ *n, 
			doublereal_ *x)
{
	/* System generated locals */
	integer_ i__1, i__2, i__3;

	/* Builtin functions */
	//double log_(doublereal_), sqrt_(doublereal_), cos_(doublereal_);

	/* Local variables */
	static integer_ i__;
	static doublereal_ u[128];
	static integer_ il, iv, il2;
	//extern /* Subroutine */ int dlaruv_(integer_ *, integer_ *, doublereal_ *);


	/*  -- LAPACK auxiliary routine (version 3.1) --   
	Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..   
	November 2006   


	Purpose   
	=======   

	DLARNV returns a vector of n random trl_real_ numbers from a uniform or   
	normal distribution.   

	Arguments   
	=========   

	IDIST   (input) INTEGER   
	Specifies the distribution of the random numbers:   
	= 1:  uniform (0,1)   
	= 2:  uniform (-1,1)   
	= 3:  normal (0,1)   

	ISEED   (input/output) INTEGER array, dimension (4)   
	On entry, the seed of the random number generator; the array   
	elements must be between 0 and 4095, and ISEED(4) must be   
	odd.   
	On exit, the seed is updated.   

	N       (input) INTEGER   
	The number of random numbers to be generated.   

	X       (output) DOUBLE PRECISION array, dimension (N)   
	The generated random numbers.   

	Further Details   
	===============   

	This routine calls the auxiliary routine DLARUV to generate random   
	real numbers from a uniform (0,1) distribution, in batches of up to   
	128 using vectorisable code. The Box-Muller method is used to   
	transform numbers from a uniform to a normal distribution.   

	=====================================================================   


	Parameter adjustments */
	--x;
	--iseed;

	/* Function Body */
	i__1 = *n;
	for (iv = 1; iv <= i__1; iv += 64) {
		/* Computing MIN */
		i__2 = 64, i__3 = *n - iv + 1;
		il = min(i__2,i__3);
		if (*idist == 3) {
			il2 = il << 1;
		} else {
			il2 = il;
		}

		/*        Call DLARUV to generate IL2 numbers from a uniform (0,1)   
		distribution (IL2 <= LV) */

		dlaruv_(&iseed[1], &il2, u);

		if (*idist == 1) {

			/*           Copy generated numbers */

			i__2 = il;
			for (i__ = 1; i__ <= i__2; ++i__) {
				x[iv + i__ - 1] = u[i__ - 1];
				/* L10: */
			}
		} else if (*idist == 2) {

			/*           Convert generated numbers to uniform (-1,1) distribution */

			i__2 = il;
			for (i__ = 1; i__ <= i__2; ++i__) {
				x[iv + i__ - 1] = u[i__ - 1] * 2. - 1.;
				/* L20: */
			}
		} else if (*idist == 3) {

			/*           Convert generated numbers to normal (0,1) distribution */

			i__2 = il;
			for (i__ = 1; i__ <= i__2; ++i__) {
				x[iv + i__ - 1] = sqrt_(log_(u[(i__ << 1) - 2]) * -2.) * cos_(u[(
					i__ << 1) - 1] * 6.2831853071795864769252867663);
				/* L30: */
			}
		}
		/* L40: */
	}
	return 0;

	/*     End of DLARNV */

} /* dlarnv_ */


int TRL::dlaruv_(integer_ *iseed, integer_ *n, doublereal_ *x)
{
	/* Initialized data */

	static integer_ mm[512]	/* was [128][4] */ = { 494,2637,255,2008,1253,
		3344,4084,1739,3143,3468,688,1657,1238,3166,1292,3422,1270,2016,
		154,2862,697,1706,491,931,1444,444,3577,3944,2184,1661,3482,657,
		3023,3618,1267,1828,164,3798,3087,2400,2870,3876,1905,1593,1797,
		1234,3460,328,2861,1950,617,2070,3331,769,1558,2412,2800,189,287,
		2045,1227,2838,209,2770,3654,3993,192,2253,3491,2889,2857,2094,
		1818,688,1407,634,3231,815,3524,1914,516,164,303,2144,3480,119,
		3357,837,2826,2332,2089,3780,1700,3712,150,2000,3375,1621,3090,
		3765,1149,3146,33,3082,2741,359,3316,1749,185,2784,2202,2199,1364,
		1244,2020,3160,2785,2772,1217,1822,1245,2252,3904,2774,997,2573,
		1148,545,322,789,1440,752,2859,123,1848,643,2405,2638,2344,46,
		3814,913,3649,339,3808,822,2832,3078,3633,2970,637,2249,2081,4019,
		1478,242,481,2075,4058,622,3376,812,234,641,4005,1122,3135,2640,
		2302,40,1832,2247,2034,2637,1287,1691,496,1597,2394,2584,1843,336,
		1472,2407,433,2096,1761,2810,566,442,41,1238,1086,603,840,3168,
		1499,1084,3438,2408,1589,2391,288,26,512,1456,171,1677,2657,2270,
		2587,2961,1970,1817,676,1410,3723,2803,3185,184,663,499,3784,1631,
		1925,3912,1398,1349,1441,2224,2411,1907,3192,2786,382,37,759,2948,
		1862,3802,2423,2051,2295,1332,1832,2405,3638,3661,327,3660,716,
		1842,3987,1368,1848,2366,2508,3754,1766,3572,2893,307,1297,3966,
		758,2598,3406,2922,1038,2934,2091,2451,1580,1958,2055,1507,1078,
		3273,17,854,2916,3971,2889,3831,2621,1541,893,736,3992,787,2125,
		2364,2460,257,1574,3912,1216,3248,3401,2124,2762,149,2245,166,466,
		4018,1399,190,2879,153,2320,18,712,2159,2318,2091,3443,1510,449,
		1956,2201,3137,3399,1321,2271,3667,2703,629,2365,2431,1113,3922,
		2554,184,2099,3228,4012,1921,3452,3901,572,3309,3171,817,3039,
		1696,1256,3715,2077,3019,1497,1101,717,51,981,1978,1813,3881,76,
		3846,3694,1682,124,1660,3997,479,1141,886,3514,1301,3604,1888,
		1836,1990,2058,692,1194,20,3285,2046,2107,3508,3525,3801,2549,
		1145,2253,305,3301,1065,3133,2913,3285,1241,1197,3729,2501,1673,
		541,2753,949,2361,1165,4081,2725,3305,3069,3617,3733,409,2157,
		1361,3973,1865,2525,1409,3445,3577,77,3761,2149,1449,3005,225,85,
		3673,3117,3089,1349,2057,413,65,1845,697,3085,3441,1573,3689,2941,
		929,533,2841,4077,721,2821,2249,2397,2817,245,1913,1997,3121,997,
		1833,2877,1633,981,2009,941,2449,197,2441,285,1473,2741,3129,909,
		2801,421,4073,2813,2337,1429,1177,1901,81,1669,2633,2269,129,1141,
		249,3917,2481,3941,2217,2749,3041,1877,345,2861,1809,3141,2825,
		157,2881,3637,1465,2829,2161,3365,361,2685,3745,2325,3609,3821,
		3537,517,3017,2141,1537 };

	/* System generated locals */
	integer_ i__1;

	/* Local variables */
	static integer_ i__, i1, i2, i3, i4, it1, it2, it3, it4;


	/*  -- LAPACK auxiliary routine (version 3.1) --   
	Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..   
	November 2006   


	Purpose   
	=======   

	DLARUV returns a vector of n random trl_real_ numbers from a uniform (0,1)   
	distribution (n <= 128).   

	This is an auxiliary routine called by DLARNV and ZLARNV.   

	Arguments   
	=========   

	ISEED   (input/output) INTEGER array, dimension (4)   
	On entry, the seed of the random number generator; the array   
	elements must be between 0 and 4095, and ISEED(4) must be   
	odd.   
	On exit, the seed is updated.   

	N       (input) INTEGER   
	The number of random numbers to be generated. N <= 128.   

	X       (output) DOUBLE PRECISION array, dimension (N)   
	The generated random numbers.   

	Further Details   
	===============   

	This routine uses a multiplicative congruential method with modulus   
	2**48 and multiplier 33952834046453 (see G.S.Fishman,   
	'Multiplicative congruential random number generators with modulus   
	2**b: an exhaustive analysis for b = 32 and a partial analysis for   
	b = 48', Math. Comp. 189, pp 331-344, 1990).   

	48-bit integer_s are stored in 4 integer_ array elements with 12 bits   
	per element. Hence the routine is portable across machines with   
	integer_s of 32 bits or more.   

	=====================================================================   

	Parameter adjustments */
	--iseed;
	--x;

	/* Function Body */

	i1 = iseed[1];
	i2 = iseed[2];
	i3 = iseed[3];
	i4 = iseed[4];

	i__1 = min(*n,128);
	for (i__ = 1; i__ <= i__1; ++i__) {

L20:

		/*        Multiply the seed by i-th power of the multiplier modulo 2**48 */

		it4 = i4 * mm[i__ + 383];
		it3 = it4 / 4096;
		it4 -= it3 << 12;
		it3 = it3 + i3 * mm[i__ + 383] + i4 * mm[i__ + 255];
		it2 = it3 / 4096;
		it3 -= it2 << 12;
		it2 = it2 + i2 * mm[i__ + 383] + i3 * mm[i__ + 255] + i4 * mm[i__ + 
			127];
		it1 = it2 / 4096;
		it2 -= it1 << 12;
		it1 = it1 + i1 * mm[i__ + 383] + i2 * mm[i__ + 255] + i3 * mm[i__ + 
			127] + i4 * mm[i__ - 1];
		it1 %= 4096;

		/*        Convert 48-bit integer_ to a trl_real_ number in the interval (0,1) */

		x[i__] = ((doublereal_) it1 + ((doublereal_) it2 + ((doublereal_) it3 + (
			doublereal_) it4 * 2.44140625e-4) * 2.44140625e-4) * 
			2.44140625e-4) * 2.44140625e-4;

		if (x[i__] == 1.) {
			/*           If a trl_real_ number has n bits of precision, and the first   
			n bits of the 48-bit integer_ above happen to be all 1 (which   
			will occur about once every 2**n calls), then X( I ) will   
			be rounded to exactly 1.0.   
			Since X( I ) is not supposed to return exactly 0.0 or 1.0,   
			the statistically correct thing to do in this situation is   
			simply to iterate again.   
			N.B. the case X( I ) = 0.0 should not be possible. */
			i1 += 2;
			i2 += 2;
			i3 += 2;
			i4 += 2;
			goto L20;
		}

		/* L10: */
	}

	/*     Return final value of seed */

	iseed[1] = it1;
	iseed[2] = it2;
	iseed[3] = it3;
	iseed[4] = it4;
	return 0;

	/*     End of DLARUV */

} /* dlaruv_ */



//////////////////////////////////////////////////
////////////////////////////////
/////////////////////////////////////////

int TRL::dsyev_(char *jobz, char *uplo, integer_ *n, doublereal_ *a,
		   integer_ *lda, doublereal_ *w, doublereal_ *work, integer_ *lwork, 
		   integer_ *info)
{
	/*  -- LAPACK driver routine (version 3.1) --   
	Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..   
	November 2006   


	Purpose   
	=======   

	DSYEV computes all eigenvalues and, optionally, eigenvectors of a   
	real symmetric matrix A.   

	Arguments   
	=========   

	JOBZ    (input) CHARACTER*1   
	= 'N':  Compute eigenvalues only;   
	= 'V':  Compute eigenvalues and eigenvectors.   

	UPLO    (input) CHARACTER*1   
	= 'U':  Upper triangle of A is stored;   
	= 'L':  Lower triangle of A is stored.   

	N       (input) INTEGER   
	The order of the matrix A.  N >= 0.   

	A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)   
	On entry, the symmetric matrix A.  If UPLO = 'U', the   
	leading N-by-N upper triangular part of A contains the   
	upper triangular part of the matrix A.  If UPLO = 'L',   
	the leading N-by-N lower triangular part of A contains   
	the lower triangular part of the matrix A.   
	On exit, if JOBZ = 'V', then if INFO = 0, A contains the   
	orthonormal eigenvectors of the matrix A.   
	If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')   
	or the upper triangle (if UPLO='U') of A, including the   
	diagonal, is destroyed.   

	LDA     (input) INTEGER   
	The leading dimension of the array A.  LDA >= max(1,N).   

	W       (output) DOUBLE PRECISION array, dimension (N)   
	If INFO = 0, the eigenvalues in ascending order.   

	WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))   
	On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

	LWORK   (input) INTEGER   
	The length of the array WORK.  LWORK >= max(1,3*N-1).   
	For optimal efficiency, LWORK >= (NB+2)*N,   
	where NB is the blocksize for DSYTRD returned by ILAENV.   

	If LWORK = -1, then a workspace query is assumed; the routine   
	only calculates the optimal size of the WORK array, returns   
	this value as the first entry of the WORK array, and no error   
	message related to LWORK is issued by XERBLA.   

	INFO    (output) INTEGER   
	= 0:  successful exit   
	< 0:  if INFO = -i, the i-th argument had an illegal value   
	> 0:  if INFO = i, the algorithm failed to converge; i   
	off-diagonal elements of an intermediate tridiagonal   
	form did not converge to zero.   

	=====================================================================   


	Test the input parameters.   

	Parameter adjustments */
	/* Table of constant values */
	static integer_ c__1 = 1;
	static integer_ c_n1 = -1;
	static integer_ c__0 = 0;
	static doublereal_ c_b17 = 1.;

	/* System generated locals */
	integer_ a_dim1, a_offset, i__1, i__2;
	doublereal_ d__1;
	/* Builtin functions */
	//double sqrt_(doublereal_);
	/* Local variables */
	static integer_ nb;
	static doublereal_ eps;
	static integer_ inde;
	static doublereal_ anrm;
	static integer_ imax;
	static doublereal_ rmin, rmax;
	//extern /* Subroutine */ int dscal_(integer_ *, doublereal_ *, doublereal_ *, 
	//	integer_ *);
	static doublereal_ sigma;
	//extern logical_ lsame_(char *, char *);
	static integer_ iinfo;
	static logical_ lower, wantz;
	//extern doublereal_ dlamch_(char *);
	static integer_ iscale;
	//extern /* Subroutine */ int dlascl_(char *, integer_ *, integer_ *, 
	//	doublereal_ *, doublereal_ *, integer_ *, integer_ *, doublereal_ *, 
	//	integer_ *, integer_ *);
	static doublereal_ safmin;
	//extern integer_ ilaenv_(integer_ *, char *, char *, integer_ *, integer_ *, 
	//	integer_ *, integer_ *, ftnlen_, ftnlen_);
	//extern /* Subroutine */ int xerbla_(char *, integer_ *);
	static doublereal_ bignum;
	static integer_ indtau;
	//extern /* Subroutine */ int dsterf_(integer_ *, doublereal_ *, doublereal_ *,
	//	integer_ *);
	//extern doublereal_ dlansy_(char *, char *, integer_ *, doublereal_ *, 
	//	integer_ *, doublereal_ *);
	static integer_ indwrk;
	//extern /* Subroutine */ int dorgtr_(char *, integer_ *, doublereal_ *, 
	//	integer_ *, doublereal_ *, doublereal_ *, integer_ *, integer_ *), dsteqr_(char *, integer_ *, doublereal_ *, doublereal_ *, 
	//	doublereal_ *, integer_ *, doublereal_ *, integer_ *), 
	//	dsytrd_(char *, integer_ *, doublereal_ *, integer_ *, doublereal_ *, 
	//	doublereal_ *, doublereal_ *, doublereal_ *, integer_ *, integer_ *);
	static integer_ llwork;
	static doublereal_ smlnum;
	static integer_ lwkopt;
	static logical_ lquery;


	a_dim1 = *lda;
	a_offset = 1 + a_dim1;
	a -= a_offset;
	--w;
	--work;

	/* Function Body */
	wantz = lsame_(jobz, "V");
	lower = lsame_(uplo, "L");
	lquery = *lwork == -1;

	*info = 0;
	if (! (wantz || lsame_(jobz, "N"))) {
		*info = -1;
	} else if (! (lower || lsame_(uplo, "U"))) {
		*info = -2;
	} else if (*n < 0) {
		*info = -3;
	} else if (*lda < max(1,*n)) {
		*info = -5;
	}

	if (*info == 0) {
		nb = ilaenv_(&c__1, "DSYTRD", uplo, n, &c_n1, &c_n1, &c_n1, (ftnlen_)6,
			(ftnlen_)1);
		/* Computing MAX */
		i__1 = 1, i__2 = (nb + 2) * *n;
		lwkopt = max(i__1,i__2);
		work[1] = (doublereal_) lwkopt;

		/* Computing MAX */
		i__1 = 1, i__2 = *n * 3 - 1;
		if (*lwork < max(i__1,i__2) && ! lquery) {
			*info = -8;
		}
	}

	if (*info != 0) {
		i__1 = -(*info);
		xerbla_("DSYEV ", &i__1);
		return 0;
	} else if (lquery) {
		return 0;
	}

	/*     Quick return if possible */

	if (*n == 0) {
		return 0;
	}

	if (*n == 1) {
		w[1] = a[a_dim1 + 1];
		work[1] = 2.;
		if (wantz) {
			a[a_dim1 + 1] = 1.;
		}
		return 0;
	}

	/*     Get machine constants. */

	safmin = dlamch_("Safe minimum");
	eps = dlamch_("Precision");
	smlnum = safmin / eps;
	bignum = 1. / smlnum;
	rmin = sqrt_(smlnum);
	rmax = sqrt_(bignum);

	/*     Scale matrix to allowable range, if necessary. */

	anrm = dlansy_("M", uplo, n, &a[a_offset], lda, &work[1]);
	iscale = 0;
	if (anrm > 0. && anrm < rmin) {
		iscale = 1;
		sigma = rmin / anrm;
	} else if (anrm > rmax) {
		iscale = 1;
		sigma = rmax / anrm;
	}
	if (iscale == 1) {
		dlascl_(uplo, &c__0, &c__0, &c_b17, &sigma, n, n, &a[a_offset], lda, 
			info);
	}

	/*     Call DSYTRD to reduce symmetric matrix to tridiagonal form. */

	inde = 1;
	indtau = inde + *n;
	indwrk = indtau + *n;
	llwork = *lwork - indwrk + 1;
	dsytrd_(uplo, n, &a[a_offset], lda, &w[1], &work[inde], &work[indtau], &
		work[indwrk], &llwork, &iinfo);

	/*     For eigenvalues only, call DSTERF.  For eigenvectors, first call   
	DORGTR to generate the orthogonal matrix, then call DSTEQR. */

	if (! wantz) {
		dsterf_(n, &w[1], &work[inde], info);
	} else {
		dorgtr_(uplo, n, &a[a_offset], lda, &work[indtau], &work[indwrk], &
			llwork, &iinfo);
		dsteqr_(jobz, n, &w[1], &work[inde], &a[a_offset], lda, &work[indtau],
			info);
	}

	/*     If matrix was scaled, then rescale eigenvalues appropriately. */

	if (iscale == 1) {
		if (*info == 0) {
			imax = *n;
		} else {
			imax = *info - 1;
		}
		d__1 = 1. / sigma;
		dscal_(&imax, &d__1, &w[1], &c__1);
	}

	/*     Set WORK(1) to optimal workspace size. */

	work[1] = (doublereal_) lwkopt;

	return 0;

	/*     End of DSYEV */

} /* dsyev_ */



int TRL::dpotrf_(char *uplo, integer_ *n, doublereal_ *a, integer_ *
			lda, integer_ *info)
{
	/*  -- LAPACK routine (version 3.1) --   
	Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..   
	November 2006   


	Purpose   
	=======   

	DPOTRF computes the Cholesky factorization of a trl_real_ symmetric   
	positive definite matrix A.   

	The factorization has the form   
	A = U**T * U,  if UPLO = 'U', or   
	A = L  * L**T,  if UPLO = 'L',   
	where U is an upper triangular matrix and L is lower triangular.   

	This is the block version of the algorithm, calling Level 3 BLAS.   

	Arguments   
	=========   

	UPLO    (input) CHARACTER*1   
	= 'U':  Upper triangle of A is stored;   
	= 'L':  Lower triangle of A is stored.   

	N       (input) INTEGER   
	The order of the matrix A.  N >= 0.   

	A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
	On entry, the symmetric matrix A.  If UPLO = 'U', the leading   
	N-by-N upper triangular part of A contains the upper   
	triangular part of the matrix A, and the strictly lower   
	triangular part of A is not referenced.  If UPLO = 'L', the   
	leading N-by-N lower triangular part of A contains the lower   
	triangular part of the matrix A, and the strictly upper   
	triangular part of A is not referenced.   

	On exit, if INFO = 0, the factor U or L from the Cholesky   
	factorization A = U**T*U or A = L*L**T.   

	LDA     (input) INTEGER   
	The leading dimension of the array A.  LDA >= max(1,N).   

	INFO    (output) INTEGER   
	= 0:  successful exit   
	< 0:  if INFO = -i, the i-th argument had an illegal value   
	> 0:  if INFO = i, the leading minor of order i is not   
	positive definite, and the factorization could not be   
	completed.   

	=====================================================================   


	Test the input parameters.   

	Parameter adjustments */
	/* Table of constant values */
	static integer_ c__1 = 1;
	static integer_ c_n1 = -1;
	static doublereal_ c_b13 = -1.;
	static doublereal_ c_b14 = 1.;

	/* System generated locals */
	integer_ a_dim1, a_offset, i__1, i__2, i__3, i__4;
	/* Local variables */
	static integer_ j, jb, nb;
	//extern /* Subroutine */ int dgemm_(char *, char *, integer_ *, integer_ *, 
	//	integer_ *, doublereal_ *, doublereal_ *, integer_ *, doublereal_ *, 
	//	integer_ *, doublereal_ *, doublereal_ *, integer_ *);
	//extern logical_ lsame_(char *, char *);
	//extern /* Subroutine */ int dtrsm_(char *, char *, char *, char *, 
	//	integer_ *, integer_ *, doublereal_ *, doublereal_ *, integer_ *, 
	//	doublereal_ *, integer_ *);
	static logical_ upper;
	//extern /* Subroutine */ int dsyrk_(char *, char *, integer_ *, integer_ *, 
	//	doublereal_ *, doublereal_ *, integer_ *, doublereal_ *, doublereal_ *, 
	//	integer_ *), dpotf2_(char *, integer_ *, 
	//	doublereal_ *, integer_ *, integer_ *), xerbla_(char *, 
	//	integer_ *);
	//extern integer_ ilaenv_(integer_ *, char *, char *, integer_ *, integer_ *, 
	//	integer_ *, integer_ *, ftnlen_, ftnlen_);


	a_dim1 = *lda;
	a_offset = 1 + a_dim1;
	a -= a_offset;

	/* Function Body */
	*info = 0;
	upper = lsame_(uplo, "U");
	if (! upper && ! lsame_(uplo, "L")) {
		*info = -1;
	} else if (*n < 0) {
		*info = -2;
	} else if (*lda < max(1,*n)) {
		*info = -4;
	}
	if (*info != 0) {
		i__1 = -(*info);
		xerbla_("DPOTRF", &i__1);
		return 0;
	}

	/*     Quick return if possible */

	if (*n == 0) {
		return 0;
	}

	/*     Determine the block size for this environment. */

	nb = ilaenv_(&c__1, "DPOTRF", uplo, n, &c_n1, &c_n1, &c_n1, (ftnlen_)6, (
		ftnlen_)1);
	if (nb <= 1 || nb >= *n) {

		/*        Use unblocked code. */

		dpotf2_(uplo, n, &a[a_offset], lda, info);
	} else {

		/*        Use blocked code. */

		if (upper) {

			/*           Compute the Cholesky factorization A = U'*U. */

			i__1 = *n;
			i__2 = nb;
			for (j = 1; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {

				/*              Update and factorize the current diagonal block and test   
				for non-positive-definiteness.   

				Computing MIN */
				i__3 = nb, i__4 = *n - j + 1;
				jb = min(i__3,i__4);
				i__3 = j - 1;
				dsyrk_("Upper", "Transpose", &jb, &i__3, &c_b13, &a[j * 
					a_dim1 + 1], lda, &c_b14, &a[j + j * a_dim1], lda);
				dpotf2_("Upper", &jb, &a[j + j * a_dim1], lda, info);
				if (*info != 0) {
					goto L30;
				}
				if (j + jb <= *n) {

					/*                 Compute the current block row. */

					i__3 = *n - j - jb + 1;
					i__4 = j - 1;
					dgemm_("Transpose", "No transpose", &jb, &i__3, &i__4, &
						c_b13, &a[j * a_dim1 + 1], lda, &a[(j + jb) * 
						a_dim1 + 1], lda, &c_b14, &a[j + (j + jb) * 
						a_dim1], lda);
					i__3 = *n - j - jb + 1;
					dtrsm_("Left", "Upper", "Transpose", "Non-unit", &jb, &
						i__3, &c_b14, &a[j + j * a_dim1], lda, &a[j + (j 
						+ jb) * a_dim1], lda);
				}
				/* L10: */
			}

		} else {

			/*           Compute the Cholesky factorization A = L*L'. */

			i__2 = *n;
			i__1 = nb;
			for (j = 1; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {

				/*              Update and factorize the current diagonal block and test   
				for non-positive-definiteness.   

				Computing MIN */
				i__3 = nb, i__4 = *n - j + 1;
				jb = min(i__3,i__4);
				i__3 = j - 1;
				dsyrk_("Lower", "No transpose", &jb, &i__3, &c_b13, &a[j + 
					a_dim1], lda, &c_b14, &a[j + j * a_dim1], lda);
				dpotf2_("Lower", &jb, &a[j + j * a_dim1], lda, info);
				if (*info != 0) {
					goto L30;
				}
				if (j + jb <= *n) {

					/*                 Compute the current block column. */

					i__3 = *n - j - jb + 1;
					i__4 = j - 1;
					dgemm_("No transpose", "Transpose", &i__3, &jb, &i__4, &
						c_b13, &a[j + jb + a_dim1], lda, &a[j + a_dim1], 
						lda, &c_b14, &a[j + jb + j * a_dim1], lda);
					i__3 = *n - j - jb + 1;
					dtrsm_("Right", "Lower", "Transpose", "Non-unit", &i__3, &
						jb, &c_b14, &a[j + j * a_dim1], lda, &a[j + jb + 
						j * a_dim1], lda);
				}
				/* L20: */
			}
		}
	}
	goto L40;

L30:
	*info = *info + j - 1;

L40:
	return 0;

	/*     End of DPOTRF */

} /* dpotrf_ */


//
int TRL::dtrtrs_(char *uplo, char *trans, char *diag, integer_ *n, 
			integer_ *nrhs, doublereal_ *a, integer_ *lda, doublereal_ *b, integer_ *
			ldb, integer_ *info  	)
{
	/*  -- LAPACK routine (version 3.1) --   
	Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..   
	November 2006   


	Purpose   
	=======   

	DTRTRS solves a triangular system of the form   

	A * X = B  or  A**T * X = B,   

	where A is a triangular matrix of order N, and B is an N-by-NRHS   
	matrix.  A check is made to verify that A is nonsingular.   

	Arguments   
	=========   

	UPLO    (input) CHARACTER*1   
	= 'U':  A is upper triangular;   
	= 'L':  A is lower triangular.   

	TRANS   (input) CHARACTER*1   
	Specifies the form of the system of equations:   
	= 'N':  A * X = B  (No transpose)   
	= 'T':  A**T * X = B  (Transpose)   
	= 'C':  A**H * X = B  (Conjugate transpose = Transpose)   

	DIAG    (input) CHARACTER*1   
	= 'N':  A is non-unit triangular;   
	= 'U':  A is unit triangular.   

	N       (input) INTEGER   
	The order of the matrix A.  N >= 0.   

	NRHS    (input) INTEGER   
	The number of right hand sides, i.e., the number of columns   
	of the matrix B.  NRHS >= 0.   

	A       (input) DOUBLE PRECISION array, dimension (LDA,N)   
	The triangular matrix A.  If UPLO = 'U', the leading N-by-N   
	upper triangular part of the array A contains the upper   
	triangular matrix, and the strictly lower triangular part of   
	A is not referenced.  If UPLO = 'L', the leading N-by-N lower   
	triangular part of the array A contains the lower triangular   
	matrix, and the strictly upper triangular part of A is not   
	referenced.  If DIAG = 'U', the diagonal elements of A are   
	also not referenced and are assumed to be 1.   

	LDA     (input) INTEGER   
	The leading dimension of the array A.  LDA >= max(1,N).   

	B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)   
	On entry, the right hand side matrix B.   
	On exit, if INFO = 0, the solution matrix X.   

	LDB     (input) INTEGER   
	The leading dimension of the array B.  LDB >= max(1,N).   

	INFO    (output) INTEGER   
	= 0:  successful exit   
	< 0: if INFO = -i, the i-th argument had an illegal value   
	> 0: if INFO = i, the i-th diagonal element of A is zero,   
	indicating that the matrix is singular and the solutions   
	X have not been computed.   

	=====================================================================   


	Test the input parameters.   

	Parameter adjustments */
	/* Table of constant values */
	static doublereal_ c_b12 = 1.;

	/* System generated locals */
	integer_ a_dim1, a_offset, b_dim1, b_offset, i__1;
	/* Local variables */
	//extern logical_ lsame_(char *, char *);
	//extern /* Subroutine */ int dtrsm_(char *, char *, char *, char *, 
	//	integer_ *, integer_ *, doublereal_ *, doublereal_ *, integer_ *, 
	//	doublereal_ *, integer_ *), xerbla_(
	//	char *, integer_ *);
	static logical_ nounit;


	a_dim1 = *lda;
	a_offset = 1 + a_dim1;
	a -= a_offset;
	b_dim1 = *ldb;
	b_offset = 1 + b_dim1;
	b -= b_offset;

	/* Function Body */
	*info = 0;
	nounit = lsame_(diag, "N");
	if (! lsame_(uplo, "U") && ! lsame_(uplo, "L")) {
		*info = -1;
	} else if (! lsame_(trans, "N") && ! lsame_(trans, 
		"T") && ! lsame_(trans, "C")) {
			*info = -2;
	} else if (! nounit && ! lsame_(diag, "U")) {
		*info = -3;
	} else if (*n < 0) {
		*info = -4;
	} else if (*nrhs < 0) {
		*info = -5;
	} else if (*lda < max(1,*n)) {
		*info = -7;
	} else if (*ldb < max(1,*n)) {
		*info = -9;
	}
	if (*info != 0) {
		i__1 = -(*info);
		xerbla_("DTRTRS", &i__1);
		return 0;
	}

	/*     Quick return if possible */

	if (*n == 0) {
		return 0;
	}

	/*     Check for singularity. */

	if (nounit) {
		i__1 = *n;
		for (*info = 1; *info <= i__1; ++(*info)) {
			if (a[*info + *info * a_dim1] == 0.) {
				return 0;
			}
			/* L10: */
		}
	}
	*info = 0;

	/*     Solve A * x = b  or  A' * x = b. */

	dtrsm_("Left", uplo, trans, diag, n, nrhs, &c_b12, &a[a_offset], lda, &b[
		b_offset], ldb);

		return 0;

		/*     End of DTRTRS */

} /* dtrtrs_ */



int TRL::dsytrd_(char *uplo, integer_ *n, doublereal_ *a, integer_ *
			lda, doublereal_ *d__, doublereal_ *e, doublereal_ *tau, doublereal_ *
			work, integer_ *lwork, integer_ *info)
{
	/*  -- LAPACK routine (version 3.1) --   
	Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..   
	November 2006   


	Purpose   
	=======   

	DSYTRD reduces a trl_real_ symmetric matrix A to trl_real_ symmetric   
	tridiagonal form T by an orthogonal similarity transformation:   
	Q**T * A * Q = T.   

	Arguments   
	=========   

	UPLO    (input) CHARACTER*1   
	= 'U':  Upper triangle of A is stored;   
	= 'L':  Lower triangle of A is stored.   

	N       (input) INTEGER   
	The order of the matrix A.  N >= 0.   

	A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
	On entry, the symmetric matrix A.  If UPLO = 'U', the leading   
	N-by-N upper triangular part of A contains the upper   
	triangular part of the matrix A, and the strictly lower   
	triangular part of A is not referenced.  If UPLO = 'L', the   
	leading N-by-N lower triangular part of A contains the lower   
	triangular part of the matrix A, and the strictly upper   
	triangular part of A is not referenced.   
	On exit, if UPLO = 'U', the diagonal and first superdiagonal   
	of A are overwritten by the corresponding elements of the   
	tridiagonal matrix T, and the elements above the first   
	superdiagonal, with the array TAU, represent the orthogonal   
	matrix Q as a product of elementary reflectors; if UPLO   
	= 'L', the diagonal and first subdiagonal of A are over-   
	written by the corresponding elements of the tridiagonal   
	matrix T, and the elements below the first subdiagonal, with   
	the array TAU, represent the orthogonal matrix Q as a product   
	of elementary reflectors. See Further Details.   

	LDA     (input) INTEGER   
	The leading dimension of the array A.  LDA >= max(1,N).   

	D       (output) DOUBLE PRECISION array, dimension (N)   
	The diagonal elements of the tridiagonal matrix T:   
	D(i) = A(i,i).   

	E       (output) DOUBLE PRECISION array, dimension (N-1)   
	The off-diagonal elements of the tridiagonal matrix T:   
	E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'.   

	TAU     (output) DOUBLE PRECISION array, dimension (N-1)   
	The scalar factors of the elementary reflectors (see Further   
	Details).   

	WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))   
	On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

	LWORK   (input) INTEGER   
	The dimension of the array WORK.  LWORK >= 1.   
	For optimum performance LWORK >= N*NB, where NB is the   
	optimal blocksize.   

	If LWORK = -1, then a workspace query is assumed; the routine   
	only calculates the optimal size of the WORK array, returns   
	this value as the first entry of the WORK array, and no error   
	message related to LWORK is issued by XERBLA.   

	INFO    (output) INTEGER   
	= 0:  successful exit   
	< 0:  if INFO = -i, the i-th argument had an illegal value   

	Further Details   
	===============   

	If UPLO = 'U', the matrix Q is represented as a product of elementary   
	reflectors   

	Q = H(n-1) . . . H(2) H(1).   

	Each H(i) has the form   

	H(i) = I - tau * v * v'   

	where tau is a trl_real_ scalar, and v is a trl_real_ vector with   
	v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in   
	A(1:i-1,i+1), and tau in TAU(i).   

	If UPLO = 'L', the matrix Q is represented as a product of elementary   
	reflectors   

	Q = H(1) H(2) . . . H(n-1).   

	Each H(i) has the form   

	H(i) = I - tau * v * v'   

	where tau is a trl_real_ scalar, and v is a trl_real_ vector with   
	v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),   
	and tau in TAU(i).   

	The contents of A on exit are illustrated by the following examples   
	with n = 5:   

	if UPLO = 'U':                       if UPLO = 'L':   

	(  d   e   v2  v3  v4 )              (  d                  )   
	(      d   e   v3  v4 )              (  e   d              )   
	(          d   e   v4 )              (  v1  e   d          )   
	(              d   e  )              (  v1  v2  e   d      )   
	(                  d  )              (  v1  v2  v3  e   d  )   

	where d and e denote diagonal and off-diagonal elements of T, and vi   
	denotes an element of the vector defining H(i).   

	=====================================================================   


	Test the input parameters   

	Parameter adjustments */
	/* Table of constant values */
	static integer_ c__1 = 1;
	static integer_ c_n1 = -1;
	static integer_ c__3 = 3;
	static integer_ c__2 = 2;
	static doublereal_ c_b22 = -1.;
	static doublereal_ c_b23 = 1.;

	/* System generated locals */
	integer_ a_dim1, a_offset, i__1, i__2, i__3;
	/* Local variables */
	static integer_ i__, j, nb, kk, nx, iws;
	//extern logical_ lsame_(char *, char *);
	static integer_ nbmin, iinfo;
	static logical_ upper;
	//extern /* Subroutine */ int dsytd2_(char *, integer_ *, doublereal_ *, 
	//	integer_ *, doublereal_ *, doublereal_ *, doublereal_ *, integer_ *), dsyr2k_(char *, char *, integer_ *, integer_ *, doublereal_ 
	//	*, doublereal_ *, integer_ *, doublereal_ *, integer_ *, doublereal_ *, 
	//	doublereal_ *, integer_ *), dlatrd_(char *, 
	//	integer_ *, integer_ *, doublereal_ *, integer_ *, doublereal_ *, 
	//	doublereal_ *, doublereal_ *, integer_ *), xerbla_(char *, 
	//	integer_ *);
	//extern integer_ ilaenv_(integer_ *, char *, char *, integer_ *, integer_ *, 
	//	integer_ *, integer_ *, ftnlen_, ftnlen_);
	static integer_ ldwork, lwkopt;
	static logical_ lquery;


	a_dim1 = *lda;
	a_offset = 1 + a_dim1;
	a -= a_offset;
	--d__;
	--e;
	--tau;
	--work;

	/* Function Body */
	*info = 0;
	upper = lsame_(uplo, "U");
	lquery = *lwork == -1;
	if (! upper && ! lsame_(uplo, "L")) {
		*info = -1;
	} else if (*n < 0) {
		*info = -2;
	} else if (*lda < max(1,*n)) {
		*info = -4;
	} else if (*lwork < 1 && ! lquery) {
		*info = -9;
	}

	if (*info == 0) {

		/*        Determine the block size. */

		nb = ilaenv_(&c__1, "DSYTRD", uplo, n, &c_n1, &c_n1, &c_n1, (ftnlen_)6,
			(ftnlen_)1);
		lwkopt = *n * nb;
		work[1] = (doublereal_) lwkopt;
	}

	if (*info != 0) {
		i__1 = -(*info);
		xerbla_("DSYTRD", &i__1);
		return 0;
	} else if (lquery) {
		return 0;
	}

	/*     Quick return if possible */

	if (*n == 0) {
		work[1] = 1.;
		return 0;
	}

	nx = *n;
	iws = 1;
	if (nb > 1 && nb < *n) {

		/*        Determine when to cross over from blocked to unblocked code   
		(last block is always handled by unblocked code).   

		Computing MAX */
		i__1 = nb, i__2 = ilaenv_(&c__3, "DSYTRD", uplo, n, &c_n1, &c_n1, &
			c_n1, (ftnlen_)6, (ftnlen_)1);
		nx = max(i__1,i__2);
		if (nx < *n) {

			/*           Determine if workspace is large enough for blocked code. */

			ldwork = *n;
			iws = ldwork * nb;
			if (*lwork < iws) {

				/*              Not enough workspace to use optimal NB:  determine the   
				minimum value of NB, and reduce NB or force use of   
				unblocked code by setting NX = N.   

				Computing MAX */
				i__1 = *lwork / ldwork;
				nb = max(i__1,1);
				nbmin = ilaenv_(&c__2, "DSYTRD", uplo, n, &c_n1, &c_n1, &c_n1,
					(ftnlen_)6, (ftnlen_)1);
				if (nb < nbmin) {
					nx = *n;
				}
			}
		} else {
			nx = *n;
		}
	} else {
		nb = 1;
	}

	if (upper) {

		/*        Reduce the upper triangle of A.   
		Columns 1:kk are handled by the unblocked method. */

		kk = *n - (*n - nx + nb - 1) / nb * nb;
		i__1 = kk + 1;
		i__2 = -nb;
		for (i__ = *n - nb + 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += 
			i__2) {

				/*           Reduce columns i:i+nb-1 to tridiagonal form and form the   
				matrix W which is needed to update the unreduced part of   
				the matrix */

				i__3 = i__ + nb - 1;
				dlatrd_(uplo, &i__3, &nb, &a[a_offset], lda, &e[1], &tau[1], &
					work[1], &ldwork);

				/*           Update the unreduced submatrix A(1:i-1,1:i-1), using an   
				update of the form:  A := A - V*W' - W*V' */

				i__3 = i__ - 1;
				dsyr2k_(uplo, "No transpose", &i__3, &nb, &c_b22, &a[i__ * a_dim1 
					+ 1], lda, &work[1], &ldwork, &c_b23, &a[a_offset], lda);

				/*           Copy superdiagonal elements back into A, and diagonal   
				elements into D */

				i__3 = i__ + nb - 1;
				for (j = i__; j <= i__3; ++j) {
					a[j - 1 + j * a_dim1] = e[j - 1];
					d__[j] = a[j + j * a_dim1];
					/* L10: */
				}
				/* L20: */
		}

		/*        Use unblocked code to reduce the last or only block */

		dsytd2_(uplo, &kk, &a[a_offset], lda, &d__[1], &e[1], &tau[1], &iinfo);
	} else {

		/*        Reduce the lower triangle of A */

		i__2 = *n - nx;
		i__1 = nb;
		for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {

			/*           Reduce columns i:i+nb-1 to tridiagonal form and form the   
			matrix W which is needed to update the unreduced part of   
			the matrix */

			i__3 = *n - i__ + 1;
			dlatrd_(uplo, &i__3, &nb, &a[i__ + i__ * a_dim1], lda, &e[i__], &
				tau[i__], &work[1], &ldwork);

			/*           Update the unreduced submatrix A(i+ib:n,i+ib:n), using   
			an update of the form:  A := A - V*W' - W*V' */

			i__3 = *n - i__ - nb + 1;
			dsyr2k_(uplo, "No transpose", &i__3, &nb, &c_b22, &a[i__ + nb + 
				i__ * a_dim1], lda, &work[nb + 1], &ldwork, &c_b23, &a[
					i__ + nb + (i__ + nb) * a_dim1], lda);

					/*           Copy subdiagonal elements back into A, and diagonal   
					elements into D */

					i__3 = i__ + nb - 1;
					for (j = i__; j <= i__3; ++j) {
						a[j + 1 + j * a_dim1] = e[j];
						d__[j] = a[j + j * a_dim1];
						/* L30: */
					}
					/* L40: */
		}

		/*        Use unblocked code to reduce the last or only block */

		i__1 = *n - i__ + 1;
		dsytd2_(uplo, &i__1, &a[i__ + i__ * a_dim1], lda, &d__[i__], &e[i__], 
			&tau[i__], &iinfo);
	}

	work[1] = (doublereal_) lwkopt;
	return 0;

	/*     End of DSYTRD */

} /* dsytrd_ */


int TRL::dorgtr_(char *uplo, integer_ *n, doublereal_ *a, integer_ *
			lda, doublereal_ *tau, doublereal_ *work, integer_ *lwork, integer_ *info)
{
	/*  -- LAPACK routine (version 3.1) --   
	Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..   
	November 2006   


	Purpose   
	=======   

	DORGTR generates a trl_real_ orthogonal matrix Q which is defined as the   
	product of n-1 elementary reflectors of order N, as returned by   
	DSYTRD:   

	if UPLO = 'U', Q = H(n-1) . . . H(2) H(1),   

	if UPLO = 'L', Q = H(1) H(2) . . . H(n-1).   

	Arguments   
	=========   

	UPLO    (input) CHARACTER*1   
	= 'U': Upper triangle of A contains elementary reflectors   
	from DSYTRD;   
	= 'L': Lower triangle of A contains elementary reflectors   
	from DSYTRD.   

	N       (input) INTEGER   
	The order of the matrix Q. N >= 0.   

	A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
	On entry, the vectors which define the elementary reflectors,   
	as returned by DSYTRD.   
	On exit, the N-by-N orthogonal matrix Q.   

	LDA     (input) INTEGER   
	The leading dimension of the array A. LDA >= max(1,N).   

	TAU     (input) DOUBLE PRECISION array, dimension (N-1)   
	TAU(i) must contain the scalar factor of the elementary   
	reflector H(i), as returned by DSYTRD.   

	WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))   
	On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

	LWORK   (input) INTEGER   
	The dimension of the array WORK. LWORK >= max(1,N-1).   
	For optimum performance LWORK >= (N-1)*NB, where NB is   
	the optimal blocksize.   

	If LWORK = -1, then a workspace query is assumed; the routine   
	only calculates the optimal size of the WORK array, returns   
	this value as the first entry of the WORK array, and no error   
	message related to LWORK is issued by XERBLA.   

	INFO    (output) INTEGER   
	= 0:  successful exit   
	< 0:  if INFO = -i, the i-th argument had an illegal value   

	=====================================================================   


	Test the input arguments   

	Parameter adjustments */
	/* Table of constant values */
	static integer_ c__1 = 1;
	static integer_ c_n1 = -1;

	/* System generated locals */
	integer_ a_dim1, a_offset, i__1, i__2, i__3;
	/* Local variables */
	static integer_ i__, j, nb;
	//extern logical_ lsame_(char *, char *);
	static integer_ iinfo;
	static logical_ upper;
	//extern /* Subroutine */ int xerbla_(char *, integer_ *);
	//extern integer_ ilaenv_(integer_ *, char *, char *, integer_ *, integer_ *, 
	//	integer_ *, integer_ *, ftnlen_, ftnlen_);
	//extern /* Subroutine */ int dorgql_(integer_ *, integer_ *, integer_ *, 
	//	doublereal_ *, integer_ *, doublereal_ *, doublereal_ *, integer_ *, 
	//	integer_ *), dorgqr_(integer_ *, integer_ *, integer_ *, doublereal_ *,
	//	integer_ *, doublereal_ *, doublereal_ *, integer_ *, integer_ *);
	static integer_ lwkopt;
	static logical_ lquery;


	a_dim1 = *lda;
	a_offset = 1 + a_dim1;
	a -= a_offset;
	--tau;
	--work;

	/* Function Body */
	*info = 0;
	lquery = *lwork == -1;
	upper = lsame_(uplo, "U");
	if (! upper && ! lsame_(uplo, "L")) {
		*info = -1;
	} else if (*n < 0) {
		*info = -2;
	} else if (*lda < max(1,*n)) {
		*info = -4;
	} else /* if(complicated condition) */ {
		/* Computing MAX */
		i__1 = 1, i__2 = *n - 1;
		if (*lwork < max(i__1,i__2) && ! lquery) {
			*info = -7;
		}
	}

	if (*info == 0) {
		if (upper) {
			i__1 = *n - 1;
			i__2 = *n - 1;
			i__3 = *n - 1;
			nb = ilaenv_(&c__1, "DORGQL", " ", &i__1, &i__2, &i__3, &c_n1, (
				ftnlen_)6, (ftnlen_)1);
		} else {
			i__1 = *n - 1;
			i__2 = *n - 1;
			i__3 = *n - 1;
			nb = ilaenv_(&c__1, "DORGQR", " ", &i__1, &i__2, &i__3, &c_n1, (
				ftnlen_)6, (ftnlen_)1);
		}
		/* Computing MAX */
		i__1 = 1, i__2 = *n - 1;
		lwkopt = max(i__1,i__2) * nb;
		work[1] = (doublereal_) lwkopt;
	}

	if (*info != 0) {
		i__1 = -(*info);
		xerbla_("DORGTR", &i__1);
		return 0;
	} else if (lquery) {
		return 0;
	}

	/*     Quick return if possible */

	if (*n == 0) {
		work[1] = 1.;
		return 0;
	}

	if (upper) {

		/*        Q was determined by a call to DSYTRD with UPLO = 'U'   

		Shift the vectors which define the elementary reflectors one   
		column to the left, and set the last row and column of Q to   
		those of the unit matrix */

		i__1 = *n - 1;
		for (j = 1; j <= i__1; ++j) {
			i__2 = j - 1;
			for (i__ = 1; i__ <= i__2; ++i__) {
				a[i__ + j * a_dim1] = a[i__ + (j + 1) * a_dim1];
				/* L10: */
			}
			a[*n + j * a_dim1] = 0.;
			/* L20: */
		}
		i__1 = *n - 1;
		for (i__ = 1; i__ <= i__1; ++i__) {
			a[i__ + *n * a_dim1] = 0.;
			/* L30: */
		}
		a[*n + *n * a_dim1] = 1.;

		/*        Generate Q(1:n-1,1:n-1) */

		i__1 = *n - 1;
		i__2 = *n - 1;
		i__3 = *n - 1;
		dorgql_(&i__1, &i__2, &i__3, &a[a_offset], lda, &tau[1], &work[1], 
			lwork, &iinfo);

	} else {

		/*        Q was determined by a call to DSYTRD with UPLO = 'L'.   

		Shift the vectors which define the elementary reflectors one   
		column to the right, and set the first row and column of Q to   
		those of the unit matrix */

		for (j = *n; j >= 2; --j) {
			a[j * a_dim1 + 1] = 0.;
			i__1 = *n;
			for (i__ = j + 1; i__ <= i__1; ++i__) {
				a[i__ + j * a_dim1] = a[i__ + (j - 1) * a_dim1];
				/* L40: */
			}
			/* L50: */
		}
		a[a_dim1 + 1] = 1.;
		i__1 = *n;
		for (i__ = 2; i__ <= i__1; ++i__) {
			a[i__ + a_dim1] = 0.;
			/* L60: */
		}
		if (*n > 1) {

			/*           Generate Q(2:n,2:n) */

			i__1 = *n - 1;
			i__2 = *n - 1;
			i__3 = *n - 1;
			dorgqr_(&i__1, &i__2, &i__3, &a[(a_dim1 << 1) + 2], lda, &tau[1], 
				&work[1], lwork, &iinfo);
		}
	}
	work[1] = (doublereal_) lwkopt;
	return 0;

	/*     End of DORGTR */

} /* dorgtr_ */


int TRL::dstein_(integer_ *n, doublereal_ *d__, doublereal_ *e, 
			integer_ *m, doublereal_ *w, integer_ *iblock, integer_ *isplit, 
			doublereal_ *z__, integer_ *ldz, doublereal_ *work, integer_ *iwork, 
			integer_ *ifail, integer_ *info)
{
	/*  -- LAPACK routine (version 3.1) --   
	Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..   
	November 2006   


	Purpose   
	=======   

	DSTEIN computes the eigenvectors of a trl_real_ symmetric tridiagonal   
	matrix T corresponding to specified eigenvalues, using inverse   
	iteration.   

	The maximum number of iterations allowed for each eigenvector is   
	specified by an internal parameter MAXITS (currently set to 5).   

	Arguments   
	=========   

	N       (input) INTEGER   
	The order of the matrix.  N >= 0.   

	D       (input) DOUBLE PRECISION array, dimension (N)   
	The n diagonal elements of the tridiagonal matrix T.   

	E       (input) DOUBLE PRECISION array, dimension (N-1)   
	The (n-1) subdiagonal elements of the tridiagonal matrix   
	T, in elements 1 to N-1.   

	M       (input) INTEGER   
	The number of eigenvectors to be found.  0 <= M <= N.   

	W       (input) DOUBLE PRECISION array, dimension (N)   
	The first M elements of W contain the eigenvalues for   
	which eigenvectors are to be computed.  The eigenvalues   
	should be grouped by split-off block and ordered from   
	smallest to largest within the block.  ( The output array   
	W from DSTEBZ with ORDER = 'B' is expected here. )   

	IBLOCK  (input) INTEGER array, dimension (N)   
	The submatrix indices associated with the corresponding   
	eigenvalues in W; IBLOCK(i)=1 if eigenvalue W(i) belongs to   
	the first submatrix from the top, =2 if W(i) belongs to   
	the second submatrix, etc.  ( The output array IBLOCK   
	from DSTEBZ is expected here. )   

	ISPLIT  (input) INTEGER array, dimension (N)   
	The splitting points, at which T breaks up into submatrices.   
	The first submatrix consists of rows/columns 1 to   
	ISPLIT( 1 ), the second of rows/columns ISPLIT( 1 )+1   
	through ISPLIT( 2 ), etc.   
	( The output array ISPLIT from DSTEBZ is expected here. )   

	Z       (output) DOUBLE PRECISION array, dimension (LDZ, M)   
	The computed eigenvectors.  The eigenvector associated   
	with the eigenvalue W(i) is stored in the i-th column of   
	Z.  Any vector which fails to converge is set to its current   
	iterate after MAXITS iterations.   

	LDZ     (input) INTEGER   
	The leading dimension of the array Z.  LDZ >= max(1,N).   

	WORK    (workspace) DOUBLE PRECISION array, dimension (5*N)   

	IWORK   (workspace) INTEGER array, dimension (N)   

	IFAIL   (output) INTEGER array, dimension (M)   
	On normal exit, all elements of IFAIL are zero.   
	If one or more eigenvectors fail to converge after   
	MAXITS iterations, then their indices are stored in   
	array IFAIL.   

	INFO    (output) INTEGER   
	= 0: successful exit.   
	< 0: if INFO = -i, the i-th argument had an illegal value   
	> 0: if INFO = i, then i eigenvectors failed to converge   
	in MAXITS iterations.  Their indices are stored in   
	array IFAIL.   

	Internal Parameters   
	===================   

	MAXITS  INTEGER, default = 5   
	The maximum number of iterations performed.   

	EXTRA   INTEGER, default = 2   
	The number of iterations performed after norm growth   
	criterion is satisfied, should be at least 1.   

	=====================================================================   


	Test the input parameters.   

	Parameter adjustments */
	/* Table of constant values */
	static integer_ c__2 = 2;
	static integer_ c__1 = 1;
	static integer_ c_n1 = -1;

	/* System generated locals */
	integer_ z_dim1, z_offset, i__1, i__2, i__3;
	doublereal_ d__1, d__2, d__3, d__4, d__5;
	/* Builtin functions */
	//double sqrt_(doublereal_);
	/* Local variables */
	static integer_ i__, j, b1, j1, bn;
	static doublereal_ xj, scl, eps, sep, nrm, tol;
	static integer_ its;
	static doublereal_ xjm, ztr, eps1;
	static integer_ jblk, nblk;
	//extern doublereal_ ddot_(integer_ *, doublereal_ *, integer_ *, doublereal_ *, 
	//	integer_ *);
	static integer_ jmax;
	//extern doublereal_ dnrm2_(integer_ *, doublereal_ *, integer_ *);
	//extern /* Subroutine */ int dscal_(integer_ *, doublereal_ *, doublereal_ *, 
	//	integer_ *);
	static integer_ iseed[4], gpind, iinfo;
	//extern doublereal_ dasum_(integer_ *, doublereal_ *, integer_ *);
	//extern /* Subroutine */ int dcopy_(integer_ *, doublereal_ *, integer_ *, 
	//	doublereal_ *, integer_ *), daxpy_(integer_ *, doublereal_ *, 
	//	doublereal_ *, integer_ *, doublereal_ *, integer_ *);
	static doublereal_ ortol;
	static integer_ indrv1, indrv2, indrv3, indrv4, indrv5;
	//extern doublereal_ dlamch_(char *);
	//extern /* Subroutine */ int dlagtf_(integer_ *, doublereal_ *, doublereal_ *,
	//	doublereal_ *, doublereal_ *, doublereal_ *, doublereal_ *, integer_ *
	//	, integer_ *);
	//extern integer_ idamax_(integer_ *, doublereal_ *, integer_ *);
	//extern /* Subroutine */ int xerbla_(char *, integer_ *), dlagts_(
	//	integer_ *, integer_ *, doublereal_ *, doublereal_ *, doublereal_ *, 
	//	doublereal_ *, integer_ *, doublereal_ *, doublereal_ *, integer_ *);
	static integer_ nrmchk;
	//extern /* Subroutine */ int dlarnv_(integer_ *, integer_ *, integer_ *, 
	//	doublereal_ *);
	static integer_ blksiz;
	static doublereal_ onenrm, dtpcrt, pertol;


	--d__;
	--e;
	--w;
	--iblock;
	--isplit;
	z_dim1 = *ldz;
	z_offset = 1 + z_dim1;
	z__ -= z_offset;
	--work;
	--iwork;
	--ifail;

	/* Function Body */
	*info = 0;
	i__1 = *m;
	for (i__ = 1; i__ <= i__1; ++i__) {
		ifail[i__] = 0;
		/* L10: */
	}

	if (*n < 0) {
		*info = -1;
	} else if (*m < 0 || *m > *n) {
		*info = -4;
	} else if (*ldz < max(1,*n)) {
		*info = -9;
	} else {
		i__1 = *m;
		for (j = 2; j <= i__1; ++j) {
			if (iblock[j] < iblock[j - 1]) {
				*info = -6;
				goto L30;
			}
			if (iblock[j] == iblock[j - 1] && w[j] < w[j - 1]) {
				*info = -5;
				goto L30;
			}
			/* L20: */
		}
L30:
		;
	}

	if (*info != 0) {
		i__1 = -(*info);
		xerbla_("DSTEIN", &i__1);
		return 0;
	}

	/*     Quick return if possible */

	if (*n == 0 || *m == 0) {
		return 0;
	} else if (*n == 1) {
		z__[z_dim1 + 1] = 1.;
		return 0;
	}

	/*     Get machine constants. */

	eps = dlamch_("Precision");

	/*     Initialize seed for random number generator DLARNV. */

	for (i__ = 1; i__ <= 4; ++i__) {
		iseed[i__ - 1] = 1;
		/* L40: */
	}

	/*     Initialize pointers. */

	indrv1 = 0;
	indrv2 = indrv1 + *n;
	indrv3 = indrv2 + *n;
	indrv4 = indrv3 + *n;
	indrv5 = indrv4 + *n;

	/*     Compute eigenvectors of matrix blocks. */

	j1 = 1;
	i__1 = iblock[*m];
	for (nblk = 1; nblk <= i__1; ++nblk) {

		/*        Find starting and ending indices of block nblk. */

		if (nblk == 1) {
			b1 = 1;
		} else {
			b1 = isplit[nblk - 1] + 1;
		}
		bn = isplit[nblk];
		blksiz = bn - b1 + 1;
		if (blksiz == 1) {
			goto L60;
		}
		gpind = b1;

		/*        Compute reorthogonalization criterion and stopping criterion. */

		onenrm = (d__1 = d__[b1], abs(d__1)) + (d__2 = e[b1], abs(d__2));
		/* Computing MAX */
		d__3 = onenrm, d__4 = (d__1 = d__[bn], abs(d__1)) + (d__2 = e[bn - 1],
			abs(d__2));
		onenrm = max(d__3,d__4);
		i__2 = bn - 1;
		for (i__ = b1 + 1; i__ <= i__2; ++i__) {
			/* Computing MAX */
			d__4 = onenrm, d__5 = (d__1 = d__[i__], abs(d__1)) + (d__2 = e[
				i__ - 1], abs(d__2)) + (d__3 = e[i__], abs(d__3));
				onenrm = max(d__4,d__5);
				/* L50: */
		}
		ortol = onenrm * .001;

		dtpcrt = sqrt_(.1 / blksiz);

		/*        Loop through eigenvalues of block nblk. */

L60:
		jblk = 0;
		i__2 = *m;
		for (j = j1; j <= i__2; ++j) {
			if (iblock[j] != nblk) {
				j1 = j;
				goto L160;
			}
			++jblk;
			xj = w[j];

			/*           Skip all the work if the block size is one. */

			if (blksiz == 1) {
				work[indrv1 + 1] = 1.;
				goto L120;
			}

			/*           If eigenvalues j and j-1 are too close, add a relatively   
			small perturbation. */

			if (jblk > 1) {
				eps1 = (d__1 = eps * xj, abs(d__1));
				pertol = eps1 * 10.;
				sep = xj - xjm;
				if (sep < pertol) {
					xj = xjm + pertol;
				}
			}

			its = 0;
			nrmchk = 0;

			/*           Get random starting vector. */

			dlarnv_(&c__2, iseed, &blksiz, &work[indrv1 + 1]);

			/*           Copy the matrix T so it won't be destroyed in factorization. */

			dcopy_(&blksiz, &d__[b1], &c__1, &work[indrv4 + 1], &c__1);
			i__3 = blksiz - 1;
			dcopy_(&i__3, &e[b1], &c__1, &work[indrv2 + 2], &c__1);
			i__3 = blksiz - 1;
			dcopy_(&i__3, &e[b1], &c__1, &work[indrv3 + 1], &c__1);

			/*           Compute LU factors with partial pivoting  ( PT = LU ) */

			tol = 0.;
			dlagtf_(&blksiz, &work[indrv4 + 1], &xj, &work[indrv2 + 2], &work[
				indrv3 + 1], &tol, &work[indrv5 + 1], &iwork[1], &iinfo);

				/*           Update iteration count. */

L70:
				++its;
				if (its > 5) {
					goto L100;
				}

				/*           Normalize and scale the righthand side vector Pb.   

				Computing MAX */
				d__2 = eps, d__3 = (d__1 = work[indrv4 + blksiz], abs(d__1));
				scl = blksiz * onenrm * max(d__2,d__3) / dasum_(&blksiz, &work[
					indrv1 + 1], &c__1);
					dscal_(&blksiz, &scl, &work[indrv1 + 1], &c__1);

					/*           Solve the system LU = Pb. */

					dlagts_(&c_n1, &blksiz, &work[indrv4 + 1], &work[indrv2 + 2], &
						work[indrv3 + 1], &work[indrv5 + 1], &iwork[1], &work[
							indrv1 + 1], &tol, &iinfo);

							/*           Reorthogonalize by modified Gram-Schmidt if eigenvalues are   
							close enough. */

							if (jblk == 1) {
								goto L90;
							}
							if ((d__1 = xj - xjm, abs(d__1)) > ortol) {
								gpind = j;
							}
							if (gpind != j) {
								i__3 = j - 1;
								for (i__ = gpind; i__ <= i__3; ++i__) {
									ztr = -ddot_(&blksiz, &work[indrv1 + 1], &c__1, &z__[b1 + 
										i__ * z_dim1], &c__1);
									daxpy_(&blksiz, &ztr, &z__[b1 + i__ * z_dim1], &c__1, &
										work[indrv1 + 1], &c__1);
									/* L80: */
								}
							}

							/*           Check the infinity norm of the iterate. */

L90:
							jmax = idamax_(&blksiz, &work[indrv1 + 1], &c__1);
							nrm = (d__1 = work[indrv1 + jmax], abs(d__1));

							/*           Continue for additional iterations after norm reaches   
							stopping criterion. */

							if (nrm < dtpcrt) {
								goto L70;
							}
							++nrmchk;
							if (nrmchk < 3) {
								goto L70;
							}

							goto L110;

							/*           If stopping criterion was not satisfied, update info and   
							store eigenvector number in array ifail. */

L100:
							++(*info);
							ifail[*info] = j;

							/*           Accept iterate as jth eigenvector. */

L110:
							scl = 1. / dnrm2_(&blksiz, &work[indrv1 + 1], &c__1);
							jmax = idamax_(&blksiz, &work[indrv1 + 1], &c__1);
							if (work[indrv1 + jmax] < 0.) {
								scl = -scl;
							}
							dscal_(&blksiz, &scl, &work[indrv1 + 1], &c__1);
L120:
							i__3 = *n;
							for (i__ = 1; i__ <= i__3; ++i__) {
								z__[i__ + j * z_dim1] = 0.;
								/* L130: */
							}
							i__3 = blksiz;
							for (i__ = 1; i__ <= i__3; ++i__) {
								z__[b1 + i__ - 1 + j * z_dim1] = work[indrv1 + i__];
								/* L140: */
							}

							/*           Save the shift to check eigenvalue spacing at next   
							iteration. */

							xjm = xj;

							/* L150: */
		}
L160:
		;
	}

	return 0;

	/*     End of DSTEIN */

} /* dstein_ */




//////////////////////////
//////////////////////////////
/////////////////////////////////////

int TRL::dsytd2_(char *uplo, integer_ *n, doublereal_ *a, integer_ *
			lda, doublereal_ *d__, doublereal_ *e, doublereal_ *tau, integer_ *info)
{
	/*  -- LAPACK routine (version 3.1) --   
	Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..   
	November 2006   


	Purpose   
	=======   

	DSYTD2 reduces a trl_real_ symmetric matrix A to symmetric tridiagonal   
	form T by an orthogonal similarity transformation: Q' * A * Q = T.   

	Arguments   
	=========   

	UPLO    (input) CHARACTER*1   
	Specifies whether the upper or lower triangular part of the   
	symmetric matrix A is stored:   
	= 'U':  Upper triangular   
	= 'L':  Lower triangular   

	N       (input) INTEGER   
	The order of the matrix A.  N >= 0.   

	A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
	On entry, the symmetric matrix A.  If UPLO = 'U', the leading   
	n-by-n upper triangular part of A contains the upper   
	triangular part of the matrix A, and the strictly lower   
	triangular part of A is not referenced.  If UPLO = 'L', the   
	leading n-by-n lower triangular part of A contains the lower   
	triangular part of the matrix A, and the strictly upper   
	triangular part of A is not referenced.   
	On exit, if UPLO = 'U', the diagonal and first superdiagonal   
	of A are overwritten by the corresponding elements of the   
	tridiagonal matrix T, and the elements above the first   
	superdiagonal, with the array TAU, represent the orthogonal   
	matrix Q as a product of elementary reflectors; if UPLO   
	= 'L', the diagonal and first subdiagonal of A are over-   
	written by the corresponding elements of the tridiagonal   
	matrix T, and the elements below the first subdiagonal, with   
	the array TAU, represent the orthogonal matrix Q as a product   
	of elementary reflectors. See Further Details.   

	LDA     (input) INTEGER   
	The leading dimension of the array A.  LDA >= max(1,N).   

	D       (output) DOUBLE PRECISION array, dimension (N)   
	The diagonal elements of the tridiagonal matrix T:   
	D(i) = A(i,i).   

	E       (output) DOUBLE PRECISION array, dimension (N-1)   
	The off-diagonal elements of the tridiagonal matrix T:   
	E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'.   

	TAU     (output) DOUBLE PRECISION array, dimension (N-1)   
	The scalar factors of the elementary reflectors (see Further   
	Details).   

	INFO    (output) INTEGER   
	= 0:  successful exit   
	< 0:  if INFO = -i, the i-th argument had an illegal value.   

	Further Details   
	===============   

	If UPLO = 'U', the matrix Q is represented as a product of elementary   
	reflectors   

	Q = H(n-1) . . . H(2) H(1).   

	Each H(i) has the form   

	H(i) = I - tau * v * v'   

	where tau is a trl_real_ scalar, and v is a trl_real_ vector with   
	v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in   
	A(1:i-1,i+1), and tau in TAU(i).   

	If UPLO = 'L', the matrix Q is represented as a product of elementary   
	reflectors   

	Q = H(1) H(2) . . . H(n-1).   

	Each H(i) has the form   

	H(i) = I - tau * v * v'   

	where tau is a trl_real_ scalar, and v is a trl_real_ vector with   
	v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),   
	and tau in TAU(i).   

	The contents of A on exit are illustrated by the following examples   
	with n = 5:   

	if UPLO = 'U':                       if UPLO = 'L':   

	(  d   e   v2  v3  v4 )              (  d                  )   
	(      d   e   v3  v4 )              (  e   d              )   
	(          d   e   v4 )              (  v1  e   d          )   
	(              d   e  )              (  v1  v2  e   d      )   
	(                  d  )              (  v1  v2  v3  e   d  )   

	where d and e denote diagonal and off-diagonal elements of T, and vi   
	denotes an element of the vector defining H(i).   

	=====================================================================   


	Test the input parameters   

	Parameter adjustments */
	/* Table of constant values */
	static integer_ c__1 = 1;
	static doublereal_ c_b8 = 0.;
	static doublereal_ c_b14 = -1.;

	/* System generated locals */
	integer_ a_dim1, a_offset, i__1, i__2, i__3;
	/* Local variables */
	static integer_ i__;
	//extern doublereal_ ddot_(integer_ *, doublereal_ *, integer_ *, doublereal_ *, 
	//	integer_ *);
	static doublereal_ taui;
	//extern /* Subroutine */ int dsyr2_(char *, integer_ *, doublereal_ *, 
	//	doublereal_ *, integer_ *, doublereal_ *, integer_ *, doublereal_ *, 
	//	integer_ *);
	static doublereal_ alpha;
	//extern logical_ lsame_(char *, char *);
	//extern /* Subroutine */ int daxpy_(integer_ *, doublereal_ *, doublereal_ *, 
	//	integer_ *, doublereal_ *, integer_ *);
	static logical_ upper;
	//extern /* Subroutine */ int dsymv_(char *, integer_ *, doublereal_ *, 
	//	doublereal_ *, integer_ *, doublereal_ *, integer_ *, doublereal_ *, 
	//	doublereal_ *, integer_ *), dlarfg_(integer_ *, doublereal_ *,
	//	doublereal_ *, integer_ *, doublereal_ *), xerbla_(char *, integer_ *
	//	);


	a_dim1 = *lda;
	a_offset = 1 + a_dim1;
	a -= a_offset;
	--d__;
	--e;
	--tau;

	/* Function Body */
	*info = 0;
	upper = lsame_(uplo, "U");
	if (! upper && ! lsame_(uplo, "L")) {
		*info = -1;
	} else if (*n < 0) {
		*info = -2;
	} else if (*lda < max(1,*n)) {
		*info = -4;
	}
	if (*info != 0) {
		i__1 = -(*info);
		xerbla_("DSYTD2", &i__1);
		return 0;
	}

	/*     Quick return if possible */

	if (*n <= 0) {
		return 0;
	}

	if (upper) {

		/*        Reduce the upper triangle of A */

		for (i__ = *n - 1; i__ >= 1; --i__) {

			/*           Generate elementary reflector H(i) = I - tau * v * v'   
			to annihilate A(1:i-1,i+1) */

			dlarfg_(&i__, &a[i__ + (i__ + 1) * a_dim1], &a[(i__ + 1) * a_dim1 
				+ 1], &c__1, &taui);
			e[i__] = a[i__ + (i__ + 1) * a_dim1];

			if (taui != 0.) {

				/*              Apply H(i) from both sides to A(1:i,1:i) */

				a[i__ + (i__ + 1) * a_dim1] = 1.;

				/*              Compute  x := tau * A * v  storing x in TAU(1:i) */

				dsymv_(uplo, &i__, &taui, &a[a_offset], lda, &a[(i__ + 1) * 
					a_dim1 + 1], &c__1, &c_b8, &tau[1], &c__1);

				/*              Compute  w := x - 1/2 * tau * (x'*v) * v */

				alpha = taui * -.5 * ddot_(&i__, &tau[1], &c__1, &a[(i__ + 1) 
					* a_dim1 + 1], &c__1);
				daxpy_(&i__, &alpha, &a[(i__ + 1) * a_dim1 + 1], &c__1, &tau[
					1], &c__1);

					/*              Apply the transformation as a rank-2 update:   
					A := A - v * w' - w * v' */

					dsyr2_(uplo, &i__, &c_b14, &a[(i__ + 1) * a_dim1 + 1], &c__1, 
						&tau[1], &c__1, &a[a_offset], lda);

					a[i__ + (i__ + 1) * a_dim1] = e[i__];
			}
			d__[i__ + 1] = a[i__ + 1 + (i__ + 1) * a_dim1];
			tau[i__] = taui;
			/* L10: */
		}
		d__[1] = a[a_dim1 + 1];
	} else {

		/*        Reduce the lower triangle of A */

		i__1 = *n - 1;
		for (i__ = 1; i__ <= i__1; ++i__) {

			/*           Generate elementary reflector H(i) = I - tau * v * v'   
			to annihilate A(i+2:n,i) */

			i__2 = *n - i__;
			/* Computing MIN */
			i__3 = i__ + 2;
			dlarfg_(&i__2, &a[i__ + 1 + i__ * a_dim1], &a[min(i__3,*n) + i__ *
				a_dim1], &c__1, &taui);
			e[i__] = a[i__ + 1 + i__ * a_dim1];

			if (taui != 0.) {

				/*              Apply H(i) from both sides to A(i+1:n,i+1:n) */

				a[i__ + 1 + i__ * a_dim1] = 1.;

				/*              Compute  x := tau * A * v  storing y in TAU(i:n-1) */

				i__2 = *n - i__;
				dsymv_(uplo, &i__2, &taui, &a[i__ + 1 + (i__ + 1) * a_dim1], 
					lda, &a[i__ + 1 + i__ * a_dim1], &c__1, &c_b8, &tau[
						i__], &c__1);

						/*              Compute  w := x - 1/2 * tau * (x'*v) * v */

						i__2 = *n - i__;
						alpha = taui * -.5 * ddot_(&i__2, &tau[i__], &c__1, &a[i__ + 
							1 + i__ * a_dim1], &c__1);
						i__2 = *n - i__;
						daxpy_(&i__2, &alpha, &a[i__ + 1 + i__ * a_dim1], &c__1, &tau[
							i__], &c__1);

							/*              Apply the transformation as a rank-2 update:   
							A := A - v * w' - w * v' */

							i__2 = *n - i__;
							dsyr2_(uplo, &i__2, &c_b14, &a[i__ + 1 + i__ * a_dim1], &c__1, 
								&tau[i__], &c__1, &a[i__ + 1 + (i__ + 1) * a_dim1], 
								lda);

							a[i__ + 1 + i__ * a_dim1] = e[i__];
			}
			d__[i__] = a[i__ + i__ * a_dim1];
			tau[i__] = taui;
			/* L20: */
		}
		d__[*n] = a[*n + *n * a_dim1];
	}

	return 0;

	/*     End of DSYTD2 */

} /* dsytd2_ */


int TRL::dsyr2k_(char *uplo, char *trans, integer_ *n, integer_ *k, 
			doublereal_ *alpha, doublereal_ *a, integer_ *lda, doublereal_ *b, 
			integer_ *ldb, doublereal_ *beta, doublereal_ *c__, integer_ *ldc)
{
	/* System generated locals */
	integer_ a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, 
		i__3;
	/* Local variables */
	static integer_ i__, j, l, info;
	static doublereal_ temp1, temp2;
	//extern logical_ lsame_(char *, char *);
	static integer_ nrowa;
	static logical_ upper;
	//extern /* Subroutine */ int xerbla_(char *, integer_ *);
	/*  Purpose   
	=======   
	DSYR2K  performs one of the symmetric rank 2k operations   
	C := alpha*A*B' + alpha*B*A' + beta*C,   
	or   
	C := alpha*A'*B + alpha*B'*A + beta*C,   
	where  alpha and beta  are scalars, C is an  n by n  symmetric matrix   
	and  A and B  are  n by k  matrices  in the  first  case  and  k by n   
	matrices in the second case.   
	Arguments   
	==========   
	UPLO   - CHARACTER*1.   
	On  entry,   UPLO  specifies  whether  the  upper  or  lower   
	triangular  part  of the  array  C  is to be  referenced  as   
	follows:   
	UPLO = 'U' or 'u'   Only the  upper triangular part of  C   
	is to be referenced.   
	UPLO = 'L' or 'l'   Only the  lower triangular part of  C   
	is to be referenced.   
	Unchanged on exit.   
	TRANS  - CHARACTER*1.   
	On entry,  TRANS  specifies the operation to be performed as   
	follows:   
	TRANS = 'N' or 'n'   C := alpha*A*B' + alpha*B*A' +   
	beta*C.   
	TRANS = 'T' or 't'   C := alpha*A'*B + alpha*B'*A +   
	beta*C.   
	TRANS = 'C' or 'c'   C := alpha*A'*B + alpha*B'*A +   
	beta*C.   
	Unchanged on exit.   
	N      - INTEGER.   
	On entry,  N specifies the order of the matrix C.  N must be   
	at least zero.   
	Unchanged on exit.   
	K      - INTEGER.   
	On entry with  TRANS = 'N' or 'n',  K  specifies  the number   
	of  columns  of the  matrices  A and B,  and on  entry  with   
	TRANS = 'T' or 't' or 'C' or 'c',  K  specifies  the  number   
	of rows of the matrices  A and B.  K must be at least  zero.   
	Unchanged on exit.   
	ALPHA  - DOUBLE PRECISION.   
	On entry, ALPHA specifies the scalar alpha.   
	Unchanged on exit.   
	A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is   
	k  when  TRANS = 'N' or 'n',  and is  n  otherwise.   
	Before entry with  TRANS = 'N' or 'n',  the  leading  n by k   
	part of the array  A  must contain the matrix  A,  otherwise   
	the leading  k by n  part of the array  A  must contain  the   
	matrix A.   
	Unchanged on exit.   
	LDA    - INTEGER.   
	On entry, LDA specifies the first dimension of A as declared   
	in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'   
	then  LDA must be at least  max( 1, n ), otherwise  LDA must   
	be at least  max( 1, k ).   
	Unchanged on exit.   
	B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is   
	k  when  TRANS = 'N' or 'n',  and is  n  otherwise.   
	Before entry with  TRANS = 'N' or 'n',  the  leading  n by k   
	part of the array  B  must contain the matrix  B,  otherwise   
	the leading  k by n  part of the array  B  must contain  the   
	matrix B.   
	Unchanged on exit.   
	LDB    - INTEGER.   
	On entry, LDB specifies the first dimension of B as declared   
	in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'   
	then  LDB must be at least  max( 1, n ), otherwise  LDB must   
	be at least  max( 1, k ).   
	Unchanged on exit.   
	BETA   - DOUBLE PRECISION.   
	On entry, BETA specifies the scalar beta.   
	Unchanged on exit.   
	C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).   
	Before entry  with  UPLO = 'U' or 'u',  the leading  n by n   
	upper triangular part of the array C must contain the upper   
	triangular part  of the  symmetric matrix  and the strictly   
	lower triangular part of C is not referenced.  On exit, the   
	upper triangular part of the array  C is overwritten by the   
	upper triangular part of the updated matrix.   
	Before entry  with  UPLO = 'L' or 'l',  the leading  n by n   
	lower triangular part of the array C must contain the lower   
	triangular part  of the  symmetric matrix  and the strictly   
	upper triangular part of C is not referenced.  On exit, the   
	lower triangular part of the array  C is overwritten by the   
	lower triangular part of the updated matrix.   
	LDC    - INTEGER.   
	On entry, LDC specifies the first dimension of C as declared   
	in  the  calling  (sub)  program.   LDC  must  be  at  least   
	max( 1, n ).   
	Unchanged on exit.   
	Level 3 Blas routine.   
	-- Written on 8-February-1989.   
	Jack Dongarra, Argonne National Laboratory.   
	Iain Duff, AERE Harwell.   
	Jeremy Du Croz, Numerical Algorithms Group Ltd.   
	Sven Hammarling, Numerical Algorithms Group Ltd.   
	Test the input parameters.   
	Parameter adjustments */
	a_dim1 = *lda;
	a_offset = 1 + a_dim1;
	a -= a_offset;
	b_dim1 = *ldb;
	b_offset = 1 + b_dim1;
	b -= b_offset;
	c_dim1 = *ldc;
	c_offset = 1 + c_dim1;
	c__ -= c_offset;
	/* Function Body */
	if (lsame_(trans, "N")) {
		nrowa = *n;
	} else {
		nrowa = *k;
	}
	upper = lsame_(uplo, "U");
	info = 0;
	if (! upper && ! lsame_(uplo, "L")) {
		info = 1;
	} else if (! lsame_(trans, "N") && ! lsame_(trans, 
		"T") && ! lsame_(trans, "C")) {
			info = 2;
	} else if (*n < 0) {
		info = 3;
	} else if (*k < 0) {
		info = 4;
	} else if (*lda < max(1,nrowa)) {
		info = 7;
	} else if (*ldb < max(1,nrowa)) {
		info = 9;
	} else if (*ldc < max(1,*n)) {
		info = 12;
	}
	if (info != 0) {
		xerbla_("DSYR2K", &info);
		return 0;
	}
	/*     Quick return if possible. */
	if (*n == 0 || (*alpha == 0. || *k == 0) && *beta == 1.) {
		return 0;
	}
	/*     And when  alpha.eq.zero. */
	if (*alpha == 0.) {
		if (upper) {
			if (*beta == 0.) {
				i__1 = *n;
				for (j = 1; j <= i__1; ++j) {
					i__2 = j;
					for (i__ = 1; i__ <= i__2; ++i__) {
						c__[i__ + j * c_dim1] = 0.;
						/* L10: */
					}
					/* L20: */
				}
			} else {
				i__1 = *n;
				for (j = 1; j <= i__1; ++j) {
					i__2 = j;
					for (i__ = 1; i__ <= i__2; ++i__) {
						c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
						/* L30: */
					}
					/* L40: */
				}
			}
		} else {
			if (*beta == 0.) {
				i__1 = *n;
				for (j = 1; j <= i__1; ++j) {
					i__2 = *n;
					for (i__ = j; i__ <= i__2; ++i__) {
						c__[i__ + j * c_dim1] = 0.;
						/* L50: */
					}
					/* L60: */
				}
			} else {
				i__1 = *n;
				for (j = 1; j <= i__1; ++j) {
					i__2 = *n;
					for (i__ = j; i__ <= i__2; ++i__) {
						c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
						/* L70: */
					}
					/* L80: */
				}
			}
		}
		return 0;
	}
	/*     Start the operations. */
	if (lsame_(trans, "N")) {
		/*        Form  C := alpha*A*B' + alpha*B*A' + C. */
		if (upper) {
			i__1 = *n;
			for (j = 1; j <= i__1; ++j) {
				if (*beta == 0.) {
					i__2 = j;
					for (i__ = 1; i__ <= i__2; ++i__) {
						c__[i__ + j * c_dim1] = 0.;
						/* L90: */
					}
				} else if (*beta != 1.) {
					i__2 = j;
					for (i__ = 1; i__ <= i__2; ++i__) {
						c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
						/* L100: */
					}
				}
				i__2 = *k;
				for (l = 1; l <= i__2; ++l) {
					if (a[j + l * a_dim1] != 0. || b[j + l * b_dim1] != 0.) {
						temp1 = *alpha * b[j + l * b_dim1];
						temp2 = *alpha * a[j + l * a_dim1];
						i__3 = j;
						for (i__ = 1; i__ <= i__3; ++i__) {
							c__[i__ + j * c_dim1] = c__[i__ + j * c_dim1] + a[
								i__ + l * a_dim1] * temp1 + b[i__ + l * 
									b_dim1] * temp2;
								/* L110: */
						}
					}
					/* L120: */
				}
				/* L130: */
			}
		} else {
			i__1 = *n;
			for (j = 1; j <= i__1; ++j) {
				if (*beta == 0.) {
					i__2 = *n;
					for (i__ = j; i__ <= i__2; ++i__) {
						c__[i__ + j * c_dim1] = 0.;
						/* L140: */
					}
				} else if (*beta != 1.) {
					i__2 = *n;
					for (i__ = j; i__ <= i__2; ++i__) {
						c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
						/* L150: */
					}
				}
				i__2 = *k;
				for (l = 1; l <= i__2; ++l) {
					if (a[j + l * a_dim1] != 0. || b[j + l * b_dim1] != 0.) {
						temp1 = *alpha * b[j + l * b_dim1];
						temp2 = *alpha * a[j + l * a_dim1];
						i__3 = *n;
						for (i__ = j; i__ <= i__3; ++i__) {
							c__[i__ + j * c_dim1] = c__[i__ + j * c_dim1] + a[
								i__ + l * a_dim1] * temp1 + b[i__ + l * 
									b_dim1] * temp2;
								/* L160: */
						}
					}
					/* L170: */
				}
				/* L180: */
			}
		}
	} else {
		/*        Form  C := alpha*A'*B + alpha*B'*A + C. */
		if (upper) {
			i__1 = *n;
			for (j = 1; j <= i__1; ++j) {
				i__2 = j;
				for (i__ = 1; i__ <= i__2; ++i__) {
					temp1 = 0.;
					temp2 = 0.;
					i__3 = *k;
					for (l = 1; l <= i__3; ++l) {
						temp1 += a[l + i__ * a_dim1] * b[l + j * b_dim1];
						temp2 += b[l + i__ * b_dim1] * a[l + j * a_dim1];
						/* L190: */
					}
					if (*beta == 0.) {
						c__[i__ + j * c_dim1] = *alpha * temp1 + *alpha * 
							temp2;
					} else {
						c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1] 
						+ *alpha * temp1 + *alpha * temp2;
					}
					/* L200: */
				}
				/* L210: */
			}
		} else {
			i__1 = *n;
			for (j = 1; j <= i__1; ++j) {
				i__2 = *n;
				for (i__ = j; i__ <= i__2; ++i__) {
					temp1 = 0.;
					temp2 = 0.;
					i__3 = *k;
					for (l = 1; l <= i__3; ++l) {
						temp1 += a[l + i__ * a_dim1] * b[l + j * b_dim1];
						temp2 += b[l + i__ * b_dim1] * a[l + j * a_dim1];
						/* L220: */
					}
					if (*beta == 0.) {
						c__[i__ + j * c_dim1] = *alpha * temp1 + *alpha * 
							temp2;
					} else {
						c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1] 
						+ *alpha * temp1 + *alpha * temp2;
					}
					/* L230: */
				}
				/* L240: */
			}
		}
	}
	return 0;
	/*     End of DSYR2K. */
} /* dsyr2k_ */





/////////
int TRL::dlatrd_(char *uplo, integer_ *n, integer_ *nb, doublereal_ *
			a, integer_ *lda, doublereal_ *e, doublereal_ *tau, doublereal_ *w, 
			integer_ *ldw)
{
	/*  -- LAPACK auxiliary routine (version 3.1) --   
	Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..   
	November 2006   


	Purpose   
	=======   

	DLATRD reduces NB rows and columns of a trl_real_ symmetric matrix A to   
	symmetric tridiagonal form by an orthogonal similarity   
	transformation Q' * A * Q, and returns the matrices V and W which are   
	needed to apply the transformation to the unreduced part of A.   

	If UPLO = 'U', DLATRD reduces the last NB rows and columns of a   
	matrix, of which the upper triangle is supplied;   
	if UPLO = 'L', DLATRD reduces the first NB rows and columns of a   
	matrix, of which the lower triangle is supplied.   

	This is an auxiliary routine called by DSYTRD.   

	Arguments   
	=========   

	UPLO    (input) CHARACTER*1   
	Specifies whether the upper or lower triangular part of the   
	symmetric matrix A is stored:   
	= 'U': Upper triangular   
	= 'L': Lower triangular   

	N       (input) INTEGER   
	The order of the matrix A.   

	NB      (input) INTEGER   
	The number of rows and columns to be reduced.   

	A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
	On entry, the symmetric matrix A.  If UPLO = 'U', the leading   
	n-by-n upper triangular part of A contains the upper   
	triangular part of the matrix A, and the strictly lower   
	triangular part of A is not referenced.  If UPLO = 'L', the   
	leading n-by-n lower triangular part of A contains the lower   
	triangular part of the matrix A, and the strictly upper   
	triangular part of A is not referenced.   
	On exit:   
	if UPLO = 'U', the last NB columns have been reduced to   
	tridiagonal form, with the diagonal elements overwriting   
	the diagonal elements of A; the elements above the diagonal   
	with the array TAU, represent the orthogonal matrix Q as a   
	product of elementary reflectors;   
	if UPLO = 'L', the first NB columns have been reduced to   
	tridiagonal form, with the diagonal elements overwriting   
	the diagonal elements of A; the elements below the diagonal   
	with the array TAU, represent the  orthogonal matrix Q as a   
	product of elementary reflectors.   
	See Further Details.   

	LDA     (input) INTEGER   
	The leading dimension of the array A.  LDA >= (1,N).   

	E       (output) DOUBLE PRECISION array, dimension (N-1)   
	If UPLO = 'U', E(n-nb:n-1) contains the superdiagonal   
	elements of the last NB columns of the reduced matrix;   
	if UPLO = 'L', E(1:nb) contains the subdiagonal elements of   
	the first NB columns of the reduced matrix.   

	TAU     (output) DOUBLE PRECISION array, dimension (N-1)   
	The scalar factors of the elementary reflectors, stored in   
	TAU(n-nb:n-1) if UPLO = 'U', and in TAU(1:nb) if UPLO = 'L'.   
	See Further Details.   

	W       (output) DOUBLE PRECISION array, dimension (LDW,NB)   
	The n-by-nb matrix W required to update the unreduced part   
	of A.   

	LDW     (input) INTEGER   
	The leading dimension of the array W. LDW >= max(1,N).   

	Further Details   
	===============   

	If UPLO = 'U', the matrix Q is represented as a product of elementary   
	reflectors   

	Q = H(n) H(n-1) . . . H(n-nb+1).   

	Each H(i) has the form   

	H(i) = I - tau * v * v'   

	where tau is a trl_real_ scalar, and v is a trl_real_ vector with   
	v(i:n) = 0 and v(i-1) = 1; v(1:i-1) is stored on exit in A(1:i-1,i),   
	and tau in TAU(i-1).   

	If UPLO = 'L', the matrix Q is represented as a product of elementary   
	reflectors   

	Q = H(1) H(2) . . . H(nb).   

	Each H(i) has the form   

	H(i) = I - tau * v * v'   

	where tau is a trl_real_ scalar, and v is a trl_real_ vector with   
	v(1:i) = 0 and v(i+1) = 1; v(i+1:n) is stored on exit in A(i+1:n,i),   
	and tau in TAU(i).   

	The elements of the vectors v together form the n-by-nb matrix V   
	which is needed, with W, to apply the transformation to the unreduced   
	part of the matrix, using a symmetric rank-2k update of the form:   
	A := A - V*W' - W*V'.   

	The contents of A on exit are illustrated by the following examples   
	with n = 5 and nb = 2:   

	if UPLO = 'U':                       if UPLO = 'L':   

	(  a   a   a   v4  v5 )              (  d                  )   
	(      a   a   v4  v5 )              (  1   d              )   
	(          a   1   v5 )              (  v1  1   a          )   
	(              d   1  )              (  v1  v2  a   a      )   
	(                  d  )              (  v1  v2  a   a   a  )   

	where d denotes a diagonal element of the reduced matrix, a denotes   
	an element of the original matrix that is unchanged, and vi denotes   
	an element of the vector defining H(i).   

	=====================================================================   


	Quick return if possible   

	Parameter adjustments */
	/* Table of constant values */
	static doublereal_ c_b5 = -1.;
	static doublereal_ c_b6 = 1.;
	static integer_ c__1 = 1;
	static doublereal_ c_b16 = 0.;

	/* System generated locals */
	integer_ a_dim1, a_offset, w_dim1, w_offset, i__1, i__2, i__3;
	/* Local variables */
	static integer_ i__, iw;
	//extern doublereal_ ddot_(integer_ *, doublereal_ *, integer_ *, doublereal_ *, 
	//	integer_ *);
	static doublereal_ alpha;
	//extern /* Subroutine */ int dscal_(integer_ *, doublereal_ *, doublereal_ *, 
	//	integer_ *);
	//extern logical_ lsame_(char *, char *);
	//extern /* Subroutine */ int dgemv_(char *, integer_ *, integer_ *, 
	//	doublereal_ *, doublereal_ *, integer_ *, doublereal_ *, integer_ *, 
	//	doublereal_ *, doublereal_ *, integer_ *), daxpy_(integer_ *, 
	//	doublereal_ *, doublereal_ *, integer_ *, doublereal_ *, integer_ *), 
	//	dsymv_(char *, integer_ *, doublereal_ *, doublereal_ *, integer_ *, 
	//	doublereal_ *, integer_ *, doublereal_ *, doublereal_ *, integer_ *), dlarfg_(integer_ *, doublereal_ *, doublereal_ *, integer_ *,
	//	doublereal_ *);


	a_dim1 = *lda;
	a_offset = 1 + a_dim1;
	a -= a_offset;
	--e;
	--tau;
	w_dim1 = *ldw;
	w_offset = 1 + w_dim1;
	w -= w_offset;

	/* Function Body */
	if (*n <= 0) {
		return 0;
	}

	if (lsame_(uplo, "U")) {

		/*        Reduce last NB columns of upper triangle */

		i__1 = *n - *nb + 1;
		for (i__ = *n; i__ >= i__1; --i__) {
			iw = i__ - *n + *nb;
			if (i__ < *n) {

				/*              Update A(1:i,i) */

				i__2 = *n - i__;
				dgemv_("No transpose", &i__, &i__2, &c_b5, &a[(i__ + 1) * 
					a_dim1 + 1], lda, &w[i__ + (iw + 1) * w_dim1], ldw, &
					c_b6, &a[i__ * a_dim1 + 1], &c__1);
				i__2 = *n - i__;
				dgemv_("No transpose", &i__, &i__2, &c_b5, &w[(iw + 1) * 
					w_dim1 + 1], ldw, &a[i__ + (i__ + 1) * a_dim1], lda, &
					c_b6, &a[i__ * a_dim1 + 1], &c__1);
			}
			if (i__ > 1) {

				/*              Generate elementary reflector H(i) to annihilate   
				A(1:i-2,i) */

				i__2 = i__ - 1;
				dlarfg_(&i__2, &a[i__ - 1 + i__ * a_dim1], &a[i__ * a_dim1 + 
					1], &c__1, &tau[i__ - 1]);
				e[i__ - 1] = a[i__ - 1 + i__ * a_dim1];
				a[i__ - 1 + i__ * a_dim1] = 1.;

				/*              Compute W(1:i-1,i) */

				i__2 = i__ - 1;
				dsymv_("Upper", &i__2, &c_b6, &a[a_offset], lda, &a[i__ * 
					a_dim1 + 1], &c__1, &c_b16, &w[iw * w_dim1 + 1], &
					c__1);
				if (i__ < *n) {
					i__2 = i__ - 1;
					i__3 = *n - i__;
					dgemv_("Transpose", &i__2, &i__3, &c_b6, &w[(iw + 1) * 
						w_dim1 + 1], ldw, &a[i__ * a_dim1 + 1], &c__1, &
						c_b16, &w[i__ + 1 + iw * w_dim1], &c__1);
					i__2 = i__ - 1;
					i__3 = *n - i__;
					dgemv_("No transpose", &i__2, &i__3, &c_b5, &a[(i__ + 1) *
						a_dim1 + 1], lda, &w[i__ + 1 + iw * w_dim1], &
						c__1, &c_b6, &w[iw * w_dim1 + 1], &c__1);
					i__2 = i__ - 1;
					i__3 = *n - i__;
					dgemv_("Transpose", &i__2, &i__3, &c_b6, &a[(i__ + 1) * 
						a_dim1 + 1], lda, &a[i__ * a_dim1 + 1], &c__1, &
						c_b16, &w[i__ + 1 + iw * w_dim1], &c__1);
					i__2 = i__ - 1;
					i__3 = *n - i__;
					dgemv_("No transpose", &i__2, &i__3, &c_b5, &w[(iw + 1) * 
						w_dim1 + 1], ldw, &w[i__ + 1 + iw * w_dim1], &
						c__1, &c_b6, &w[iw * w_dim1 + 1], &c__1);
				}
				i__2 = i__ - 1;
				dscal_(&i__2, &tau[i__ - 1], &w[iw * w_dim1 + 1], &c__1);
				i__2 = i__ - 1;
				alpha = tau[i__ - 1] * -.5 * ddot_(&i__2, &w[iw * w_dim1 + 1], 
					&c__1, &a[i__ * a_dim1 + 1], &c__1);
				i__2 = i__ - 1;
				daxpy_(&i__2, &alpha, &a[i__ * a_dim1 + 1], &c__1, &w[iw * 
					w_dim1 + 1], &c__1);
			}

			/* L10: */
		}
	} else {

		/*        Reduce first NB columns of lower triangle */

		i__1 = *nb;
		for (i__ = 1; i__ <= i__1; ++i__) {

			/*           Update A(i:n,i) */

			i__2 = *n - i__ + 1;
			i__3 = i__ - 1;
			dgemv_("No transpose", &i__2, &i__3, &c_b5, &a[i__ + a_dim1], lda, 
				&w[i__ + w_dim1], ldw, &c_b6, &a[i__ + i__ * a_dim1], &
				c__1);
			i__2 = *n - i__ + 1;
			i__3 = i__ - 1;
			dgemv_("No transpose", &i__2, &i__3, &c_b5, &w[i__ + w_dim1], ldw, 
				&a[i__ + a_dim1], lda, &c_b6, &a[i__ + i__ * a_dim1], &
				c__1);
			if (i__ < *n) {

				/*              Generate elementary reflector H(i) to annihilate   
				A(i+2:n,i) */

				i__2 = *n - i__;
				/* Computing MIN */
				i__3 = i__ + 2;
				dlarfg_(&i__2, &a[i__ + 1 + i__ * a_dim1], &a[min(i__3,*n) + 
					i__ * a_dim1], &c__1, &tau[i__]);
				e[i__] = a[i__ + 1 + i__ * a_dim1];
				a[i__ + 1 + i__ * a_dim1] = 1.;

				/*              Compute W(i+1:n,i) */

				i__2 = *n - i__;
				dsymv_("Lower", &i__2, &c_b6, &a[i__ + 1 + (i__ + 1) * a_dim1]
				, lda, &a[i__ + 1 + i__ * a_dim1], &c__1, &c_b16, &w[
					i__ + 1 + i__ * w_dim1], &c__1);
					i__2 = *n - i__;
					i__3 = i__ - 1;
					dgemv_("Transpose", &i__2, &i__3, &c_b6, &w[i__ + 1 + w_dim1], 
						ldw, &a[i__ + 1 + i__ * a_dim1], &c__1, &c_b16, &w[
							i__ * w_dim1 + 1], &c__1);
							i__2 = *n - i__;
							i__3 = i__ - 1;
							dgemv_("No transpose", &i__2, &i__3, &c_b5, &a[i__ + 1 + 
								a_dim1], lda, &w[i__ * w_dim1 + 1], &c__1, &c_b6, &w[
									i__ + 1 + i__ * w_dim1], &c__1);
									i__2 = *n - i__;
									i__3 = i__ - 1;
									dgemv_("Transpose", &i__2, &i__3, &c_b6, &a[i__ + 1 + a_dim1], 
										lda, &a[i__ + 1 + i__ * a_dim1], &c__1, &c_b16, &w[
											i__ * w_dim1 + 1], &c__1);
											i__2 = *n - i__;
											i__3 = i__ - 1;
											dgemv_("No transpose", &i__2, &i__3, &c_b5, &w[i__ + 1 + 
												w_dim1], ldw, &w[i__ * w_dim1 + 1], &c__1, &c_b6, &w[
													i__ + 1 + i__ * w_dim1], &c__1);
													i__2 = *n - i__;
													dscal_(&i__2, &tau[i__], &w[i__ + 1 + i__ * w_dim1], &c__1);
													i__2 = *n - i__;
													alpha = tau[i__] * -.5 * ddot_(&i__2, &w[i__ + 1 + i__ * 
														w_dim1], &c__1, &a[i__ + 1 + i__ * a_dim1], &c__1);
													i__2 = *n - i__;
													daxpy_(&i__2, &alpha, &a[i__ + 1 + i__ * a_dim1], &c__1, &w[
														i__ + 1 + i__ * w_dim1], &c__1);
			}

			/* L20: */
		}
	}

	return 0;

	/*     End of DLATRD */

} /* dlatrd_ */



int TRL::dsyr2_(char *uplo, integer_ *n, doublereal_ *alpha, 
		   doublereal_ *x, integer_ *incx, doublereal_ *y, integer_ *incy, 
		   doublereal_ *a, integer_ *lda)
{
	/* System generated locals */
	integer_ a_dim1, a_offset, i__1, i__2;
	/* Local variables */
	static integer_ i__, j, ix, iy, jx, jy, kx, ky, info;
	static doublereal_ temp1, temp2;
	//extern logical_ lsame_(char *, char *);
	//extern /* Subroutine */ int xerbla_(char *, integer_ *);
	/*  Purpose   
	=======   
	DSYR2  performs the symmetric rank 2 operation   
	A := alpha*x*y' + alpha*y*x' + A,   
	where alpha is a scalar, x and y are n element vectors and A is an n   
	by n symmetric matrix.   
	Arguments   
	==========   
	UPLO   - CHARACTER*1.   
	On entry, UPLO specifies whether the upper or lower   
	triangular part of the array A is to be referenced as   
	follows:   
	UPLO = 'U' or 'u'   Only the upper triangular part of A   
	is to be referenced.   
	UPLO = 'L' or 'l'   Only the lower triangular part of A   
	is to be referenced.   
	Unchanged on exit.   
	N      - INTEGER.   
	On entry, N specifies the order of the matrix A.   
	N must be at least zero.   
	Unchanged on exit.   
	ALPHA  - DOUBLE PRECISION.   
	On entry, ALPHA specifies the scalar alpha.   
	Unchanged on exit.   
	X      - DOUBLE PRECISION array of dimension at least   
	( 1 + ( n - 1 )*abs( INCX ) ).   
	Before entry, the incremented array X must contain the n   
	element vector x.   
	Unchanged on exit.   
	INCX   - INTEGER.   
	On entry, INCX specifies the increment for the elements of   
	X. INCX must not be zero.   
	Unchanged on exit.   
	Y      - DOUBLE PRECISION array of dimension at least   
	( 1 + ( n - 1 )*abs( INCY ) ).   
	Before entry, the incremented array Y must contain the n   
	element vector y.   
	Unchanged on exit.   
	INCY   - INTEGER.   
	On entry, INCY specifies the increment for the elements of   
	Y. INCY must not be zero.   
	Unchanged on exit.   
	A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).   
	Before entry with  UPLO = 'U' or 'u', the leading n by n   
	upper triangular part of the array A must contain the upper   
	triangular part of the symmetric matrix and the strictly   
	lower triangular part of A is not referenced. On exit, the   
	upper triangular part of the array A is overwritten by the   
	upper triangular part of the updated matrix.   
	Before entry with UPLO = 'L' or 'l', the leading n by n   
	lower triangular part of the array A must contain the lower   
	triangular part of the symmetric matrix and the strictly   
	upper triangular part of A is not referenced. On exit, the   
	lower triangular part of the array A is overwritten by the   
	lower triangular part of the updated matrix.   
	LDA    - INTEGER.   
	On entry, LDA specifies the first dimension of A as declared   
	in the calling (sub) program. LDA must be at least   
	max( 1, n ).   
	Unchanged on exit.   
	Level 2 Blas routine.   
	-- Written on 22-October-1986.   
	Jack Dongarra, Argonne National Lab.   
	Jeremy Du Croz, Nag Central Office.   
	Sven Hammarling, Nag Central Office.   
	Richard Hanson, Sandia National Labs.   
	Test the input parameters.   
	Parameter adjustments */
	--x;
	--y;
	a_dim1 = *lda;
	a_offset = 1 + a_dim1;
	a -= a_offset;
	/* Function Body */
	info = 0;
	if (! lsame_(uplo, "U") && ! lsame_(uplo, "L")) {
		info = 1;
	} else if (*n < 0) {
		info = 2;
	} else if (*incx == 0) {
		info = 5;
	} else if (*incy == 0) {
		info = 7;
	} else if (*lda < max(1,*n)) {
		info = 9;
	}
	if (info != 0) {
		xerbla_("DSYR2 ", &info);
		return 0;
	}
	/*     Quick return if possible. */
	if (*n == 0 || *alpha == 0.) {
		return 0;
	}
	/*     Set up the start points in X and Y if the increments are not both   
	unity. */
	if (*incx != 1 || *incy != 1) {
		if (*incx > 0) {
			kx = 1;
		} else {
			kx = 1 - (*n - 1) * *incx;
		}
		if (*incy > 0) {
			ky = 1;
		} else {
			ky = 1 - (*n - 1) * *incy;
		}
		jx = kx;
		jy = ky;
	}
	/*     Start the operations. In this version the elements of A are   
	accessed sequentially with one pass through the triangular part   
	of A. */
	if (lsame_(uplo, "U")) {
		/*        Form  A  when A is stored in the upper triangle. */
		if (*incx == 1 && *incy == 1) {
			i__1 = *n;
			for (j = 1; j <= i__1; ++j) {
				if (x[j] != 0. || y[j] != 0.) {
					temp1 = *alpha * y[j];
					temp2 = *alpha * x[j];
					i__2 = j;
					for (i__ = 1; i__ <= i__2; ++i__) {
						a[i__ + j * a_dim1] = a[i__ + j * a_dim1] + x[i__] * 
							temp1 + y[i__] * temp2;
						/* L10: */
					}
				}
				/* L20: */
			}
		} else {
			i__1 = *n;
			for (j = 1; j <= i__1; ++j) {
				if (x[jx] != 0. || y[jy] != 0.) {
					temp1 = *alpha * y[jy];
					temp2 = *alpha * x[jx];
					ix = kx;
					iy = ky;
					i__2 = j;
					for (i__ = 1; i__ <= i__2; ++i__) {
						a[i__ + j * a_dim1] = a[i__ + j * a_dim1] + x[ix] * 
							temp1 + y[iy] * temp2;
						ix += *incx;
						iy += *incy;
						/* L30: */
					}
				}
				jx += *incx;
				jy += *incy;
				/* L40: */
			}
		}
	} else {
		/*        Form  A  when A is stored in the lower triangle. */
		if (*incx == 1 && *incy == 1) {
			i__1 = *n;
			for (j = 1; j <= i__1; ++j) {
				if (x[j] != 0. || y[j] != 0.) {
					temp1 = *alpha * y[j];
					temp2 = *alpha * x[j];
					i__2 = *n;
					for (i__ = j; i__ <= i__2; ++i__) {
						a[i__ + j * a_dim1] = a[i__ + j * a_dim1] + x[i__] * 
							temp1 + y[i__] * temp2;
						/* L50: */
					}
				}
				/* L60: */
			}
		} else {
			i__1 = *n;
			for (j = 1; j <= i__1; ++j) {
				if (x[jx] != 0. || y[jy] != 0.) {
					temp1 = *alpha * y[jy];
					temp2 = *alpha * x[jx];
					ix = jx;
					iy = jy;
					i__2 = *n;
					for (i__ = j; i__ <= i__2; ++i__) {
						a[i__ + j * a_dim1] = a[i__ + j * a_dim1] + x[ix] * 
							temp1 + y[iy] * temp2;
						ix += *incx;
						iy += *incy;
						/* L70: */
					}
				}
				jx += *incx;
				jy += *incy;
				/* L80: */
			}
		}
	}
	return 0;
	/*     End of DSYR2 . */
} /* dsyr2_ */



//

int TRL::dsymv_(char *uplo, integer_ *n, doublereal_ *alpha, 
		   doublereal_ *a, integer_ *lda, doublereal_ *x, integer_ *incx, doublereal_ 
		   *beta, doublereal_ *y, integer_ *incy)
{
	/* System generated locals */
	integer_ a_dim1, a_offset, i__1, i__2;
	/* Local variables */
	static integer_ i__, j, ix, iy, jx, jy, kx, ky, info;
	static doublereal_ temp1, temp2;
	//extern logical_ lsame_(char *, char *);
	//extern /* Subroutine */ int xerbla_(char *, integer_ *);
	/*  Purpose   
	=======   
	DSYMV  performs the matrix-vector  operation   
	y := alpha*A*x + beta*y,   
	where alpha and beta are scalars, x and y are n element vectors and   
	A is an n by n symmetric matrix.   
	Arguments   
	==========   
	UPLO   - CHARACTER*1.   
	On entry, UPLO specifies whether the upper or lower   
	triangular part of the array A is to be referenced as   
	follows:   
	UPLO = 'U' or 'u'   Only the upper triangular part of A   
	is to be referenced.   
	UPLO = 'L' or 'l'   Only the lower triangular part of A   
	is to be referenced.   
	Unchanged on exit.   
	N      - INTEGER.   
	On entry, N specifies the order of the matrix A.   
	N must be at least zero.   
	Unchanged on exit.   
	ALPHA  - DOUBLE PRECISION.   
	On entry, ALPHA specifies the scalar alpha.   
	Unchanged on exit.   
	A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).   
	Before entry with  UPLO = 'U' or 'u', the leading n by n   
	upper triangular part of the array A must contain the upper   
	triangular part of the symmetric matrix and the strictly   
	lower triangular part of A is not referenced.   
	Before entry with UPLO = 'L' or 'l', the leading n by n   
	lower triangular part of the array A must contain the lower   
	triangular part of the symmetric matrix and the strictly   
	upper triangular part of A is not referenced.   
	Unchanged on exit.   
	LDA    - INTEGER.   
	On entry, LDA specifies the first dimension of A as declared   
	in the calling (sub) program. LDA must be at least   
	max( 1, n ).   
	Unchanged on exit.   
	X      - DOUBLE PRECISION array of dimension at least   
	( 1 + ( n - 1 )*abs( INCX ) ).   
	Before entry, the incremented array X must contain the n   
	element vector x.   
	Unchanged on exit.   
	INCX   - INTEGER.   
	On entry, INCX specifies the increment for the elements of   
	X. INCX must not be zero.   
	Unchanged on exit.   
	BETA   - DOUBLE PRECISION.   
	On entry, BETA specifies the scalar beta. When BETA is   
	supplied as zero then Y need not be set on input.   
	Unchanged on exit.   
	Y      - DOUBLE PRECISION array of dimension at least   
	( 1 + ( n - 1 )*abs( INCY ) ).   
	Before entry, the incremented array Y must contain the n   
	element vector y. On exit, Y is overwritten by the updated   
	vector y.   
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
	a_offset = 1 + a_dim1;
	a -= a_offset;
	--x;
	--y;
	/* Function Body */
	info = 0;
	if (! lsame_(uplo, "U") && ! lsame_(uplo, "L")) {
		info = 1;
	} else if (*n < 0) {
		info = 2;
	} else if (*lda < max(1,*n)) {
		info = 5;
	} else if (*incx == 0) {
		info = 7;
	} else if (*incy == 0) {
		info = 10;
	}
	if (info != 0) {
		xerbla_("DSYMV ", &info);
		return 0;
	}
	/*     Quick return if possible. */
	if (*n == 0 || *alpha == 0. && *beta == 1.) {
		return 0;
	}
	/*     Set up the start points in  X  and  Y. */
	if (*incx > 0) {
		kx = 1;
	} else {
		kx = 1 - (*n - 1) * *incx;
	}
	if (*incy > 0) {
		ky = 1;
	} else {
		ky = 1 - (*n - 1) * *incy;
	}
	/*     Start the operations. In this version the elements of A are   
	accessed sequentially with one pass through the triangular part   
	of A.   
	First form  y := beta*y. */
	if (*beta != 1.) {
		if (*incy == 1) {
			if (*beta == 0.) {
				i__1 = *n;
				for (i__ = 1; i__ <= i__1; ++i__) {
					y[i__] = 0.;
					/* L10: */
				}
			} else {
				i__1 = *n;
				for (i__ = 1; i__ <= i__1; ++i__) {
					y[i__] = *beta * y[i__];
					/* L20: */
				}
			}
		} else {
			iy = ky;
			if (*beta == 0.) {
				i__1 = *n;
				for (i__ = 1; i__ <= i__1; ++i__) {
					y[iy] = 0.;
					iy += *incy;
					/* L30: */
				}
			} else {
				i__1 = *n;
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
	if (lsame_(uplo, "U")) {
		/*        Form  y  when A is stored in upper triangle. */
		if (*incx == 1 && *incy == 1) {
			i__1 = *n;
			for (j = 1; j <= i__1; ++j) {
				temp1 = *alpha * x[j];
				temp2 = 0.;
				i__2 = j - 1;
				for (i__ = 1; i__ <= i__2; ++i__) {
					y[i__] += temp1 * a[i__ + j * a_dim1];
					temp2 += a[i__ + j * a_dim1] * x[i__];
					/* L50: */
				}
				y[j] = y[j] + temp1 * a[j + j * a_dim1] + *alpha * temp2;
				/* L60: */
			}
		} else {
			jx = kx;
			jy = ky;
			i__1 = *n;
			for (j = 1; j <= i__1; ++j) {
				temp1 = *alpha * x[jx];
				temp2 = 0.;
				ix = kx;
				iy = ky;
				i__2 = j - 1;
				for (i__ = 1; i__ <= i__2; ++i__) {
					y[iy] += temp1 * a[i__ + j * a_dim1];
					temp2 += a[i__ + j * a_dim1] * x[ix];
					ix += *incx;
					iy += *incy;
					/* L70: */
				}
				y[jy] = y[jy] + temp1 * a[j + j * a_dim1] + *alpha * temp2;
				jx += *incx;
				jy += *incy;
				/* L80: */
			}
		}
	} else {
		/*        Form  y  when A is stored in lower triangle. */
		if (*incx == 1 && *incy == 1) {
			i__1 = *n;
			for (j = 1; j <= i__1; ++j) {
				temp1 = *alpha * x[j];
				temp2 = 0.;
				y[j] += temp1 * a[j + j * a_dim1];
				i__2 = *n;
				for (i__ = j + 1; i__ <= i__2; ++i__) {
					y[i__] += temp1 * a[i__ + j * a_dim1];
					temp2 += a[i__ + j * a_dim1] * x[i__];
					/* L90: */
				}
				y[j] += *alpha * temp2;
				/* L100: */
			}
		} else {
			jx = kx;
			jy = ky;
			i__1 = *n;
			for (j = 1; j <= i__1; ++j) {
				temp1 = *alpha * x[jx];
				temp2 = 0.;
				y[jy] += temp1 * a[j + j * a_dim1];
				ix = jx;
				iy = jy;
				i__2 = *n;
				for (i__ = j + 1; i__ <= i__2; ++i__) {
					ix += *incx;
					iy += *incy;
					y[iy] += temp1 * a[i__ + j * a_dim1];
					temp2 += a[i__ + j * a_dim1] * x[ix];
					/* L110: */
				}
				y[jy] += *alpha * temp2;
				jx += *incx;
				jy += *incy;
				/* L120: */
			}
		}
	}
	return 0;
	/*     End of DSYMV . */
} /* dsymv_ */




//

int TRL::dlarfg_(integer_ *n, doublereal_ *alpha, doublereal_ *x, 
			integer_ *incx, doublereal_ *tau)
{
	/*  -- LAPACK auxiliary routine (version 3.1) --   
	Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..   
	November 2006   


	Purpose   
	=======   

	DLARFG generates a trl_real_ elementary reflector H of order n, such   
	that   

	H * ( alpha ) = ( beta ),   H' * H = I.   
	(   x   )   (   0  )   

	where alpha and beta are scalars, and x is an (n-1)-element trl_real_   
	vector. H is represented in the form   

	H = I - tau * ( 1 ) * ( 1 v' ) ,   
	( v )   

	where tau is a trl_real_ scalar and v is a trl_real_ (n-1)-element   
	vector.   

	If the elements of x are all zero, then tau = 0 and H is taken to be   
	the unit matrix.   

	Otherwise  1 <= tau <= 2.   

	Arguments   
	=========   

	N       (input) INTEGER   
	The order of the elementary reflector.   

	ALPHA   (input/output) DOUBLE PRECISION   
	On entry, the value alpha.   
	On exit, it is overwritten with the value beta.   

	X       (input/output) DOUBLE PRECISION array, dimension   
	(1+(N-2)*abs(INCX))   
	On entry, the vector x.   
	On exit, it is overwritten with the vector v.   

	INCX    (input) INTEGER   
	The increment between elements of X. INCX > 0.   

	TAU     (output) DOUBLE PRECISION   
	The value tau.   

	=====================================================================   


	Parameter adjustments */
	/* System generated locals */
	integer_ i__1;
	doublereal_ d__1;
	/* Builtin functions */
	//double d_sign(doublereal_ *, doublereal_ *);
	/* Local variables */
	static integer_ j, knt;
	static doublereal_ beta;
	//extern doublereal_ dnrm2_(integer_ *, doublereal_ *, integer_ *);
	//extern /* Subroutine */ int dscal_(integer_ *, doublereal_ *, doublereal_ *, 
	//	integer_ *);
	static doublereal_ xnorm;
	//extern doublereal_ dlapy2_(doublereal_ *, doublereal_ *), dlamch_(char *);
	static doublereal_ safmin, rsafmn;

	--x;

	/* Function Body */
	if (*n <= 1) {
		*tau = 0.;
		return 0;
	}

	i__1 = *n - 1;
	xnorm = dnrm2_(&i__1, &x[1], incx);

	if (xnorm == 0.) {

		/*        H  =  I */

		*tau = 0.;
	} else {

		/*        general case */

		d__1 = dlapy2_(alpha, &xnorm);
		beta = -d_sign(&d__1, alpha);
		safmin = dlamch_("S") / dlamch_("E");
		if (abs(beta) < safmin) {

			/*           XNORM, BETA may be inaccurate; scale X and recompute them */

			rsafmn = 1. / safmin;
			knt = 0;
L10:
			++knt;
			i__1 = *n - 1;
			dscal_(&i__1, &rsafmn, &x[1], incx);
			beta *= rsafmn;
			*alpha *= rsafmn;
			if (abs(beta) < safmin) {
				goto L10;
			}

			/*           New BETA is at most 1, at least SAFMIN */

			i__1 = *n - 1;
			xnorm = dnrm2_(&i__1, &x[1], incx);
			d__1 = dlapy2_(alpha, &xnorm);
			beta = -d_sign(&d__1, alpha);
			*tau = (beta - *alpha) / beta;
			i__1 = *n - 1;
			d__1 = 1. / (*alpha - beta);
			dscal_(&i__1, &d__1, &x[1], incx);

			/*           If ALPHA is subnormal, it may lose relative accuracy */

			*alpha = beta;
			i__1 = knt;
			for (j = 1; j <= i__1; ++j) {
				*alpha *= safmin;
				/* L20: */
			}
		} else {
			*tau = (beta - *alpha) / beta;
			i__1 = *n - 1;
			d__1 = 1. / (*alpha - beta);
			dscal_(&i__1, &d__1, &x[1], incx);
			*alpha = beta;
		}
	}

	return 0;

	/*     End of DLARFG */

} /* dlarfg_ */




//




















