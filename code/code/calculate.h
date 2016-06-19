#ifndef _CALCULATE_H_
#define _CALCULATE_H_

#include <cmath>
#include <limits>
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

#include <Eigen/Dense> //for matrix calculation
#include <boost/math/distributions/non_central_chi_squared.hpp> //for non-central chi squared
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/thread.hpp>   
#include <boost/date_time.hpp>  
#include <boost/bind.hpp>
#include <boost/ref.hpp>

//#define _DEBUG_HWU_WEI_

typedef Eigen::MatrixXd GTmat_;

typedef Eigen::VectorXd GTvec_;

typedef Eigen::ArrayXXd GTarr_;//two dims array

typedef Eigen::MatrixXf GTmatF_;

typedef Eigen::VectorXf GTvecF_;

typedef Eigen::ArrayXXf GTarrF_;//two dims array

typedef Eigen::MatrixXi GTimat_;



namespace Stat_fuc
{
	template<class T>
	inline T SQR(const T a) {return a*a;}

	double auc_frm_LR(vector<double> & LR, vector<int> & n_patient, vector<int> & n_health);
	double auc_frm_LR(vector<double> & LR, vector<int> & n_patient, vector<int> & n_health, double & var);
	double auc_frm_LR(vector<double> & indi_LR, vector<bool> & if_disease);
	double auc_frm_LR(vector<double> & indi_LR, vector<bool> & if_disease, double & var);
	double auc_frm_LR(vector<double> & indi_LR, vector< vector<int> > & compare_struct);
	double auc_frm_LR(vector<double> & indi_LR, vector< vector<int> > & compare_struct, double & var);

	void score_frm_LR(vector<double> & indi_LR, vector< vector<int> > & compare_struct, vector<double> & fam_score);


	void score_frm_LR(vector<double> & LR, vector< vector<int> > & patient_idx, vector< vector<int> > & health_idx, int tnp, int tnh, vector<double> & score);
	void score_frm_LR(vector<double> & LR, vector<bool> & if_disease, vector<double> & score);
	double auc_frm_score(vector<double> & score, vector<bool> & if_disease);
	double auc_frm_score(vector<double> & score, vector<bool> & if_disease, double & var);
	void aucVarCov_frm_score(vector<double> & pre_score, double pre_auc, vector<double> & score, vector<bool> & if_disease, double & auc, double & var, double & cov);

	int max_index(vector<double> &arr);
	int first_peak(vector<double> &arr);
	void btstrp_sampling(int n, vector<int> & rst, int & seed);

	double chi_sqrq(double k_sqr, double df);//return q_value, tail value(lower tail)
	double chi_sqrp(double k_sqr, double df);//return p_value
	double std_norm_p(double z);//return p_value, when z is positive, it's Cumulative distribution function
	double std_norm_p1(double z);//return p_value, higher one tail
	double std_norm_q(double p);//return q value
	double P_gamma(double x, double shape, double scale);
	double gammq(const double a, const double x);
	double gammp(const double a, const double x);
	void gser(double &gamser, const double a, const double x, double &gln);
	double gammln(const double xx);
	void gcf(double &gammcf, const double a, const double x, double &gln);

	void sampling(vector<int> & ori, int size, vector<int> & rst, int & seed);
	void ran_devide(vector<int> & ori, int fold, vector< vector<int> > & rst, int & seed, vector< vector<int> > & rst2);//random devide to equal size
	double ran1(int &idum);//seed should be negative value for the first time
	double ran2(int &idum);//seed should be negative value for the first time

	bool emprical_CI(vector<double> & arr, double & Upper, double & Lower, double pct);
	double median(vector<double> & arr);
	double mean(vector<double> & arr);
	double sum(vector<double> & arr);
	void indexx(vector<double> &arr, vector<int> & indx);


	inline void SWAP(int &a, int &b)
	{
		int temp=a;
		a=b;
		b=temp;
	}

	//get rank
	void crank(vector<double> &w);
	void normalQuantile(vector<double> & ori, vector<double> & rst);
	void unifQuantile(vector<double> & ori, vector<double> & rst);
	void RanNormal(double mu, double sig, int size, int & seed, vector<double> & rst);
	

	struct Erf {
		static const int ncof=28;
		static const double cof[28];

		inline double SQR(double x){
			return x*x;
		}

		inline double erf(double x) {
			if (x >=0.) return 1.0 - erfccheb(x);
			else return erfccheb(-x) - 1.0;
		}

		inline double erfc(double x) {
			if (x >= 0.) return erfccheb(x);
			else return 2.0 - erfccheb(-x);
		}

		double erfccheb(double z){
			int j;
			double t,ty,tmp,d=0.,dd=0.;
			if (z < 0.) cerr<<"erfccheb requires nonnegative argument";
			t = 2./(2.+z);
			ty = 4.*t - 2.;
			for (j=ncof-1;j>0;j--) {
				tmp = d;
				d = ty*d - dd + cof[j];
				dd = tmp;
			}
			return t*exp(-z*z + 0.5*(cof[0] + ty*d) - dd);
		}

		double inverfc(double p) {
			double x,err,t,pp;
			if (p >= 2.0) return -100.;
			if (p <= 0.0) return 100.;
			pp = (p < 1.0)? p : 2. - p;
			t = sqrt(-2.*log(pp/2.));
			x = -0.70711*((2.30753+t*0.27061)/(1.+t*(0.99229+t*0.04481)) - t);
			for (int j=0;j<2;j++) {
				err = erfc(x) - pp;
				x += err/(1.12837916709551257*exp(-(x*x))-x*err);
			}
			return (p < 1.0? x : -x);
		}

		inline double inverf(double p) {return inverfc(1.-p);}

		double erfcc(const double x)
		{
			double t,z=fabs(x),ans;
			t=2./(2.+z);
			ans=t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
				t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
				t*(-0.82215223+t*0.17087277)))))))));
			return (x >= 0.0 ? ans : 2.0-ans);
		}


	};


	struct Normaldist : Erf {
		double mu, sig;
		Normaldist(double mmu = 0., double ssig = 1.) : mu(mmu), sig(ssig) {
			if (sig <= 0.) cerr<<"bad sig in Normaldist";
		}
		double p(double x) {
			return (0.398942280401432678/sig)*exp(-0.5*SQR((x-mu)/sig));
		}
		double cdf(double x) {
			return 0.5*erfc(-0.707106781186547524*(x-mu)/sig);
		}
		double invcdf(double p) {
			if (p <= 0. || p >= 1.) cerr<<"bad p in Normaldist";
			return -1.41421356237309505*sig*inverfc(2.*p)+mu;
		}
	};







	///////////////
	struct Gauleg18 {
		static const int ngau = 18;
		static const double y[18];
		static const double w[18];

		inline double SIGN(double a, double b)
		{
			return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
		}

		inline double MAX(double a, double b)
		{
			return b > a ? (b) : (a);
		}

		inline double MIN(double a, double b) 
		{
			return b < a ? (b) : (a);
		}

		inline double SQR(double a) 
		{
			return a*a;
		}

		inline int MAX(int a, int b)
		{
			return b > a ? (b) : (a);
		}

		inline int MIN(int a, int b) 
		{
			return b < a ? (b) : (a);
		}

	};
	
	struct Gamma : Gauleg18 {
		static const int ASWITCH=100;
		static const double EPS;
		static const double FPMIN;
		double gln;

		double gammp(const double a, const double x) {
			if (x < 0.0 || a <= 0.0) throw("bad args in gammp");
			if (x == 0.0) return 0.0;
			else if ((int)a >= ASWITCH) return gammpapprox(a,x,1);
			else if (x < a+1.0) return gser(a,x);
			else return 1.0-gcf(a,x);
		}

		double gammq(const double a, const double x) {
			if (x < 0.0 || a <= 0.0) throw("bad args in gammq");
			if (x == 0.0) return 1.0;
			else if ((int)a >= ASWITCH) return gammpapprox(a,x,0);
			else if (x < a+1.0) return 1.0-gser(a,x);
			else return gcf(a,x);
		}

		double gser(const double a, const double x) {
			double sum,del,ap;
			gln=gammln(a);
			ap=a;
			del=sum=1.0/a;
			for (;;) {
				++ap;
				del *= x/ap;
				sum += del;
				if (fabs(del) < fabs(sum)*EPS) {
					return sum*exp(-x+a*log(x)-gln);
				}
			}
		}

		double gcf(const double a, const double x) {
			int i;
			double an,b,c,d,del,h;
			gln=gammln(a);
			b=x+1.0-a;
			c=1.0/FPMIN;
			d=1.0/b;
			h=d;
			for (i=1;;i++) {
				an = -i*(i-a);
				b += 2.0;
				d=an*d+b;
				if (fabs(d) < FPMIN) d=FPMIN;
				c=b+an/c;
				if (fabs(c) < FPMIN) c=FPMIN;
				d=1.0/d;
				del=d*c;
				h *= del;
				if (fabs(del-1.0) <= EPS) break;
			}
			return exp(-x+a*log(x)-gln)*h;
		}

		double gammpapprox(double a, double x, int psig) {
			int j;
			double xu,t,sum,ans;
			double a1 = a-1.0, lna1 = log(a1), sqrta1 = sqrt(a1);
			gln = gammln(a);
			if (x > a1) xu = MAX(a1 + 11.5*sqrta1, x + 6.0*sqrta1);
			else xu = MAX(0.,MIN(a1 - 7.5*sqrta1, x - 5.0*sqrta1));
			sum = 0;
			for (j=0;j<ngau;j++) {
				t = x + (xu-x)*y[j];
				sum += w[j]*exp(-(t-a1)+a1*(log(t)-lna1));
			}
			ans = sum*(xu-x)*exp(a1*(lna1-1.)-gln);
			return (psig?(ans>0.0? 1.0-ans:-ans):(ans>=0.0? ans:1.0+ans));
		}

		double invgammp(double p, double a);

	};


	
	struct Beta : Gauleg18 {
		static const int SWITCH=3000;
		static const double EPS, FPMIN;

		double betai(const double a, const double b, const double x) {
			double bt;
			if (a <= 0.0 || b <= 0.0) throw("Bad a or b in routine betai");
			if (x < 0.0 || x > 1.0) throw("Bad x in routine betai");
			if (x == 0.0 || x == 1.0) return x;
			if (a > SWITCH && b > SWITCH) return betaiapprox(a,b,x);
			bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x));
			if (x < (a+1.0)/(a+b+2.0)) return bt*betacf(a,b,x)/a;
			else return 1.0-bt*betacf(b,a,1.0-x)/b;
		}

		double betacf(const double a, const double b, const double x) {
			int m,m2;
			double aa,c,d,del,h,qab,qam,qap;
			qab=a+b;
			qap=a+1.0;
			qam=a-1.0;
			c=1.0;
			d=1.0-qab*x/qap;
			if (fabs(d) < FPMIN) d=FPMIN;
			d=1.0/d;
			h=d;
			for (m=1;m<10000;m++) {
				m2=2*m;
				aa=m*(b-m)*x/((qam+m2)*(a+m2));
				d=1.0+aa*d;
				if (fabs(d) < FPMIN) d=FPMIN;
				c=1.0+aa/c;
				if (fabs(c) < FPMIN) c=FPMIN;
				d=1.0/d;
				h *= d*c;
				aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
				d=1.0+aa*d;
				if (fabs(d) < FPMIN) d=FPMIN;
				c=1.0+aa/c;
				if (fabs(c) < FPMIN) c=FPMIN;
				d=1.0/d;
				del=d*c;
				h *= del;
				if (fabs(del-1.0) <= EPS) break;
			}
			return h;
		}

		double betaiapprox(double a, double b, double x) {
			int j;
			double xu,t,sum,ans;
			double a1 = a-1.0, b1 = b-1.0, mu = a/(a+b);
			double lnmu=log(mu),lnmuc=log(1.-mu);
			t = sqrt(a*b/(SQR(a+b)*(a+b+1.0)));
			if (x > a/(a+b)) {
				if (x >= 1.0) return 1.0;
				xu = MIN(1.,MAX(mu + 10.*t, x + 5.0*t));
			} else {
				if (x <= 0.0) return 0.0;
				xu = MAX(0.,MIN(mu - 10.*t, x - 5.0*t));
			}
			sum = 0;
			for (j=0;j<18;j++) {
				t = x + (xu-x)*y[j];
				sum += w[j]*exp(a1*(log(t)-lnmu)+b1*(log(1-t)-lnmuc));
			}
			ans = sum*(xu-x)*exp(a1*lnmu-gammln(a)+b1*lnmuc-gammln(b)+gammln(a+b));
			return ans>0.0? 1.0-ans : -ans;
		}

		double invbetai(double p, double a, double b) {
			const double EPS = 1.e-8;
			double pp,t,u,err,x,al,h,w,afac,a1=a-1.,b1=b-1.;
			int j;
			if (p <= 0.) return 0.;
			else if (p >= 1.) return 1.;
			else if (a >= 1. && b >= 1.) {
				pp = (p < 0.5)? p : 1. - p;
				t = sqrt(-2.*log(pp));
				x = (2.30753+t*0.27061)/(1.+t*(0.99229+t*0.04481)) - t;
				if (p < 0.5) x = -x;
				al = (SQR(x)-3.)/6.;
				h = 2./(1./(2.*a-1.)+1./(2.*b-1.));
				w = (x*sqrt(al+h)/h)-(1./(2.*b-1)-1./(2.*a-1.))*(al+5./6.-2./(3.*h));
				x = a/(a+b*exp(2.*w));
			} else {
				double lna = log(a/(a+b)), lnb = log(b/(a+b));
				t = exp(a*lna)/a;
				u = exp(b*lnb)/b;
				w = t + u;
				if (p < t/w) x = pow(a*w*p,1./a);
				else x = 1. - pow(b*w*(1.-p),1./b);
			}
			afac = -gammln(a)-gammln(b)+gammln(a+b);
			for (j=0;j<10;j++) {
				if (x == 0. || x == 1.) return x;
				err = betai(a,b,x) - p;
				t = exp(a1*log(x)+b1*log(1.-x) + afac);
				u = err/t;
				x -= (t = u/(1.-0.5*MIN(1.,u*(a1/x - b1/(1.-x)))));
				if (x <= 0.) x = 0.5*(x + t);
				if (x >= 1.) x = 0.5*(x + t + 1.);
				if (fabs(t) < EPS*x && j > 0) break;
			}
			return x;
		}

	};

	struct Gammadist : Gamma {
		double alph, bet, fac;
		Gammadist(double aalph, double bbet = 1.) : alph(aalph), bet(bbet) {
			if (alph <= 0. || bet <= 0.) throw("bad alph,bet in Gammadist");
			fac = alph*log(bet)-gammln(alph);
		}
		double p(double x) {
			if (x <= 0.) throw("bad x in Gammadist");
			return exp(-bet*x+(alph-1.)*log(x)+fac);
		}
		double cdf(double x) {
			if (x < 0.) throw("bad x in Gammadist");
			return gammp(alph,bet*x);
		}
		double invcdf(double p) {
			if (p < 0. || p >= 1.) throw("bad p in Gammadist");
			return invgammp(p,alph)/bet;
		}
	};
	struct Betadist : Beta {
		double alph, bet, fac;
		Betadist(double aalph, double bbet) : alph(aalph), bet(bbet) {
			if (alph <= 0. || bet <= 0.) throw("bad alph,bet in Betadist");
			fac = gammln(alph+bet)-gammln(alph)-gammln(bet);
		}
		double p(double x) {
			if (x <= 0. || x >= 1.) throw("bad x in Betadist");
			return exp((alph-1.)*log(x)+(bet-1.)*log(1.-x)+fac);
		}
		double cdf(double x) {
			if (x < 0. || x > 1.) throw("bad x in Betadist");
			return betai(alph,bet,x);
		}
		double invcdf(double p) {
			if (p < 0. || p > 1.) throw("bad p in Betadist");
			return invbetai(p,alph,bet);
		}
	};
	struct Studenttdist : Beta {
		double nu, mu, sig, np, fac;
		Studenttdist(double nnu, double mmu = 0., double ssig = 1.)
			: nu(nnu), mu(mmu), sig(ssig) {
				if (sig <= 0. || nu <= 0.) throw("bad sig,nu in Studentdist");
				np = 0.5*(nu + 1.);
				fac = gammln(np)-gammln(0.5*nu);
		}
		double p(double t) {
			return exp(-np*log(1.+SQR((t-mu)/sig)/nu)+fac)
				/(sqrt(3.14159265358979324*nu)*sig);
		}
		double cdf(double t) {
			double p = 0.5*betai(0.5*nu, 0.5, nu/(nu+SQR((t-mu)/sig)));
			if (t >= mu) return 1. - p;
			else return p;
		}
		double invcdf(double p) {
			if (p <= 0. || p >= 1.) throw("bad p in Studentdist");
			double x = invbetai(2.*MIN(p,1.-p), 0.5*nu, 0.5);
			x = sig*sqrt(nu*(1.-x)/x);
			return (p >= 0.5? mu+x : mu-x);
		}
		double aa(double t) {
			if (t < 0.) throw("bad t in Studentdist");
			return 1.-betai(0.5*nu, 0.5, nu/(nu+SQR(t)));
		}
		double invaa(double p) {
			if (p < 0. || p >= 1.) throw("bad p in Studentdist");
			double x = invbetai(1.-p, 0.5*nu, 0.5);
			return sqrt(nu*(1.-x)/x);
		}
	};
	struct Poissondist : Gamma {
		double lam;
		Poissondist(double llam) : lam(llam) {
			if (lam <= 0.) throw("bad lam in Poissondist");	
		}
		double p(int n) {
			if (n < 0) throw("bad n in Poissondist");
			return exp(-lam + n*log(lam) - gammln(n+1.));
		}
		double cdf(int n) {
			if (n < 0) throw("bad n in Poissondist");
			if (n == 0) return 0.;
			return gammq((double)n,lam);
		}
		int invcdf(double p) {
			int n,nl,nu,inc=1;
			if (p <= 0. || p >= 1.) throw("bad p in Poissondist");
			if (p < exp(-lam)) return 0;
			n = (int)MAX(sqrt(lam),5.);
			if (p < cdf(n)) {
				do {
					n = MAX(n-inc,0);
					inc *= 2;
				} while (p < cdf(n));
				nl = n; nu = n + inc/2;
			} else {
				do {
					n += inc;
					inc *= 2;
				} while (p > cdf(n));
				nu = n; nl = n - inc/2;
			}
			while (nu-nl>1) {
				n = (nl+nu)/2;
				if (p < cdf(n)) nu = n;
				else nl = n;
			}
			return nl;
		}
	};
	struct Binomialdist : Beta {
		int n;
		double pe, fac;
		Binomialdist(int nn, double ppe) : n(nn), pe(ppe) {
			if (n <= 0 || pe <= 0. || pe >= 1.) throw("bad args in Binomialdist");
			fac = gammln(n+1.);
		}
		double p(int k) {
			if (k < 0) throw("bad k in Binomialdist");
			if (k > n) return 0.;
			return exp(k*log(pe)+(n-k)*log(1.-pe)
				+fac-gammln(k+1.)-gammln(n-k+1.));
		}
		double cdf(int k) {
			if (k < 0) throw("bad k in Binomialdist");
			if (k == 0) return 0.;
			if (k > n) return 1.;
			return 1. - betai((double)k,n-k+1.,pe);
		}
		int invcdf(double p) {
			int k,kl,ku,inc=1;
			if (p <= 0. || p >= 1.) throw("bad p in Binomialdist");
			k = MAX(0,MIN(n,(int)(n*pe)));
			if (p < cdf(k)) {
				do {
					k = MAX(k-inc,0);
					inc *= 2;
				} while (p < cdf(k));
				kl = k; ku = k + inc/2;
			} else {
				do {
					k = MIN(k+inc,n+1);
					inc *= 2;
				} while (p > cdf(k));
				ku = k; kl = k - inc/2;
			}
			while (ku-kl>1) {
				k = (kl+ku)/2;
				if (p < cdf(k)) ku = k;
				else kl = k;
			}
			return kl;
		}
	};
	struct Chisqdist : Gamma {
		double nu,fac;
		Chisqdist(double nnu) : nu(nnu) {
			if (nu <= 0.) throw("bad nu in Chisqdist");
			fac = 0.693147180559945309*(0.5*nu)+gammln(0.5*nu);
		}
		double p(double x2) {
			if (x2 <= 0.) throw("bad x2 in Chisqdist");
			return exp(-0.5*(x2-(nu-2.)*log(x2))-fac);
		}
		double cdf(double x2) {
			if (x2 < 0.) throw("bad x2 in Chisqdist");
			return gammp(0.5*nu,0.5*x2);
		}
		double invcdf(double p) {
			if (p < 0. || p >= 1.) throw("bad p in Chisqdist");
			return 2.*invgammp(p,0.5*nu);
		}
	};
	struct Fdist : Beta {
		double nu1,nu2;
		double fac;
		Fdist(double nnu1, double nnu2) : nu1(nnu1), nu2(nnu2) {
			if (nu1 <= 0. || nu2 <= 0.) throw("bad nu1,nu2 in Fdist");
			fac = 0.5*(nu1*log(nu1)+nu2*log(nu2))+gammln(0.5*(nu1+nu2))
				-gammln(0.5*nu1)-gammln(0.5*nu2);
		}
		double p(double f) {
			if (f <= 0.) throw("bad f in Fdist");
			return exp((0.5*nu1-1.)*log(f)-0.5*(nu1+nu2)*log(nu2+nu1*f)+fac);
		}
		double cdf(double f) {
			if (f < 0.) throw("bad f in Fdist");
			return betai(0.5*nu1,0.5*nu2,nu1*f/(nu2+nu1*f));
		}
		double invcdf(double p) {
			if (p <= 0. || p >= 1.) throw("bad p in Fdist");
			double x = invbetai(p,0.5*nu1,0.5*nu2);
			return nu2*x/(nu1*(1.-x));
		}
	};




	static Gamma cal_GAMMA;

	struct Ran {
		unsigned long long int u,v,w;
		Ran(unsigned long long int j) : v(4101842887655102017LL), w(1) {
			u = j ^ v; int64();
			v = u; int64();
			w = v; int64();
		}
		inline unsigned long long int int64() {
			u = u * 2862933555777941757LL + 7046029254386353087LL;
			v ^= v >> 17; v ^= v << 31; v ^= v >> 8;
			w = 4294957665U*(w & 0xffffffff) + (w >> 32);
			unsigned long long int x = u ^ (u << 21); x ^= x >> 35; x ^= x << 4;
			return (x + v) ^ w;
		}
		inline double doub() { return 5.42101086242752217E-20 * int64(); }
		inline unsigned int int32() { return (unsigned int)int64(); }
	};
	struct Ranq1 {
		unsigned long long int v;
		Ranq1(unsigned long long int j) : v(4101842887655102017LL) {
			v ^= j;
			v = int64();
		}
		inline unsigned long long int int64() {
			v ^= v >> 21; v ^= v << 35; v ^= v >> 4;
			return v * 2685821657736338717LL;
		}
		inline double doub() { return 5.42101086242752217E-20 * int64(); }
		inline unsigned int int32() { return (unsigned int)int64(); }
	};
	struct Ranq2 {
		unsigned long long int v,w;
		Ranq2(unsigned long long int j) : v(4101842887655102017LL), w(1) {
			v ^= j;
			w = int64();
			v = int64();
		}
		inline unsigned long long int int64() {
			v ^= v >> 17; v ^= v << 31; v ^= v >> 8;
			w = 4294957665U*(w & 0xffffffff) + (w >> 32);
			return v ^ w;
		}
		inline double doub() { return 5.42101086242752217E-20 * int64(); }
		inline unsigned int int32() { return (unsigned int)int64(); }
	};
	struct Ranhash {
		inline unsigned long long int int64(unsigned long long int u) {
			unsigned long long int v = u * 3935559000370003845LL + 2691343689449507681LL;
			v ^= v >> 21; v ^= v << 37; v ^= v >> 4;
			v *= 4768777513237032717LL;
			v ^= v << 20; v ^= v >> 41; v ^= v << 5;
			return  v;
		}
		inline unsigned int int32(unsigned long long int u)
		{ return (unsigned int)(int64(u) & 0xffffffff) ; }
		inline double doub(unsigned long long int u)
		{ return 5.42101086242752217E-20 * int64(u); }
	};
	struct Ranbyte {
		int s[256],i,j,ss;
		unsigned int v;
		Ranbyte(int u) {
			v = 2244614371U ^ u;
			for (i=0; i<256; i++) {s[i] = i;}
			for (j=0, i=0; i<256; i++) {
				ss = s[i];
				j = (j + ss + (v >> 24)) & 0xff;
				s[i] = s[j]; s[j] = ss;
				v = (v << 24) | (v >> 8);
			}
			i = j = 0;
			for (int k=0; k<256; k++) int8();
		}
		inline unsigned char int8() {
			i = (i+1) & 0xff;
			ss = s[i];
			j = (j+ss) & 0xff;
			s[i] = s[j]; s[j] = ss;
			return (unsigned char)(s[(s[i]+s[j]) & 0xff]);
		}
		unsigned int int32() {
			v = 0;
			for (int k=0; k<4; k++) {
				i = (i+1) & 0xff;
				ss = s[i];
				j = (j+ss) & 0xff;
				s[i] = s[j]; s[j] = ss;
				v = (v << 8) | s[(s[i]+s[j]) & 0xff];
			}
			return v;
		}
		double doub() {
			return 2.32830643653869629E-10 * ( int32() +
				2.32830643653869629E-10 * int32() );
		}
	};
	struct Ranfib {
		double dtab[55], dd;
		int inext, inextp;
		Ranfib(unsigned long long int j) : inext(0), inextp(31) {
			Ranq1 init(j);
			for (int k=0; k<55; k++) dtab[k] = init.doub();
		}
		double doub() {
			if (++inext == 55) inext = 0;
			if (++inextp == 55) inextp = 0;
			dd = dtab[inext] - dtab[inextp];
			if (dd < 0) dd += 1.0;
			return (dtab[inext] = dd);
		}
		inline unsigned long int32()
		{ return (unsigned long)(doub() * 4294967295.0);}
	};
	struct Ranlim32 {
		unsigned int u,v,w1,w2;
		Ranlim32(unsigned int j) : v(2244614371U), w1(521288629U), w2(362436069U) {
			u = j ^ v; int32();
			v = u; int32();
		}
		inline unsigned int int32() {
			u = u * 2891336453U + 1640531513U;
			v ^= v >> 13; v ^= v << 17; v ^= v >> 5;
			w1 = 33378 * (w1 & 0xffff) + (w1 >> 16);
			w2 = 57225 * (w2 & 0xffff) + (w2 >> 16);
			unsigned int x = u ^ (u << 9); x ^= x >> 17; x ^= x << 6;
			unsigned int y = w1 ^ (w1 << 17); y ^= y >> 15; y ^= y << 5;
			return (x + v) ^ (y + w2);
		}
		inline double doub() { return 2.32830643653869629E-10 * int32(); }
		inline double truedoub() {
			return 2.32830643653869629E-10 * ( int32() +
				2.32830643653869629E-10 * int32() );
		}
	};

	struct Expondev : Ran {
		double beta;
		Expondev(double bbeta, unsigned long long int i) : Ran(i), beta(bbeta) {}
		double dev() {
			double u;
			do u = doub(); while (u == 0.);
			return -log(u)/beta;
		}
	};


	
	struct Logisticdev : Ran {
		double mu,sig;
		Logisticdev(double mmu, double ssig, unsigned long long int i) : Ran(i), mu(mmu), sig(ssig) {}
		double dev() {
			double u;
			do u = doub(); while (u*(1.-u) == 0.);
			return mu + 0.551328895421792050*sig*log(u/(1.-u));
		}
	};
	struct Normaldev_BM : Ran {
		double mu,sig;
		double storedval;
		Normaldev_BM(double mmu, double ssig, unsigned long long int i)
			: Ran(i), mu(mmu), sig(ssig), storedval(0.) {}
		double dev() {
			double v1,v2,rsq,fac;
			if (storedval == 0.) {
				do {
					v1=2.0*doub()-1.0;
					v2=2.0*doub()-1.0;
					rsq=v1*v1+v2*v2;
				} while (rsq >= 1.0 || rsq == 0.0);
				fac=sqrt(-2.0*log(rsq)/rsq);
				storedval = v1*fac;
				return mu + sig*v2*fac;
			} else {
				fac = storedval;
				storedval = 0.;
				return mu + sig*fac;
			}
		}
	};
	struct Cauchydev : Ran {
		double mu,sig;
		Cauchydev(double mmu, double ssig, unsigned long long int i) : Ran(i), mu(mmu), sig(ssig) {}
		double dev() {
			double v1,v2;
			do {
				v1=2.0*doub()-1.0;
				v2=doub();
			} while (SQR(v1)+SQR(v2) >= 1. || v2 == 0.);
			return mu + sig*v1/v2;
		}
	};

	struct Normaldev : Ran {
		double mu,sig;
		Normaldev(double mmu, double ssig, unsigned long long int i)
			: Ran(i), mu(mmu), sig(ssig){}
		double dev() {
			double u,v,x,y,q;
			do {
				u = doub();
				v = 1.7156*(doub()-0.5);
				x = u - 0.449871;
				y = abs(v) + 0.386595;
				q = SQR(x) + y*(0.19600*y-0.25472*x);
			} while (q > 0.27597
				&& (q > 0.27846 || SQR(v) > -4.*log(u)*SQR(u)));
			return mu + sig*v/u;
		}
	};
	struct Gammadev : Normaldev {
		double alph, oalph, bet;
		double a1,a2;
		Gammadev(double aalph, double bbet, unsigned long long int i)
			: Normaldev(0.,1.,i), alph(aalph), oalph(aalph), bet(bbet) {
				if (alph <= 0.) throw("bad alph in Gammadev");
				if (alph < 1.) alph += 1.;
				a1 = alph-1./3.;
				a2 = 1./sqrt(9.*a1);
		}
		double dev() {
			double u,v,x;
			do {
				do {
					x = Normaldev::dev();
					v = 1. + a2*x;
				} while (v <= 0.);
				v = v*v*v;
				u = doub();
			} while (u > 1. - 0.331*SQR(SQR(x)) &&
				log(u) > 0.5*SQR(x) + a1*(1.-v+log(v)));
			if (alph == oalph) return a1*v/bet;
			else {
				do u=doub(); while (u == 0.);
				return pow(u,1./oalph)*a1*v/bet;
			}
		}
	};
}


namespace Mat_fuc
{
	void svdcmp(vector< vector<double> > &a, vector<double> &w, vector< vector<double> > &v);//a=u*w*v(T),a is changed to u
	double SVD_inverse1(vector< vector<double> > & mat, vector< vector<double> > & inverse);//return rank, normal index(same as used in math)
	double SVD_inverse2(vector< vector<double> > & mat, vector< vector<double> > & inverse);//transfered index
	double pythag(const double a, const double b);
	void svbksb(vector< vector<double> > &u, vector<double> &w, vector< vector<double> > &v, vector<double> &b, vector<double> &x);
	void solve(vector< vector<double> > &X, vector<double> &b, vector<double> &Y);

	void printVec(string wt, vector<double> & vec);
	void printMat(string wt, vector< vector<double> > & mat);


	inline double SIGN(double a, double b);
	inline double MAX(double a, double b);
	inline double MIN(double a, double b);
	inline double SQR(double a);

	inline int MAX(int a, int b);
	inline int MIN(int a, int b);

	struct Symmeig {
		int n;
		vector< vector<double> > z;
		vector<double> d,e;
		bool yesvecs;

		Symmeig(vector< vector<double> > &a, bool yesvec=true) : n(a.size()), z(a), d(n),
			e(n), yesvecs(yesvec)
		{
			tred2();
			tqli();
			sort();
		}
		void sort() {
			if (yesvecs)
				eigsrt(d,&z);
			else
				eigsrt(d);
		}
		void tred2();
		void tqli();
		double pythag(const double a, const double b);
		void eigsrt(vector<double> &d, vector< vector<double> > *v=NULL);
	};

}

inline double Mat_fuc::SIGN(double a, double b)
{
	return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
}

inline double Mat_fuc::MAX(double a, double b)
{
	return b > a ? (b) : (a);
}

inline double Mat_fuc::MIN(double a, double b) 
{
	return b < a ? (b) : (a);
}

inline double Mat_fuc::SQR(double a) 
{
	return a*a;
}

inline int Mat_fuc::MAX(int a, int b)
{
	return b > a ? (b) : (a);
}

inline int Mat_fuc::MIN(int a, int b) 
{
	return b < a ? (b) : (a);
}



#endif