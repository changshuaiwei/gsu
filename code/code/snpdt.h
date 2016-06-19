#ifndef _SNPDT_H_
#define _SNPDT_H_

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <functional>
#include <new>
#include <bitset>
#include <algorithm>

#include "parset.h"

class Individual;
class Locus;
class CSNP;

class StrPair;

using namespace std;

typedef vector<double> dvec_t;
typedef vector<bool> bvec_t;
typedef vector< vector<bool> > bmatrix_t;
typedef vector<int> intvec_t;
typedef vector< vector<double> > dmatrix_t;
typedef vector< vector<int> > intmatrix_t;

typedef vector<Individual*>::iterator iIndividual;
typedef vector<Locus*>::iterator iLocus;
typedef vector<CSNP*>::iterator iSNP;

class DmatBlock{
public:
	DmatBlock(){
		nrow=0;
		ncol=0;
		M=0;
	}
	DmatBlock(int nr, int nc){
		nrow=nr;
		ncol=nc;
		double * tmp=new double[nrow*ncol];
		M=tmp;
	}

	DmatBlock(vector< vector<double> > mt){
		nrow=0; ncol=0;
		nrow=mt.size();
		if(nrow>0){
			ncol=mt[0].size();
			if(ncol>0){
				M=new double[nrow*ncol];

				for(int i=0; i<mt.size(); i++){
					for(int j=0; j<mt[i].size(); j++){
						M[i*nrow+j]=mt[i][j];
					}
				}
			}
		}

	}

	~DmatBlock(){
		delete [] M;
	}

	void clear(){
		double* tmp=NULL;
		delete [] M;
		M=tmp;
	}

	void copy(vector< vector<double> > mt){
		if(mt.size()>0){
			if(mt[0].size()>0){
				nrow=mt.size();
				ncol=mt[0].size();

				double * tmp=new double[nrow*ncol];
				delete [] M;
				M=tmp;

				for(int i=0; i<mt.size(); i++){
					for(int j=0; j<mt[i].size(); j++){
						M[i*nrow+j]=mt[i][j];
					}
				}
			}
		}
	}

	inline void assign(int i_row, int j_col, double t){
		//without checking to speed up
		M[i_row*nrow + j_col]=t;
	}

	inline double get(int i_row, int j_col){
		//without checking to speed up
		return M[i_row*nrow + j_col];
	}

	double * getPointer(){
		return M;
	}

private:
	int nrow;
	int ncol;
	double * M;
};

class Dmat{
public:
	Dmat(){
		nrow=0;
		ncol=0;
		M=0;
	}
	Dmat(int nr, int nc){
		nrow=nr;
		ncol=nc;
		double ** tmp=new double* [nrow];
		for(int i=0; i<nrow; i++){
			double * tmptmp=new double [ncol];
			tmp[i]=tmptmp;
		}
		M=tmp;
	}

	Dmat(vector< vector<double> > mt){
		nrow=0; ncol=0;
		nrow=mt.size();
		if(nrow>0){
			ncol=mt[0].size();
			if(ncol>0){

				double ** tmp=new double* [nrow];
				for(int i=0; i<nrow; i++){
					double * tmptmp=new double [ncol];
					tmp[i]=tmptmp;
				}

				M=tmp;

				for(int i=0; i<mt.size(); i++){
					for(int j=0; j<mt[i].size(); j++){
						//M[i*nrow+j]=mt[i][j];
						M[i][j]=mt[i][j];
					}
				}
			}
		}

	}

	~Dmat(){
		for(int i=0; i<nrow; i++){
			delete [] M[i];
		}
		delete [] M;
	}

	void clear(){
		double** tmp=NULL;
		for(int i=0; i<nrow; i++){
			delete [] M[i];
		}
		delete [] M;
		M=tmp;
	}

	void copy(vector< vector<double> > mt){
		if(mt.size()>0){

			if(mt[0].size()>0){
				nrow=mt.size();
				ncol=mt[0].size();

				double ** tmp=new double* [nrow];
				for(int i=0; i<nrow; i++){
					double * tmptmp=new double [ncol];
					tmp[i]=tmptmp;
				}

				for(int i=0; i<nrow; i++){
					delete [] M[i];
				}
				delete [] M;
				M=tmp;

				for(int i=0; i<mt.size(); i++){
					for(int j=0; j<mt[i].size(); j++){
						//M[i*nrow+j]=mt[i][j];
						M[i][j]=mt[i][j];
					}
				}
			}
		}
	}

	inline void assign(int i_row, int j_col, double t){
		//without checking to speed up
		//M[i_row*nrow + j_col]=t;
		M[i_row][j_col]=t;
	}

	inline double get(int i_row, int j_col){
		//without checking to speed up
		//return M[i_row*nrow + j_col];
		return M[i_row][j_col];
	}

	double ** getPointer(){
		return M;
	}

private:
	int nrow;
	int ncol;
	double ** M;
};

class Dvec{
public:
	Dvec(){
		size=0;
		V=0;
	}
	Dvec(int s){
		size=s;
		V=new double[size];
	}
	Dvec(vector<double> vec){
		size=vec.size();
		if(size>0){
			V=new double[size];
			for(int i=0; i<size; i++){
				V[i]=vec[i];
			}
		}
	}

	~Dvec(){
		delete [] V;
	}

	void clear(){
		double* tmp=NULL;
		delete [] V;
		V=tmp;
	}

	void copy(vector<double> vec){
		size=vec.size();
		if(size>0){
			//first creat then delete, important!!!
			double* tmp=new double[size];
			delete [] V;
			V=tmp;

			for(int i=0; i<size; i++){
				V[i]=vec[i];
			}
		}else{
			double* tmp=NULL;
			delete [] V;
			V=tmp;
		}
	}

	inline void assign(int idx, double x){
		V[idx]=x;
	}
	inline double get(int idx){
		return V[idx];
	}

	double * getPointer(){
		return V;
	}

	void getVec(vector<double> & vec){
		vec.clear();
		vec.resize(size);
		for(int i=0; i<vec.size(); i++){
			vec[i]=V[i];
		}
	}
private:
	int size;
	double * V;
};


//class for individual
class Individual
{
public:
	Individual() { 
		fid=iid=pat=mat=""; 
		sex=-9; phenotype=-9;
		pheno_str="";
		sexcode="";
		aff=-9;
		nfid=-9;
		clist.resize(0);
		clistMissing.resize(0);
		plist.resize(0);
		plistMissing.resize(0);

		M_addr=0;
		F_addr=0;
		tmp_child=0;
		children.clear();
		inner_ID=-9;
	}

	~Individual(){
		fid=iid=pat=mat=""; 
		sex=-9; phenotype=-9;
		pheno_str="";
		sexcode="";
		aff=-9;
		nfid=-9;
		clist.resize(0);
		clistMissing.resize(0);
		plist.resize(0);
		plistMissing.resize(0);

		M_addr=0;
		F_addr=0;
		tmp_child=0;
		children.clear();
		inner_ID=-9;
	}

	string fid;
	string iid;

	int nfid;

	// Parental codes
	string pat;      
	string mat;

	int sex;//sex, 0 for female, 1 for male, -9 for missing
	string sexcode;
	double phenotype;//quantitative
	string pheno_str;
	int aff;//binary phenotype, 0 for unaff, 1 for aff, -9 for missing

	dvec_t clist; // multiple covariates
	vector<bool> clistMissing;

	dvec_t plist; // multiple phenotypes
	vector<bool> plistMissing;

	vector<double> cov;
	vector<double> sndcov;
	vector<double> phe_vec;

	Individual * M_addr;
	Individual * F_addr;
	vector<Individual*> children;
	Individual * tmp_child;
	int inner_ID;

};





//a toy example
class StrPair {
public:
	StrPair(){strKey=""; strEntry="";}

	string strKey;
	string strEntry;

	// Copy constructor
	StrPair(const StrPair& h1) { copy(h1); }
	StrPair & operator= (const StrPair & h1) { copy(h1); return *this; }

	void copy(const StrPair &h1)
	{
		strKey = h1.strKey;
		strEntry = h1.strEntry;
	}

	bool operator< (const StrPair & p2) const
	{
		return (strKey < p2.strKey) ;
	}

	bool operator== (const StrPair & p2) const
	{
		return (strKey == p2.strKey && strEntry == p2.strEntry  );
	}
};

//class for genetic locus
class Locus {
public:
	Locus() { chr=0; name=""; allele1=""; allele2=""; freq=0; pos=0; bp=0; nm=0; }

	int chr;
	string name;
	string allele1;
	string allele2;

	double freq;     // of allele1
	double pos;      // cM map positions
	int bp;          // base-pair position
	int nm;          // number of non-missing alleles

	// Copy constructor
	Locus(const Locus& h1) { copy(h1); }
	Locus & operator= (const Locus & h1) { copy(h1); return *this; }

	void copy(const Locus &h1)
	{
		chr = h1.chr;
		name = h1.name;
		allele1 = h1.allele1;
		allele2 = h1.allele2;
		freq = h1.freq;
		pos = h1.pos;
		bp = h1.bp;
		nm = h1.nm;
	}

	bool operator< (const Locus & p2) const
	{
		return (chr < p2.chr || (chr == p2.chr && bp < p2.bp) );
	}

	bool operator== (const Locus & p2) const
	{
		return ( name == p2.name );
	}


};

// Main genotype storage, ordered by SNP
class CSNP
{
public:

	vector<bool> one; // SNP-major mode genotypes
	vector<bool> two; 

};

//define less<Locus*> operator
namespace std
{
	template<>
	class less<Locus*> {
	public:
		bool operator()(Locus const* p1, Locus const* p2)
		{

			// Locus comparison based first on distance, 
			// but then pointers in case we have a degenerate map 
			// file (i.e. so we can sort on position, but so that 
			// set<Locus*> still works

			if(!p1)
				return true;
			if(!p2)
				return false;

			if (p1->chr < p2->chr)
				return true;

			if (p1->chr > p2->chr)
				return false;

			if (p1->bp < p2->bp) 
				return true;

			return false;

		}
	};

};






//the data class snpdt
class snpdt
{
public:
	snpdt();

	~snpdt();

	void test();

	void readBinData();

	bool openBinaryFile(string, ifstream &);

	void readFamFile(string filename);
	void getFamInfo(vector<Individual*> & fam);
	void readPopuFile(string filename);
	void getPopuInfo(vector<Individual*> & popu);

	void readCovFile(string filename);
	void read2ndCovFile(string filename);

	void readPheMat(string filename);

	void readPheFile(string filename);
	void readWtFile(string filename);

	void readSnpSet();
	void readSnpSet2Cols();

	vector<int> readExtSnpF(string filename);

	void clear();

	void write_BITFILE();

	void write_BITFILE(string file_root);

	void subset(vector<int> sample, vector<int> locus);
	void subset(vector<int> locus);
	void sampling(vector<int> sample);
	void reset();//use the backup original data, and clear the subset data

	string genotypeToStr(int, int);

	string intToGenotype(int, int);

	int genoStrToInt(int snp, string A1, string A2);

	void getCovMat(vector< vector<double> > & cov);
	void get2ndCovMat(vector< vector<double> > & cov);
	void getPheMat(vector< vector<double> > & phe);
	void getWtVec(vector<double> & weight);
	void getPheVec(vector< double > & phe);
	void getSnpSet(vector< vector<int> > & snpset);
	inline string getSnpSetNameI(int idx);

	inline int genotypeToInt(int indi, int snp);
	inline int allel1ToInt(int indi, int snp);
	inline int allel2ToInt(int indi, int snp);

	void writeIntToGenotype(int, int, int);

	void writeGenMapPed();
	void writeGenMapPed(string file_root);

	void readGenMapPed();
	void readMtxFile();

	vector<int> matchName(vector<string> names);

	void matchNameAllel(vector<string> names, vector<string> allels, vector<int> & colidx, vector<int> & geno_type);

	

	////function for extract information
	inline int totalIndi();
	inline int totalLocus();
	inline int nHealth();
	inline int nDisease();
	inline int bPheno(int indi);//return binary phenotype
	inline double qPheno(int indi);//return quant phenotype
	inline string snpName(int snp); 
	inline Locus * getLocus(int snp);
	inline Individual * getIndividual(int indi);
	inline int fID(int indi);

private:
	// Genotype/phenotype per individual file
	vector<Individual*>	_sample; 
	vector<string> _cov_name;
	vector<string> _sndcov_name;
	vector<string> _phs_name;
	vector<double> _phs_weight;

	//whole population, for calculation of kinship
	vector<Individual*> _population;

	// SNP information (ordered by SNP/individual)
	vector<CSNP*>	_SNP;

	// Locus information
	vector<Locus*>	_locus; 

	vector<double>	_phenotype;

	// Genotype/phenotype per individual file
	vector<Individual*>	_sample_bkup; 

	// SNP information (ordered by SNP/individual)
	vector<CSNP*>	_SNP_bkup;

	// Locus information
	vector<Locus*>	_locus_bkup; 

	//snpset information
	vector<string> _setname;
	vector<vector<int> > _snpset;

	bool _resample;
	bool _relocus;
};

inline int snpdt::totalIndi()
{
	return _sample.size();
}

inline string snpdt::getSnpSetNameI(int idx)
{
	return _setname[idx];
}

inline int snpdt::totalLocus()
{
	return _locus.size();
}

inline int snpdt::nHealth()
{
	int count=0;
	for(int i=0; i<_sample.size();i++)
		if (_sample[i]->aff == 0) count++;

	return count;
	
}

inline int snpdt::nDisease()
{
	int count=0;
	for(int i=0; i<_sample.size();i++)
		if (_sample[i]->aff == 1) count++;

	return count;
}

inline int snpdt::genotypeToInt(int indi, int snp)
{

	// Return a genotype in suitable text format to be 
	// written to most output files


	bool s1 = _SNP[snp]->one[indi];
	bool s2 = _SNP[snp]->two[indi];

	if((s1) && (!s2)) return -9;
	else return ( int(s1)+int(s2) );

}

inline int snpdt::allel1ToInt(int indi, int snp)
{

	// Return a genotype in suitable text format to be 
	// written to most output files


	bool s1 = _SNP[snp]->one[indi];
	bool s2 = _SNP[snp]->two[indi];

	if(s1){
		if(s2) return 0;
		else return -9;
	}else{
		return 1;
	}
}


inline int snpdt::allel2ToInt(int indi, int snp)
{

	// Return a genotype in suitable text format to be 
	// written to most output files


	bool s1 = _SNP[snp]->one[indi];
	bool s2 = _SNP[snp]->two[indi];

	if(s2){
		return 1;
	}else{
		if(s1) return -9;
		else return 0;
	}
}

inline int snpdt::bPheno(int indi)
{
	return _sample[indi]->aff;
}

inline double snpdt::qPheno(int indi)
{
	return _sample[indi]->phenotype;
}

inline string snpdt::snpName(int snp)
{
	return _locus[snp]->name;
}


inline Locus * snpdt::getLocus(int snp)
{
	return _locus[snp];
}

inline Individual * snpdt::getIndividual(int indi)
{
	return _sample[indi];
}

inline int snpdt::fID(int indi)
{
	return _sample[indi]->nfid;
}
#endif
