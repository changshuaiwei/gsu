#ifndef _PAMU_H
#define _PAMU_H

#include <Eigen/Dense> //for matrix calculation
#include <boost/math/distributions/non_central_chi_squared.hpp> //for non-central chi squared

//need to first include these two, otherwise, there will be some strange bug under linux

#include <functional>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <list>
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <new>

#include "parset.h"
#include "snpdt.h"
#include "lmw.h"
#include "kfold.h"
#include "tamw.h"
#include "hwu.h"
#include "gtu.h"
#include "format.h"


class Individual;
class Locus;
class CSNP;
class snpdt;

using namespace std;

class Pamu
{
public:
	Pamu(){
		_data=0;
	}
	~Pamu(){
		_data=0;
	}
	
	void initialize(snpdt *);
	void clear();
	void printLOG(string);
	void readBinData();
	void readPedData();
	void readCovData();
	void readPheData();
	void readPheMat();
	void readPheWt();
	void readPopuInfo();
	void readSnpSet();
	void write_BITFILE();
	void write_PEDFILE();
	void assocGTU();
	void assocFGTU();
	void assocHWU();

	void calKP();
	void assocLMW();
	void cvLMW();
	void aplLMW();
	void assocFLMW();
	void cvFLMW();
	void aplFLMW();
	void scanTAMW();
	void aplTAMW();
	void scanFTAMW();
	void aplFTAMW();
	void fmtFMT();
	snpdt* dataPointer();

private:
	snpdt * _data;

};


#endif