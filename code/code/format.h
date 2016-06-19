#ifndef _FORMAT_H_
#define _FORMAT_H_

#include <cmath>
#include <limits>
#include <iostream>
#include <fstream>
#include <vector>
#include "calculate.h"

#include "parset.h"
#include "snpdt.h"

using namespace std;

namespace FMT{
	void csv2tped(string infile, string outfile);
}

#endif