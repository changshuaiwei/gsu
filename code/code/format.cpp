#include "format.h"

void FMT::csv2tped(string infile, string outfile)
{
	gfun::printLOG("Reading csv file in [ " 
		+ infile + " ] and writing tped file in [ " 
		+ outfile + " ]\n");

	gfun::checkFileExists( infile );

	ifstream COV;
	COV.open(infile.c_str());
	COV.clear();

	ofstream OUTF;
	OUTF.open(outfile.c_str());
	OUTF.clear();

	string str_buf;
	vector<string> indi_name;
	getline(COV,str_buf);
	int Ncov=gfun::split_string(str_buf,indi_name);

	string outfile2=outfile + ".indi.txt";
	ofstream OUTF2;
	OUTF2.open(outfile2.c_str());
	OUTF.clear();
	for(int i=1; i<indi_name.size(); i++){
		OUTF2<<i<<"\t"
			<<indi_name[i]<<"\t"
			<<"0\t"
			<<"0\t"
			<<"0\t"
			<<"0\t"
		<<"\n";
	}
	OUTF2.close();

	int c=0;

	while(!COV.eof())
	{
		cout<<c;
		string tmpstr0, tmpstr;
		vector<string> tmpstrvec0, tmpstrvec;

		COV >> tmpstr0;
		int Ncov_tmp=gfun::split_string(tmpstr0,tmpstrvec0,",");
		if(Ncov_tmp!=Ncov){
			stringstream ss;
			ss << c;
			string str = ss.str();
			gfun::printLOG("\nsome line " + str + "  have wrong number of SNPs!!\n");
			break;
		}
		for(int i=0; i<Ncov; i++){
			tmpstr=tmpstrvec0[i];
			if(i==0){
				tmpstrvec.clear();
				gfun::split_string(tmpstr,tmpstrvec,"_");
				OUTF<<tmpstrvec[0]<<"\t"
					<<tmpstr<<"\t"
					<<"0\t"
					<<tmpstrvec[1]<<"\t";
			}else{
				OUTF<<tmpstr[0]<<" "<<tmpstr[1]<<"\t";
			}
		}

		OUTF<<"\n";
		c++;
		cout<<"\r";
	}

	COV.clear();
	COV.close();

	OUTF.close();
}