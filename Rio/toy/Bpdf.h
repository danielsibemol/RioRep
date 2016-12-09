#include <stdio.h>
#include <TMath.h>
#include <TROOT.h>
#include <math.h>
#include "BkgNorm.h"

using namespace std;

double Bpdf(vector<TComplex> &bkg, vector<double> &bkg_coefs, vector<double> &bkg_normalization, double &s12, double &s13, double &s23){

	// This function calculates the background PDF
	
	double bpdf;
	bpdf = 0;
	for (unsigned int i=0; i<bkg.size(); i++) {
		bpdf = bpdf + bkg_coefs.at(i)*bkg.at(i).Re()/bkg_normalization.at(i);
	}
	return bpdf;
}
