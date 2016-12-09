#include <stdio.h>
#include <TMath.h>
#include <TROOT.h>
#include <math.h>

using namespace std;

int number_of_iterations = 0; 

double Spdf(vector<TComplex> &sig, vector< vector<TComplex> > &coefs_products, double &s12, double &s13){
	double spdf, sig_spdf_product;
	int number_of_resonances = sig.size();
	spdf = 0;

	// This function calculates the Signal PDF
	for (int i = 0; i < number_of_resonances; i++){
		for (int j = i; j < number_of_resonances; j++){
			if (i == j) sig_spdf_product = coefs_products[i][j].Rho()*sig[i].Rho()*sig[j].Rho();
			else sig_spdf_product = 2*(coefs_products[i][j]*sig[i]*sig[j].Conjugate(sig[j])).Re();
		spdf += sig_spdf_product;
		}
	}
	number_of_iterations++;
	double AccInt= UseAcceptance?acceptance(s12, s13):1;
	return spdf*AccInt;
}

double Bpdf(vector<TComplex> &bkg, double par[], vector<double> &bkg_normalization, int &number_of_resonances, int &number_of_bkg_components){
	double temp_params, bpdf;
	vector<double> bkg_coefs; 

	// This function calculates the background PDF
	for (int i = 0; i < number_of_bkg_components; i++) {
		temp_params = par[(i+8*number_of_resonances)];
		bkg_coefs.push_back(temp_params);
	}

	bpdf = 0;
	for (unsigned int i=0; i<bkg.size(); i++) {
		bpdf = bpdf + bkg_coefs.at(i)*bkg.at(i).Re()/bkg_normalization.at(i);
	}
	return bpdf;
}

double TotPdf(double &bkg_frac, vector<TComplex> &sig, vector< vector<TComplex> > &coefs_products, double &s12, double &s13, double &nls, vector<TComplex> &bkg, double par[], vector<double> &bkg_normalization, int &number_of_resonances, int &number_of_bkg_components ){

	double totFunc = 0;
	totFunc =	(1-bkg_frac)*Spdf(sig,coefs_products, s12, s13)/nls + bkg_frac*Bpdf(bkg, par,bkg_normalization,number_of_resonances, number_of_bkg_components);

	return totFunc;

}



