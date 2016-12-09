#include <iostream>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TFile.h>
#include <vector>
#include "Spdf.h"
#include "Bpdf.h"

using namespace std;

vector<int> GLOBAL_resonances, GLOBAL_bkg_components;
vector<TComplex> GLOBAL_coefs, GLOBAL_pwa_coefs, GLOBAL_pwa_coefs_prime;
vector<double> GLOBAL_bkg_coefs, GLOBAL_bkg_normalization, GLOBAL_KK_bin_limits; 
int GLOBAL_final_state, GLOBAL_number_of_resonances, GLOBAL_number_of_bkg_components;
bool GLOBAL_is_bkg = 0;

Long64_t nentries;


void FCN(int& nDim, double* gout, double& result, double par[], int flg){

	double PDF, Mass, S12, S13, S23, cos13_12, m1, m2, m3, m1sq, m2sq, m3sq, S, S12min, S12max, S13min, S13max, S23min, S23max;
	vector<TComplex> sig_amp, bkg_amp;

	int dummy_int;
	double * dummy_double;
	dummy_int = nDim;
	dummy_int = flg;
	dummy_double = gout;

	// This function calculates the PDF to be maximized. It is called every in every iteration of minuit, until the PDF is maximized.

	if (GLOBAL_final_state == 0) {
		m1 = mK;
		m2 = mK;
		m3 = mpi;
	}
	else if (GLOBAL_final_state == 1) {
		m1 = mK;
		m2 = mpi;
		m3 = mpi;
	}
	else if (GLOBAL_final_state == 2) {
		m1 = mpi;
		m2 = mpi;
		m3 = mpi;
	}
	else if (GLOBAL_final_state == 3) {
		m1 = mK;
		m2 = mK;
		m3 = mK;
	}
	else{
		cout << "FCN - ERROR: Incorrect final state inserted";
		return;
	}

	Mass = par[0];
	S12 = par[1];
	S13 = par[2];
	m1sq = m1*m1;
	m2sq = m2*m2;
	m3sq = m3*m3;
	S = Mass*Mass;    
	S23 = S + m1sq + m2sq + m3sq - S12 - S13;

	S12min = (m1 + m2)*(m1 + m2);
	S12max = (Mass - m3)*(Mass - m3);
	S13min = (m1 + m3)*(m1 + m3);
	S13max = (Mass - m2)*(Mass - m2);
	S23min = (m2 + m3)*(m2 + m3);
	S23max = (Mass - m1)*(Mass - m1);

	if (S12 < S12min || S12 > S12max || S13 < S13min || S13 > S13max || S23 < S23min || S23 > S23max){
		result =  0;
		return;
	}


	// Defines vectors with apropriate sizes for Amplitudes function.

	sig_amp.resize(GLOBAL_number_of_resonances);
	bkg_amp.resize(GLOBAL_number_of_bkg_components);

	// Sets initial value of the pdf.

	PDF = 0;

	// Calculates the pdf

	cos13_12 = ((S - S12 - m3sq)*(S12 + m1sq - m2sq) + 2*S12*(m1sq + m3sq - S13))/sqrt(lambda(S, S12, m3sq)*lambda(S12, m1sq, m2sq));

	if (fabs(cos13_12)>1){
		result =  0;
		return;
	}
	Amplitudes(GLOBAL_final_state, Mass, S12, S13, S23, sig_amp, bkg_amp, GLOBAL_resonances, GLOBAL_bkg_components, GLOBAL_KK_bin_limits, GLOBAL_pwa_coefs, GLOBAL_pwa_coefs_prime);
	if (GLOBAL_is_bkg) PDF = Bpdf(bkg_amp, GLOBAL_bkg_coefs, GLOBAL_bkg_normalization, S12, S13, S23);
	else PDF = Spdf(sig_amp, GLOBAL_coefs, S12, S13, S23);

	// Gives result to minuit.
    
	result = -PDF;
}
