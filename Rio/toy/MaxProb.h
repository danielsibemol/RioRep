#include <stdio.h>
#include <TMath.h>
#include <TRandom1.h>
#include <TROOT.h>
#include <TComplex.h>
#include <TFitter.h>
#include "FCN.h"

using namespace std;

void MaxProb(int final_state, bool is_gaussian, double Mass_min, double Mass_max, vector<TComplex> coefs, vector<double> bkg_coefs, vector<double> bkg_normalization, vector<int> resonances, vector<int> bkg_components, double pdfmax[], TRandom1 * Random, vector<double> KK_bin_limits, vector<TComplex> pwa_coefs, vector<TComplex> pwa_coefs_prime ){
	double M, s12, s13, s23, max_sig_M, max_sig_s12, max_sig_s13, max_sig_s23, max_bkg_M, max_bkg_s12, max_bkg_s13, max_bkg_s23, m1, m2, m3, m1sq, m2sq, m3sq,
	       cos13_12, cos12_13, cos31_12, Bkg_par1 = 0, Bkg_par2 = 0;
	int ngenevents = 1000000;
	int number_of_resonances, number_of_bkg_components;
	bool is_bkg = 0;
	vector<TComplex> sig_amp;
	vector<TComplex> bkg_amp;
	double *points = new double[number_of_points];
	double *weights = new double[number_of_points];
	TF1 tf1;
        double spdfTemp, bpdfTemp;

	// Defines integration weights and points to be used in generation, in case of gaussian mass.

	tf1.CalcGaussLegendreSamplingPoints(number_of_points,points,weights,rel_tolerance);

	TFitter* minimizer = new TFitter(4);
	double p1 = -1;
	minimizer->ExecuteCommand("SET PRINTOUT",&p1, 1);

	number_of_resonances = coefs.size();
	number_of_bkg_components = bkg_coefs.size();

	pdfmax[0] = 0;
	pdfmax[1] = 0;
	if (final_state == 0) {
		m1 = mK;
		m2 = mK;
		m3 = mpi;
	}
	else if (final_state == 1) {
		m1 = mK;
		m2 = mpi;
		m3 = mpi;
	}
	else if (final_state == 2) {
		m1 = mpi;
		m2 = mpi;
		m3 = mpi;
	}
	else if (final_state == 3) {
		m1 = mK;
		m2 = mK;
		m3 = mK;
	}
	else{
		cout << "MaxProb - ERROR: Incorrect final state inserted";
		exit(-1);
	}

	m1sq = m1*m1;
	m2sq = m2*m2;
	m3sq = m3*m3;

	// Creates the vectors to be inserted in the amplitudes function

	sig_amp.resize(number_of_resonances);
	bkg_amp.resize(number_of_bkg_components);


	// Generates several point in the phase space and chooses the one with higher PDF value as its maximum, for both signal and background.

	for(int i = 0; i < ngenevents; i++){
		Generator(final_state, is_gaussian, Mass_min, Mass_max, is_bkg, Bkg_par1, Bkg_par2, M, s12, s13, s23, cos13_12, cos12_13, cos31_12, Random,
				number_of_points, points, weights);
		Amplitudes(final_state, M, s12, s13, s23, sig_amp, bkg_amp, resonances, bkg_components, KK_bin_limits, pwa_coefs, pwa_coefs_prime);

		spdfTemp =   Spdf(sig_amp, coefs, s12, s13, s23); 
		if (pdfmax[0] < spdfTemp) {
			pdfmax[0] = spdfTemp;	
			max_sig_M = M;
			max_sig_s12 = s12;
			max_sig_s13 = s13;
			max_sig_s23 = s23;
		}

		bpdfTemp = Bpdf(bkg_amp, bkg_coefs, bkg_normalization, s12, s13, s23); 
		if (pdfmax[1] < bpdfTemp) {
			pdfmax[1] = bpdfTemp;
			max_bkg_M = M;
			max_bkg_s12 = s12;
			max_bkg_s13 = s13;
			max_bkg_s23 = s23;
		}		
		spdfTemp=0; 
		bpdfTemp=0;
	}

	// Fills the global variables to be used by the minimizer

	GLOBAL_final_state = final_state;
	GLOBAL_resonances = resonances;
	GLOBAL_bkg_components = bkg_components;
	GLOBAL_coefs = coefs;
	GLOBAL_bkg_coefs = bkg_coefs;
	GLOBAL_bkg_normalization = bkg_normalization;

	GLOBAL_number_of_resonances = number_of_resonances;
	GLOBAL_number_of_bkg_components = number_of_bkg_components;

	GLOBAL_KK_bin_limits = KK_bin_limits;
	GLOBAL_pwa_coefs = pwa_coefs;
	GLOBAL_pwa_coefs_prime = pwa_coefs_prime;
	GLOBAL_is_bkg = 0;

	// Uses the point found in the MC method as an input to a maximiztion (minimization of -Spdf and -Bpdf)

	if (coefs.size() > 0){
		if (is_gaussian) {
			minimizer->SetParameter(0,"M",max_sig_M,0.00001,Mass_min,Mass_max);
		}else{
			minimizer->SetParameter(0,"M",M,0.00001,M-0.1*M,M+0.1*M);
			minimizer->FixParameter(0);
		}
		minimizer->SetParameter(1,"s12",max_sig_s12,0.00000001,max_sig_s12-0.2*max_sig_s12,max_sig_s12+0.2*max_sig_s12);
		minimizer->SetParameter(2,"s13",max_sig_s13,0.00000001,max_sig_s13-0.2*max_sig_s13,max_sig_s13+0.2*max_sig_s13);

		minimizer->SetFCN(FCN);

		minimizer->ExecuteCommand("MIGRAD",0,0);

		max_sig_M = minimizer->GetParameter(0);
		max_sig_s12 = minimizer->GetParameter(1);
		max_sig_s13 = minimizer->GetParameter(2);
		max_sig_s23 = max_sig_M*max_sig_M + m1sq + m2sq + m3sq - max_sig_s12 - max_sig_s13;



		Amplitudes(final_state, max_sig_M, max_sig_s12, max_sig_s13, max_sig_s23, sig_amp, bkg_amp, resonances, bkg_components, KK_bin_limits, pwa_coefs, pwa_coefs_prime);

		pdfmax[0] = 1.3*Spdf(sig_amp, coefs, max_sig_s12, max_sig_s13, max_sig_s23);
	}

	GLOBAL_is_bkg = 1;

	if (bkg_coefs.size() > 0){
		if (is_gaussian) {
			minimizer->SetParameter(0,"M",max_bkg_M,0.00001,Mass_min,Mass_max);
		}

		minimizer->SetParameter(1,"s12",max_bkg_s12,0.00000001,max_bkg_s12-0.2*max_bkg_s12,max_bkg_s12+0.2*max_bkg_s12);
		minimizer->SetParameter(2,"s13",max_bkg_s13,0.00000001,max_bkg_s13-0.2*max_bkg_s13,max_bkg_s13+0.2*max_bkg_s13);

		minimizer->SetFCN(FCN);

		minimizer->ExecuteCommand("MIGRAD",0,0);

		max_bkg_M = minimizer->GetParameter(0);
		max_bkg_s12 = minimizer->GetParameter(1);
		max_bkg_s13 = minimizer->GetParameter(2);
		max_bkg_s23 = max_bkg_M*max_bkg_M + m1sq + m2sq + m3sq - max_bkg_s12 - max_bkg_s13;

		Amplitudes(final_state, max_bkg_M, max_bkg_s12, max_bkg_s13, max_bkg_s23, sig_amp, bkg_amp, resonances, bkg_components, KK_bin_limits, pwa_coefs, pwa_coefs_prime);

		pdfmax[1] = 1.3*Bpdf(bkg_amp, bkg_coefs, bkg_normalization, max_bkg_s12, max_bkg_s13, max_bkg_s23);
	}

	sig_amp.clear();
	bkg_amp.clear();
	cout << "------------------------------------------------------------------------" << endl;
	cout << endl;
	cout << "Maximum of signal PDF: " << pdfmax[0] << endl;
	cout << "Maximum of background PDF: " << pdfmax[1] << endl;
	cout << endl;

	GLOBAL_bkg_coefs.clear();
	GLOBAL_bkg_components.clear();
	GLOBAL_bkg_normalization.clear();
	GLOBAL_coefs.clear();
	GLOBAL_KK_bin_limits.clear();
	GLOBAL_pwa_coefs.clear();
	GLOBAL_pwa_coefs_prime.clear();
	GLOBAL_resonances.clear();
	}
