#include <stdio.h>
#include <TMath.h>
#include <TROOT.h>
#include "../src/Amplitudes.h"
#include "IntegrationFunctions.h"

using namespace std;

void ComplexSigNorm(int final_state, bool is_gaussian, double Mass_min, double Mass_max, double Bkg_par1, double Bkg_par2, int number_of_resonances,
		int number_of_bkg_components, vector<int> resonances, vector<int> bkg_components, vector<double> res_masses, vector<double> res_widths,
		vector< vector<double> > res_extra_pars,bool resonances_calc_norm[], vector< vector<TComplex> > &sig_normalization,
		vector<double> &bkg_normalization, vector<double> KK_bin_limits, vector<TComplex> pwa_coefs, vector<TComplex> pwa_coefs_prime){
	double M, m1, m2, m3, s13min, s13max, thin_res_border_max;
	int sig_index1 = 0, sig_index2 = 0, bkg_index = 0;
	string SigRe = "SigRe", SigIm = "SigIm", Bkg = "Bkg";
	double *points = new double[number_of_points];
	double *weights = new double[number_of_points];
	vector<double> integral_masses(2), integral_widths(2);
	vector< vector<double> > integral_res_extra_pars(2);
	TwoDimIntFunctor TwoDimSigReIntObject, TwoDimSigImIntObject, TwoDimBkgIntObject;
	OneDimIntFunctor OneDimSigIntObject, OneDimBkgIntObject;
	MassIntFunctor MassSigReIntObject, MassSigImIntObject, MassBkgIntObject;

	// This funtcion calculates the mean value of the amplitude of each resonant channel and background component across the Dalitz plot

	// Defines masses

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
		cout << "ComplexSigNorm - ERROR: Incorrect final state inserted";
		exit(-1);
	}

	M = D_Mass;

	// Calculates integration limits

	s13min = (m1 + m3)*(m1 + m3);
	s13max = (M - m2)*(M - m2);

	// Defines vectors with apropriate sizes, to insert in amplitudes function
	vector< vector< TComplex> > sig_vector(number_of_resonances, vector< TComplex>(number_of_resonances));
	if (sig_normalization.size() == 0) sig_normalization = sig_vector;

	if (bkg_normalization.size() == 0) bkg_normalization.resize(number_of_bkg_components);

	// Calculates integration points and weights

	TF1 function;
	function.CalcGaussLegendreSamplingPoints(number_of_points,points,weights,rel_tolerance);
	for (int j = 0; j < AVAILABLE_RESONANCES; j++){
		sig_index2 = 0;
		if (resonances[j] == 1) {
			for (int k = 0; k < AVAILABLE_RESONANCES; k++){
				if (resonances[k] == 1 ){
					if (resonances_calc_norm[j] == 1 || resonances_calc_norm[k] == 1) {
						// Defines which masses and widths will be used in this components integration
						if (sig_index2 < sig_index1){
							integral_masses[0] = res_masses[sig_index2];
							integral_masses[1] = res_masses[sig_index1];
							integral_widths[0] = res_widths[sig_index2];
							integral_widths[1] = res_widths[sig_index1];
							integral_res_extra_pars[0] = res_extra_pars[sig_index2];
							integral_res_extra_pars[1] = res_extra_pars[sig_index1];
						}else {
							integral_masses[0] = res_masses[sig_index1];
							integral_masses[1] = res_masses[sig_index2];
							integral_widths[0] = res_widths[sig_index1];
							integral_widths[1] = res_widths[sig_index2];
							integral_res_extra_pars[0] = res_extra_pars[sig_index1];
							integral_res_extra_pars[1] = res_extra_pars[sig_index2];
						}

						if (!is_gaussian) {
							// Defines integration functions

							TwoDimSigReIntObject.SetValues(OneDimSigIntObject, number_of_points, points, weights, SigRe, final_state, M, integral_masses, integral_widths,
									integral_res_extra_pars, KK_bin_limits, pwa_coefs, pwa_coefs_prime, j, k);

							TwoDimSigImIntObject.SetValues(OneDimSigIntObject, number_of_points, points, weights, SigIm, final_state, M, integral_masses, integral_widths,
									integral_res_extra_pars, KK_bin_limits, pwa_coefs, pwa_coefs_prime, j, k);

							if( final_state != 0 && fabs(integral_widths[0]) < 0.02){
								thin_res_border_max = pow(integral_masses[0] + 3*integral_widths[0],2);
							}else if (final_state != 0 && fabs(integral_widths[1]) < 0.02){
								thin_res_border_max = pow(integral_masses[1] + 3*integral_widths[1],2); 
							}
							TF1 TwoDimSigReIntFunction("TwoDimSigReIntFunction", TwoDimSigReIntObject, s13min, s13max, 0);
							TF1 TwoDimSigImIntFunction("TwoDimSigImIntFunction", TwoDimSigImIntObject, s13min, s13max, 0);

							// Calculates the integral. The line below calculates the integral along the s12 axis for each s13 point and then integrates along the s13 axis.
							if (sig_index1 > sig_index2) {
								sig_normalization[sig_index1][sig_index2] = sig_normalization[sig_index2][sig_index1].Conjugate(sig_normalization[sig_index2][sig_index1]);
							}else{
								if (final_state != 0 && (fabs(integral_widths[0]) < 0.02 || fabs(integral_widths[1]) < 0.02)){							
									sig_normalization[sig_index1][sig_index2](TwoDimSigReIntFunction.IntegralFast(number_of_points, points, weights, s13min, thin_res_border_max) + 
											TwoDimSigReIntFunction.IntegralFast(number_of_points, points, weights, thin_res_border_max, s13max),
											TwoDimSigImIntFunction.IntegralFast(number_of_points, points, weights, s13min, thin_res_border_max) + 
											TwoDimSigImIntFunction.IntegralFast(number_of_points, points, weights, thin_res_border_max, s13max));
								} else {
									sig_normalization[sig_index1][sig_index2](TwoDimSigReIntFunction.IntegralFast(number_of_points, points, weights, s13min,s13max),
											TwoDimSigImIntFunction.IntegralFast(number_of_points, points, weights, s13min,s13max));
								}	
							}


						}else{
							// Defines integration functions

							MassSigReIntObject.SetValues(TwoDimSigReIntObject, OneDimSigIntObject, number_of_points, points, weights, SigRe, final_state, Mass_min, Mass_max,
									Bkg_par1, Bkg_par2, integral_masses, integral_widths, integral_res_extra_pars, KK_bin_limits, pwa_coefs, pwa_coefs_prime, j, k);
							TF1 MassSigReIntFunction("MassSigReIntFunction", MassSigReIntObject, Mass_min, Mass_max, 0);

							MassSigImIntObject.SetValues(TwoDimSigImIntObject, OneDimSigIntObject, number_of_points, points, weights, SigIm, final_state, Mass_min, Mass_max,
									Bkg_par1, Bkg_par2, integral_masses, integral_widths, integral_res_extra_pars, KK_bin_limits, pwa_coefs, pwa_coefs_prime, j, k);
							TF1 MassSigImIntFunction("MassSigImIntFunction", MassSigImIntObject, Mass_min, Mass_max, 0);

							// Calculates the integral. The line below calculates the integral along the s12 axis for each s13 point and then integrates along the s13 axis.

							if (sig_index1 > sig_index2) {
								sig_normalization[sig_index1][sig_index2] = sig_normalization[sig_index2][sig_index1].Conjugate(sig_normalization[sig_index2][sig_index1]);
							}else {
								sig_normalization[sig_index1][sig_index2](MassSigReIntFunction.IntegralFast(number_of_points, points, weights, Mass_min, Mass_max),
										MassSigImIntFunction.IntegralFast(number_of_points, points, weights, Mass_min, Mass_max));
							}
						}
					}
					sig_index2++;
				}
			}
			sig_index1++;
		}
	}

	for (int j = 0; j < AVAILABLE_BKG_COMPONENTS; j++){
		if (bkg_components[j] == 1){
			if (!is_gaussian) {
				// Defines integration function

				TwoDimBkgIntObject.SetValues(OneDimBkgIntObject, number_of_points, points, weights, Bkg, final_state, M, integral_masses, integral_widths, res_extra_pars, KK_bin_limits, 
						pwa_coefs, pwa_coefs_prime, j, 0);
				TF1 TwoDimBkgIntFunction("TwoDimBkgIntFunction", TwoDimBkgIntObject, s13min, s13max, 0);

				// Calculates the integral. The line below calculates the integral along the s12 axis for each s13 point and then integrates along the s13 axis.

				bkg_normalization.at(bkg_index) = TwoDimBkgIntFunction.IntegralFast(number_of_points, points, weights, s13min,s13max);
				bkg_index++;
			} else{
				// Defines integration function

				MassBkgIntObject.SetValues(TwoDimBkgIntObject, OneDimBkgIntObject, number_of_points, points, weights, Bkg, final_state, Mass_min, Mass_max, Bkg_par1, Bkg_par2, 
						integral_masses, integral_widths, res_extra_pars, KK_bin_limits, pwa_coefs,pwa_coefs_prime, j, 0);
				TF1 MassBkgIntFunction("MassBkgIntFunction", MassBkgIntObject, Mass_min, Mass_max, 0);

				// Calculates the integral. The line below calculates the integral along the s12 axis for each s13 point and then integrates along the s13 axis.

				bkg_normalization.at(bkg_index) = MassBkgIntFunction.IntegralFast(number_of_points, points, weights, Mass_min,Mass_max);
				bkg_index++;
			}
		}
	}
}

double SigNorm(vector< vector<TComplex> > complex_sig_normalization, int number_of_resonances, double par[], bool real_and_imaginary){
	double sig_normalization;
	TComplex tmp;
	vector<TComplex> params;
	TComplex product;

	// This function calculates the normalization factor for the signal

	// Gets fit parameters
	for (int i = 0; i < number_of_resonances; i++) {
		if (real_and_imaginary) tmp(par[8*i], par[8*i + 1]);
		else tmp(par[8*i], par[8*i + 1], 1);
		params.push_back(tmp);
	}

	sig_normalization = 0;

	// Calculates the signal normalization using the mean value of each resonant channel's amplitude

	for (int j = 0; j < number_of_resonances; j++){
		for (int k = 0; k < number_of_resonances; k++){
			product = params.at(j)*params.at(k).Conjugate(params.at(k));
			product = product*complex_sig_normalization[j][k];
			sig_normalization = sig_normalization + product.Re();
		}
	}
	return sig_normalization;
}

