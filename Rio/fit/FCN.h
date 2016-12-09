#include <iostream>
#include <TTree.h>
#include <TFile.h>
#include <vector>
#include <TStopwatch.h>
#include "SigNorm.h"
#include "TotPdf.h"

using namespace std;

vector<double> S12;
vector<double> S13;
vector<double> S23;
vector<double> Mass;
vector<double> KK_bin_limit;
vector< vector<TComplex> > Events_Signal_Amplitudes;
vector< vector<TComplex> > Events_Bkg_Amplitudes;
vector<double> sWeight;
double sFactor;

double last_fcn;

vector<double> Acceptance;

double bkg_fraction, mass_min, mass_max, bkg_par1, bkg_par2, FCN_value;
vector<int> resonances, bkg_components;
vector<TComplex> coefs, bkg_coefs;
int final_state, number_of_resonances, number_of_bkg_components, number_of_pwa_bins, iteration_number = 0;

Long64_t nentries;

vector< vector<TComplex> > sig_normalization;
vector<double> bkg_normalization, last_iteration_pwa_parameters, last_iteration_masses, last_iteration_widths;
vector< vector<double> > last_iteration_res_extra_pars;
bool is_gaussian, real_and_imaginary;

void FCN(int& nDim, double* gout, double& result, double par[], int flg){

	double NL_s, log_likelihood, L_s;
	vector<TComplex> sig_amp, bkg_amp, pwa_coefs, pwa_coefs_prime;
	vector< vector<TComplex> > coefs_product(number_of_resonances,vector<TComplex>(number_of_resonances));
	TComplex temp_pwa_coef, tmp;
	vector<double> res_masses, res_widths, tmp_extra_pars(4), recalc_masses, recalc_widths;
        vector< vector<double> > res_extra_pars, recalc_res_extra_pars;
	vector<int> recalc_res_index, recalc_resonances;
	bool calc_norm = 0, pwa_calc_norm = 0, resonances_calc_norm[AVAILABLE_RESONANCES] = {}, extra_pars_calc_norm = 0;
	int npt, res_number = 0;
	UNUSED(nDim);
	UNUSED(flg);
	UNUSED(gout);
	cout.precision(12);

	// This function calculates the log likelihhod to be maximized. It is called every in every iteration of minuit, until the likelihood is maximized.
	// Checks if any of the parameters used in the pwa Amplitudes have changed
	for (int i = 0; i<number_of_pwa_bins; i++) {
		if (last_iteration_pwa_parameters[2*i] == par[8*number_of_resonances + number_of_bkg_components + 2*i] &&
				last_iteration_pwa_parameters[2*i+1] == par[8*number_of_resonances + number_of_bkg_components + 2*i+1]){
			continue;
		}
		pwa_calc_norm = 1;
		break;
	}

	// Checks if any of the resonances masses or widths have changed

	for (int i = 0; i<AVAILABLE_RESONANCES; i++){
		if (resonances[i] == 1) {
			res_masses.push_back(par[8*res_number+2]);
			res_widths.push_back(par[8*res_number+3]);
			for (int extra_par = 0; extra_par < 4; extra_par++){
                               tmp_extra_pars[extra_par] = par[8*res_number+4+extra_par];
			       if (tmp_extra_pars[extra_par] == last_iteration_res_extra_pars[res_number][extra_par] && extra_pars_calc_norm == 0){
				       extra_pars_calc_norm = 0;
			       }else{
				       extra_pars_calc_norm = 1;
			       }
			       last_iteration_res_extra_pars[res_number][extra_par] = tmp_extra_pars[extra_par];
			}
			res_extra_pars.push_back(tmp_extra_pars);
			res_number++;
			if ((res_masses[res_number - 1] == last_iteration_masses[res_number - 1] && res_widths[res_number - 1] == last_iteration_widths[res_number - 1]) && !(i == PWA_INDEX && pwa_calc_norm) && extra_pars_calc_norm == 0){
				recalc_resonances.push_back(resonances_calc_norm[i]);
				continue;
			}  // Comenta para fazer a massa e largura livres
			calc_norm = 1;
			recalc_res_index.push_back(res_number - 1);
			resonances_calc_norm[i] = 1;
			last_iteration_masses[res_number - 1] = res_masses[res_number - 1];
			last_iteration_widths[res_number - 1] = res_widths[res_number - 1];
			recalc_masses.push_back(res_masses[res_number - 1]);
			recalc_widths.push_back(res_widths[res_number - 1]);
			recalc_res_extra_pars.push_back(tmp_extra_pars);
		}
		recalc_resonances.push_back(resonances_calc_norm[i]);
	}

	coefs.resize(number_of_resonances);

	// Converts the double pars of the minimizer to a TComplex vector, while checking if they changed since last iteration
	for (int i = 0; i < number_of_resonances; i++) {
		if (real_and_imaginary) tmp(par[8*i], par[8*i + 1]);
		else tmp(par[8*i], par[8*i + 1], 1);
		coefs[i] = tmp;
		//cout << tmp << endl;
	}

	// Caclculates the products of the coefficients
	for (int res_i = 0; res_i < number_of_resonances; res_i++) {
		for (int res_j = 0; res_j < number_of_resonances; res_j++) {
			coefs_product[res_i][res_j] = coefs[res_i]*coefs[res_j].Conjugate(coefs[res_j]);
		}
	}

	// Fills the vectors with the current pwa parameters

	for(unsigned long pwa_par_number = 0; pwa_par_number < KK_bin_limit.size(); pwa_par_number++){
		if (real_and_imaginary){
			temp_pwa_coef(par[8*number_of_resonances+number_of_bkg_components+2*pwa_par_number], par[8*number_of_resonances+number_of_bkg_components+2*pwa_par_number+1]);
		}else{
			temp_pwa_coef(par[8*number_of_resonances+number_of_bkg_components+2*pwa_par_number], par[8*number_of_resonances+number_of_bkg_components+2*pwa_par_number+1], 1);
		}
		last_iteration_pwa_parameters[2*pwa_par_number] = par[8*number_of_resonances+number_of_bkg_components+2*pwa_par_number];
		last_iteration_pwa_parameters[2*pwa_par_number+1] = par[8*number_of_resonances+number_of_bkg_components+2*pwa_par_number+1];
		pwa_coefs.push_back(temp_pwa_coef);
	}

	npt = KK_bin_limit.size();

	// Calculates derivatives

	Complex_Derivative(KK_bin_limit, pwa_coefs, npt, pwa_coefs_prime);

	// Calculates integrals, if necessary

	if (calc_norm) {
		cout << "------------------------------------------------------------------------" << endl;
		cout << endl;
		cout << "Calculating amplitudes integrals for signal and background" << endl;
		cout << endl;

		ComplexSigNorm(final_state, is_gaussian, mass_min, mass_max, bkg_par1, bkg_par2, number_of_resonances, number_of_bkg_components, resonances, bkg_components, res_masses, 
				res_widths, res_extra_pars, resonances_calc_norm, sig_normalization, bkg_normalization, KK_bin_limit, pwa_coefs, pwa_coefs_prime);

		/*		for (int res_i = 0; res_i < number_of_resonances; res_i++) {
				for (int res_j = 0; res_j < number_of_resonances; res_j++) {
				cout << "sig_normalization[" << res_i << "][" << res_j << "] = " << sig_normalization[res_i][res_j] << endl;
				}
				}
		 */	}

		// Calculates the Dalitz plot normalization.

		NL_s = SigNorm(sig_normalization, number_of_resonances, par, real_and_imaginary);
		cout << "NL_s = " << NL_s << endl;

		// Defines vectors with apropriate sizes for Amplitudes function.

		sig_amp.resize(number_of_resonances);
		bkg_amp.resize(number_of_bkg_components);

		// Sets initial value of the likelihood.

		log_likelihood = 0;

		// Calculates the log likelihood for every event and sum it up.

		cout << "------------------------------------------------------------------------" << endl;
		cout << endl;
		cout << "Calculating the likelihood" << endl; 
		for (unsigned int jentry=0; jentry<Mass.size();jentry++) {

			if (calc_norm){
				/*			cout << "iteration number = " << iteration_number << endl;
				//cout << "final_state = " << final_state << endl;
				//cout << "Mass = " << Mass[jentry] << endl;
				cout << "s12 = " << S12[jentry] << endl;
				cout << "s13 = " << S13[jentry] << endl;
				//cout << "s23 = " << S23[jentry] << endl;
				//cout << "sig_amp = " << sig_amp.size() << endl;
				//cout << "bkg_amp = " << bkg_amp.size() << endl;
				cout << "recalc_resonances = " << recalc_resonances.size() << endl;
				//cout << "bkg_components = " << bkg_components.size() << endl;
				cout << "res_masses = " << res_masses.size() << endl;
				cout << "res_widths = " << res_widths.size() << endl;

				for (int i = 0; i<number_of_resonances;  ++i){
				cout << "par[4*i+2] = " << par[4*i+2] << ", par[4*i+3] = " << par[4*i+3]  << endl; 
				cout << "res_masses[" << i << "] = " << res_masses[i] << ", res_widths[" << i << "] = " << res_widths[i] << endl;
				}
				 */		Amplitudes(final_state, Mass.at(jentry), S12.at(jentry), S13.at(jentry), S23.at(jentry), sig_amp, bkg_amp, recalc_resonances, bkg_components, recalc_masses, recalc_widths, recalc_res_extra_pars, KK_bin_limit, pwa_coefs,
						 pwa_coefs_prime);
				 for (unsigned int recalc_res_i = 0; recalc_res_i<recalc_res_index.size(); recalc_res_i++) {
					 Events_Signal_Amplitudes[jentry][recalc_res_index[recalc_res_i]] = sig_amp[recalc_res_i];
				 }
				 Events_Bkg_Amplitudes[jentry] = bkg_amp;
			}
			L_s = 	TotPdf(bkg_fraction,Events_Signal_Amplitudes[jentry], coefs_product, S12.at(jentry), S13.at(jentry), NL_s, Events_Bkg_Amplitudes[jentry], par, bkg_normalization, number_of_resonances, number_of_bkg_components);
			//                if (sqrt(S12[jentry]) < 1.02 && iteration_number == 0) cout << "S12 = " << S12[jentry] << ", S13 = " << S13[jentry] << ", L_s = " << L_s << ", jentry = " << jentry << endl; 	
			if(L_s<=0) continue;        		 
			log_likelihood += log(L_s);

		}		

		// Gives result to minuit.

		result = -2*log_likelihood;
		FCN_value = result;

		cout << "FCN = " << result << endl;

		iteration_number++;

		last_fcn = FCN_value;

}
