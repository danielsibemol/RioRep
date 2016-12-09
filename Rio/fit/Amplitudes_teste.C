#include <stdio.h>
#include <TMath.h>
#include <TComplex.h>
#include <TRandom.h>
#include <TFitter.h>
#include <TH2F.h>
#include <vector>
#include <fstream>
#include <iostream>
#include "GenericFunctions.h"
#include "Amplitudes.h"

using namespace std;

void teste(){

	double s12 = 1.5;
	double s13 = 1.5;
	double s23, m1, m2, m3;
	double M = D_Mass;
	vector<double> a;   
	vector<TComplex> b, c;

	vector<double> last_iteration_pwa_parameters, last_iteration_masses(3), last_iteration_widths(3);
	vector<double> KK_bin_limit, temp_extra_pars(4);
	vector< vector<double> > last_iteration_res_extra_pars(3, temp_extra_pars);
	vector<TComplex> pwa_coefs, pwa_coefs_prime;

	last_iteration_masses[0] = 0.00; 
	last_iteration_widths[0] = 0.00; 
	last_iteration_masses[1] = 0.965; 
	last_iteration_widths[1] = 0.00; 
	last_iteration_masses[2] = 1.019461; 
	last_iteration_widths[2] = 0.004266; 

	cout << "teste 0" << endl;

	cout << "teste 1" << endl;
	int final_state = 1;

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
                cout << "Generator - ERROR: Incorrect final state inserted";
                exit(-1);
        }

	s23 = M*M + m1*m1 + m2*m2 + m3*m3 - s12 - s13;

	cout << "s12 = " << s12 << ", s13 = " << s13 << ", s23 = " << s23 << endl;

	cout << "teste 2" << endl;
	vector<TComplex> sig, bkg;
	sig.resize(3);
	bkg.resize(0);
	vector<int> resonances;
	resonances.resize(37);
	resonances[5] = 1;
	resonances[31] = 1;
	resonances[34] = 1;
	vector<int> bkg_components;
	bkg_components.resize(4);

	cout << "teste 3" << endl;

	Amplitudes(final_state, M, s12, s13, s23, sig, bkg, resonances, bkg_components, last_iteration_masses, last_iteration_widths, last_iteration_res_extra_pars, KK_bin_limit, pwa_coefs, pwa_coefs_prime);
	cout << sig[0] << endl;
	cout << sig[1] << endl;
	cout << sig[2] << endl;

}
