#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <TROOT.h>
#include <TRandom1.h>
#include <TTree.h>
#include <TFile.h>
#include <TComplex.h>
#include <sstream>
#include "Generator.h"
#include "../src/Amplitudes.h"
#include "MaxProb.h"
#include "Read_Parameters.h"


using namespace std;

void ToyMCGenerator(string input_file_name){

	double M, s12, s13, s23, s_low, s_high, cos13_12, cos12_13, cos31_12, pdfmax[2], random_number_1, random_number_2, spdf, bkg_fraction, bpdf, Mass_min, Mass_max,
	       Bkg_par1, Bkg_par2;
	int number_of_bkg_components, number_of_available_bkg_components, number_of_available_resonances, number_of_resonances, final_state, number_of_events, number_of_samples, number_of_generated_signal_events, number_of_generated_bkg_events, seed, number_of_pwa_bins, npt;
	bool is_gaussian, is_bkg, last_generation_succeeded,real_and_imaginary;
	vector<double> bkg_coefs, bkg_normalization, KK_bin_limits;
	vector<int> resonances, bkg_components;
	vector<TComplex> coefs, sig_amp, bkg_amp, pwa_coefs, pwa_coefs_prime;
	string output_file_name;
	TRandom1 * Random = new TRandom1(time(NULL));
	double *points = new double[number_of_points];
	double *weights = new double[number_of_points];

	TF1 tf1;

	tf1.CalcGaussLegendreSamplingPoints(number_of_points,points,weights,rel_tolerance);


	// Read generation parameters from input file

	Read_Parameters(input_file_name, final_state, is_gaussian, Mass_min, Mass_max, Bkg_par1, Bkg_par2, number_of_events, number_of_samples, resonances, bkg_components, coefs, bkg_coefs, bkg_fraction, output_file_name, seed, number_of_pwa_bins, pwa_coefs, KK_bin_limits, real_and_imaginary, UseAcceptance, Acceptance_Ntuple_Name, Acceptance_Histo_Name, UseBackHisto, Back_Ntuple_Name, Back_Histo_Name);

	// If the seed is not required to be random, fix the seed

	if (seed != 0) {
		Random->SetSeed(seed);
	}


	// compute derivatives for binned S-wave amplitude (necessary for spline interpolation.
	npt = number_of_pwa_bins;

	Complex_Derivative(KK_bin_limits, pwa_coefs, npt, pwa_coefs_prime);

	// Counts number of resonant channels available in the model, of resonant channels to be generated, of background components available in the model and number of background components the user chose to generate 

	number_of_available_resonances = resonances.size();
	number_of_resonances = coefs.size();
	number_of_available_bkg_components = bkg_components.size();
	number_of_bkg_components = bkg_coefs.size();

	// Generates "number_of_samples" samples

	for (int sample = 0; sample<number_of_samples; sample++){
		ostringstream oss;
		oss << output_file_name.c_str() << "_" << sample << ".root";
		// Define ntuple and its variables

		TFile *ntuple = new TFile(oss.str().c_str(), "RECREATE");
		TTree *DecayTree = new TTree("DecayTree","signal");
		DecayTree->Branch("s12",&s12,"s12/D");
		DecayTree->Branch("s23",&s23,"s23/D");
		DecayTree->Branch("s13",&s13,"s13/D");
		DecayTree->Branch("s_low",&s_low,"s_low/D");
		DecayTree->Branch("s_high",&s_high,"s_high/D");
		DecayTree->Branch("M",&M,"M/D");
                
		// Calculates the background normalization
		cout << "------------------------------------------------------------------------" << endl;
		cout << endl;
		cout << "Calculating Bkg normalization" << endl;
		cout << endl;
                
		BkgNorm(final_state, is_gaussian, Mass_min, Mass_max, Bkg_par1, Bkg_par2, number_of_bkg_components, bkg_components, bkg_normalization);

		// Calculates the Maximum of signal and background PDFs 
		cout << "------------------------------------------------------------------------" << endl;
		cout << endl;
		cout << "Calculating the signal PDF maximum" << endl;
		cout << endl;
		

		MaxProb(final_state, is_gaussian, Mass_min, Mass_max, coefs, bkg_coefs, bkg_normalization, resonances, bkg_components, pdfmax, Random, KK_bin_limits, pwa_coefs, pwa_coefs_prime);

		sig_amp.resize(number_of_resonances);

		bkg_amp.resize(number_of_bkg_components);

		// Loop until all events are generated
		number_of_generated_signal_events = 0;
		number_of_generated_bkg_events = 0;
		last_generation_succeeded = 0;
		random_number_1 = Random->Rndm();
		for(int i = 0; i < number_of_events; i++){

			// Generates 2 random numbers, one to assign the event to signal or background and the other to compare with the respective normalized PDFs

			if (last_generation_succeeded) random_number_1 = Random->Rndm();
			random_number_2 = Random->Rndm();

			// Assigns to event signal or background 

			if (random_number_1 < bkg_fraction) {

				is_bkg = 1;

				// Generates a point in phase space 

				Generator(final_state, is_gaussian, Mass_min, Mass_max, is_bkg, Bkg_par1, Bkg_par2, M, s12, s13, s23, cos13_12, cos12_13, cos31_12, Random, number_of_points, points, weights);

				// Calculates the amplitudes values at the generated point of the phase space

				Amplitudes(final_state, M, s12, s13, s23, sig_amp, bkg_amp, resonances, bkg_components, KK_bin_limits, pwa_coefs, pwa_coefs_prime);

				// Calculates background PDF

				bpdf = Bpdf(bkg_amp, bkg_coefs, bkg_normalization, s12, s13, s23)/pdfmax[1];

				// Compares random number with background PDF

				if (random_number_2 > bpdf) {
					i--;
					last_generation_succeeded = 0;
					continue;
				}

				// Assigns s_low and s_high values

				if (s12>=s13){
					s_high = s12;
					s_low = s13;
				}else{
					s_high = s13;
					s_low = s12;
				}

				// Takes account of the generated event

				last_generation_succeeded = 1;
				number_of_generated_bkg_events++;
			}
			else {

				is_bkg = 0;

				// Generates a point in phase space 

				Generator(final_state, is_gaussian, Mass_min, Mass_max, is_bkg, Bkg_par1, Bkg_par2, M, s12, s13, s23, cos13_12, cos12_13, cos31_12, Random, number_of_points, points, weights);

				// Calculates the amplitudes values at the generated point of the phase space

				Amplitudes(final_state, M, s12, s13, s23, sig_amp, bkg_amp, resonances, bkg_components, KK_bin_limits, pwa_coefs, pwa_coefs_prime);

				// Calculates signal PDF

				spdf = Spdf(sig_amp,coefs, s12, s13, s23)/pdfmax[0];

				// Compares random number with signal PDF

				if (random_number_2 > spdf) {
					i--;
					last_generation_succeeded = 0;
					continue;
				}

				// Assigns s_low and s_high values

				if (s12>=s13){
					s_high = s12;
					s_low = s13;
				}else{
					s_high = s13;
					s_low = s12;
				}

				// Takes account of the generated event

				last_generation_succeeded = 1;
				number_of_generated_signal_events++;
			}

			// Saves event in ntuple if event number < normalized pdf

			if(i%1000==0 && (i) != 0)cout << "Generating sample number " << sample <<  " - " << double(i)*100.0/double(number_of_events) << "%    \r" << flush;

			DecayTree->Fill();

		}
		cout << "Generated  sample number " << sample <<  " - 100%        \r" << flush;
		cout << endl;
		cout << endl;
		cout << "Number of signal events generated in this sample: " << number_of_generated_signal_events << endl;
		cout << "Number of background events generated in this sample: " << number_of_generated_bkg_events << endl;
		cout << endl;

		// Cleans vectors

		sig_amp.clear();
		bkg_amp.clear();

		// Saves and closes ntuple

		ntuple->cd();
		DecayTree->Write("", TObject::kOverwrite);
		ntuple->Close();
	} 

}


