#include <TFitter.h>
#include <TTree.h>
#include <TFile.h>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <ostream>
#include "TGraph.h"
#include "TCut.h"
#include "Read_Parameters.h"
#include "FCN.h"
#include "Fractions.h"
#include "PlotFCN.h"

#include <TRandom1.h>


double rnd(double valPar) {
	return (long long) (valPar * 1000 + (valPar > 0 ? 0.5 : -0.5)) / 1000.0;
}

void Log_Likelihood_Fitter(string input_txt_file_name, int sample_number)
{
	//float M, s12, s13, s23, dummyS12, dummyS13, dummyS23;
	double M, s12, s13, s23, dummyS12, dummyS13, dummyS23;


	string input_ntuple_name, output_pwa_txt_file;

	vector<string> included_resonant_channel_string;
	vector<double> tmp_res_extra_pars(4);
	bool is_gaussian;
	int seed, resonance_number, bkg_component_number, npt, number_of_samples, event_number=0;
	TRandom1 * Random = new TRandom1(time(NULL));
	double best_coef1, error_coef1, best_coef2, error_coef2, best_mass, best_width, error_mass, error_width, rndNumber;
	vector<int> fix_parameter_index;
	vector<TComplex> pwa_coefs, pwa_coefs_prime, sig_amp, bkg_amp;
	bool resonances_calc_norm[AVAILABLE_RESONANCES];
	TStopwatch clock;
	clock.Start();

	// Define a minuit Fitter
	TFitter * minimizer = new TFitter(200);

	number_of_resonances = 0;
	number_of_bkg_components = 0;
	resonances.clear();
	fix_parameter_index.clear();
	bkg_components.clear();
	KK_bin_limit.clear();
	pwa_coefs.clear();
	iteration_number = 0;

	// Read the input parameters of the fit
	Read_Parameters(input_txt_file_name, minimizer, real_and_imaginary, final_state, is_gaussian, bkg_par1, bkg_par2, number_of_samples, resonances, bkg_components, number_of_resonances,
			number_of_bkg_components, bkg_fraction, input_ntuple_name, seed, fix_parameter_index, number_of_pwa_bins, output_pwa_txt_file, pwa_coefs, KK_bin_limit, UsesWeights, UseAcceptance, Acceptance_Ntuple_Name, Acceptance_Histo_Name,UseBackHisto, Back_Ntuple_Name, Back_Histo_Name );

	// Checks if the seed is required to be random or not
	/*
	   if (seed != 0) {
	   Random->SetSeed(seed);
	   }
	 */
	// Defines amplitudes vectors sizes
	sig_amp.resize(number_of_resonances);
	bkg_amp.resize(number_of_bkg_components);

	// Define input ntuple
	ostringstream oss;

	if (number_of_samples == 1) oss << input_ntuple_name.c_str();
	else oss << input_ntuple_name.c_str() << "_" << sample_number << ".root";

	TFile *f = new TFile(oss.str().c_str());

	if (f->IsZombie()) {
		std::cout << "LogLikelihoodFitter - Error opening file" << input_ntuple_name << std::endl;
		exit(-1);
	}             

	TTree *t = (TTree*)f->Get( "DecayTree" );
	t->SetBranchAddress( "s12",&s12 );
	t->SetBranchAddress( "s13",&s13 );
	t->SetBranchAddress( "s23",&s23 );
	//  t->SetBranchAddress( "M",&M );
	nentries = t->GetEntries( );
	cout<< " NENTRIES = " << nentries << endl;
	npt = KK_bin_limit.size( );

	Complex_Derivative(KK_bin_limit, pwa_coefs, npt, pwa_coefs_prime);

	// fills the first pwa parameters vector

	for (int i = 0; i<number_of_pwa_bins; i++) {
		last_iteration_pwa_parameters.push_back(minimizer->GetParameter(8*number_of_resonances + number_of_bkg_components + 2*i));
		last_iteration_pwa_parameters.push_back(minimizer->GetParameter(8*number_of_resonances + number_of_bkg_components + 2*i + 1));
	}

	// fills the masses and widths vectors
	for (int res_number = 0; res_number < number_of_resonances; res_number++) {
		last_iteration_masses.push_back(minimizer->GetParameter(8*res_number+2));
		last_iteration_widths.push_back(minimizer->GetParameter(8*res_number+3));
		for (int extra_par = 0;  extra_par < 4; extra_par++) {
			tmp_res_extra_pars[extra_par] = minimizer->GetParameter(8*res_number+4+extra_par);
		}
		last_iteration_res_extra_pars.push_back(tmp_res_extra_pars);
	}

	// Read the ntuple
	M = D_Mass;
	mass_min = 0;
	mass_max = 9999999999999;

	cout << "------------------------------------------------------------------------" << endl;
	cout << endl;
	cout << "Reading ntuple" << endl;
	cout << endl;
	//WeightSum=0, SqrWeightSum=0;
	reading_ntuple = 0;

	int nentries0=0;
	Long64_t jentry = 0;
	for ( Long64_t jentry0 = 0; jentry0 < nentries; jentry0++ ) {
		t->GetEntry( jentry0 );
		if( jentry0 < 10 )  cout<<" jentry0 = "<<jentry<<  endl;
		if ( jentry0 == 0)
		{
			mass_max = M;
			mass_min = M;
		}
		else
		{
			if ( mass_max < M ) mass_max = M;
			if ( mass_min > M ) mass_min = M;
		}

		if ( final_state == 0 )
		{

			S12.push_back( s12 );
			S13.push_back( s13 );
			S23.push_back( s23 );
			Mass.push_back( D_Mass );	

		}

		else if ( final_state == 2 )
		{
//			if(!veto( s12 )) continue;
//			if(!veto( s13 )) continue;

			S12.push_back(s12);
			S13.push_back(s13);
			S23.push_back(s23);
			Mass.push_back(D_Mass);

		}
		else if ( final_state == 3 ){

			Random->SetSeed( seed*jentry );
			rndNumber = Random->Rndm();
			dummyS12  =  double( rndNumber>0.5?s12:s13 );
			dummyS13  =  double( rndNumber>0.5?s13:s12 );				
			dummyS23  =  double( s23 );

			S12.push_back( dummyS12 );
			S13.push_back( dummyS13 );	
			S23.push_back( dummyS23 );
			Mass.push_back( D_Mass );	

		}else{
			cout << "Log_Likelihood_Fitter - ERROR: Unavailable final state" << endl;
			exit(-1);
		}

		Amplitudes(final_state, Mass.at(event_number), S12.at(event_number), S13.at(event_number), S23.at(event_number), sig_amp, bkg_amp, resonances, bkg_components, last_iteration_masses, 
				last_iteration_widths, last_iteration_res_extra_pars, KK_bin_limit, pwa_coefs, pwa_coefs_prime);
		Events_Signal_Amplitudes.push_back(sig_amp);
		Events_Bkg_Amplitudes.push_back(bkg_amp);
		event_number++;

		if( jentry0 < 10 )  	cout<<" eoloop "<< jentry0 <<" s12 "<<s12<<" s13 "<<s13<<" s23 "<<s23<< endl;
		++nentries0;
	}

	nentries = nentries0;
	cout<<" nentries0 = "<<nentries0<<" nentries "<<nentries<<  endl;

	cout<<"Loglike: Mass size : "<<Mass.size()<<" S12 size : "<<S12.size()<<endl;

	reading_ntuple = 0;
	cout << "S12.size() = " << S12.size() << endl;

	for (int res_number = 0; res_number < AVAILABLE_RESONANCES; res_number++) {
		resonances_calc_norm[res_number] = resonances[res_number];
	}
	// Calculate Signal and Background normalization

	cout << "------------------------------------------------------------------------" << endl;
	cout << endl;
	cout << "Calculating Dalitz plot normalization factors" << endl;
	cout << endl;

	ComplexSigNorm(final_state, is_gaussian, mass_min, mass_max, bkg_par1, bkg_par2, number_of_resonances, number_of_bkg_components, resonances, bkg_components, last_iteration_masses, last_iteration_widths, last_iteration_res_extra_pars, resonances_calc_norm, sig_normalization, bkg_normalization, KK_bin_limit, pwa_coefs, pwa_coefs_prime);

	for (int bkg_comp = 0; bkg_comp < number_of_bkg_components; bkg_comp++){
		cout << "bkg_normalization[" << bkg_comp << "} = " << bkg_normalization[bkg_comp] << endl;
	}

	// Set the FCN function of the minimizer
	minimizer->SetFCN(FCN);

	cout << "------------------------------------------------------------------------" << endl;
	cout << endl;
	cout << "Minimizing" << endl;
	cout << endl;

	int status = 0;

	// Run the simplex minimizer to get close to the minimum 
	minimizer->ExecuteCommand("SIMPLEX",0,0);

	// Run the migrad minimizer (an extended Powell's method) to improve the fit. 
	minimizer->ExecuteCommand("MIGRAD",0,0);

	// Run the migrad minimizer (an extended Powell's method) to improve the fit. 
	status = minimizer->ExecuteCommand("HESSE",0,0);

	ostringstream oss2;
	oss2 << "fitparameters.txt";  
	ofstream myfile;  //file to output the fit result - ofstream : Stream class to write on files
	myfile.open (oss2.str().c_str(),ios::app); //ios:app :  append to the end of the file

	// Prints results, with phases in degrees and errors.

	cout << "------------------------------------------------------------------------" << endl;
	cout << endl;
	cout << "FIT RESULTS" << endl;
	cout << endl;

	cout << "########################################################################################################################################################" << endl;
	if (real_and_imaginary) cout << left << setw(20) << "Resonant channel"  << setw(15) << "Real part" << setw(20) << "Real part error" << setw(15)
		<< "Imaginary part" << setw(25) << "Imaginary part error" << setw(10)  << "Mass" << setw(15) << "Mass error" << setw(10)  << "Width" << setw(15)
			<< "Width error" << endl;
	else  cout << left << setw(20) << "Resonant channel"  << setw(15) << "Magnitude" << setw(20) << "Magnitude error" << setw(20) << "Phase (degrees)"
		<< setw(25) << "Phase error (Degrees)" << setw(10) << "Mass" << setw(15) << "Mass error" << setw(10) << "Width" << setw(15) << "Width error" << endl;
	cout << "########################################################################################################################################################" << endl;

	included_resonant_channel_string.resize(number_of_resonances);
	resonance_number = 0;
	double texpars[ 8*AVAILABLE_RESONANCES ];
	for (int i=0; i<AVAILABLE_RESONANCES; i++) {
		if (resonances[i] == 1) {

			best_coef1 = minimizer->GetParameter(8*resonance_number);
			error_coef1 = minimizer->GetParError(8*resonance_number);
			best_coef2 = minimizer->GetParameter(8*resonance_number+1);
			error_coef2 = minimizer->GetParError(8*resonance_number+1);
			best_mass = minimizer->GetParameter(8*resonance_number+2);
			error_mass = minimizer->GetParError(8*resonance_number+2);
			best_width = minimizer->GetParameter(8*resonance_number+3);
			error_width = minimizer->GetParError(8*resonance_number+3);
			cout << left << setw(20) << resonant_channel_string[i] << setw(15) << rnd(best_coef1) << setw(20) << rnd(error_coef1);
			if (real_and_imaginary) cout << setw(15) << rnd(best_coef2) << setw(25) << rnd(error_coef2);
			else cout << setw(20) << rnd(best_coef2*180/PI) << setw(25) << rnd(error_coef2*180/PI);
			cout << setw(10) << best_mass << setw(15) << error_mass << setw(10) << best_width << setw(15) << error_width << endl;

			// writing fit parameters
			myfile << resonant_channel_string[i] << " ";
			myfile << best_coef1 << " " << error_coef1;
			if (real_and_imaginary) myfile << " "  << best_coef2 << " "  << error_coef2;
			else myfile << " "  << best_coef2*180/PI << " "  << error_coef2*180/PI;
			myfile << " "  << best_mass << " "  << error_mass << " "  << best_width << " "  << error_width << " ";
			myfile << endl;


			included_resonant_channel_string[resonance_number] = resonant_channel_string_tex[i];
			texpars[8*resonance_number]=best_coef1;
			texpars[8*resonance_number+1]=error_coef1;
			texpars[8*resonance_number+2]=best_coef2;
			texpars[8*resonance_number+3]=error_coef2;

			resonance_number++;
		}
	}
	cout << "#######################################################################################################################################################" << endl;
	cout << endl;
	cout << "################################################################" << endl;
	cout << left << setw(25) << "Background component"  << setw(15) << "Fraction" << setw(20) << "Fraction error" << endl;
	cout << "################################################################" << endl;

 	myfile << "status " << status << endl;
	myfile << "fcn "    << last_fcn << endl;
	myfile.close( );

	bkg_component_number = 0;
	for (int i=0; i<AVAILABLE_BKG_COMPONENTS; i++) {
		if (bkg_components[i] == 1) {

			best_coef1 = minimizer->GetParameter(8*number_of_resonances+bkg_component_number);
			error_coef1 = minimizer->GetParError(8*number_of_resonances+bkg_component_number);
			cout << left << setw(25) << bkg_component_string[i] << setw(15) << best_coef1 << setw(20) << error_coef1 << endl;
			bkg_component_number++;
		}
	}
	cout << "################################################################" << endl;
	cout << endl;
	cout << "number of iterations: " << iteration_number << endl;
	cout << endl;

	// Define PWA output txt file name
	ostringstream out_string;

	out_string << output_pwa_txt_file.c_str() << "_" << sample_number << ".txt";

	ofstream PWA_output_file(out_string.str().c_str());

	for (int PWA_bin = 0; PWA_bin < npt; PWA_bin++) {
		PWA_output_file << minimizer->GetParameter(8*resonance_number+bkg_component_number+2*PWA_bin) << " "
			<< minimizer->GetParError(8*resonance_number+bkg_component_number+2*PWA_bin) << " "
			<< minimizer->GetParameter(8*resonance_number+bkg_component_number+2*PWA_bin+1) << " "
			<< minimizer->GetParError(8*resonance_number+bkg_component_number+2*PWA_bin+1) << endl;
	}
	PWA_output_file.close();



	UseAcceptance=kFALSE;
	vector<vector<TComplex> > NormForFractions;

	ComplexSigNorm(final_state, is_gaussian, mass_min, mass_max, bkg_par1, bkg_par2, number_of_resonances, number_of_bkg_components, resonances, bkg_components, last_iteration_masses, last_iteration_widths, last_iteration_res_extra_pars, resonances_calc_norm, NormForFractions, bkg_normalization, KK_bin_limit, pwa_coefs, pwa_coefs_prime);





	Fractions(minimizer, real_and_imaginary,NormForFractions,number_of_resonances, fix_parameter_index, included_resonant_channel_string, input_ntuple_name, oss2.str() ,texpars);


	cout<<" UseAcceptance before "<<UseAcceptance<<endl;

	UseAcceptance=kTRUE;


	bool drawPlot = kTRUE;
	cout << "S12.size() = " << S12.size() << endl;

	cout << " (no plot) fitting time of sample " << sample_number <<" = "<< endl;
	clock.Print("u");
	if( drawPlot ) 
		PlotFCN(	minimizer, 
				final_state, 
				is_gaussian, 
				mass_min, 
				mass_max, 
				bkg_par1, 
				bkg_par2, 
				number_of_resonances, 
				number_of_bkg_components, 
				number_of_pwa_bins, 
				resonances, 
				bkg_components, 
				S12,
				S13, 
				S23, 
				input_txt_file_name,
				oss2.str() );

	clock.Stop();
	cout << "full fitting time of sample " << sample_number << endl;
	clock.Print("u");
	S12.clear();
	S13.clear();
	S23.clear();
	Mass.clear();
	last_iteration_masses.clear();
	last_iteration_pwa_parameters.clear();
	last_iteration_widths.clear();
	Events_Bkg_Amplitudes.clear();
	Events_Signal_Amplitudes.clear();
	delete f;
}

