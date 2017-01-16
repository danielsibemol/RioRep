#include <stdio.h>
#include <math.h>
#include <string>
#include <iomanip>
#include <TEnv.h>

using namespace std;

void Read_Parameters(string input_txt_file_name, int &final_state, bool &is_gaussian, double &Mass_min, double &Mass_max, double &Bkg_par1, double &Bkg_par2, 
		int &number_of_events, 
		int &number_of_samples, vector<int> &resonances, vector<int> &bkg_components, vector<TComplex> &coefs, vector<double> &bkg_coefs, double &bkg_fraction, 
		string &output_file_name, int &seed, int &number_of_pwa_bins, vector<TComplex> &pwa_coefs, vector<double> &KK_bin_limits, bool &real_and_imaginary, 
		bool &UseAcceptance,  string &Acceptance_Ntuple_Name, string &Acceptance_Histo_Name, bool &UseBackHisto,  string &Back_Ntuple_Name, string &Back_Histo_Name){

	int available_resonances = AVAILABLE_RESONANCES;
	int available_bkg_components = AVAILABLE_BKG_COMPONENTS;
	int temp_pwa_fix;
	string str_coef1, str_coef2, str_fraction, pwa_txt_file;
	double temp_coef1, temp_coef2, temp_bkg_coef, temp_mKK;
	TComplex temp_coef;

	TEnv *input_txt_file = new TEnv(input_txt_file_name.c_str());

	// This function reads the parameters to be used in the generation 

	// Checks if the input .txt file was loaded correctly

	if (input_txt_file->IsZombie()) {
		std::cout << "ReadParameters - ERROR: File not found - " << input_txt_file_name << std::endl;
		exit(-1);
	}

	// Reads the name of the output file

	output_file_name = input_txt_file->GetValue("output_file_name", "ntuple");

	cout << "Output file name: " << output_file_name << endl;
	cout << endl;
	cout << "------------------------------------------------------------------------" << endl;
	cout << endl;

	// Reads an integer associated with the final state of the sample

	final_state = input_txt_file->GetValue("final_state",-1);

	// Reads the number of events to be generated

	number_of_events = input_txt_file->GetValue("number_of_events",0);

	cout << "Number of events to be generated: " << number_of_events << endl;
	cout << endl;
	cout << "------------------------------------------------------------------------" << endl;
	cout << endl;

	// Reads the number of samples to be generated

	number_of_samples = input_txt_file->GetValue("number_of_samples",0);

	cout << "Number of samples to be generated: " << number_of_samples << endl;
	cout << endl;
	cout << "------------------------------------------------------------------------" << endl;
	cout << endl;

	// Reads the choice for the seed: random or not

	seed = input_txt_file->GetValue("seed",0);

	// Reads the choice for the mass distribution: delta or gaussian

	is_gaussian = input_txt_file->GetValue("is_gaussian",0);

	if(is_gaussian){
		cout << "Gaussian mass distribution"; 
		cout << endl;
		// Reads the Mass window

		Mass_min = input_txt_file->GetValue("Mass_min",0.);
		Mass_max = input_txt_file->GetValue("Mass_max",0.);

		cout << "Mass window: from " << Mass_min << " to " << Mass_max << endl;
		cout << endl;
		cout << "------------------------------------------------------------------------" << endl;
		cout << endl;
		// Reads the Background parameters for mass distribution

		Bkg_par1 = input_txt_file->GetValue("Bkg_par1",0.);;
		Bkg_par2 = input_txt_file->GetValue("Bkg_par2",0.);;

		cout << "Background mass distribution parameters: " << Bkg_par1 << " and " << Bkg_par2 << endl;
		cout << endl;
		cout << "------------------------------------------------------------------------" << endl;
		cout << endl;
	}else{
		cout << "Delta Mass Distribution";
		cout << endl;
		cout << "------------------------------------------------------------------------" << endl;
		cout << endl;
	}
	// Checks if the final state is kkpi

	if (final_state == 0) {
		cout << "The chosen final state was kkpi" << endl;
		cout << endl;
		cout << "------------------------------------------------------------------------" << endl;
		cout << endl;
	}	
	else if (final_state == 1){
		cout << "The chosen final state was kpipi" << endl;
		cout << endl;
		cout << "------------------------------------------------------------------------" << endl;
		cout << endl;
	} 
	else if (final_state == 2){
		cout << "The chosen final state was pipipi" << endl;
		cout << endl;
		cout << "------------------------------------------------------------------------" << endl;
		cout << endl;
	} 
	else if (final_state == 3){
		cout << "The chosen final state was kkk" << endl;
		cout << endl;
		cout << "------------------------------------------------------------------------" << endl;
		cout << endl;
	} 
	else {
		cout << "ReadParameters - ERROR - Invalid Final State" << endl;
		exit(-1);
	}	

	// Reads the choice for the complex coefficients, if given as magn and phase or real and imaginary parts

	real_and_imaginary = input_txt_file->GetValue("real_and_imaginary",0);

	// Reads the resonant real and imaginary parts

	cout << "Signal generation parameters: " << endl;
	cout << endl;
	cout << "##########################################################" << endl;
	if (real_and_imaginary) cout << left << setw(20) << "Resonant Channel"  << setw(15) << "Real Part" << setw(15) << "Imaginary Part" << endl;
	else cout << left << setw(20) << "Resonant Channel"  << setw(15) << "Magnitude" << setw(15) << "Phase (degrees)" << endl;
	cout << "##########################################################" << endl;

	for (int i=0; i<available_resonances; i++) {

		// Defines the identification strings associated with each parameter 

		if (real_and_imaginary) {
			str_coef1 = resonant_channel_string[i];
			str_coef1 += "_re";
			str_coef2 = resonant_channel_string[i];
			str_coef2 += "_im";
		}else{
			str_coef1 = resonant_channel_string[i];
			str_coef1 += "_amp";
			str_coef2 = resonant_channel_string[i];
			str_coef2 += "_phs";
		}
		// Reads the magnitude and phase

		temp_coef1 = input_txt_file->GetValue(str_coef1.c_str(),0.);
		temp_coef2 = input_txt_file->GetValue(str_coef2.c_str(),0.);

		// Checks if the coefficients are not zero. If they are, the resonant channel will be excluded from the model

		if (temp_coef1 == 0 && temp_coef2 == 0) resonances.push_back(0);
		else resonances.push_back(1);

		// Fills the amplitudes and phases that will be used in the model 

		if(resonances.at(i) == 1){
			if (real_and_imaginary){
				temp_coef(temp_coef1, temp_coef2);
				cout << left << setw(20) << resonant_channel_string[i] << setw(15) << temp_coef.Re() << setw(15) << temp_coef.Im() << endl;
			} else {
				temp_coef(temp_coef1,temp_coef2*PI/180,1);
				cout << left << setw(20) << resonant_channel_string[i] << setw(15) << temp_coef.Rho() << setw(15) << temp_coef.Theta()*180/PI << endl;
			}
			coefs.push_back(temp_coef);
		}
	}
	cout << "##########################################################" << endl;

	cout << endl;
	cout << coefs.size() << " resonances in total" << endl;
	cout << endl;
	cout << "------------------------------------------------------------------------" << endl;
	cout << endl;

	cout << "Background generation parameters: " << endl;
	cout << endl;
	cout << "########################################" << endl;
	cout << left << setw(25) << "Background Component"  << setw(15) << "Fraction" << endl;
	cout << "########################################" << endl;

	// Reads the background coefficients

	for (int i=0; i<available_bkg_components; i++) {

		// Defines the identification strings associated with each background parameter 

		str_fraction = bkg_component_string[i];
		str_fraction += "_fraction"; 

		// Checks if a fraction is defined for the current background component. If it is, the program will include the background component in the model

		bkg_components.push_back(input_txt_file->Defined(str_fraction.c_str()));

		// Reads the fraction

		temp_bkg_coef = input_txt_file->GetValue(str_fraction.c_str(),0.);

		// Checks if fraction is not zero. If it is, the background component will be excluded from the model

		if (temp_bkg_coef == 0) bkg_components[i] = 0;

		// Fills the fractions that will be used in the model 

		if(bkg_components.at(i) == 1){
			bkg_coefs.push_back(temp_bkg_coef);
			cout << left << setw(25) << bkg_component_string[i] << setw(15) << temp_bkg_coef << endl;
		}
	}
	cout << "########################################" << endl;
	cout << endl;
	cout << bkg_coefs.size() << " background components in total" << endl;
	cout << endl;
	cout << "------------------------------------------------------------------------" << endl;
	cout << endl;

	// Reads the background fraction

	bkg_fraction = input_txt_file->GetValue("bkg_fraction",0.);;
	cout << "Background fraction: "<< bkg_fraction << endl;
	cout << endl;
	cout << "------------------------------------------------------------------------" << endl;
	cout << endl;
	cout << endl;

	//  Reads the name of input and output files for pwa analisys

	number_of_pwa_bins = input_txt_file->GetValue("number_of_pwa_bins",0);
	pwa_txt_file = input_txt_file->GetValue("pwa_txt_file","pwa_coefs.txt");

	// Reads the input file for pwa analisys

	ifstream infile(pwa_txt_file.c_str(),ios::in);

	if (!infile.good()) {
		std::cout << "ReadParameters - ERROR: File not found - " << pwa_txt_file << std::endl;
		exit(-1);
	}

	for (int i=0; i<number_of_pwa_bins; i++) {
		infile >> temp_mKK;
		infile >> temp_coef1;
		infile >> temp_coef2;
		infile >> temp_pwa_fix;
		if (real_and_imaginary) temp_coef(temp_coef1,temp_coef2);
		else temp_coef(temp_coef1,temp_coef2*PI/180,1);
		pwa_coefs.push_back(temp_coef);
		KK_bin_limits.push_back(temp_mKK);

	}
	/////////////////////////////////////////////////////////////////////////////////////////////
	// Reads the acceptance files
	UseAcceptance = input_txt_file->GetValue("UseAcceptance",0);
	UseAcceptance?cout<<"Fiting with Acceptance"<<endl:cout<<" Fiting without Acceptance"<<endl;
	cout << endl;
	cout << "------------------------------------------------------------------------" << endl;
	cout << endl;
	if(UseAcceptance){
		Acceptance_Ntuple_Name = input_txt_file->GetValue("Acceptance_Ntuple_Name","fAcceptance.root");
		Acceptance_Histo_Name    = input_txt_file->GetValue("Acceptance_Histo_Name","hAcceptance");
		cout<<" Reading acceptance histo "<<Acceptance_Histo_Name<<" from file "<<Acceptance_Ntuple_Name<<" ---- oK! "<< endl;
		///////////////////////////////////////////voy a leer la acept desde aqui//////////////////////////////////////////////////
		cout<<"s12MinForAcc "<<s12MinForAcc << "s12MaxForAcc  "<< s12MaxForAcc<<endl;

		AccNtpFile = new TFile(Acceptance_Ntuple_Name.c_str());
		if (AccNtpFile->IsZombie()) {std::cout << "Acceptance - Error opening file " << Acceptance_Ntuple_Name << std::endl;exit(-1);}    			 	
		AccHist = (TH2D*)AccNtpFile->Get(Acceptance_Histo_Name.c_str());
		TAxis* s12AxisForAcc = AccHist->GetYaxis();
		TAxis* s13AxisForAcc = AccHist->GetXaxis();
		s12MinForAcc = s12AxisForAcc->GetXmin();
		s12MaxForAcc = s12AxisForAcc->GetXmax();
		s13MinForAcc = s13AxisForAcc->GetXmin();
		s13MaxForAcc = s13AxisForAcc->GetXmax();
		NBins12ForAcc = AccHist->GetNbinsY();
		NBins13ForAcc = AccHist->GetNbinsX();
		s12ExtForAcc = s12MaxForAcc - s12MinForAcc;
		s13ExtForAcc = s13MaxForAcc - s13MinForAcc;

		cout<<"s12MinForAcc "<<s12MinForAcc << " s12MaxForAcc  "<< s12MaxForAcc <<" NBins12ForAcc "<<NBins12ForAcc<<" s12ExtForAcc " << s12ExtForAcc<<endl;

		cout<<"s13MinForAcc "<<s13MinForAcc << " s13MaxForAcc  "<< s13MaxForAcc <<" NBins13ForAcc "<<NBins13ForAcc<<" s13ExtForAcc " << s13ExtForAcc<<endl;

	}
	/////////////////////////////////////////////////////////////////////////////////////////////
	// Reads the back files
	/////////////////////////////////////////////////////////////////////////////////////////////

	UseBackHisto = input_txt_file->GetValue("UseBackHisto",0);
	UseBackHisto?cout<<"Fiting background with histo"<<endl:cout<<" Fiting background without histo"<<endl;
	cout << endl;
	cout << "------------------------------------------------------------------------" << endl;
	cout << endl;
	if(UseBackHisto){
		Back_Ntuple_Name = input_txt_file->GetValue("Back_Ntuple_Name","fBack.root");
		Back_Histo_Name    = input_txt_file->GetValue("Back_Histo_Name","hBack");
		cout<<" Reading acceptance histo "<<Back_Histo_Name<<" from file "<<Back_Ntuple_Name<<" ---- oK! "<< endl;
		///////////////////////////////////////////voy a leer la acept desde aqui//////////////////////////////////////////////////
		cout<<"s12MinForBack "<<s12MinForBack << " s12MaxForBack  "<< s12MaxForBack<<endl;

		BackNtpFile = new TFile(Back_Ntuple_Name.c_str());

		cout<<" ntuple bien"<<endl;
		if (BackNtpFile->IsZombie()) {std::cout << "Background - Error opening file " << Back_Ntuple_Name << std::endl;exit(-1);}
		cout<<" después de zombie"<<endl;

		BackHist = (TH2D*)BackNtpFile->Get(Back_Histo_Name.c_str());
		cout<<" después de bachist"<<endl;
		TAxis* s12AxisForBack = BackHist->GetYaxis();
		TAxis* s13AxisForBack = BackHist->GetXaxis();
		cout<<" después de axis"<<endl;
		s12MinForBack = s12AxisForBack->GetXmin();
		cout<<" después de minsS12"<<endl;
		s12MaxForBack = s12AxisForBack->GetXmax();
		cout<<" después de maxS12"<<endl;
		s13MinForBack = s13AxisForBack->GetXmin();
		cout<<" después de minsS13"<<endl;
		s13MaxForBack = s13AxisForBack->GetXmax();
		cout<<" después de maxS13"<<endl;
		NBins12ForBack = BackHist->GetNbinsY();
		NBins13ForBack = BackHist->GetNbinsX();
		s12ExtForBack = s12MaxForBack - s12MinForBack;
		s13ExtForBack = s13MaxForBack - s13MinForBack;
		cout<<" después de todo back"<<endl;
		cout<<"s12MinForBack "<<s12MinForBack << " s12MaxForBack  "<< s12MaxForBack <<" NBins12ForBack "<<NBins12ForBack<<" s12ExtForBack " << s12ExtForBack<<endl;

		cout<<"s13MinForBack "<<s13MinForBack << " s13MaxForBack  "<< s13MaxForBack <<" NBins13ForBack "<<NBins13ForBack<<" s13ExtForAcc " << s13ExtForBack<<endl;

	}
	/////////////////////////////////////////////////////////////////////////////////////////////



}
