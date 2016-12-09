#include <stdio.h>
#include <math.h>
#include <string>
#include <TComplex.h>
#include <TEnv.h>
#include "../src/GenericFunctions.h"

using namespace std;

void Read_Parameters(string input_txt_file_name, TFitter* minimizer, bool &real_and_imaginary, int &final_state, bool &is_gaussian, double &Bkg_par1, double &Bkg_par2, int &number_of_samples, 
		vector<int> &resonances, vector<int> &bkg_components, int &number_of_resonances, int &number_of_bkg_components, double &bkg_fraction, string &input_ntuple_name, int &seed, 
		vector<int> &fix_parameter_index, int &number_of_pwa_bins, string &output_pwa_txt_file, vector<TComplex> &pwa_coefs, vector<double> &KK_bin_limits, bool &UsesWeights, 
		bool &UseAcceptance,  string &Acceptance_Ntuple_Name, string &Acceptance_Histo_Name, bool &UseBackHisto,  string &Back_Ntuple_Name, string &Back_Histo_Name){

	int available_resonances = AVAILABLE_RESONANCES;
	int available_bkg_components = AVAILABLE_BKG_COMPONENTS;
	string str_coef1, str_coef1_lower_limit, str_coef1_upper_limit, str_coef1_fix, str_coef2, str_coef2_lower_limit, str_coef2_upper_limit, str_coef2_fix, str_fraction,
	       str_fraction_lower_limit, str_fraction_upper_limit, str_fraction_fix, pwa_txt_file, pwa_coef1, pwa_coef2, str_res_mass, str_res_mass_lower_limit,
	       str_res_mass_upper_limit, str_res_mass_fix, str_res_width, str_res_width_lower_limit, str_res_width_upper_limit, str_res_width_fix, str_res_extra_par_1, 
	       str_res_extra_par_1_lower_limit, str_res_extra_par_1_upper_limit, str_res_extra_par_1_fix, str_res_extra_par_2, str_res_extra_par_2_lower_limit, 
	       str_res_extra_par_2_upper_limit, str_res_extra_par_2_fix, str_res_extra_par_3, str_res_extra_par_3_lower_limit, str_res_extra_par_3_upper_limit, 
	       str_res_extra_par_3_fix, str_res_extra_par_4, str_res_extra_par_4_lower_limit, str_res_extra_par_4_upper_limit, str_res_extra_par_4_fix;
	int temp_coef1_fix, temp_coef2_fix, temp_mass_fix, temp_width_fix, temp_bkg_fix, temp_pwa_fix, temp_res_extra_par_1_fix, temp_res_extra_par_2_fix, temp_res_extra_par_3_fix, 
	    temp_res_extra_par_4_fix;
	double temp_coef1, temp_coef2, temp_mass, temp_width, temp_bkg_coef, temp_coef1_lower_limit, temp_coef2_lower_limit, temp_coef1_upper_limit, temp_coef2_upper_limit, temp_mass_lower_limit,
	       temp_mass_upper_limit, temp_width_lower_limit, temp_width_upper_limit, temp_mKK, temp_res_extra_par_1, temp_res_extra_par_1_upper_limit, temp_res_extra_par_1_lower_limit, 
	       temp_res_extra_par_2, temp_res_extra_par_2_upper_limit, temp_res_extra_par_2_lower_limit, temp_res_extra_par_3, temp_res_extra_par_3_upper_limit, temp_res_extra_par_3_lower_limit, 
	       temp_res_extra_par_4, temp_res_extra_par_4_upper_limit, temp_res_extra_par_4_lower_limit;
	TComplex temp_coef;

	TEnv *input_txt_file = new TEnv(input_txt_file_name.c_str());

	// This function reads the parameters to be used in the fit 

	// Checks if the input .txt file was loaded correctly

	if (input_txt_file->IsZombie()) {
		std::cout << "ReadParameters - ERROR: File not found - " << input_txt_file_name << std::endl;
		exit(-1);
	}

	// Reads the name of the input file

	input_ntuple_name = input_txt_file->GetValue("input_ntuple_name", "ntuple");

	cout << "------------------------------------------------------------------------" << endl;
	cout << endl;
	cout << "input file name: " << input_ntuple_name << endl;
	cout << endl;

	// Reads the number of samples

	number_of_samples = input_txt_file->GetValue("number_of_samples", 1);

	cout << "number of samples: " << number_of_samples << endl;
	cout << endl;
	cout << "------------------------------------------------------------------------" << endl;
	cout << endl;

	// Reads an integer associated with the final state of the sample

	final_state = input_txt_file->GetValue("final_state",-1);

	// Reads the seed

	seed = input_txt_file->GetValue("seed",0);

	// Reads the choice for the mass distribution: delta or gaussian

	is_gaussian = input_txt_file->GetValue("is_gaussian",0);

	if(is_gaussian){
		cout << "Gaussian mass distribution";
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
	// Checks if the final state is kkpi or pipipi

	if (final_state == 0) {
		cout << "The chosen final state was kkpi" << endl;
		cout << endl;
		cout << "------------------------------------------------------------------------" << endl;
		cout << endl;
	}	
	else if (final_state == 1) {
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
	else if (final_state == 3) {
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

	// Reads if  sWeights are used

	UsesWeights = input_txt_file->GetValue("UsesWeights",0);
	UsesWeights?cout<<"Fiting with sWeights"<<endl:cout<<" Fiting without sWeights"<<endl;
	cout << endl;
	cout << "------------------------------------------------------------------------" << endl;
	cout << endl;

	// Reads the fit parameters

	cout << "Signal fit parameters: " << endl;

	for (int i=0; i<available_resonances; i++) {

		// Defines the identification strings associated with each parameter 

		if (real_and_imaginary) {
			str_coef1 = resonant_channel_string[i];
			str_coef1 += "_re";
			str_coef1_lower_limit = resonant_channel_string[i];
			str_coef1_lower_limit += "_re_lower_limit";
			str_coef1_upper_limit = resonant_channel_string[i];
			str_coef1_upper_limit += "_re_upper_limit";
			str_coef1_fix = resonant_channel_string[i];
			str_coef1_fix += "_re_fix" ;
			str_coef2 = resonant_channel_string[i];
			str_coef2 += "_im";
			str_coef2_lower_limit = resonant_channel_string[i];
			str_coef2_lower_limit += "_im_lower_limit";
			str_coef2_upper_limit = resonant_channel_string[i];
			str_coef2_upper_limit += "_im_upper_limit";
			str_coef2_fix = resonant_channel_string[i];
			str_coef2_fix += "_im_fix"; 
		} else {
			str_coef1 = resonant_channel_string[i];
			str_coef1 += "_amp";
			str_coef1_lower_limit = resonant_channel_string[i];
			str_coef1_lower_limit += "_amp_lower_limit";
			str_coef1_upper_limit = resonant_channel_string[i];
			str_coef1_upper_limit += "_amp_upper_limit";
			str_coef1_fix = resonant_channel_string[i];
			str_coef1_fix += "_amp_fix" ;
			str_coef2 = resonant_channel_string[i];
			str_coef2 += "_phs";
			str_coef2_lower_limit = resonant_channel_string[i];
			str_coef2_lower_limit += "_phs_lower_limit";
			str_coef2_upper_limit = resonant_channel_string[i];
			str_coef2_upper_limit += "_phs_upper_limit";
			str_coef2_fix = resonant_channel_string[i];
			str_coef2_fix += "_phs_fix"; 

		}
		str_res_mass = resonant_channel_string[i];
		str_res_mass += "_mass";
		str_res_mass_lower_limit = resonant_channel_string[i];
		str_res_mass_lower_limit += "_mass_lower_limit";
		str_res_mass_upper_limit = resonant_channel_string[i];
		str_res_mass_upper_limit += "_mass_upper_limit";
		str_res_mass_fix = resonant_channel_string[i];
		str_res_mass_fix += "_mass_fix" ;
		str_res_width = resonant_channel_string[i];
		str_res_width += "_width";
		str_res_width_lower_limit = resonant_channel_string[i];
		str_res_width_lower_limit += "_width_lower_limit";
		str_res_width_upper_limit = resonant_channel_string[i];
		str_res_width_upper_limit += "_width_upper_limit";
		str_res_width_fix = resonant_channel_string[i];
		str_res_width_fix += "_width_fix" ;

		str_res_extra_par_1 = resonant_channel_string[i];
		str_res_extra_par_1 += "_res_extra_par_1";
		str_res_extra_par_1_lower_limit = resonant_channel_string[i];
		str_res_extra_par_1_lower_limit += "_res_extra_par_1_lower_limit";
		str_res_extra_par_1_upper_limit = resonant_channel_string[i];
		str_res_extra_par_1_upper_limit += "_res_extra_par_1_upper_limit";
		str_res_extra_par_1_fix = resonant_channel_string[i];
		str_res_extra_par_1_fix += "_res_extra_par_1_fix" ;
		str_res_extra_par_2 = resonant_channel_string[i];
		str_res_extra_par_2 += "_res_extra_par_2";
		str_res_extra_par_2_lower_limit = resonant_channel_string[i];
		str_res_extra_par_2_lower_limit += "_res_extra_par_2_lower_limit";
		str_res_extra_par_2_upper_limit = resonant_channel_string[i];
		str_res_extra_par_2_upper_limit += "_res_extra_par_2_upper_limit";
		str_res_extra_par_2_fix = resonant_channel_string[i];
		str_res_extra_par_2_fix += "_res_extra_par_2_fix" ;
		str_res_extra_par_3 = resonant_channel_string[i];
		str_res_extra_par_3 += "_res_extra_par_3";
		str_res_extra_par_3_lower_limit = resonant_channel_string[i];
		str_res_extra_par_3_lower_limit += "_res_extra_par_3_lower_limit";
		str_res_extra_par_3_upper_limit = resonant_channel_string[i];
		str_res_extra_par_3_upper_limit += "_res_extra_par_3_upper_limit";
		str_res_extra_par_3_fix = resonant_channel_string[i];
		str_res_extra_par_3_fix += "_res_extra_par_3_fix" ;
		str_res_extra_par_4 = resonant_channel_string[i];
		str_res_extra_par_4 += "_res_extra_par_4";
		str_res_extra_par_4_lower_limit = resonant_channel_string[i];
		str_res_extra_par_4_lower_limit += "_res_extra_par_4_lower_limit";
		str_res_extra_par_4_upper_limit = resonant_channel_string[i];
		str_res_extra_par_4_upper_limit += "_res_extra_par_4_upper_limit";
		str_res_extra_par_4_fix = resonant_channel_string[i];
		str_res_extra_par_4_fix += "_res_extra_par_4_fix" ;

		// Reads the magnitude, phase, limits and if it is fix or not

		temp_coef1 = input_txt_file->GetValue(str_coef1.c_str(),0.);
        	if (real_and_imaginary) {
		temp_coef1_lower_limit = input_txt_file->GetValue(str_coef1_lower_limit.c_str(),-1000);
		temp_coef1_upper_limit = input_txt_file->GetValue(str_coef1_upper_limit.c_str(),1000);
            }
        
            else {
                temp_coef1_lower_limit = input_txt_file->GetValue(str_coef1_lower_limit.c_str(),0.);
                temp_coef1_upper_limit = input_txt_file->GetValue(str_coef1_upper_limit.c_str(),1000);
                
            }
		temp_coef1_fix = input_txt_file->GetValue(str_coef1_fix.c_str(),1);
		temp_coef2 = input_txt_file->GetValue(str_coef2.c_str(),0.);
        if (real_and_imaginary) {
		temp_coef2_lower_limit = input_txt_file->GetValue(str_coef2_lower_limit.c_str(),-1000.);
		temp_coef2_upper_limit = input_txt_file->GetValue(str_coef2_upper_limit.c_str(),1000.);
        }
        else {
            temp_coef2_lower_limit = input_txt_file->GetValue(str_coef2_lower_limit.c_str(),-180.);
            temp_coef2_upper_limit = input_txt_file->GetValue(str_coef2_upper_limit.c_str(),180.);
        }
		temp_coef2_fix = input_txt_file->GetValue(str_coef2_fix.c_str(),1);
		temp_mass = input_txt_file->GetValue(str_res_mass.c_str(),resonant_channel_mass[i]);
		temp_mass_lower_limit = input_txt_file->GetValue(str_res_mass_lower_limit.c_str(), 0.0);
		temp_mass_upper_limit = input_txt_file->GetValue(str_res_mass_upper_limit.c_str(),resonant_channel_mass[i] + 1);
		temp_mass_fix = input_txt_file->GetValue(str_res_mass_fix.c_str(),1);
		temp_width = input_txt_file->GetValue(str_res_width.c_str(),resonant_channel_width[i]);
		temp_width_lower_limit = input_txt_file->GetValue(str_res_width_lower_limit.c_str(), 0);
		temp_width_upper_limit = input_txt_file->GetValue(str_res_width_upper_limit.c_str(),resonant_channel_width[i] + 1);
		temp_width_fix = input_txt_file->GetValue(str_res_width_fix.c_str(),1);

		temp_res_extra_par_1 = input_txt_file->GetValue(str_res_extra_par_1.c_str(),-99999999.);
		temp_res_extra_par_1_lower_limit = input_txt_file->GetValue(str_res_extra_par_1_lower_limit.c_str(),-100000000);
		temp_res_extra_par_1_upper_limit = input_txt_file->GetValue(str_res_extra_par_1_upper_limit.c_str(),-99999998);
		temp_res_extra_par_1_fix = input_txt_file->GetValue(str_res_extra_par_1_fix.c_str(),1);
		temp_res_extra_par_2 = input_txt_file->GetValue(str_res_extra_par_2.c_str(),-99999999.);
		temp_res_extra_par_2_lower_limit = input_txt_file->GetValue(str_res_extra_par_2_lower_limit.c_str(),-100000000);
		temp_res_extra_par_2_upper_limit = input_txt_file->GetValue(str_res_extra_par_2_upper_limit.c_str(),-99999998);
		temp_res_extra_par_2_fix = input_txt_file->GetValue(str_res_extra_par_2_fix.c_str(),1);
		temp_res_extra_par_3 = input_txt_file->GetValue(str_res_extra_par_3.c_str(),-99999999.);
		temp_res_extra_par_3_lower_limit = input_txt_file->GetValue(str_res_extra_par_3_lower_limit.c_str(),-100000000);
		temp_res_extra_par_3_upper_limit = input_txt_file->GetValue(str_res_extra_par_3_upper_limit.c_str(),-99999998);
		temp_res_extra_par_3_fix = input_txt_file->GetValue(str_res_extra_par_3_fix.c_str(),1);
		temp_res_extra_par_4 = input_txt_file->GetValue(str_res_extra_par_4.c_str(),-99999999.);
		temp_res_extra_par_4_lower_limit = input_txt_file->GetValue(str_res_extra_par_4_lower_limit.c_str(),-100000000);
		temp_res_extra_par_4_upper_limit = input_txt_file->GetValue(str_res_extra_par_4_upper_limit.c_str(),-99999998);
		temp_res_extra_par_4_fix = input_txt_file->GetValue(str_res_extra_par_4_fix.c_str(),1);

		// Checks if real and imaginary parts are fixed at zero. If they are, the resonant channel will be excluded from the model

		if (temp_coef1 == 0 && temp_coef1_fix == 1 && temp_coef2 == 0 && temp_coef2_fix == 1) resonances.push_back(0);
		else resonances.push_back(1);

		// Fills the amplitudes, phases and associated parameters that will be used in the model 

		if(resonances.at(i) == 1){
			if (real_and_imaginary) {
				temp_coef(temp_coef1, temp_coef2);
				minimizer->SetParameter(8*number_of_resonances,str_coef1.c_str(),temp_coef.Re(), 0.00001,temp_coef1_lower_limit,temp_coef1_upper_limit);
				minimizer->SetParameter(8*number_of_resonances + 1,str_coef2.c_str(),temp_coef.Im(), 0.00001,temp_coef2_lower_limit,temp_coef2_upper_limit);
			} else{
				temp_coef(temp_coef1, temp_coef2*PI/180, 1);
				minimizer->SetParameter(8*number_of_resonances,str_coef1.c_str(),temp_coef.Rho(), 0.00001,temp_coef1_lower_limit,temp_coef1_upper_limit);
				minimizer->SetParameter(8*number_of_resonances + 1,str_coef2.c_str(),temp_coef.Theta(), 0.00001,temp_coef2_lower_limit*PI/180,temp_coef2_upper_limit*PI/180);
			}
			minimizer->SetParameter(8*number_of_resonances + 2,str_res_mass.c_str(),temp_mass, 0.00001,temp_mass_lower_limit,temp_mass_upper_limit);
			minimizer->SetParameter(8*number_of_resonances + 3,str_res_width.c_str(),temp_width, 0.00001,temp_width_lower_limit,temp_width_upper_limit);
			cout << "temp_res_extra_par_1 = " << temp_res_extra_par_1 << endl;
			cout << "temp_res_extra_par_2 = " << temp_res_extra_par_2 << endl;
			cout << "temp_res_extra_par_3 = " << temp_res_extra_par_3 << endl;
			cout << "temp_res_extra_par_4 = " << temp_res_extra_par_4 << endl;
			minimizer->SetParameter(8*number_of_resonances + 4,str_res_extra_par_1.c_str(),temp_res_extra_par_1, 0.00001,temp_res_extra_par_1_lower_limit,temp_res_extra_par_1_upper_limit);
			minimizer->SetParameter(8*number_of_resonances + 5,str_res_extra_par_2.c_str(),temp_res_extra_par_2, 0.00001,temp_res_extra_par_2_lower_limit,temp_res_extra_par_2_upper_limit);
			minimizer->SetParameter(8*number_of_resonances + 6,str_res_extra_par_3.c_str(),temp_res_extra_par_3, 0.00001,temp_res_extra_par_3_lower_limit,temp_res_extra_par_3_upper_limit);
			minimizer->SetParameter(8*number_of_resonances + 7,str_res_extra_par_4.c_str(),temp_res_extra_par_4, 0.00001,temp_res_extra_par_4_lower_limit,temp_res_extra_par_4_upper_limit);
			if (temp_coef1_fix == 1){
				minimizer->FixParameter(8*number_of_resonances);
				fix_parameter_index.push_back(8*number_of_resonances);
			}
			if (temp_coef2_fix == 1){
				minimizer->FixParameter(8*number_of_resonances + 1);
				fix_parameter_index.push_back(8*number_of_resonances + 1);
			}
			if (temp_mass_fix == 1){
				minimizer->FixParameter(8*number_of_resonances + 2);
				fix_parameter_index.push_back(8*number_of_resonances + 2);
			}
			if (temp_width_fix == 1){
				minimizer->FixParameter(8*number_of_resonances + 3);
				fix_parameter_index.push_back(8*number_of_resonances + 3);
			}
			if (temp_res_extra_par_1_fix == 1){
				minimizer->FixParameter(8*number_of_resonances + 4);
				fix_parameter_index.push_back(8*number_of_resonances + 4);
			}
			if (temp_res_extra_par_2_fix == 1){
				minimizer->FixParameter(8*number_of_resonances + 5);
				fix_parameter_index.push_back(8*number_of_resonances + 5);
			}
			if (temp_res_extra_par_3_fix == 1){
				minimizer->FixParameter(8*number_of_resonances + 6);
				fix_parameter_index.push_back(8*number_of_resonances + 6);
			}
			if (temp_res_extra_par_4_fix == 1){
				minimizer->FixParameter(8*number_of_resonances + 7);
				fix_parameter_index.push_back(8*number_of_resonances + 7);
			}
			number_of_resonances++;
		}
	}
	cout << endl;
	cout << number_of_resonances << " resonances in total" << endl;
	cout << endl;
	cout << "------------------------------------------------------------------------" << endl;
	cout << endl;

	cout << "Background fit parameters: " << endl;

	// Reads the background coefficients

	for (int i=0; i<available_bkg_components; i++) {

		// Defines the identification strings associated with each background parameter 

		str_fraction = bkg_component_string[i];
		str_fraction += "_fraction"; 
		str_fraction_fix = bkg_component_string[i];
		str_fraction_fix += "_fraction_fix"; 

		// Checks if a fraction is defined for the current background component. If it is, the program will include the background component in the model

		bkg_components.push_back(input_txt_file->Defined(str_fraction.c_str()));				

		// Reads the fraction and if it is fix or not

		temp_bkg_coef = input_txt_file->GetValue(str_fraction.c_str(),0.);
		temp_bkg_fix = input_txt_file->GetValue(str_fraction_fix.c_str(),0);

		// Checks if fraction is fixed at zero. If it is, the background component will be excluded from the model

		if (temp_bkg_coef == 0 && temp_bkg_fix == 1) bkg_components[i] = 0;

		// Fills the fractions that will be used in the model 

		if(bkg_components.at(i) == 1){
			minimizer->SetParameter(number_of_bkg_components + 8*number_of_resonances,str_fraction.c_str(),temp_bkg_coef, 0.00,0,1);
			if (temp_bkg_fix == 1){
				minimizer->FixParameter(number_of_bkg_components + 8*number_of_resonances);
			}
			number_of_bkg_components++;
		}
	}
	cout << endl;
	cout << number_of_bkg_components << " background components in total" << endl;
	cout << endl;
	cout << "------------------------------------------------------------------------" << endl;
	cout << endl;

	// Reads the background fractions

	bkg_fraction = input_txt_file->GetValue("bkg_fraction",0.);;
	cout << "Background fraction: "<< bkg_fraction << endl;
	cout << endl;
	cout << "------------------------------------------------------------------------" << endl;
	cout << endl;
	cout << endl;

	// Reads the input file for pwa analisys

	pwa_txt_file = input_txt_file->GetValue("pwa_txt_file","pwa_coefs.txt");
	output_pwa_txt_file = input_txt_file->GetValue("output_pwa_txt_file","pwa_output.txt");

	ifstream infile(pwa_txt_file.c_str(),ios::in);

	if (!infile.good()) {
		std::cout << "ReadParameters - ERROR: File not found - " << pwa_txt_file << std::endl;
		exit(-1);
	}

	number_of_pwa_bins = input_txt_file->GetValue("number_of_pwa_bins",0);

	for (int i=0; i<number_of_pwa_bins; i++) {
		infile >> temp_mKK;
		infile >> temp_coef1;
		infile >> temp_coef2;
		infile >> temp_pwa_fix;
		KK_bin_limits.push_back(temp_mKK);
		char buff[128];
		sprintf(buff, "%i", i);
		if (real_and_imaginary){
			temp_coef(temp_coef1, temp_coef2);
			pwa_coefs.push_back(temp_coef);
			pwa_coef1 = "pwa_re ";
			pwa_coef1 += buff;
			pwa_coef2 = "pwa_im ";
			pwa_coef2 += buff;
			minimizer->SetParameter(2*i + number_of_bkg_components + 8*number_of_resonances,pwa_coef1.c_str(),temp_coef.Re(), 0.00001,-100.,100.);
			minimizer->SetParameter(2*i + 1 + number_of_bkg_components + 8*number_of_resonances,pwa_coef2.c_str(),temp_coef.Im(), 0.00001,-100,100);
		}else {
			temp_coef(temp_coef1, temp_coef2*PI/180, 1);
			pwa_coefs.push_back(temp_coef);
			pwa_coef1 = "pwa_amp ";
			pwa_coef1 += buff;
			pwa_coef2 = "pwa_phs ";
			pwa_coef2 += buff;
			minimizer->SetParameter(2*i + number_of_bkg_components + 8*number_of_resonances,pwa_coef1.c_str(),temp_coef.Rho(), 0.00001,0.,1000.);
			minimizer->SetParameter(2*i + 1 + number_of_bkg_components + 8*number_of_resonances,pwa_coef2.c_str(),temp_coef.Theta(), 0.00001,-180,180);			
		}


		if (temp_pwa_fix == 1) {
			minimizer->FixParameter(2*i + number_of_bkg_components + 8*number_of_resonances);
			minimizer->FixParameter(2*i + 1 + number_of_bkg_components + 8*number_of_resonances);
		}
	}

	cout << "------------------------------------------------------------------------" << endl;
	cout << endl;
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
		cout<<"s12MinForAcc "<<s12MinForAcc << ", s12MaxForAcc  "<< s12MaxForAcc<<endl;

		AccNtpFile = new TFile(Acceptance_Ntuple_Name.c_str());
		cout << "teste 0 " << endl;
		if (AccNtpFile->IsZombie()) {std::cout << "Acceptance - Error opening file " << Acceptance_Ntuple_Name << std::endl;exit(-1);}   
		cout << "teste 1 " << endl;
		AccHist = (TH2D*)AccNtpFile->Get(Acceptance_Histo_Name.c_str());			
		cout << AccHist->GetMaximum();	
		cout << "teste 2 " << endl;
		TAxis* s12AxisForAcc = AccHist->GetYaxis();
		TAxis* s13AxisForAcc = AccHist->GetXaxis();
		cout << "teste 3 " << endl;
		s12MinForAcc = s12AxisForAcc->GetXmin();
		s12MaxForAcc = s12AxisForAcc->GetXmax();
		s13MinForAcc = s13AxisForAcc->GetXmin();
		s13MaxForAcc = s13AxisForAcc->GetXmax();
		cout << "teste 4 " << endl;
		NBins12ForAcc = AccHist->GetNbinsY();
		NBins13ForAcc = AccHist->GetNbinsX();
		cout << "teste 5 " << endl;
		s12ExtForAcc = s12MaxForAcc - s12MinForAcc;
		s13ExtForAcc = s13MaxForAcc - s13MinForAcc;				
		cout << "teste 6 " << endl;

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
		cout<<"s12MinForBack "<<s12MinForBack << "s12MaxForBack  "<< s12MaxForBack<<endl;

		BackNtpFile = new TFile(Back_Ntuple_Name.c_str());
		if (BackNtpFile->IsZombie()) {std::cout << "Background - Error opening file " << Back_Ntuple_Name << std::endl;exit(-1);}

		BackHist = (TH2D*)BackNtpFile->Get(Back_Histo_Name.c_str());
		TAxis* s12AxisForBack = BackHist->GetYaxis();
		TAxis* s13AxisForBack = BackHist->GetXaxis();
		s12MinForBack = s12AxisForBack->GetXmin();
		s12MaxForBack = s12AxisForBack->GetXmax();
		s13MinForBack = s13AxisForBack->GetXmin();
		s13MaxForBack = s13AxisForBack->GetXmax();
		NBins12ForBack = BackHist->GetNbinsY();
		NBins13ForBack = BackHist->GetNbinsX();
		s12ExtForBack = s12MaxForBack - s12MinForBack;
		s13ExtForBack = s13MaxForBack - s13MinForBack;

		cout<<"s12MinForBack "<<s12MinForBack << " s12MaxForBack  "<< s12MaxForBack <<" NBins12ForBack "<<NBins12ForBack<<" s12ExtForBack " << s12ExtForBack<<endl;

		cout<<"s13MinForBack "<<s13MinForBack << " s13MaxForBack  "<< s13MaxForBack <<" NBins13ForBack "<<NBins13ForBack<<" s13ExtForAcc " << s13ExtForBack<<endl;

	}
	/////////////////////////////////////////////////////////////////////////////////////////////

}
