#include <stdlib.h>
#include <iostream>   // std::cout
#include <string>  

void makePlots(){

int NN = 500; // number of files
int n_res = 3; // number of resonaces

double b_status,b_edm,b_chi2,b_fcn;
double b_res1_amp,b_res1_amp_err,b_res1_phs,b_res1_phs_err,b_res1_mass,b_res1_mass_err,b_res1_width,b_res1_width_err,b_res1_frac,b_res1_frac_err;
double b_res2_amp,b_res2_amp_err,b_res2_phs,b_res2_phs_err,b_res2_mass,b_res2_mass_err,b_res2_width,b_res2_width_err,b_res2_frac,b_res2_frac_err;
double b_res3_amp,b_res3_amp_err,b_res3_phs,b_res3_phs_err,b_res3_mass,b_res3_mass_err,b_res3_width,b_res3_width_err,b_res3_frac,b_res3_frac_err;

TTree tree("tree","tree");
 tree.Branch("status",&b_status,"status/D");
 tree.Branch("edm"   ,&b_edm,"edm/D");
 tree.Branch("chi2"  ,&b_chi2,"chi2/D");
 tree.Branch("fcn"   ,&b_fcn,"fcn/D");

 tree.Branch("res1_amp",       &b_res1_amp       ,"res1_amp/D");
 tree.Branch("res1_amp_err",   &b_res1_amp_err   ,"res1_amp_err/D");
 tree.Branch("res1_phs",       &b_res1_phs       ,"res1_phs/D");
 tree.Branch("res1_phs_err",   &b_res1_phs_err   ,"res1_phs_err/D");
 tree.Branch("res1_mass",     &b_res1_mass     ,"res1_mass/D");
 tree.Branch("res1_mass_err", &b_res1_mass_err ,"res1_mass_err/D");
 tree.Branch("res1_width",    &b_res1_width    ,"res1_width/D");
 tree.Branch("res1_width_err",&b_res1_width_err,"res1_width_err/D");
 tree.Branch("res1_frac",     &b_res1_frac     ,"res1_frac/D");
 tree.Branch("res1_frac_err", &b_res1_frac_err ,"res1_frac_err/D");

 tree.Branch("res2_amp",       &b_res2_amp       ,"res2_amp/D");
 tree.Branch("res2_amp_err",   &b_res2_amp_err   ,"res2_amp_err/D");
 tree.Branch("res2_phs",       &b_res2_phs       ,"res2_phs/D");
 tree.Branch("res2_phs_err",   &b_res2_phs_err   ,"res2_phs_err/D");
 tree.Branch("res2_mass",     &b_res2_mass     ,"res2_mass/D");
 tree.Branch("res2_mass_err", &b_res2_mass_err ,"res2_mass_err/D");
 tree.Branch("res2_width",    &b_res2_width    ,"res2_width/D");
 tree.Branch("res2_width_err",&b_res2_width_err,"res2_width_err/D");
 tree.Branch("res2_frac",     &b_res2_frac     ,"res2_frac/D");
 tree.Branch("res2_frac_err", &b_res2_frac_err ,"res2_frac_err/D");

 tree.Branch("res3_amp",       &b_res3_amp       ,"res3_amp/D");
 tree.Branch("res3_amp_err",   &b_res3_amp_err   ,"res3_amp_err/D");
 tree.Branch("res3_phs",       &b_res3_phs       ,"res3_phs/D");
 tree.Branch("res3_phs_err",   &b_res3_phs_err   ,"res3_phs_err/D");
 tree.Branch("res3_mass",     &b_res3_mass     ,"res3_mass/D");
 tree.Branch("res3_mass_err", &b_res3_mass_err ,"res3_mass_err/D");
 tree.Branch("res3_width",    &b_res3_width    ,"res3_width/D");
 tree.Branch("res3_width_err",&b_res3_width_err,"res3_width_err/D");
 tree.Branch("res3_frac",     &b_res3_frac     ,"res3_frac/D");
 tree.Branch("res3_frac_err", &b_res3_frac_err ,"res3_frac_err/D");

for(int i=1 ; i<=NN; i++){ 

	TEnv *input_txt_file = new TEnv( Form("FitPars/fitparameters%i.txt",i) );
	//TEnv *input_txt_file = new TEnv( Form("check/check2_model3_notrandom/FitPars/fitparameters%i.txt",i) );

	if (input_txt_file->IsZombie()) {
		std::cout << "ERROR: File not found - " << input_txt_file_name << std::endl;
		exit(-1);
	}

	string reader;
	string status;
	string edm;
	string chi2;
	string fcn;

	string res1_amp,  res1_amp_error; 
	string res1_phs,  res1_phs_error; 
	string res1_mass, res1_mass_error; 
	string res1_width,res1_width_error; 
	string res1_frac, res1_frac_error; 

	string res2_amp,  res2_amp_error; 
	string res2_phs,  res2_phs_error; 
	string res2_mass, res2_mass_error; 
	string res2_width,res2_width_error; 
	string res2_frac, res2_frac_error; 

	string res3_amp,  res3_amp_error; 
	string res3_phs,  res3_phs_error; 
	string res3_mass, res3_mass_error; 
	string res3_width,res3_width_error; 
	string res3_frac, res3_frac_error; 

	//reader = input_txt_file->GetValue("a0+pi", "ERROR");
	reader = input_txt_file->GetValue("NR", "ERROR");
    	stringstream ss(reader); // Insert the string into a stream
	ss >> res1_amp >> res1_amp_error >> res1_phs >> res1_phs_error >> res1_mass >> res1_mass_error >> res1_width >> res1_width_error;
	cout << reader << endl;


	//reader = input_txt_file->GetValue("f0X+pi", "ERROR");
	reader = input_txt_file->GetValue("s_f0980+K", "ERROR");
    	stringstream ss(reader); // Insert the string into a stream
	ss >> res2_amp >> res2_amp_error >> res2_phs >> res2_phs_error >> res2_mass >> res2_mass_error >> res2_width >> res2_width_error;
	cout << reader << endl;


	reader = input_txt_file->GetValue("s_phi+K", "ERROR");
    	stringstream ss(reader); // Insert the string into a stream
	ss >> res3_amp >> res3_amp_error >> res3_phs >> res3_phs_error >> res3_mass >> res3_mass_error >> res3_width >> res3_width_error;
	cout << reader << endl;


	reader = input_txt_file->GetValue("status", "ERROR");
    	stringstream ss(reader); // Insert the string into a stream
	ss >> status;
	cout << reader << endl;

	reader = input_txt_file->GetValue("fcn", "ERROR");
    	stringstream ss(reader); // Insert the string into a stream
	ss >> fcn;
	cout << reader << endl;

	reader = input_txt_file->GetValue("frac_NR", "ERROR");
    	stringstream ss(reader); // Insert the string into a stream
	ss >> res1_frac >> res1_frac_error;
	cout << reader << endl;

	reader = input_txt_file->GetValue("frac_f_\{0\}\(980\)+K$", "ERROR");
    	stringstream ss(reader); // Insert the string into a stream
	ss >> res2_frac >> res2_frac_error;
	cout << reader << endl;

	reader = input_txt_file->GetValue("frac_#phi\_+K", "ERROR");
    	stringstream ss(reader); // Insert the string into a stream
	ss >> res3_frac >> res3_frac_error;
	cout << reader << endl;

	//reader = input_txt_file->GetValue("edm", "ERROR");
    	//stringstream ss(reader); // Insert the string into a stream
	//ss >> edm;
	//cout << reader << endl;

	//reader = input_txt_file->GetValue("chi2", "ERROR");
    	//stringstream ss(reader); // Insert the string into a stream
	//ss >> chi2;
	//cout << reader << endl;

	b_res1_amp        =  atof(res1_amp.c_str());
	b_res1_amp_err    =  atof(res1_amp_error.c_str());
	b_res1_phs        =  atof(res1_phs.c_str());
	b_res1_phs_err    =  atof(res1_phs_error.c_str());
	b_res1_mass      =  atof(res1_mass.c_str());
	b_res1_mass_err  =  atof(res1_mass_error.c_str());
	b_res1_width     =  atof(res1_width.c_str());
	b_res1_width_err =  atof(res1_width_error.c_str());
	b_res1_frac      =  atof(res1_frac.c_str());
	b_res1_frac_err  =  atof(res1_frac_error.c_str());
                         =  
	b_res2_amp        =  atof(res2_amp.c_str());
	b_res2_amp_err    =  atof(res2_amp_error.c_str());
	b_res2_phs        =  atof(res2_phs.c_str());
	b_res2_phs_err    =  atof(res2_phs_error.c_str());
	b_res2_mass      =  atof(res2_mass.c_str());
	b_res2_mass_err  =  atof(res2_mass_error.c_str());
	b_res2_width     =  atof(res2_width.c_str());
	b_res2_width_err =  atof(res2_width_error.c_str());
	b_res2_frac      =  atof(res2_frac.c_str());
	b_res2_frac_err  =  atof(res2_frac_error.c_str());
                        
        b_res3_amp        = atof(res3_amp.c_str());
	b_res3_amp_err    = atof(res3_amp_error.c_str());
	b_res3_phs        = atof(res3_phs.c_str());
	b_res3_phs_err    = atof(res3_phs_error.c_str());
	b_res3_mass      = atof(res3_mass.c_str());
	b_res3_mass_err  = atof(res3_mass_error.c_str());
	b_res3_width     = atof(res3_width.c_str());
	b_res3_width_err = atof(res3_width_error.c_str());
	b_res3_frac      = atof(res3_frac.c_str());
	b_res3_frac_err  = atof(res3_frac_error.c_str());
	

	b_status = atof(status.c_str());
	b_edm = atof(edm.c_str());
	b_chi2 = atof(chi2.c_str());
	b_fcn = atof(fcn.c_str());

	if(b_fcn != 0){
		tree.Fill();
	}

}
cout << "done" << endl;
TFile f("tree.root","recreate");
tree.Write();
f.Close();

}
