// basic file operations
#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

void makeRandomInput() {


  string resonances[] = {"a0+pi","f0980+pi","phi+K"};
  //string resonances[] = {"f0980+pi","phi+K"};
  //string resonances[] = {"f0980+pi","phi+K","f0X+pi"};
  const int NN = 3;
  bool ampphs = 1;

  int res_mass[] = {0,0,0};
  int res_fix[] = {0,0,1};
  double res_pars[NN][5];

 gRandom = new TRandom3(time(NULL));
  for(int kk=0;kk<NN;kk++){
	if(ampphs){
		res_pars[kk][0] = gRandom->Uniform(0,20);
		res_pars[kk][1] = gRandom->Uniform(-90,90);
	}else{
		res_pars[kk][0] = gRandom->Uniform(-10,10);
		res_pars[kk][1] = gRandom->Uniform(-10,10);
	}
	//res_pars[kk][2] = gRandom->Uniform(1.0,2.0);
	//res_pars[kk][3] = gRandom->Uniform(0.05,0.90);
 }
  // fix phi
 res_pars[2][0] = 1.0;
 res_pars[2][1] = 0.0;

  ofstream myfile;
  myfile.open ("Input.txt");
  myfile << "input_ntuple_name                                  Toy.root\n";
  myfile << "number_of_samples                                  1\n";
  myfile << "final_state                                        3\n";
  myfile << "is_gaussian                                        0\n";
  myfile << "Bkg_par1                                           3000\n";
  myfile << "Bkg_par2                                           2000\n";
  myfile << "real_and_imaginary                                 "<<!ampphs<<"\n";
  myfile << "UsesWeights                                        0\n";
  myfile << endl;

 for(int ii=0;ii<NN;ii++){
	if( ampphs ){ 
		myfile << resonances[ii] << "_amp                             " << res_pars[ii][0] << "\n";	
		myfile << resonances[ii] << "_amp_lower_limit                 0 \n";	
		myfile << resonances[ii] << "_amp_upper_limit                 100 \n";	
		myfile << resonances[ii] << "_amp_fix                         " << res_fix[ii] <<" \n";	
		myfile << resonances[ii] << "_phs                             " << res_pars[ii][1] << "\n";	
		myfile << resonances[ii] << "_phs_lower_limit                 -180 \n";	
		myfile << resonances[ii] << "_phs_upper_limit                  180 \n";	
		myfile << resonances[ii] << "_phs_fix                         " << res_fix[ii] <<" \n";	
	}
	else{ 
		myfile << resonances[ii] << "_re                             " << res_pars[ii][0] << "\n";	
		myfile << resonances[ii] << "_re_lower_limit                 -100 \n";	
		myfile << resonances[ii] << "_re_upper_limit                 100 \n";	
		myfile << resonances[ii] << "_re_fix                         " << res_fix[ii] <<" \n";	
		myfile << resonances[ii] << "_im                             " << res_pars[ii][1] << "\n";	
		myfile << resonances[ii] << "_im_lower_limit                 -100 \n";	
		myfile << resonances[ii] << "_im_upper_limit                 100 \n";	
		myfile << resonances[ii] << "_im_fix                         " << res_fix[ii] <<" \n";	
	}
	if(res_mass[ii]){
		myfile << resonances[ii] << "_mass                         " << res_pars[ii][2] << "\n";	
		myfile << resonances[ii] << "_mass_lower_limit              1.0\n";	
		myfile << resonances[ii] << "_mass_upper_limit              2.0 \n";	
		myfile << resonances[ii] << "_mass_fix                     " << res_fix[ii] <<" \n";	
		myfile << resonances[ii] << "_width                        " << res_pars[ii][3] << "\n";	
		myfile << resonances[ii] << "_width_lower_limit            0.05 \n";	
		myfile << resonances[ii] << "_width_upper_limit            0.90 \n";	
		myfile << resonances[ii] << "_width_fix                     " << res_fix[ii] <<" \n";	
	}
  }
 
  myfile << endl;
  myfile << "UseBackHisto                                      1\n";
  myfile << "HISTOGRAM_fraction                                1\n";
  myfile << "HISTOGRAM_fix                                     1\n";
  myfile << "Back_Ntuple_Name                                  ../Fit_Input/config_v2r4/bkg_histo_new_300bins.root\n";
  myfile << "Back_Histo_Name                                   bkgHist_acc\n"; 
  myfile << endl;

  myfile << "bkg_fraction                                           0.0954584\n";

  myfile << endl;
  myfile << "number_of_pwa_bins                                30\n";
  myfile << "pwa_txt_file                                      fix_PWA_COEFFS_30.txt\n";
  myfile << "output_pwa_txt_file                               foo_output_pwa.txt\n";
  myfile << "UseAcceptance                                     1\n";
  myfile << "Acceptance_Ntuple_Name                            ../Fit_Input/config_v2r4/effspline10_300.root\n";
  myfile << "Acceptance_Histo_Name                             eff_spline\n"; 


  myfile.close();
} 
