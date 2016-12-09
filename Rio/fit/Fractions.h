#include <stdio.h>
#include <TMath.h>
#include <TROOT.h>
#include <TComplex.h>
#include <math.h>
#include <fstream>
#include <sstream>

void zz(){
	for(int i=0;i<64;++i)cout<<"#";cout<<endl;
}

double rnd2(double valPar) {
	return (long long) (valPar * 1000 + (valPar > 0 ? 0.5 : -0.5)) / 1000.0;
}

void Fractions(TFitter* minimizer, bool real_and_imaginary, vector<vector<TComplex> > &nn, int n_res, vector<int> fix_parameter_index, vector<string> channel, string rootfile, string oss2, double texpars[]){
	double pars[8*n_res],m_err[8*n_res][8*n_res],f[n_res],dnf_dpar[8*n_res],df_dpar[n_res][8*n_res],err_f[n_res], number_of_fix_parameters;
	double nf, amp, amp_err, phs, phs_err, totalFracSum=0;
	double interFrac[n_res][n_res];
	vector<double> dummy_bkg_normalization, dummy_binned_spdf, dummy_binned_bpdf;
	vector< vector< vector<TComplex> > > dummy_bin_sig_normalization;
	vector< vector<double> > dummy_bin_bkg_normalization;

	// Read the fit parameters, define an array for the covariance matrix
	for(int i=0; i<8*n_res; i++){
		pars[i]=minimizer->GetParameter(i);
		for(int j=0; j<8*n_res; j++) {
			m_err[i][j]=0;
		}
		dnf_dpar[i]=0;
	}

	number_of_fix_parameters = fix_parameter_index.size();

	int c=0, m=0;
	for(int i=0;i<8*n_res;i++){
		for(int k=0;k<number_of_fix_parameters;k++){
			if(i==fix_parameter_index[k]){
				c++;
				goto  a;
			}
		}
		for(int j=0;j<8*n_res;j++){
			for(int k=0;k<number_of_fix_parameters;k++){
				if(j==fix_parameter_index[k]){
					m++;
					goto b;
				}
			}
			(i>j)?(m_err[i][j]=minimizer->GetCovarianceMatrixElement(i-c,j-m)):(m_err[i][j]=minimizer->GetCovarianceMatrixElement(j-m,i-c));
b:continue;
		}m=0;
a: continue;
	}

    // calculate normalization N_f
	nf= SigNorm(nn, n_res,  pars, real_and_imaginary);

  // calculate diagonal fit fractions (DFF)
	for(int i=0; i<n_res; i++) {

		f[i]=(pars[8*i]*pars[8*i]+pars[8*i+1]*pars[8*i+1]*(real_and_imaginary?1:0))*nn[i][i].Rho()/nf;


//		FALTA IMPLEMENTAR ERROS DE PARAMETROS EXTRA
		for(int j=0; j<n_res; j++) { 
			if(real_and_imaginary){

				dnf_dpar[8*i]   += 2*(   pars[8*j]*nn[i][j].Re() ); //   + pars[8*j+1]*nn[i][j].Im() );
				dnf_dpar[8*i+1] += 2*( pars[8*j+1]*nn[i][j].Re() ); // - pars[8*j]*nn[i][j].Im() );
				//dnf_dpar[8*i+2] += 2*( pars[8*j+2]*nn[i][j].Re(); // - pars[8*j]*nn[i][j].Im() );
				//dnf_dpar[8*i+3] += 2*( pars[8*j+3]*nn[i][j].Re(); // - pars[8*j]*nn[i][j].Im() );
				//dnf_dpar[8*i+4] += 2*( pars[8*j+4]*nn[i][j].Re();
                //dnf_dpar[8*i+5] += 2*( pars[8*j+5]*nn[i][j].Re();
                //dnf_dpar[8*i+6] += 2*( pars[8*j+6]*nn[i][j].Re();
                //dnf_dpar[8*i+7] += 2*( pars[8*j+7]*nn[i][j].Re();

			}         

			else{

				dnf_dpar[8*i]   += 2*pars[8*j]*((nn[i][j]*TComplex( cos(pars[8*i+1]-pars[8*j+1]),sin(pars[8*i+1]-pars[8*j+1]))).Re());
				dnf_dpar[8*i+1] -= 2*pars[8*i]*pars[8*j]*((nn[i][j]*TComplex(cos(pars[8*i+1]-pars[8*j+1] ),sin(pars[8*i+1]-pars[8*j+1] ))).Im());
				//dnf_dpar[8*i+2] -= 2*pars[8*i]*pars[8*j]*((nn[i][j]*TComplex(cos(pars[8*i+1]-pars[8*j+1] ),sin(pars[8*i+1]-pars[8*j+1] ))).Im());
				//dnf_dpar[8*i+3] -= 2*pars[8*i]*pars[8*j]*((nn[i][j]*TComplex(cos(pars[8*i+1]-pars[8*j+1] ),sin(pars[8*i+1]-pars[8*j+1] ))).Im());
				//dnf_dpar[8*i+4] -= 2*pars[8*i]*pars[8*j]*((nn[i][j]*TComplex(cos(pars[8*i+1]-pars[8*j+1] ),sin(pars[8*i+1]-pars[8*j+1] ))).Im());
				//dnf_dpar[8*i+5] -= 2*pars[8*i]*pars[8*j]*((nn[i][j]*TComplex(cos(pars[8*i+1]-pars[8*j+1] ),sin(pars[8*i+1]-pars[8*j+1] ))).Im());
				//dnf_dpar[8*i+6] -= 2*pars[8*i]*pars[8*j]*((nn[i][j]*TComplex(cos(pars[8*i+1]-pars[8*j+1] ),sin(pars[8*i+1]-pars[8*j+1] ))).Im());
				//dnf_dpar[8*i+7] -= 2*pars[8*i]*pars[8*j]*((nn[i][j]*TComplex(cos(pars[8*i+1]-pars[8*j+1] ),sin(pars[8*i+1]-pars[8*j+1] ))).Im());

			}          
		}
	}
    
 // calculate interf fracs (just the upper diagonal)
    for (int i=0; i<n_res; i++) {
        for (int j=i+1; j<n_res; j++) {
            if(!real_and_imaginary){
                // mag phase
                //TComplex c(pars[8*i]*pars[8*j], pars[8*i+1]-pars[8*j+1],kTRUE);
                TComplex c(pars[8*i]*pars[8*j], pars[8*i+1]-pars[8*j+1],kTRUE);
                TComplex d = c*nn[i][j];
                interFrac[i][j] = 2*d.Re()/nf;
            }else{
			//re-im (a+ib)*(c-id) = (ac+bd) + i(bc-ab)
                TComplex c((pars[8*i]*pars[8*j]+pars[8*i+1]*pars[8*j+1]), (pars[8*i+1]*pars[8*j]-pars[8*i]*pars[8*j+1]),kFALSE);
                TComplex d = c*nn[i][j];
                interFrac[i][j] = 2*d.Re()/nf;
            }
            totalFracSum+=interFrac[i][j];
        }}
    
   // complete the interf table with zeros and the DFF
   for(int i=0; i<n_res; i++){
   	for(int j=0; j<n_res; j++){
   		if(i>j)interFrac[i][j]=0;
   		if(i==j) interFrac[i][j]=f[i];
   		
   	}		
   }


	for(int i=0; i<n_res; i++){ 

		for(int j=0; j<8*n_res; j++){

			if(real_and_imaginary)  df_dpar[i][j]=2*(pars[8*i]*((8*i==j)?1:0)+pars[8*i+1]*(((8*i+1)==j)?1:0))*nn[i][i].Rho()/nf-(pars[8*i]*pars[8*i]+pars[8*i+1]*pars[8*i+1])*dnf_dpar[j]*nn[i][i].Rho()/(nf*nf);  
			else    df_dpar[i][j]=2*pars[8*i]*((8*i==j)?1:0)*nn[i][i].Rho()/nf-(pars[8*i]*pars[8*i]*dnf_dpar[j]*nn[i][i].Rho())/(nf*nf);
		
		}
	}


    /////////////////////// PRINT ///////////////////////////////
	ofstream vals;
	string frag =rootfile.substr(0,rootfile.length()-4);
	string outtxt=frag+"txt";
	vals.open(outtxt.c_str());     
	(real_and_imaginary)?cout<<" Calculating fractions for real-imag form"<<endl:cout<<" Calculating fractions for amp-phase form  (UseAcceptance =) "<<UseAcceptance<<endl;   
	zz();
	cout<<left <<setw(15)<<"Channel"<< setw(15) <<" "<<" Fraction (%)"<< setw(15)<<" "<<" Fraction error (%)"<< endl;
	zz();

	ofstream myfile;  //file to output the fit result - ofstream : Stream class to write on files
	myfile.open (oss2.c_str(),ios::app); //ios:app :  append to the end of the file

	for(int i=0; i<n_res; i++){
		err_f[i]=0;
		for(int j=0; j<8*n_res; j++){
			for(int k=0;k<8*n_res; k++){
				err_f[i]+=df_dpar[i][j]*df_dpar[i][k]*m_err[j][k]; 
			}
		}
		err_f[i]=sqrt(err_f[i]);
		cout << left << setw(15) << channel[i] << setw(15) <<std::setprecision(5) <<" "<< f[i]*100 << setw(15)<<" "<<err_f[i]*100<< endl;//setw(15) <<std::setprecision(5) << interFrac[i][i+1]*100<<   endl;
		vals << setw(15) <<std::setprecision(5) << f[i] <<setw(15) <<err_f[i]<<setw(15) << interFrac[i][i+1]<<    endl;
        	myfile << "frac_" << channel[i] << " " << f[i] << " " << err_f[i] << " " << endl ;
	}
	myfile << endl;
	myfile.close();
	zz();
	vals.close();


	for(int i=0; i<n_res; i++){

		for(int j=0; j<n_res; j++){

			if(i>j)interFrac[i][j]=0;
			else 	if(i==j) interFrac[i][j]=f[i];
		}	
	}

// Print as interferencias
for(int i=0; i<n_res; i++){
	cout<< left << setw(15) << channel[i]<<" ";
	for(int j=0; j<n_res; j++){
		cout<<setw(15) <<std::setprecision(5) << interFrac[i][j]*100;	
		} 
	cout<<endl;
}

ofstream myfiletex;  //file to output the fit result - ofstream : Stream class to write on files
string outex= frag+"tex";

myfiletex.open(outex.c_str());

string legend = "";
double ffsum = 0; 
//////////////////////////////////////////////////
myfiletex <<"\\begin{tabular}{@{}p{3.0cm}p{3cm}p{3cm}p{2.5cm}@{}}"<<endl
<<"\\bf Resonance  & \\bf ";
if(real_and_imaginary)myfiletex <<" Real"<<" & \\bf Imag.";
else myfiletex<<" Amplitude"<<" & \\bf Phase ";
myfiletex<<" & \\bf Fraction (\\%)"<<" \\\\"<<endl
<<"\\multicolumn{4}{c}{\\vspace{0.02cm}}"<<" \\"<<"\\"<<endl
<<"\\hline"<<endl;

for(int i=0; i<n_res; i++){

	myfiletex << "$" <<channel[i] << "$" <<" & "<< rnd2(texpars[8*i]) << " $\\pm$ " << rnd2(texpars[8*i+1])<<" & ";
	if(real_and_imaginary)myfiletex<< rnd2(texpars[8*i+2])<< " $\\pm$ "  <<rnd2(texpars[8*i+3])<<" & ";
	else myfiletex<< rnd2(texpars[8*i+2]*180/PI)<< " $\\pm$ "  <<rnd2(texpars[8*i+3]*180/PI)<<" & ";
	myfiletex<<rnd2(f[i]*100) << " $\\pm$ "<< rnd2(err_f[i]*100) <<" \\"<<"\\"<<endl;
	ffsum += f[i]*100;

}

myfiletex <<"\\hline"<<endl
<< "Total Fit Fraction & & & " << rnd2(ffsum) << " \\\\ \\hline" << endl
<<"\\multicolumn{4}{c}{\\vspace{0.06cm}}"<<" \\"<<"\\"<<endl
<<"\\end{tabular}"<<endl;

myfiletex.close();
cout<<"Created tex table "<<outex<<endl;



//Sum of fraction
double f_sum =0;
double f_sum_err =0;
for(int i=0; i<n_res; i++){
    f_sum += f[i]; 	
    f_sum_err += err_f[i]*err_f[i]; 
}
f_sum_err = sqrt(f_sum_err);
cout << " Frac Sum = " << f_sum << " +- " << f_sum_err << endl;


ofstream myfiletex2;
string frag2 ="Results_PDF_fit/"+rootfile.substr(0,rootfile.length()-4);
string outex2= frag+"_amp_phs.tex";
myfiletex2.open(outex2.c_str());

myfiletex2 <<"\\begin{tabular}{@{}p{3.0cm}p{2cm}p{2cm}p{2.5cm}@{}}"<<endl
//   <<"\\bf Resonance  & \\bf Amplitude & \\bf Phase & \\bf Fraction"<<" \\"<<"\\"<<"\\toprule"<<endl
<<"\\bf Resonance  & \\bf ";
if(real_and_imaginary)myfiletex2 <<" Magnitude"<<" & \\bf Phase";
else myfiletex2<<" Amplitude"<<" & \\bf Phase ";
myfiletex2<<" & \\bf Fraction (\%)"<<" \\\\"<<endl
<<"\\multicolumn{4}{c}{\\vspace{0.02cm}}"<<" \\"<<"\\"<<endl
<<"\\hline"<<endl;

for(int i=0; i<n_res; i++){
	amp = sqrt(pow(texpars[8*i],2) + pow(texpars[8*i+2],2));
	amp_err = sqrt(pow(2*texpars[8*i]*texpars[8*i+1],2) + pow(2*texpars[8*i+2]*texpars[8*i+3],2));
	phs = atan2(texpars[8*i],texpars[8*i+2]);
	phs_err = sqrt(pow((1/texpars[8*i+2])*(1/(1+pow(texpars[8*i]/texpars[8*i+2],2)))*texpars[8*i+1],2) 
			+ pow((texpars[8*i]/pow(texpars[8*i+2],2))*(1/(1+pow(texpars[8*i]/texpars[8*i+2],2)))*texpars[8*i+3],2));
	//cout << "amp = " << amp << ", amp_err = " << amp_err << ", phs = " << phs << ", phs_err = " << phs_err << endl; 
	myfiletex2 << "$" <<channel[i] << "$" <<" & "<< rnd2(amp) << " $\\pm$ " << rnd2(amp_err)<<" & ";
	if(real_and_imaginary)myfiletex2<< rnd2(phs*180/PI)<< " $\\pm$ "  <<rnd2(phs_err*180/PI)<<" & ";
	else myfiletex2<< rnd2(texpars[8*i+2]*180/PI)<< " $\\pm$ "  <<rnd2(texpars[8*i+3]*180/PI)<<" & ";
	myfiletex2<<rnd2(f[i]*100) << " $\\pm$ "<< rnd2(err_f[i]*100) <<" \\"<<"\\"<<endl;

}

myfiletex2 <<"\\hline"<<endl
<< "Total Fit Fraction & & & " << rnd2(ffsum) << " \\\\ \\hline" << endl
<<"\\multicolumn{4}{c}{\\vspace{0.06cm}}"<<" \\"<<"\\"<<endl
<<"\\end{tabular}"<<endl;

myfiletex2.close();
cout<<"Created tex 2 table "<<outex2<<endl;
cout<<"totalFracSum = "<<totalFracSum<<endl;
}

