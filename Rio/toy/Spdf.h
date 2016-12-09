#include <stdio.h>
#include <TMath.h>
#include <TROOT.h>
#include <math.h>

using namespace std;

double Spdf(vector<TComplex> &sig, vector<TComplex> &coefs, double &s12, double &s13, double &s23){
	double direct, interference, spdf;
	TComplex temp_Ai, temp_Aj;
	direct = 0;
	interference = 0;

	// This function calculates the Signal PDF
	
	for (unsigned int i = 0; i < sig.size(); i++){
		for (unsigned int j = i; j < sig.size(); j++){
			if (i == j){
				direct = direct + coefs.at(i).Rho()*coefs.at(i).Rho()*sig.at(i).Rho()*sig.at(i).Rho();
			}
			else{
				temp_Ai = coefs.at(i)*sig.at(i);
				temp_Aj = coefs.at(j)*sig.at(j);
				temp_Aj = temp_Aj.Conjugate(temp_Aj);
				interference = interference + 2*(temp_Ai*temp_Aj).Re();
			}
		}
	}
	spdf = direct + interference;
	double AccInt=UseAcceptance?acceptance(s12, s13):1;
	return spdf*AccInt;
}
