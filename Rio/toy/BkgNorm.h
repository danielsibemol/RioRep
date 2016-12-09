#include <stdio.h>
#include <TMath.h>
#include <TROOT.h>
#include <TComplex.h>
#include <ctime>
#include "IntegrationFunctions.h"

using namespace std;

void BkgNorm(int final_state, bool is_gaussian, double Mass_min, double Mass_max, double Bkg_par1, double Bkg_par2, int number_of_bkg_components, vector<int> bkg_components, vector<double> &bkg_normalization){
    double M, s13min, s13max, m1, m2, m3;
    int bkg_index = 0;
    string SigRe = "SigRe", SigIm = "SigIm", Bkg = "Bkg";
    double *points = new double[number_of_points];
    double *weights = new double[number_of_points];
 	vector<double> KK_bin_limit;
	vector<TComplex> pwa_coefs, pwa_coefs_prime;
    TwoDimIntFunctor TwoDimBkgIntObject;
    OneDimIntFunctor OneDimBkgIntObject;
    MassIntFunctor MassBkgIntObject;
    
    clock_t begin = clock();
    
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
        cout << "BkgNorm - ERROR: Incorrect final state inserted";
        exit(-1);
    }
    
    M = D_Mass;
	
    s13min = (m1 + m3)*(m1 + m3);
	s13max = (M - m2)*(M - m2);
    
    bkg_normalization.resize(number_of_bkg_components);
	
    // Calculates integration points and weights
    
    TF1 Function;
    Function.CalcGaussLegendreSamplingPoints(number_of_points,points,weights,rel_tolerance);
 
    // Calculates the normalization of each component

    for (int j = 0; j < AVAILABLE_BKG_COMPONENTS; j++){
        if (bkg_components[j] == 1){
            if (!is_gaussian) {
            
                // Defines integration function

                TwoDimBkgIntObject.SetValues(OneDimBkgIntObject, number_of_points, points, weights, Bkg, final_state, M, KK_bin_limit, pwa_coefs, pwa_coefs_prime, j, 0);
                TF1 TwoDimBkgIntFunction("TwoDimBkgIntFunction", TwoDimBkgIntObject, s13min, s13max, 0);
                
                // Calculates the integral. The line below calculates the integral along the s12 axis for each s13 point and then integrates along the s13 axis.
                
                bkg_normalization.at(bkg_index) = TwoDimBkgIntFunction.IntegralFast(number_of_points, points, weights, s13min,s13max);
                bkg_index++;
             }else{
                // Defines integration function
                
                MassBkgIntObject.SetValues(TwoDimBkgIntObject, OneDimBkgIntObject, number_of_points, points, weights, Bkg, final_state, Mass_min, Mass_max, Bkg_par1, Bkg_par2,
                                                                       KK_bin_limit, pwa_coefs, pwa_coefs_prime, j, 0);
                TF1 MassBkgIntFunction("MassBkgIntFunction", MassBkgIntObject, Mass_min, Mass_min, 0);
                
                // Calculates the integral. The line below calculates the integral along the s12 axis for each s13 point and then integrates along the s13 axis.
                
                bkg_normalization.at(bkg_index) = MassBkgIntFunction.IntegralFast(number_of_points, points, weights, Mass_min,Mass_max);
                cout << "bkg_normalization[" << bkg_index << "] = " << bkg_normalization[bkg_index] << endl;
                bkg_index++;
            }
        }
    }
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    cout << "BkgNorm took " << elapsed_secs << "s" << endl;
}
