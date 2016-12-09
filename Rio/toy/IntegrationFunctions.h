#include <iostream>
#include "TF1.h"
#include "TH2F.h"
#include "TMath.h"
#include "Math/WrappedTF1.h"
#include "Math/GaussLegendreIntegrator.h"

using namespace std;

int number_of_points = 250;
double rel_tolerance = 0.0000001;

class OneDimIntFunctor
{
	public:
		OneDimIntFunctor ()
		{}
		void SetValues (string &integration_type, int &final_state, double &M, double &s13, vector<double> &KK_bin_limits,
				vector<TComplex> &pwa_coefs, vector<TComplex> &pwa_coefs_prime, int &index1, int &index2){

			_integration_type = integration_type;
			_final_state = final_state;
			_M = M;
			_s13 = s13;
			_KK_bin_limits = KK_bin_limits;
			_pwa_coefs = pwa_coefs;
			_pwa_coefs_prime = pwa_coefs_prime;
			_index1 = index1;
			_index2 = index2;

			// Defines masses of daughters

			if (_final_state == 0) {
				m1 = mK;
				m2 = mK;
				m3 = mpi;
			}
			else if (_final_state == 1) {
				m1 = mK;
				m2 = mpi;
				m3 = mpi;
			}
			else if (_final_state == 2) {
				m1 = mpi;
				m2 = mpi;
				m3 = mpi;
			}
			else if (_final_state == 3) {
				m1 = mK;
				m2 = mK;
				m3 = mK;
			}
			else{
				cout << "OneDimBkgIntFunctor - ERROR: Incorrect final state inserted";
				exit(-1);
			}

			// Calculates s23 from s12 and s13

			resonances.resize(AVAILABLE_RESONANCES);
			bkg_components.resize(AVAILABLE_BKG_COMPONENTS);

		}
		double operator() (double *s12, double *p) {

			UNUSED(p);

			s23 = _M*_M + m1*m1 + m2*m2 + m3*m3 - s12[0] - _s13;

			// Defines vectors with appropriate sizes

			// Includes in the model only the element of the current integral calculation.

			if (_integration_type.compare("Bkg") == 0 || _integration_type.compare("BinBkg") == 0){
				bkg_amp.resize(1);
				bkg_components[_index1] = 1;
			} else if (_integration_type.compare("SigRe") == 0){
				sig_amp.resize(2);
				resonances[_index1] = 1;
				resonances[_index2] = 1;
			} else if (_integration_type.compare("SigIm") == 0){
				sig_amp.resize(2);
				resonances[_index1] = 1;
				resonances[_index2] = 1;
			} else if (_integration_type.compare("BinSigRe") == 0){
				sig_amp.resize(2);
				resonances[_index1] = 1;
				resonances[_index2] = 1;
			} else if (_integration_type.compare("BinSigIm") == 0){
				sig_amp.resize(2);
				resonances[_index1] = 1;
				resonances[_index2] = 1;
			} else {
				cout << "OneDimBkgIntFunctor - ERROR: Invalid integration type: " << _integration_type << endl;
				exit(-1);
			}


			// Calculates the element value in the current phase space point

			Amplitudes(_final_state, _M, s12[0], _s13, s23, sig_amp, bkg_amp, resonances, bkg_components, _KK_bin_limits, _pwa_coefs, _pwa_coefs_prime);


			if (_integration_type.compare("Bkg") == 0 || _integration_type.compare("BinBkg") == 0){
				return bkg_amp[0].Re();
			} else if (_integration_type.compare("SigRe") == 0){
				if (_index1 < _index2) temp_sig = sig_amp[0]*sig_amp[1].Conjugate(sig_amp[1]);
				else if(_index1 > _index2) temp_sig = sig_amp[1]*sig_amp[0].Conjugate(sig_amp[0]);
				else temp_sig = sig_amp[0]*sig_amp[0].Conjugate(sig_amp[0]);
				return temp_sig.Re();
			} else if (_integration_type.compare("SigIm") == 0){
				if (_index1 < _index2) temp_sig = sig_amp[0]*sig_amp[1].Conjugate(sig_amp[1]);
				else if(_index1 > _index2) temp_sig = sig_amp[1]*sig_amp[0].Conjugate(sig_amp[0]);
				else temp_sig = sig_amp[0]*sig_amp[0].Conjugate(sig_amp[0]);
				return temp_sig.Im();
			} else if(_integration_type.compare("BinSigRe") == 0){
				if (_index1 < _index2) temp_sig = sig_amp[0]*sig_amp[1].Conjugate(sig_amp[1]);
				else if(_index1 > _index2) temp_sig = sig_amp[1]*sig_amp[0].Conjugate(sig_amp[0]);
				else temp_sig = sig_amp[0]*sig_amp[0].Conjugate(sig_amp[0]);
				return temp_sig.Re();
			} else if (_integration_type.compare("BinSigIm") == 0){
				if (_index1 < _index2) temp_sig = sig_amp[0]*sig_amp[1].Conjugate(sig_amp[1]);
				else if(_index1 > _index2) temp_sig = sig_amp[1]*sig_amp[0].Conjugate(sig_amp[0]);
				else temp_sig = sig_amp[0]*sig_amp[0].Conjugate(sig_amp[0]);
				return temp_sig.Im();
			} else{
				cout << "OneDimBkgIntFunctor - ERROR: Invalid integration type: " << _integration_type << endl;
				exit(-1);
			}
		}
	private:
		string _integration_type;
		int _final_state;
		double _M, _s13;
		vector<double> _KK_bin_limits; 
		vector<TComplex> _pwa_coefs, _pwa_coefs_prime;
		int _index1, _index2;
		double s23, m1, m2, m3;
		TComplex temp_sig;
		vector<int> resonances, bkg_components;
		vector<TComplex> sig_amp, bkg_amp;
};

class TwoDimIntFunctor
{
	public:
		TwoDimIntFunctor ()
		{}
		void SetValues (OneDimIntFunctor  &OneDimIntObject, int &npoints, double *points, double *weights, string &integration_type, int &final_state, double &M,
				vector<double> KK_bin_limits, vector<TComplex> pwa_coefs, vector<TComplex> pwa_coefs_prime,
				int &index1, int index2 = 0, double Bins12Min = 0, double Bins12Max = 0){
			_OneDimIntObject = OneDimIntObject;
			_npoints = npoints;
			_points = points;
			_weights = weights;
			_integration_type = integration_type;
			_final_state = final_state;
			_M = M;
			_KK_bin_limits = KK_bin_limits;
			_pwa_coefs = pwa_coefs;
			_pwa_coefs_prime = pwa_coefs_prime;
			_index1 = index1;
			_index2 = index2;
			_Bins12Min = Bins12Min;
			_Bins12Max = Bins12Max;

			// Defines masses of daughters

			if (_final_state == 0) {
				m1 = mK;
				m2 = mK;
				m3 = mpi;
			}
			else if (_final_state == 1) {
				m1 = mK;
				m2 = mpi;
				m3 = mpi;
			}
			else if (_final_state == 2) {
				m1 = mpi;
				m2 = mpi;
				m3 = mpi;
			}
			else if (_final_state == 3) {
				m1 = mK;
				m2 = mK;
				m3 = mK;
			}
			else{
				cout << "TwoDimBkgFunctor - ERROR: Incorrect final state inserted";
				exit(-1);
			}

			s = _M*_M;
			m3sq = m3*m3;
			m2sq = m2*m2;
			m1sq = m1*m1;

		}
		double operator() (double *s13 , double *p){

			UNUSED(p);

			// Calculates the 12 integration limits for the current s13 value
			c1 = _M*_M - s13[0] - m2*m2;
			c2 = s13[0] + m1*m1 - m3*m3;
			c3 = sqrt(lambda (s,s13[0],m2sq));
			c4 = sqrt(lambda (s13[0],m3sq,m1sq));
			s12min = m1*m1 + m2*m2 + (c1*c2 - c3*c4) / (2*s13[0]);
			s12max = m1*m1 + m2*m2 + (c1*c2 + c3*c4) / (2*s13[0]);
			if (_integration_type == "BinSigRe" || _integration_type == "BinSigIm" || _integration_type == "BinBkg"){
				if (s12min < _Bins12Min) s12_lower_limit = _Bins12Min;
				else s12_lower_limit = s12min;
				if (s12max > _Bins12Max) s12_upper_limit = _Bins12Max;
				else s12_upper_limit = s12max;
			}else{
				s12_lower_limit = s12min;
				s12_upper_limit = s12max;
			}

			// Defines integration function

			_OneDimIntObject.SetValues(_integration_type, _final_state, _M, s13[0], _KK_bin_limits, _pwa_coefs, _pwa_coefs_prime, _index1, _index2);
			TF1 OneDimIntFunction("OneDimIntFunction", _OneDimIntObject, s12_lower_limit, s12_upper_limit, 0);

			// Integrates in the s12 axis
			if (_integration_type == "BinSigRe" || _integration_type == "BinSigIm" || _integration_type == "BinBkg"){
				if (_Bins12Max > s12min && _Bins12Min < s12max){
					if (_final_state == 2 || _final_state == 3) {
						if (s13[0] < s12_upper_limit) {
							if (s13[0] > s12_lower_limit) {
								s12_lower_limit = s13[0];
							}
							Integral = OneDimIntFunction.IntegralFast(_npoints, _points, _weights, s12_lower_limit, s12_upper_limit);
						}else Integral = 0;
					} else Integral = OneDimIntFunction.IntegralFast(_npoints, _points, _weights, s12_lower_limit, s12_upper_limit);

				}
				else{
					Integral = 0;
				}
			}else Integral =  OneDimIntFunction.IntegralFast(_npoints, _points, _weights, s12min, s12max);
			return Integral;
		}
	private:
		OneDimIntFunctor _OneDimIntObject;
		int _npoints;
		double *_points, *_weights;
		string _integration_type;
		int _final_state;
		double _M;
		vector<double> _KK_bin_limits; 
		vector<TComplex> _pwa_coefs, _pwa_coefs_prime;
		int _index1, _index2;
		double _Bins12Min, _Bins12Max;
		double m1, m2, m3, s12max, s12min, s12_lower_limit, s12_upper_limit, c1, c2, c3, c4, Integral, s, m1sq, m2sq, m3sq;
};

class MassIntFunctor
{
	public:
		MassIntFunctor ()
		{}
		void SetValues (TwoDimIntFunctor  &TwoDimIntObject, OneDimIntFunctor  &OneDimIntObject, int &npoints, double *points, double *weights, string &integration_type,
				int &final_state, double &Mass_min, double &Mass_max, double &Bkg_par1, double &Bkg_par2,
				vector<double> &KK_bin_limits, vector<TComplex> &pwa_coefs, vector<TComplex> &pwa_coefs_prime, int &index1,
				int index2 = 0, double Bins13Min = 0, double Bins13Max = 0, double Bins12Min = 0, double Bins12Max = 0){
			_TwoDimIntObject = TwoDimIntObject;
			_OneDimIntObject = OneDimIntObject;
			_npoints = npoints;
			_points = points;
			_weights = weights;
			_integration_type = integration_type;
			_final_state = final_state;
			_Mass_min = Mass_min;
			_Mass_max = Mass_max;
			_Bkg_par1 = Bkg_par1;
			_Bkg_par2 = Bkg_par2;
			_KK_bin_limits = KK_bin_limits;
			_pwa_coefs = pwa_coefs;
			_pwa_coefs_prime = pwa_coefs_prime;
			_index1 = index1;
			_index2 = index2;
			_Bins13Min = Bins13Min;
			_Bins13Max = Bins13Max;
			_Bins12Min = Bins12Min;
			_Bins12Max = Bins12Max;
		}
		double operator() (double *M , double *p){
			double m1, m2, m3, s13max, s13min, s13_lower_limit, s13_upper_limit, Mass_pdf, Integral, Mother_Mass, Mother_Width;

			UNUSED(p);

			// Defines masses of daughters

			if (_final_state == 0) {
				m1 = mK;
				m2 = mK;
				m3 = mpi;
			}
			else if (_final_state == 1) {
				m1 = mK;
				m2 = mpi;
				m3 = mpi;
			}
			else if (_final_state == 2) {
				m1 = mpi;
				m2 = mpi;
				m3 = mpi;
			}
			else if (_final_state == 3) {
				m1 = mK;
				m2 = mK;
				m3 = mK;
			}
			else{
				cout << "TwoDimBkgFunctor - ERROR: Incorrect final state inserted";
				exit(-1);
			}

			Mother_Mass = D_Mass;
			Mother_Width = D_width;
            Mass_pdf = 0;

			if (_integration_type == "Bkg" || _integration_type == "BinBkg"){
				Mass_pdf = Bkg_Mass(M[0], _Bkg_par1, _Bkg_par2, _Mass_min, _Mass_max);
			}else if(_integration_type == "SigRe"|| _integration_type == "SigIm" || _integration_type == "BinSigRe"|| _integration_type == "BinSigIm"){
				Mass_pdf = Gaussian_Mass (M[0], Mother_Mass, Mother_Width, _Mass_min, _Mass_max, _npoints, _points, _weights);
			}

			// Calculates the 13 integration limits for the current s13 value

			s13min = (m1 + m3)*(m1 + m3);
			s13max = (M[0] - m2)*(M[0] - m2);

			if (_integration_type == "BinSigRe" || _integration_type == "BinSigIm" || _integration_type == "BinBkg"){
				if (s13min < _Bins13Min) s13_lower_limit = _Bins13Min;
				else s13_lower_limit = s13min;
				if (s13max > _Bins13Max) s13_upper_limit = _Bins13Max;
				else s13_upper_limit = s13max;
			}else{
				s13_lower_limit = s13min;
				s13_upper_limit = s13max;
			}


			_TwoDimIntObject.SetValues(_OneDimIntObject, _npoints, _points, _weights, _integration_type, _final_state, M[0], _KK_bin_limits,
					_pwa_coefs, _pwa_coefs_prime, _index1, _index2, _Bins12Min, _Bins12Max);
			TF1 TwoDimIntFunction("TwoDimIntFunction", _TwoDimIntObject, s13_lower_limit, s13_upper_limit, 0);

			// Integrates in the s12 axis
			if (_integration_type == "BinSigRe"|| _integration_type == "BinSigIm" || _integration_type == "BinBkg"){
				if (_Bins13Max > s13min && _Bins13Min < s13max) Integral = Mass_pdf*TwoDimIntFunction.IntegralFast(_npoints, _points, _weights, s13_lower_limit, s13_upper_limit);
				else{
					Integral = 0;
				}
			}else Integral = Mass_pdf*TwoDimIntFunction.IntegralFast(_npoints, _points, _weights, s13min, s13max);
			return Integral;
		}
	private:
		TwoDimIntFunctor _TwoDimIntObject;
		OneDimIntFunctor _OneDimIntObject;
		int _npoints;
		double *_points, *_weights;
		string _integration_type;
		int _final_state;
		double _Mass_min, _Mass_max, _Bkg_par1, _Bkg_par2;
		vector<double> _KK_bin_limits; 
		vector<TComplex> _pwa_coefs, _pwa_coefs_prime;
		int _index1, _index2;
		double _Bins13Min, _Bins13Max, _Bins12Min, _Bins12Max;
};
