#include <stdio.h>
#include <TMath.h>
#include <TROOT.h>
#include <TComplex.h>
#include "Functions-pwa.h"
#include "Magalhaes.h"

using namespace std;

bool reading_ntuple = 0;

///// for the Log
void Amplitudes(int &final_state, 
		double &M, 
		double &s12, 
		double &s13, 
		double &s23, 
		vector<TComplex> &sig, 
		vector<TComplex> &bkg, 
		vector<int> &resonances, 
		vector<int> &bkg_components, 
		vector<double> &res_masses, 
		vector<double> &res_widths,
		vector< vector<double> > &res_extra_pars,
		vector<double> &KK_bin_limits, 
		vector<TComplex> &pwa_coefs, 
		vector<TComplex> &pwa_coefs_prime) 

{

	double s, m_0, w_0, m12, m13, m23, m1sq, m2sq, m3sq, gamma, f_theta, fR, fD, m1, m2, m3, aa, bb, ak, bk, aa3, bb3, dmKK, real, imaginary, Cd, Cm, real_O, imag_O, M_O, W_O, g_K, g_pi;
	int khi,klo,npt, size_resonances = 0, size_bkg = 0, spin, khi12, klo12, khi13, klo13;
	TComplex pole, Real_factor, BWbkg, temp_sig, temp_bkg, sig12, sig13, gamma_pi, gamma_K, gamma_f0;


	// Define masses and square masses

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
		cout << "Amplitudes - ERROR: Incorrect final state inserted";
		return;
	}

	m12 = sqrt(s12);
	m13 = sqrt(s13);
	m23 = sqrt(s23);

	m1sq = m1*m1;
	m2sq = m2*m2;
	m3sq = m3*m3;
	s = M*M;

	////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////
	//  K K PI AMPLITUDES
	////////////////////////////////////////////////////////////     
	////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////
	//   CALCULATING THE  D -> K*(892) K COMPLEX AMPLITUDE
	//////////////////////////////////////////////////////////////     
	if (resonances.at(0) == 1) {

		//if (reading_ntuple == 1) cout << "Begin - K*(892)" << endl;
		m_0 = res_masses[size_resonances];
		w_0 = res_widths[size_resonances];
		spin = 1;

		// Blatt-Weisskopf form factors

		fD = Form_Factor_Mother_Decay(spin, M, s13, m2sq, m_0);
		fR = Form_Factor_Resonance_Decay(spin, m_0, s13, m1sq, m3sq);

		//  Mass-dependent width

		gamma = Gamma(spin, m_0, w_0, m13, m1sq, m3sq);

		// Angular distribution

		f_theta = Angular_Distribution(spin, M, m3, m1, m2, s13, s12, s23);

		// Complex amplitude for D -> K*(892) K

		Real_factor(fD * fR * f_theta,0);

		// cout << "fD = " << fD << ", fR = " << fR << ", gamma = " << gamma << ", f_theta = " << f_theta << endl;

		temp_sig = Real_factor*BW(s13, m_0, gamma);
		sig.at(size_resonances) = temp_sig;
		size_resonances++;
	}	
	////////////////////////////////////////////////////////////
	//   CALCULATING THE  D -> K*(1430) K COMPLEX AMPLITUDE
	////////////////////////////////////////////////////////////     
	if (resonances.at(1) == 1) {

		//if (reading_ntuple == 1) cout << "Begin - K*(1430)" << endl;
		m_0 = res_masses[size_resonances];
		w_0 = res_widths[size_resonances];
		spin = 0;

		gamma = Gamma(spin, m_0, w_0, m13, m1sq, m3sq);

		// Complex amplitude for D -> K*(1430) K
		sig.at(size_resonances) = BW(s13, m_0, gamma);
		size_resonances++;
	}
	////////////////////////////////////////////////////////////
	//   CALCULATING THE  D -> PHI PI+ COMPLEX AMPLITUDE
	////////////////////////////////////////////////////////////     
	if (resonances.at(2) == 1) {

		m_0 = res_masses[size_resonances];
		w_0 = res_widths[size_resonances];
		spin = 1;

		// Blatt-Weisskopf form factors

		fD = Form_Factor_Mother_Decay(spin, M, s12, m3sq, m_0);
		fR = Form_Factor_Resonance_Decay(spin, m_0, s12, m1sq, m2sq);

		// Mass-dependent width

		gamma = Gamma(spin, m_0, w_0, m12, m1sq, m2sq);

		// Angular distribution

		f_theta = Angular_Distribution(spin, M, m1, m2, m3, s12, s23, s13);

		// Complex amplitude for D -> PHI PI+

		Real_factor(fD * fR * f_theta,0);

		sig.at(size_resonances) = Real_factor*BW (s12, m_0, gamma);
		size_resonances++;
	}
	////////////////////////////////////////////////////////////
	//   CALCULATING THE  D -> a0(1450) pi COMPLEX AMPLITUDE
	////////////////////////////////////////////////////////////     
	if (resonances.at(3) == 1) {

		m_0 = res_masses[size_resonances];
		w_0 = res_widths[size_resonances];
		spin = 0;

		// Mass-dependent width

		gamma = Gamma(spin, m_0, w_0, m12, m1sq, m2sq);

		sig.at(size_resonances) = BW(s12, m_0, gamma);
		size_resonances++;

	}
	////////////////////////////////////////////////////////////
	//   CALCULATING THE  D -> kappa K COMPLEX AMPLITUDE
	////////////////////////////////////////////////////////////     
	if (resonances.at(4) == 1) {
		real = res_masses[size_resonances];
		imaginary = res_widths[size_resonances];		

		sig.at(size_resonances) = Pole(s13, real, imaginary);
		size_resonances++;
	}
	////////////////////////////////////////////////////////////
	//   D -> pi(K) pi(K) pi(K) (NR)
	////////////////////////////////////////////////////////////     
	if (resonances.at(5) == 1) {
		sig.at(size_resonances) = sig.at(size_resonances).One();
		size_resonances++;
	}
	////////////////////////////////////////////////////////////
	//   CALCULATING THE D -> (KK)_S wave pi AMPLITUDE
	//////////////////////////////////////////////////////////// 

	if (resonances.at(6) == 1) {

		npt = KK_bin_limits.size();

		find_bin_KKP(m12,KK_bin_limits,npt,khi,klo);
		dmKK = KK_bin_limits[khi] - KK_bin_limits[klo];
		aa = (KK_bin_limits[khi] - m12)/dmKK;
		bb = (m12 - KK_bin_limits[klo])/dmKK;
		aa3 = aa * aa * aa;
		bb3 = bb * bb * bb;
		ak = aa * pwa_coefs[klo].Re() + bb * pwa_coefs[khi].Re() + ((aa3 - aa)*pwa_coefs_prime[klo].Re() + (bb3 - bb) * pwa_coefs_prime[khi].Re()) * (dmKK*dmKK)/6.0;
		bk = aa * pwa_coefs[klo].Im() + bb * pwa_coefs[khi].Im() + ((aa3 - aa)*pwa_coefs_prime[klo].Im() + (bb3 - bb) * pwa_coefs_prime[khi].Im()) * (dmKK*dmKK)/6.0;	

		temp_sig(ak,bk);

		sig.at(size_resonances) = temp_sig;
		size_resonances++;
	}				
	////////////////////////////////////////////////////////////
	//   D -> f0(980)(->KK) pi
	////////////////////////////////////////////////////////////     
	if (resonances.at(7) == 1) {
		m_0 = res_masses[size_resonances];
		w_0 = res_widths[size_resonances];
		g_K = res_extra_pars[size_resonances][0];
		g_pi = res_extra_pars[size_resonances][1];
		spin = 0;


		// Mass-dependent width
		sig12 = flatte(s12, m_0,g_K,g_pi);

		sig.at(size_resonances) = sig12;

		size_resonances++;			
	}
	////////////////////////////////////////////////////////////
	//   D -> f0(X)(->KK) pi
	////////////////////////////////////////////////////////////     
	if (resonances.at(8) == 1) {
		m_0 = res_masses[size_resonances];
		w_0 = res_widths[size_resonances];
		spin = 0;

		// Mass-dependent width

		gamma = Gamma(spin, m_0, w_0, m12, m1sq, m2sq);

		sig12 = BW(s12, m_0, gamma);

		sig.at(size_resonances) = sig12; 
		size_resonances++;			
	}	
	////////////////////////////////////////////////////////////
	//   D -> K2(1430) K
	////////////////////////////////////////////////////////////     
	if (resonances.at(9) == 1) {
		m_0 = res_masses[size_resonances];
		w_0 = res_widths[size_resonances];
		spin = 2;
		// Blatt-Weisskopf form factors

		fD = K2_Form_Factor_Mother_Decay(spin, M, s13, m2sq, m_0);
		fR = Form_Factor_Resonance_Decay(spin, m_0, s13, m1sq, m3sq);

		// Mass-dependent width

		gamma = Gamma(spin, m_0, w_0, m13, m1sq, m3sq);

		// Angular distribution

		f_theta = Angular_Distribution(spin, M, m3, m1, m2, s13, s12, s23);

		// Complex amplitude for D -> rho(770) PI+

		Real_factor(fD*fR*f_theta, 0);

		sig13 = Real_factor*BW(s13, m_0, gamma);

		sig.at(size_resonances) = sig13;
		size_resonances++;			
	}
	////////////////////////////////////////////////////////////
	//   CALCULATING THE  D -> PHI(1680) PI+ COMPLEX AMPLITUDE
	////////////////////////////////////////////////////////////     
	if (resonances.at(10) == 1) {

		m_0 = res_masses[size_resonances];
		w_0 = res_widths[size_resonances];
		spin = 1;

		// Blatt-Weisskopf form factors

		fD = Form_Factor_Mother_Decay(spin, M, s12, m3sq, m_0);
		fR = Form_Factor_Resonance_Decay(spin, m_0, s12, m1sq, m2sq);

		// Mass-dependent width

		gamma = Gamma(spin, m_0, w_0, m12, m1sq, m2sq);

		// Angular distribution

		f_theta = Angular_Distribution(spin, M, m1, m2, m3, s12, s23, s13);

		// Complex amplitude for D -> PHI PI+

		Real_factor(fD * fR * f_theta,0);

		sig.at(size_resonances) = Real_factor*BW (s12, m_0, gamma);
		size_resonances++;
	}
	////////////////////////////////////////////////////////////
	//   CALCULATING THE  D -> kappa (BW) K COMPLEX AMPLITUDE
	////////////////////////////////////////////////////////////     
	if (resonances.at(11) == 1) {
		m_0 = res_masses[size_resonances];
		w_0 = res_widths[size_resonances];
		spin = 0;

		// Complex amplitude for D -> K*(1430) K
		sig.at(size_resonances) = BW(s13, m_0, w_0);
		size_resonances++;
	}
	////////////////////////////////////////////////////////////
	//   CALCULATING THE  D -> f0(1500) pi COMPLEX AMPLITUDE
	////////////////////////////////////////////////////////////     
	if (resonances.at(12) == 1) {
		m_0 = res_masses[size_resonances];
		w_0 = res_widths[size_resonances];
		spin = 0;

		gamma = Gamma(spin, m_0, w_0, m12, m1sq, m2sq);

		sig.at(size_resonances) = BW(s12, m_0, gamma);
		size_resonances++;
	}
	////////////////////////////////////////////////////////////
	//   D -> K1(1410) K
	////////////////////////////////////////////////////////////     
	if (resonances.at(13) == 1) {
		m_0 = res_masses[size_resonances];
		w_0 = res_widths[size_resonances];
		spin = 1;
		// Blatt-Weisskopf form factors

		fD = K1_Form_Factor_Mother_Decay(spin, M, s13, m2sq, m_0);
		fR = Form_Factor_Resonance_Decay(spin, m_0, s13, m1sq, m3sq);

		// Mass-dependent width

		gamma = Gamma(spin, m_0, w_0, m13, m1sq, m3sq);

		// Angular distribution

		f_theta = Angular_Distribution(spin, M, m3, m1, m2, s13, s12, s23);

		// Complex amplitude for D -> rho(770) PI+

		Real_factor(fD*fR*f_theta, 0);

		sig13 = Real_factor*BW(s13, m_0, gamma);

		sig.at(size_resonances) =  sig13;
		size_resonances++;			
	}
	////////////////////////////////////////////////////////////
	//   CALCULATING THE  D -> a2(1320) PI+ COMPLEX AMPLITUDE
	////////////////////////////////////////////////////////////     
	if (resonances.at(14) == 1) {

		m_0 = res_masses[size_resonances];
		w_0 = res_widths[size_resonances];
		spin = 2;

		// Blatt-Weisskopf form factors

		fD = Form_Factor_Mother_Decay(spin, M, s12, m3sq, m_0);
		fR = Form_Factor_Resonance_Decay(spin, m_0, s12, m1sq, m2sq);

		// Mass-dependent width

		gamma = Gamma(spin, m_0, w_0, m12, m1sq, m2sq);

		// Angular distribution

		f_theta = Angular_Distribution(spin, M, m1, m2, m3, s12, s23, s13);

		// Complex amplitude for D -> PHI PI+

		Real_factor(fD * fR * f_theta,0);

		sig.at(size_resonances) = Real_factor*BW (s12, m_0, gamma);
		size_resonances++;
	}
	////////////////////////////////////////////////////////////
	//   CALCULATING THE  D -> f2(1270) PI+ COMPLEX AMPLITUDE
	////////////////////////////////////////////////////////////     
	if (resonances.at(15) == 1) {

		m_0 = res_masses[size_resonances];
		w_0 = res_widths[size_resonances];
		spin = 2;

		// Blatt-Weisskopf form factors

		fD = Form_Factor_Mother_Decay(spin, M, s12, m3sq, m_0);
		fR = Form_Factor_Resonance_Decay(spin, m_0, s12, m1sq, m2sq);

		// Mass-dependent width

		gamma = Gamma(spin, m_0, w_0, m12, m1sq, m2sq);

		// Angular distribution

		f_theta = Angular_Distribution(spin, M, m1, m2, m3, s12, s23, s13);

		// Complex amplitude for D -> PHI PI+

		Real_factor(fD * fR * f_theta,0);

		sig.at(size_resonances) = Real_factor*BW (s12, m_0, gamma);
		size_resonances++;
	}
	////////////////////////////////////////////////////////////
	//   CALCULATING THE  D -> K*(1680) K COMPLEX AMPLITUDE
	////////////////////////////////////////////////////////////     
	if (resonances.at(16) == 1) {

		m_0 = res_masses[size_resonances];
		w_0 = res_widths[size_resonances];
		spin = 0;

		gamma = Gamma(spin, m_0, w_0, m13, m1sq, m3sq);

		// Complex amplitude for D -> K*(1430) K
		sig.at(size_resonances) = BW(s13, m_0, gamma);
		size_resonances++;
	}
	////////////////////////////////////////////////////////////
	//   CALCULATING THE  D -> f0(1710) pi COMPLEX AMPLITUDE
	////////////////////////////////////////////////////////////     
	if (resonances.at(17) == 1) {
		m_0 = res_masses[size_resonances];
		w_0 = res_widths[size_resonances];
		spin = 0;

		gamma = Gamma(spin, m_0, w_0, m12, m1sq, m2sq);

		sig.at(size_resonances) = BW(s12, m_0, gamma);
		size_resonances++;
	}
	////////////////////////////////////////////////////////////
	//   CALCULATING THE  D -> f2p(1525) PI+ COMPLEX AMPLITUDE
	////////////////////////////////////////////////////////////     
	if (resonances.at(18) == 1) {

		m_0 = res_masses[size_resonances];
		w_0 = res_widths[size_resonances];
		spin = 2;

		// Blatt-Weisskopf form factors

		fD = Form_Factor_Mother_Decay(spin, M, s12, m3sq, m_0);
		fR = Form_Factor_Resonance_Decay(spin, m_0, s12, m1sq, m2sq);

		// Mass-dependent width

		gamma = Gamma(spin, m_0, w_0, m12, m1sq, m2sq);

		// Angular distribution

		f_theta = Angular_Distribution(spin, M, m1, m2, m3, s12, s23, s13);

		// Complex amplitude for D -> PHI PI+

		Real_factor(fD * fR * f_theta,0);

		sig.at(size_resonances) = Real_factor*BW (s12, m_0, gamma);
		size_resonances++;
	}
	////////////////////////////////////////////////////////////
	//   CALCULATING THE  Magalhaes Kpi S-wave COMPLEX AMPLITUDE
	////////////////////////////////////////////////////////////     
	if (resonances.at(19) == 1) {
		m_0 = res_masses[size_resonances];
		w_0 = res_widths[size_resonances];
		Cd = res_extra_pars[size_resonances][0];
		Cm = res_extra_pars[size_resonances][1];
		spin = 0;

		sig.at(size_resonances) = Magalhaes(s13, m_0, Cd, Cm);
		size_resonances++;
	}

	////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////
	//  PI PI PI AMPLITUDES
	////////////////////////////////////////////////////////////     
	////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////
	//   D -> rho omega
	////////////////////////////////////////////////////////////     

	if (resonances.at(20) == 1) {

		m_0 = res_masses[size_resonances];
		w_0 = res_widths[size_resonances];

		real_O = res_extra_pars[size_resonances][0];
		imag_O = res_extra_pars[size_resonances][1];
		spin = 1;

		M_O = m0_omega; // = 0.78265;  //0.939;
		W_O = w0_omega; // = 0.00849;

		fD = Form_Factor_Mother_Decay(spin, M, s12, m3sq, m_0);
		fR = Form_Factor_Resonance_Decay(spin, m_0, s12, m1sq, m2sq);

		// Mass-dependent width

		gamma = Gamma(spin, m_0, w_0, m12, m1sq, m2sq);

		// Angular distribution

		f_theta = Angular_Distribution(spin, M, m1, m2, m3, s12, s23, s13);

		// Complex amplitude for D -> rho(770) PI+

		Real_factor(fD*fR*f_theta, 0);

		sig12 = Real_factor*GS(s12, m_0, w_0, gamma, m1sq, m2sq )*( TComplex( 1.0, 0.0) + TComplex(real_O, imag_O)*BW(s12, M_O, W_O) );

		// same for s13

		// Blatt-Weisskopf form factors

		fD = Form_Factor_Mother_Decay(spin, M, s13, m2sq, m_0);
		fR = Form_Factor_Resonance_Decay(spin, m_0, s13, m1sq, m3sq);

		// Mass-dependent width

		gamma = Gamma(spin, m_0, w_0, m13, m1sq, m3sq);

		// Angular distribution

		f_theta = Angular_Distribution(spin, M, m1, m3, m2, s13, s23, s12);

		// Complex amplitude for D -> rho(770) PI+

		Real_factor(fD*fR*f_theta, 0);

		sig13 = Real_factor*GS(s13,  m_0, w_0, gamma, m1sq, m3sq )*( TComplex( 1.0, 0.0) + TComplex(real_O, imag_O)*BW(s13, M_O, W_O) );

		sig.at(size_resonances) = sig12 + sig13;
		size_resonances++;
	}

	////////////////////////////////////////////////////////////
	//   D -> f0(980) pi
	////////////////////////////////////////////////////////////     

	if (resonances.at(21) == 1) {
		m_0 = res_masses[size_resonances];
		w_0 = res_widths[size_resonances];
		g_K = res_extra_pars[size_resonances][0];
		g_pi = res_extra_pars[size_resonances][1];

		sig12 = flatte(s12, m_0, g_K, g_pi);
		sig13 = flatte(s13, m_0, g_K, g_pi);

		sig.at(size_resonances) = sig12 + sig13;
		size_resonances++;						
	}	

	////////////////////////////////////////////////////////////
	//   D -> f2(1270) pi
	////////////////////////////////////////////////////////////     
	if (resonances.at(22) == 1) {
		m_0 = res_masses[size_resonances];
		w_0 = res_widths[size_resonances];
		spin = 2;

		// Blatt-Weisskopf form factors

		fD = Form_Factor_Mother_Decay(spin, M, s12, m3sq, m_0);
		fR = Form_Factor_Resonance_Decay(spin, m_0, s12, m1sq, m2sq);

		// Mass-dependent width

		gamma = Gamma(spin, m_0, w_0, m12, m1sq, m2sq);

		// Angular distribution

		f_theta = Angular_Distribution(spin, M, m1, m2, m3, s12, s23, s13);

		// Complex amplitude for D -> rho(770) PI+

		Real_factor(fD*fR*f_theta, 0);

		sig12 = Real_factor*BW(s12, m_0, gamma);

		// same for s13

		// Blatt-Weisskopf form factors

		fD = Form_Factor_Mother_Decay(spin, M, s13, m2sq, m_0);
		fR = Form_Factor_Resonance_Decay(spin, m_0, s13, m1sq, m3sq);

		// Mass-dependent width

		gamma = Gamma(spin, m_0, w_0, m13, m1sq, m3sq);

		// Angular distribution

		f_theta = Angular_Distribution(spin, M, m1, m3, m2, s13, s23, s12);

		// Complex amplitude for D -> rho(770) PI+

		Real_factor(fD*fR*f_theta, 0);

		sig13 = Real_factor*BW(s13, m_0, gamma);

		sig.at(size_resonances) = sig12 + sig13;
		size_resonances++;			
	}
	////////////////////////////////////////////////////////////
	//   D -> rho(1450) pi
	////////////////////////////////////////////////////////////     
	if (resonances.at(23) == 1) {
		m_0 = res_masses[size_resonances];
		w_0 = res_widths[size_resonances];
		spin = 1;

		// Blatt-Weisskopf form factors

		fD = Form_Factor_Mother_Decay(spin, M, s12, m3sq, m_0);
		fR = Form_Factor_Resonance_Decay(spin, m_0, s12, m1sq, m2sq);

		// Mass-dependent width

		gamma = Gamma(spin, m_0, w_0, m12, m1sq, m2sq);

		// Angular distribution

		f_theta = Angular_Distribution(spin, M, m1, m2, m3, s12, s23, s13);

		// Complex amplitude for D -> rho(770) PI+

		Real_factor(fD*fR*f_theta, 0);

		sig12 = Real_factor*GS(s12, m_0, w_0, gamma, m1sq, m2sq );

		// same for s13

		// Blatt-Weisskopf form factors

		fD = Form_Factor_Mother_Decay(spin, M, s13, m2sq, m_0);
		fR = Form_Factor_Resonance_Decay(spin, m_0, s13, m1sq, m3sq);

		// Mass-dependent width

		gamma = Gamma(spin, m_0, w_0, m13, m1sq, m3sq);

		// Angular distribution

		f_theta = Angular_Distribution(spin, M, m1, m3, m2, s13, s23, s12);

		// Complex amplitude for D -> rho(770) PI+

		Real_factor(fD*fR*f_theta, 0);

		sig13 = Real_factor*GS(s13, m_0, w_0, gamma, m1sq, m3sq );

		sig.at(size_resonances) = sig12 + sig13;
		size_resonances++;
	}
	////////////////////////////////////////////////////////////
	//   D -> f0(X) pi
	////////////////////////////////////////////////////////////     
	if (resonances.at(24) == 1) {
		m_0 = res_masses[size_resonances];
		w_0 = res_widths[size_resonances];
		spin = 0;


		// Mass-dependent width

		gamma = Gamma(spin, m_0, w_0, m12, m1sq, m2sq);

		sig12 = BW(s12, m_0, gamma);

		// same for s13

		gamma = Gamma(spin, m_0, w_0, m13, m1sq, m3sq);

		sig13 = BW(s13, m_0, gamma);

		sig.at(size_resonances) = sig12 + sig13; 
		size_resonances++;			
	}	
	////////////////////////////////////////////////////////////
	//   D -> sigma pi
	////////////////////////////////////////////////////////////     
	if (resonances.at(25) == 1) {
		///usando BW0.560,    0.500
		m_0 = res_masses[size_resonances];
		w_0 = res_widths[size_resonances]; 
		spin = 0;

		sig12 = sigma(s12);
		sig13 = sigma(s13); 
		sig.at(size_resonances) = sig12 + sig13;
		size_resonances++;
	}	

	////////////////////////////////////////////////////////////
	//   CALCULATING THE D -> (KK)_S wave pi AMPLITUDE
	//////////////////////////////////////////////////////////// 

	if (resonances.at(26) == 1) {

		npt = KK_bin_limits.size();

		find_bin( m12, m13, KK_bin_limits, npt, khi12, klo12, khi13, klo13 );

		dmKK = KK_bin_limits[khi12] - KK_bin_limits[klo12];
		aa = ( KK_bin_limits[khi12] - m12 )/dmKK;
		bb = ( m12 - KK_bin_limits[klo12])/dmKK;
		aa3 = aa * aa * aa;
		bb3 = bb * bb * bb;

		ak = aa * pwa_coefs[klo12].Re() + bb * pwa_coefs[khi12].Re() + ((aa3 - aa)*pwa_coefs_prime[klo12].Re() + (bb3 - bb) * pwa_coefs_prime[khi12].Re()) * (dmKK*dmKK)/6.0;
		bk = aa * pwa_coefs[klo12].Im() + bb * pwa_coefs[khi12].Im() + ((aa3 - aa)*pwa_coefs_prime[klo12].Im() + (bb3 - bb) * pwa_coefs_prime[khi12].Im()) * (dmKK*dmKK)/6.0;	

		sig12 = temp_sig( ak, bk );

		dmKK = KK_bin_limits[khi13] - KK_bin_limits[klo13];
		aa = ( KK_bin_limits[khi13] - m13 )/dmKK;
		bb = ( m13 - KK_bin_limits[klo13] )/dmKK;
		aa3 = aa * aa * aa;
		bb3 = bb * bb * bb;

		ak = aa * pwa_coefs[klo13].Re() + bb * pwa_coefs[khi13].Re() + ((aa3 - aa)*pwa_coefs_prime[klo13].Re() + (bb3 - bb) * pwa_coefs_prime[khi13].Re()) * (dmKK*dmKK)/6.0;
		bk = aa * pwa_coefs[klo13].Im() + bb * pwa_coefs[khi13].Im() + ((aa3 - aa)*pwa_coefs_prime[klo13].Im() + (bb3 - bb) * pwa_coefs_prime[khi13].Im()) * (dmKK*dmKK)/6.0;	

		sig13 = temp_sig( ak, bk );

		sig.at(size_resonances) = sig12 + sig13;;
		size_resonances++;
	}				

	////////////////////////////////////////////////////////////
	//   CALCULATING THE  D -> f0(1500) pi COMPLEX AMPLITUDE
	////////////////////////////////////////////////////////////     
	if (resonances.at(27) == 1) {
		m_0 = res_masses[size_resonances];
		w_0 = res_widths[size_resonances];
		spin = 0;

		gamma = Gamma(spin, m_0, w_0, m12, m1sq, m2sq);

		sig12 = BW(s12, m_0, gamma);

		// same for s13

		gamma = Gamma(spin, m_0, w_0, m13, m1sq, m3sq);

		sig13 = BW(s13, m_0, gamma);
		sig.at(size_resonances) = sig12+sig13;
		size_resonances++;
	}

	////////////////////////////////////////////////////////////
	//   SIGMA BW
	////////////////////////////////////////////////////////////  
	if (resonances.at(28) == 1) {
		m_0 = res_masses[size_resonances];
		w_0 = res_widths[size_resonances];
		spin = 0;


		// Mass-dependent width

		gamma = Gamma(spin, m_0, w_0, m12, m1sq, m2sq);

		sig12 = BW(s12, m_0, gamma);

		// same for s13

		gamma = Gamma(spin, m_0, w_0, m13, m1sq, m3sq);

		sig13 = BW(s13, m_0, gamma);

		sig.at(size_resonances) = sig12 + sig13; 
		size_resonances++;			
	}	

	////////////////////////////////////////////////////////////
	//   F0(980) BW
	////////////////////////////////////////////////////////////  
	if (resonances.at(29) == 1) {
		m_0 = res_masses[size_resonances];
		w_0 = res_widths[size_resonances];
		spin = 0;


		// Mass-dependent width

		gamma = Gamma(spin, m_0, w_0, m12, m1sq, m2sq);

		sig12 = BW(s12, m_0, gamma);

		// same for s13

		gamma = Gamma(spin, m_0, w_0, m13, m1sq, m3sq);

		sig13 = BW(s13, m_0, gamma);

		sig.at(size_resonances) = sig12 + sig13; 
		size_resonances++;			
	}	



	////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////
	//  K K K AMPLITUDES
	////////////////////////////////////////////////////////////     
	////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////
	//   CALCULATING THE  D -> a0(1450) K COMPLEX AMPLITUDE
	////////////////////////////////////////////////////////////     
	if (resonances.at(30) == 1) {

		m_0 = GLOBAL_TOY? resonant_channel_mass[30]:res_masses[size_resonances];
		w_0 = GLOBAL_TOY? resonant_channel_width[30]:res_widths[size_resonances];
		spin = 0;

		//cout << "a0 pi - " << "m0 = " << m_0 << ", w_0 = " << w_0 << endl;
		// Mass-dependent width

		gamma = Gamma(spin, m_0, w_0, m12, m1sq, m2sq);
		sig12 = BW(s12, m_0, gamma);

		gamma = Gamma(spin, m_0, w_0, m13, m1sq, m3sq);
		sig13 = BW(s13, m_0, gamma);

		sig.at(size_resonances) = sig12 + sig13;
		size_resonances++;

	}

	////////////////////////////////////////////////////////////
	//   D -> f0(980) pi
	////////////////////////////////////////////////////////////     
	if (resonances.at(31) == 1) {

		m_0 = GLOBAL_TOY? resonant_channel_mass[31]:res_masses[size_resonances];
		w_0 = GLOBAL_TOY? resonant_channel_width[31]:res_widths[size_resonances];

		g_K = GLOBAL_TOY? gK:res_extra_pars[size_resonances][0];
		g_pi = GLOBAL_TOY? gpi:res_extra_pars[size_resonances][1];
		spin = 0;

		sig12 = flatte(s12, m_0, g_K, g_pi);
		sig13 = flatte(s13, m_0, g_K, g_pi);

		sig.at(size_resonances) = sig12 + sig13;
		size_resonances++;
	}	
	////////////////////////////////////////////////////////////
	//   D -> f2(1270) pi
	////////////////////////////////////////////////////////////     
	if (resonances.at(32) == 1) {

		m_0 = GLOBAL_TOY? resonant_channel_mass[32]:res_masses[size_resonances];
		w_0 = GLOBAL_TOY? resonant_channel_width[32]:res_widths[size_resonances];
		spin = 2;

		//cout << "f2(1270) pi - " << "m0 = " << m_0 << ", w_0 = " << w_0 << endl;

		// Blatt-Weisskopf form factors

		fD = Form_Factor_Mother_Decay(spin, M, s12, m3sq, m_0);
		fR = Form_Factor_Resonance_Decay(spin, m_0, s12, m1sq, m2sq);

		// Mass-dependent width

		gamma = Gamma(spin, m_0, w_0, m12, m1sq, m2sq);

		// Angular distribution

		f_theta = Angular_Distribution(spin, M, m1, m2, m3, s12, s23, s13);

		// Complex amplitude for D -> rho(770) PI+

		Real_factor(fD*fR*f_theta, 0);

		sig12 = Real_factor*BW(s12, m_0, gamma);

		// same for s13

		// Blatt-Weisskopf form factors

		fD = Form_Factor_Mother_Decay(spin, M, s13, m2sq, m_0);
		fR = Form_Factor_Resonance_Decay(spin, m_0, s13, m1sq, m3sq);

		// Mass-dependent width

		gamma = Gamma(spin, m_0, w_0, m13, m1sq, m3sq);

		// Angular distribution

		f_theta = Angular_Distribution(spin, M, m1, m3, m2, s13, s23, s12);

		// Complex amplitude for D -> rho(770) PI+

		Real_factor(fD*fR*f_theta, 0);

		sig13 = Real_factor*BW(s13, m_0, gamma);

		sig.at(size_resonances) = sig12 + sig13;
		size_resonances++;			
	}

	////////////////////////////////////////////////////////////
	//   D -> f0(X) K
	////////////////////////////////////////////////////////////     
	if (resonances.at(33) == 1) {

		m_0 = GLOBAL_TOY? resonant_channel_mass[33]:res_masses[size_resonances];
		w_0 = GLOBAL_TOY? resonant_channel_width[33]:res_widths[size_resonances];
		spin = 0;

		// Mass-dependent width

		gamma = Gamma(spin, m_0, w_0, m12, m1sq, m2sq);

		sig12 = BW(s12, m_0, gamma);

		// same for s13

		gamma = Gamma(spin, m_0, w_0, m13, m1sq, m3sq);

		sig13 = BW(s13, m_0, gamma);

		sig.at(size_resonances) = sig12 + sig13; 
		size_resonances++;			
	}	
	////////////////////////////////////////////////////////////
	//   CALCULATING THE  D -> PHI K+ COMPLEX AMPLITUDE
	////////////////////////////////////////////////////////////
	if (resonances.at(34) == 1) {

		m_0 = GLOBAL_TOY? resonant_channel_mass[34]:res_masses[size_resonances];
		w_0 = GLOBAL_TOY? resonant_channel_width[34]:res_widths[size_resonances];
		spin = 1;

		// Blatt-Weisskopf form factors

		fD = Form_Factor_Mother_Decay(spin, M, s12, m3sq, m_0);
		fR = Form_Factor_Resonance_Decay(spin, m_0, s12, m1sq, m2sq);

		// Mass-dependent width

		gamma = Gamma(spin, m_0, w_0, m12, m1sq, m2sq);

		// Angular distribution

		f_theta = Angular_Distribution(spin, M, m1, m2, m3, s12, s23, s13);

		// Complex amplitude for D -> PHI PI+

		Real_factor( fD * fR * f_theta, 0 );

		sig12 = Real_factor*BW (s12, m_0, gamma);

		// same for s13

		// Blatt-Weisskopf form factors

		fD = Form_Factor_Mother_Decay(spin, M, s13, m2sq, m_0);
		fR = Form_Factor_Resonance_Decay(spin, m_0, s13, m1sq, m3sq);

		// Mass-dependent width

		gamma = Gamma(spin, m_0, w_0, m13, m1sq, m3sq);

		// Angular distribution

		f_theta = Angular_Distribution(spin, M, m1, m3, m2, s13, s23, s12);

		// Complex amplitude for D -> PHI PI+

		Real_factor(fD * fR * f_theta,0);

		sig13 = Real_factor*BW (s13, m_0, gamma);

		sig.at(size_resonances) = sig12 + sig13;
		size_resonances++;
	}

	////////////////////////////////////////////////////////////
	//   D -> f0(1500) pi
	////////////////////////////////////////////////////////////     
	if (resonances.at(35) == 1) {

		m_0 = GLOBAL_TOY? resonant_channel_mass[35]:res_masses[size_resonances];
		w_0 = GLOBAL_TOY? resonant_channel_width[35]:res_widths[size_resonances];
		spin = 0;


		// Mass-dependent width

		gamma = Gamma(spin, m_0, w_0, m12, m1sq, m2sq);

		sig12 = BW(s12, m_0, gamma);

		// same for s13

		gamma = Gamma(spin, m_0, w_0, m13, m1sq, m3sq);

		sig13 = BW(s13, m_0, gamma);

		sig.at(size_resonances) = sig12 + sig13; 
		size_resonances++;			
	}	

	////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////
	//  END OF SIGNAL AMPLITUDES
	////////////////////////////////////////////////////////////     
	////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////
	//   RANDOM PHI
	////////////////////////////////////////////////////////////     
	if (bkg_components.at(0) == 1) {

		m_0 = m0_phi;
		w_0 = w0_phi;

		BWbkg = BW(s12, m_0, w_0);
		temp_bkg(BWbkg.Rho(),0);
		bkg.at(size_bkg) = temp_bkg;
		size_bkg++;
	}
	////////////////////////////////////////////////////////////
	//   RANDOM K*
	////////////////////////////////////////////////////////////     
	if (bkg_components.at(1) == 1) {

		m_0 = m0_K892;
		w_0 = w0_K892;

		BWbkg = BW(s13, m_0, w_0);
		temp_bkg(BWbkg.Rho(),0);
		bkg.at(size_bkg) = temp_bkg;
		size_bkg++;
	}
	////////////////////////////////////////////////////////////
	//   COMBINATORIAL
	////////////////////////////////////////////////////////////     
	if (bkg_components.at(2) == 1) {

		temp_bkg(1,0);
		bkg.at(size_bkg) = temp_bkg;
		size_bkg++;
	}
	////////////////////////////////////////////////////////////
	//   HISTOGRAM
	////////////////////////////////////////////////////////////
	if (bkg_components.at(3) == 1) {

		temp_bkg(backgroundF(s12, s13),0);
		bkg.at(size_bkg) = temp_bkg;
		//bkg.at(size_bkg) = backgroundF(s12, s13);

		size_bkg++;
	}


}
//// for the Generator 

void Amplitudes( int &final_state, 
		double &M, 
		double &s12, 
		double &s13, 
		double &s23, 
		vector<TComplex> &sig, 
		vector<TComplex> &bkg, 
		vector<int> &resonances, 
		vector<int> &bkg_components, 
		vector<double> &KK_bin_limits, 
		vector<TComplex> &pwa_coefs, 
		vector<TComplex> &pwa_coefs_prime){ 

	vector<double> res_masses; 
	vector<double> res_widths;
	vector< vector<double> > res_extra_pars;

        Amplitudes(final_state,M,s12,s13,s23,sig,bkg,resonances,bkg_components,res_masses,res_widths,res_extra_pars,KK_bin_limits,pwa_coefs,pwa_coefs_prime);
}
