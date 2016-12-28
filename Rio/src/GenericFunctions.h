#include <stdio.h>
#include <iostream>
#include <math.h>
#include <TComplex.h>
#include <TRandom1.h>
#include <TF1.h>
#include "Constants.h"
#include <TH2F.h>
#include <TH2Poly.h>
#include <TTree.h>
#include <TFile.h>
#include <fstream>
#include <sstream>
#include <iosfwd>

using namespace std;

// This file defines a few functions and constant values.

//This function provides the Gaussian shape to the signal mass distribution
////////////////////////////////////////////////////////////////////////
double Gaussian_Mass (double &M, double &Peak_Mass, double &width, double &Mass_min, double &Mass_max, int &npoints, double * points, double * weights){

	double Mass_pdf;
	TF1 *gauss = new TF1("gauss","gaus", Mass_min, Mass_max);
	gauss->SetParameters(1/(sqrt(2*PI)*width), Peak_Mass, width);
	Mass_pdf = gauss->Eval(M)/gauss->IntegralFast(npoints, points, weights, Mass_min,Mass_max);
	delete gauss;

	return Mass_pdf;

}

//This function provides the Linear shape to the background mass distribution
////////////////////////////////////////////////////////////////////////
double Bkg_Mass (double &M, double &Bkg_par1, double &Bkg_par2, double &Mass_min, double &Mass_max){

	double Bkg;
	Bkg = (Bkg_par1 + (Bkg_par2-Bkg_par1)*(M - Mass_min)/(Mass_max - Mass_min))/(Bkg_par1 + (Bkg_par2-Bkg_par1)*((Mass_max - Mass_min)/2 - Mass_min)/(Mass_max - Mass_min));

	return Bkg;

}
//This function gets a mass value within the Gaussian mass distribution for signal

double Get_Gaussian_Mass(TRandom1 * Random, double &Peak_Mass, double &width, double &Mass_min, double &Mass_max, int &npoints, double * points, double * weights){
	double M, Random_number1, Random_number2, Mass_window, Gaussian;
	bool Ok;
	Ok = 0;
	while (!Ok) {
		Random_number1 = Random->Rndm();
		Random_number2 = Random->Rndm();
		Mass_window = Mass_max - Mass_min;
		M = Mass_min + Random_number1*Mass_window;
		Gaussian = Gaussian_Mass(M,Peak_Mass,width, Mass_min, Mass_max, npoints, points, weights)/Gaussian_Mass(Peak_Mass,Peak_Mass,width, Mass_min, Mass_max,
				npoints, points, weights);
		if (Random_number2 < Gaussian) Ok = 1; 
	}
	return M;
}

//This function gets a mass value within the linear shape mass distribution for background

double Get_Bkg_Mass(TRandom1 * Random, double &Bkg_par1, double &Bkg_par2, double &Mass_min, double &Mass_max){
	double M, Random_number1, Random_number2, Mass_window, Bkg;
	bool Ok;
	Ok = 0;
	while (!Ok) {
		Random_number1 = Random->Rndm();
		Random_number2 = Random->Rndm();
		Mass_window = Mass_max - Mass_min;
		M = Mass_min + Random_number1*Mass_window;
		if (Bkg_par1 > Bkg_par2) Bkg = Bkg_Mass(M,Bkg_par1,Bkg_par2, Mass_min, Mass_max)/Bkg_Mass(Mass_min,Bkg_par1,Bkg_par2, Mass_min, Mass_max);
		else Bkg = Bkg_Mass(M,Bkg_par1,Bkg_par2, Mass_min, Mass_max)/Bkg_Mass(Mass_max,Bkg_par1,Bkg_par2, Mass_min, Mass_max);
		if (Random_number2 < Bkg) Ok = 1; 
	}
	return M;
}


// lambda function - used in generator and amplitudes
////////////////////////////////////////////////////////////////////////
double  lambda (double &x, double &y, double &z){

	double l;
	l = (x - y - z)*(x - y - z) - 4*y*z;
	if (l<0){
		cout << "lambda - Error - lambda < 0" << endl;
		exit(-1);
	}

	return l;

}

// lambda function - used in K2 resonance
////////////////////////////////////////////////////////////////////////
double  K2lambda (double &x, double &y, double &z){

	double l = -(x - y - z)*(x - y - z) + 4*y*z;
	if (l<0){ 
		cout << "K2lambda - Error - lambda > 0" << endl;
		exit(-1);
	}

	return -(x - y - z)*(x - y - z) + 4*y*z;

}

////////////////////////////////////////////////////////////////////////
double  BinAreasAndIntlambda (double &x, double &y, double &z){

	return (x - y - z)*(x - y - z) - 4*y*z;

}


// Blatt-Weisskopf Form Factors - used in amplitudes
////////////////////////////////////////////////////////////////////////
double Form_Factor_Mother_Decay(int &spin, double &M, double &sab, double &mcsq, double &mR){

	double s = M*M, mRsq = mR*mR;
	double fD, fD0, pstr, pstr0, q2;

	if (spin == 0) return 1;
	else if (spin == 1) {
		pstr0 = sqrt(lambda(s,mRsq,mcsq))/(2*M);
		q2 = rD2*pstr0*pstr0;
		fD0 = sqrt(1 + q2);

		pstr = sqrt(lambda(s,sab,mcsq))/(2*M);
		q2 = rD2*pstr*pstr;
		fD = fD0/sqrt(1 + q2);
		return fD;
	}
	else if(spin == 2){
		pstr0 = sqrt(lambda(s,mRsq,mcsq))/(2*M);
		q2 = rD2*pstr0*pstr0;
		fD0 = sqrt(9 + 3*q2 + q2*q2);

		pstr = sqrt(lambda(s,sab,mcsq))/(2*M);
		q2 = rD2*pstr*pstr;
		fD = fD0/sqrt(9 + 3*q2 + q2*q2);
		return fD;
	}else{
		cout << "Form_Factor_Mother_Decay - ERROR - invalid spin" << endl;
		return -1;
	}

}

////////////////////////////////////////////////////////////////////////
double K2_Form_Factor_Mother_Decay(int &spin, double &M, double &sab, double &mcsq, double &mR){

	double s = M*M, mRsq = mR*mR, pstr, q2, fD, fD0;
	TComplex pstr0;
	UNUSED(spin);
	//pstr0 = sqrt(lambda(s,mRsq,mcsq))/(2*mR); // Lembrar que Rio usa M (não mR)
	pstr0(0,sqrt(K2lambda(s,mRsq,mcsq))/(2*M)); 
	q2 = rD2*pstr0.Rho2();
	fD0 = sqrt(9. + 3.*q2 + q2*q2);

	//pstr = sqrt(lambda(s,sab,mcsq))/(2*sqrt(sab)); // Lembrar que Rio usa M (não mab) (de acordo com eq. 46.20b , p. 509 chinesse pdg)
	pstr = sqrt(lambda(s,sab,mcsq))/(2*M); 
	q2 = rD2*pstr*pstr;
	fD = fD0/sqrt(9. + 3.*q2 + q2*q2);

	return fD;

}

////////////////////////////////////////////////////////////////////////
double K1_Form_Factor_Mother_Decay(int &spin, double &M, double &sab, double &mcsq, double &mR){

	double s = M*M, mRsq = mR*mR, pstr, q2, fD, fD0;
	UNUSED(spin);
	TComplex pstr0;
	//pstr0 = sqrt(lambda(s,mRsq,mcsq))/(2*mR); // Lembrar que Rio usa M (não mR)
	pstr0(0,sqrt(K2lambda(s,mRsq,mcsq))/(2*M)); 
	q2 = rD2*pstr0.Rho2();
	fD0 = sqrt(1 + q2);

	//pstr = sqrt(lambda(s,sab,mcsq))/(2*sqrt(sab)); // Lembrar que Rio usa M (não mab) (de acordo com eq. 46.20b , p. 509 chinesse pdg)
	pstr = sqrt(lambda(s,sab,mcsq))/(2*M); 
	q2 = rD2*pstr*pstr;
	fD = fD0/sqrt(1 + q2);

	return fD;

}

////////////////////////////////////////////////////////////////////////
double Form_Factor_Resonance_Decay(int &spin, double &mR, double &sab, double &masq, double &mbsq){

	double mRsq = mR*mR;
	double fR, fR0, pstr, pstr0, q2;

	if (spin == 0) return 1;
	else if (spin == 1) {

		pstr0 = sqrt(lambda(mRsq,masq,mbsq))/(2*mR);
		q2 = rR2*pstr0*pstr0;
		fR0 = sqrt(1 + q2);

		pstr = sqrt(lambda(sab,masq,mbsq))/(2*sqrt(sab));
		q2 = rR2*pstr*pstr;
		fR = fR0/sqrt(1 + q2);

		return fR;

	}
	else if(spin == 2){

		pstr0 = sqrt((mRsq - masq - mbsq)*(mRsq - masq - mbsq) - 4*masq*mbsq)/(2*mR);
		q2 = rR2*pstr0*pstr0;
		fR0 = sqrt(9 + 3*q2 + q2*q2);

		//pstr = sqrt(lambda(sab,masq,mbsq))/(2*sqrt(sab));
		pstr = sqrt((sab - masq + mbsq)*(sab - masq + mbsq) - 4*sab*mbsq)/(2*sqrt(sab));
		q2 = rR2*pstr*pstr;
		fR = fR0/sqrt(9 + 3*q2 + q2*q2);

		return fR;

	}else{
		cout << "Form_Factor_Resonance_Decay - ERROR - invalid spin" << endl;

		return -1;

	}
}

// Angular distribution
////////////////////////////////////////////////////////////////////////
double Angular_Distribution(int &l, double &M, double &ma, double &mb, double &mc, double &sab, double &sbc, double &sac){

	double spin, spin1, spin2, A, B, C, D;

	spin1 = sbc - sac + (M*M-mc*mc)*(ma*ma-mb*mb)/sab;
	A = sab-2*M*M-2*mc*mc;
	B = ((M*M-mc*mc)*(M*M-mc*mc))/sab;
	C = sab-2*ma*ma-2*mb*mb;
	D = (ma*ma-mb*mb)*(ma*ma-mb*mb)/sab;
	spin2 = spin1*spin1-((A+B)*(C+D))/3.;
	if (l==0) spin=1.0;
	else if (l==1) spin = -0.5*spin1;
	else if (l==2) spin = 3./8*spin2;
	else {
		cout << "Angular_Distribution - ERROR - invalid spin" << endl;
		spin = 0;
		exit(-1);
	}

	return spin;

}

// Mass-dependent width
////////////////////////////////////////////////////////////////////////
double Gamma(int &spin, double &mR, double &width, double &mab, double &masq, double &mbsq){

	double pstr, pstr0,fR, mRsq = mR*mR,sab = mab*mab;

	pstr0 = sqrt(lambda(mRsq,masq,mbsq))/(2*mR);
	pstr = sqrt(lambda(sab,masq,mbsq))/(2*mab);
	if (spin == 0) return width*(pstr/pstr0)*(mR/mab);
	else if (spin == 1){
		fR = Form_Factor_Resonance_Decay(spin, mR, sab, masq, mbsq);
		return width*pow((pstr/pstr0),3)*(mR/mab)*fR*fR;
	}else if (spin == 2){
		fR = Form_Factor_Resonance_Decay(spin, mR, sab, masq, mbsq);
		return width*pow((pstr/pstr0),5)*(mR/mab)*fR*fR;
	}
	else {
		cout << "Gamma - ERROR - invalid spin" << endl;

		return -1;

	}

}

// Mass dependent width for kappa resonance
////////////////////////////////////////////////////////////////////////
double Gammakappa(double &sab, double &mR, double &masq, double &mbsq){

	double pstr, rho, f;
	pstr = sqrt(lambda(sab,masq,mbsq))/(2*sqrt(sab));
	rho = 2*pstr/sqrt(sab);
	f = b2*sab + b1; 

	return rho*((sab - sA)/(mR*mR - sA))*f*exp(-(sab-mR*mR)/A);

}

// Pole function - used for D-> kappa + K resonant channel
////////////////////////////////////////////////////////////////////////
TComplex Pole(double &sab, double Re_pole, double Im_pole){

	TComplex csab(sab,0);
	TComplex pole(Re_pole,Im_pole);

	return csab.One()/(pole*pole - csab);
}

// Breit-Wigner - used in amplitudes
////////////////////////////////////////////////////////////////////////
TComplex BW(double &sij, double &mR, double &gamma){

	TComplex csij(sij,0), cmR(mR,0), cgamma(gamma,0), unity;
	return unity.One()/(cmR*cmR - csij - unity.I()*cmR*cgamma);

}


// Breit-Wigner - used in amplitudes for kappa, with Bugg (2009)
////////////////////////////////////////////////////////////////////////
TComplex BWBugg(double &sab, double &s){

	double  p = sqrt((s-mK*mK-mpi*mpi)*(s-mK*mK-mpi*mpi)-4*mK*mK*mpi*mpi )/(2*sqrt(sab));
	double rho = 2.*p/sqrt(sab);
	TComplex rho2(0,0);
	double F = exp(-alpha*p*p);
	double F2 = 0;

	return 1./(1.-A*(s-sth)-TComplex(0.,1.)*(g2Kpi*(s-sA)*F*rho) - TComplex(0.,1.)*(g2Ketap*(s-sA2)*F2*rho2));

}

// Flatté - used in the f0 resonance amplitude
////////////////////////////////////////////////////////////////////////
TComplex Flattef0(double &sab, double &mR){

	double re, im;
	TComplex gamma_pi, gamma_K, gamma_f0, csab(sab,0), cmR(mR,0), unity;

	gamma_pi(sqrt(0.25*sab - mpi*mpi), 0.);

	im = 0.;
	re = 0.;
	if (sab > 4.*mK0sq){
		im = 0;
		re = sqrt(0.25*sab - mKsq) + sqrt(0.25*sab - mK0sq);
	}
	else if (sab < 4.*mKsq){
		im = sqrt(mKsq - 0.25*sab) + sqrt(mK0sq - 0.25*sab);
		re = 0;
	}
	else if (sab < 4.*mK0sq && sab > 4.*mKsq){
		re = sqrt (0.25*sab - mKsq);
		im = sqrt (mK0sq - 0.25*sab);
	}

	gamma_K(re,im);

	gamma_f0 = gpi*gamma_pi + 0.5*gK*gamma_K;

	return unity.One()/(cmR*cmR - csab - unity.I()*cmR*gamma_f0);

}

// Flatté - used in the k*1430 resonance amplitude
////////////////////////////////////////////////////////////////////////
TComplex FlatteKstar(double &sab, double &mR){

	TComplex rho_Kpi, rho_Ketap, gamma, csab(sab,0), cmR(mR,0), unity;
	rho_Kpi(sqrt(sab - (mK+mpi)*(mK+mpi)),0);
	rho_Ketap(0,sqrt((mK+m0_etap)*(mK+m0_etap) - sab));

	gamma = ((sab - sA)/(mR*mR - sA))*(g_Kpi*rho_Kpi + g_Ketap*rho_Ketap);

	return unity.One()/(cmR*cmR - csab - unity.I()*cmR*gamma);

}

//flatte lineshape (used in f0(980)pi)
////////////////////////////////////////////////////////////////////////
TComplex flatte(Double_t &sij, Double_t &resMass, Double_t gK, Double_t gpi ){

	// constant factors from BES data
	// resMass should be 0.965 +/- 0.008 +/- 0.006 GeV/c^2
	// g1Val = 0.165;       // +/- 0.010 +/- 0.015 GeV/c^2
	// g2Val = g1Val*4.21;  // +/- 0.25 +/- 0.21

	// or from E791
	// g1Val = 0.09;
	// g2Val = 0.02;

	// or from CERN/WA76
	// g1Val = 0.28;
	// g2Val = 0.56;

	//o lo que usamos en Rio+
	// Double_t g1 = 0.22;  //gpi = 0.22;
	// Double_t g2 = 0.76; //gK = 0.76;

	// Double_t resMass = this->getMass();

	Double_t g1        = 0.165,
		 g2        = 0.695, //g1*4.21,
		 mpi0      = 0.1349766,
		 mK0       = 0.497648,
		 resMassSq = resMass*resMass,
		 mSumSq0_   = ( mpi0 + mpi0 ) * ( mpi0 + mpi0 ),
		 mSumSq1_  = ( mpi + mpi ) * ( mpi + mpi ),
		 mSumSq2_  = ( mK + mK ) * ( mK + mK ),
		 mSumSq3_   = ( mK0 + mK0 ) * ( mK0 + mK0 ),
		 dMSq      = resMassSq - sij,
		 rho1(0.0), 
		 rho2(0.0);

	if (gK != -99999999) g2 = gK;              
	if (gpi != -99999999) g1 = gpi;              

	if (sij > mSumSq0_) {
		rho1 = TMath::Sqrt(1.0 - mSumSq0_/sij)/3.0;
		if (sij > mSumSq1_) {
			rho1 += 2.0*TMath::Sqrt(1.0 - mSumSq1_/sij)/3.0;
			if (sij > mSumSq2_) {
				rho2 = 0.5*TMath::Sqrt(1.0 - mSumSq2_/sij);
				if (sij > mSumSq3_) {
					rho2 += 0.5*TMath::Sqrt(1.0 - mSumSq3_/sij );
				} else {
					// Continue analytically below higher channel thresholds
					// This contributes to the real part of the amplitude denominator
					dMSq += g2*resMass*0.5*TMath::Sqrt(mSumSq3_/sij - 1.0);
				}
			} else {
				// Continue analytically below higher channel thresholds
				// This contributes to the real part of the amplitude denominator
				rho2 = 0.0;
				dMSq += g2*resMass*(0.5*TMath::Sqrt(mSumSq2_/sij - 1.0) + 0.5*TMath::Sqrt(mSumSq3_/sij - 1.0));
			}
		} else {
			// Continue analytically below higher channel thresholds
			// This contributes to the real part of the amplitude denominator
			dMSq += g1*resMass*2.0*TMath::Sqrt(mSumSq1_/sij - 1.0)/3.0;
		}
	}

	else {

		cout<<"ERROR -- GenericFunctions -> flatte lineshape: out of phase space "<<sij<<endl;	     
	}

	Double_t width1 = g1*rho1*resMass,
		 width2 = g2*rho2*resMass,
		 widthTerm = width1 + width2;

	TComplex resAmplitude( dMSq, widthTerm );

	Double_t denomFactor = dMSq*dMSq + widthTerm*widthTerm,
		 invDenomFactor = 0.0;

	if ( denomFactor > 1e-10 ) { invDenomFactor = 1.0/denomFactor; }

	resAmplitude *= invDenomFactor;

	return resAmplitude;	

}

// sigma from BES
////////////////////////////////////////////////////////////////////////
TComplex sigma(double sij ){

	Double_t     mPiSq	= mpi*mpi;
	Double_t 	mPiSq4 = 4.0*mPiSq;
	Double_t 	sAdler   = mPiSq*0.5; // Adler zero at 0.5*(mpi)^2

	// constant factors from BES data
	Double_t b1 = 0.5843;
	Double_t b2 = 1.6663;
	Double_t A = 1.082;
	Double_t m0 = 0.9264;

	Double_t	m0Sq = m0*m0;
	Double_t	denom = m0Sq - sAdler;

	//	Double_t s = mass*mass; // Invariant mass squared combination for the system
	Double_t rho(0.0); // Phase-space factor
	if (sij > mPiSq4) {rho = TMath::Sqrt(1.0 - mPiSq4/sij);}

	Double_t f = b2*sij + b1; // f(s) function
	Double_t numerator = sij - sAdler;
	Double_t gamma(0.0);
	if (TMath::Abs(denom) > 1e-10 && TMath::Abs(A) > 1e-10) {   
		// Decay width of the system
		gamma = rho*(numerator/denom)*f*TMath::Exp(-(sij - m0Sq)/A);
	}

	// Now form the complex amplitude - use relativistic BW form (without barrier factors)    
	// Note that the M factor in the denominator is not the "pole" at ~500 MeV, but is 
	// m0 = 0.9264, the mass when the phase shift goes through 90 degrees.

	Double_t dMSq = m0Sq - sij;
	Double_t widthTerm = gamma*m0;
	TComplex resAmplitude(dMSq, widthTerm);

	Double_t denomFactor = dMSq*dMSq + widthTerm*widthTerm;

	Double_t invDenomFactor = 0.0;
	if (denomFactor > 1e-10) {invDenomFactor = 1.0/denomFactor;}

	resAmplitude*=invDenomFactor;

	return resAmplitude;

}

//Function Gounaris  - Sakurai
////////////////////////////////////////////////////////////////////////
TComplex GS( double &sij, double &mR, double &wR, double &gamma, double masq, double mbsq ){

	TComplex resAmp(0.0, 0.0), unity;
	Double_t mRSq = mR*mR,
		 mass = TMath::Sqrt( sij );

	const Double_t massSqTerm = mRSq - sij;
	const Double_t q0_ = sqrt( lambda( mRSq, masq, mbsq ) )/( 2*mR );
	const Double_t q   = sqrt( lambda( sij, masq, mbsq ) )/( 2*mass );

	if (mass < 1e-10) {
		cout << "WARNING in GounarisSakuraiRes::amplitude : mass < 1e-10." << std::endl;
		return TComplex(0.0, 0.0);
	} else if (q0_ < 1e-30) {

		cout << "WARNING in GounarisSakuraiRes::q0_ : mass < 1e-10." << std::endl;		
		return TComplex(0.0, 0.0);
	}

	const Double_t h0_ = ( 2.0/PI ) * q0_/mR * TMath::Log( ( mR + 2.0*q0_ )/( 2.0*mpi ) );
	const Double_t dhdm0_ = h0_ * ( 1.0/( 8.0*q0_*q0_ ) - 1.0/( 2.0*mRSq ) ) + 1.0/( 2*PI*mRSq );
	const Double_t d_ = ( 3.0/PI ) * mpi*mpi/( q0_*q0_ ) * TMath::Log( ( mR + 2.0*q0_ )/( 2.0*mpi ) ) + mR/( 2*PI*q0_ ) - mpi*mpi*mR/( PI*q0_*q0_*q0_ );

	const Double_t h = (2.0/PI) * q/( mass )* TMath::Log(( mass + 2.0*q)/(2.0*mpi));
	const Double_t f = gamma * mRSq/(q0_*q0_*q0_) * (q*q * (h - h0_) + massSqTerm * q0_*q0_ * dhdm0_);

	resAmp = (1 + d_ * wR/mR )*unity.One()/( massSqTerm + f - unity.I()*mR*gamma );

	return resAmp;

}

////function to test if a point is inside the DP limits
////////////////////////////////////////////////////////////////////////
bool DPLims ( double s12,  double s13, double M, double m1, double m2, double m3 ) {

	if ( s12 < pow( m1 + m2, 2 ) ) return false; 

	if ( s12 > pow( M - m3, 2 ) )  return false;  

	double p1str = 0.5 * ( s12 - m2*m2 + m1*m1 ) / sqrt( s12 ),
	       p3str = 0.5 * ( M*M - s12 - m3*m3 ) / sqrt( s12 ),
	       minD  = pow( p1str + p3str, 2 ) - pow( sqrt( p1str*p1str - m1*m1 ) + sqrt( p3str*p3str - m3*m3 ), 2 );

	if ( s13 < minD ) return false;

	double maxD = pow( p1str + p3str, 2 ) - pow( sqrt( p1str*p1str - m1*m1 ) - sqrt( p3str*p3str - m3*m3 ), 2 );

	if ( s13 > maxD ) return false;

	return true; 

}

////////////////////////////////////////////////////////////////////////
double thetaPrime( double MD, double sab, double sac, double ma, double mb, double mc ) {

	double term_ma= ma*ma;
	double term_mb= mb*mb;
	double term_mc= mc*mc;
	double term_md= MD*MD;

	double pa       =    sqrt( lambda( sab, term_ma, term_mb ) )/( 2*sqrt( sab ) ), //qi
	       pc       =    sqrt( lambda( term_md, sab, term_mc ) )/( 2*sqrt( sab ) ), /// qk
	       easab    =  ( sab - term_mb + term_ma )/( 2*sqrt( sab ) ),
	       ecsac    =  ( MD*MD - sab - term_mc )/( 2*sqrt( sab ) ),
	       thetaHel = -( sac - term_ma - term_mc - 2*easab*ecsac )/( 2*pa*pc );

	return ( 1/PI )*acos( thetaHel );

}

////////////////////////////////////////////////////////////////////////
double mPrime( double sab, double MD, double ma, double mb, double mc ) {

	double m = 2*( sqrt( sab ) - ( ma + mb ) )/( MD - mc - ( ma + mb ) ) - 1.0;

	return ( 1/PI )*acos( m );

}

////////////////////////////////////////////////////////////////////////
double CalcJacob ( double sab, double sac, double MD, double ma, double mb, double mc ) {

	if( !DPLims ( sab, sac, MD, ma, mb, mc ) )  return 0;

	double term_ma= ma*ma;
	double term_mb= mb*mb;
	double term_mc= mc*mc;
	double term_md= MD*MD;


	double pa       = sqrt( lambda( sab, term_ma, term_mb ) )/( 2*sqrt( sab ) ),
	       pc       = sqrt( lambda( term_md, sab, term_mc ) )/( 2*sqrt( sab ) ),      
	       mp       = mPrime( sab, MD, ma, mb, mc ),
	       tp       = thetaPrime( MD, sab, sac, ma, mb, mc ),
	       JacobDet = 2*pa*pc*sqrt( sab )*PI*PI*sin( PI*mp )*( MD - 3*ma )*sin( PI*tp );

	return JacobDet;

}


// Insert here the acceptance function of your sample - used in Spdf
string Acceptance_Ntuple_Name, Acceptance_Histo_Name;
bool UsesWeights, UseAcceptance;
TFile *AccNtpFile;
TH2D *AccHist;
double s12MinForAcc, s12MaxForAcc, s13MinForAcc, s13MaxForAcc, s12ExtForAcc, s13ExtForAcc;
int NBins12ForAcc, NBins13ForAcc;


////////////////////////////////////////////////////////////////////////
double acceptance_sqDp( double s12, double s13 )
{
	double M  = D_Mass,
	       ma = mpi,
	       mb = mpi,
	       mc = mpi;

	double  tprime    =   thetaPrime( M, s12,  s13, ma, mb, mc );

	double massprime =   mPrime(  s12, M, ma, mb, mc );  

	int s12Flaw = int( ( massprime - s12MinForAcc )*NBins12ForAcc / s12ExtForAcc ),
	    s13Flaw = int( ( tprime - s13MinForAcc )*NBins13ForAcc / s13ExtForAcc );	 

	double Acc = AccHist->GetBinContent( s13Flaw + 1, s12Flaw + 1 );

	return Acc;

}

////////////////////////////////////////////////////////////////////////
double acceptance_Sij( double s12, double s13 ){

	int s12Flaw = int( ( s12 - s12MinForAcc )*NBins12ForAcc / s12ExtForAcc ),
	    s13Flaw = int( ( s13 - s13MinForAcc )*NBins13ForAcc / s13ExtForAcc );
	double  Acc = AccHist->GetBinContent( s13Flaw + 1, s12Flaw + 1 );

	return Acc;

}

////////////////////////////////////////////////////////////////////////
double acceptance( double s12, double s13 )
{

	return acceptance_Sij( s12, s13 );

}


// Insert here the background histo
////////////////////////////////////////////////////////////////////////
string Back_Ntuple_Name, Back_Histo_Name;
bool UseBackHisto;
TFile *BackNtpFile;
TH2D *BackHist;
double s12MinForBack, s12MaxForBack, s13MinForBack, s13MaxForBack, s12ExtForBack, s13ExtForBack;
int NBins12ForBack, NBins13ForBack;

////////////////////////////////////////////////////////////////////////
double backgroundF_sQ ( double s12, double s13 ){

	double M  = D_Mass, ma = mpi, mb = mpi, mc = mpi;

	double  tprime    =   thetaPrime( M, s12, s13, ma, mb, mc );

	double massprime =   mPrime(  s12, M, ma, mb, mc );  


	int s12BFlaw = int(( massprime - s12MinForBack )*NBins12ForBack / s12ExtForBack );
	int s13BFlaw = int(( tprime - s13MinForBack )*NBins13ForBack / s13ExtForBack );
	double  Back = BackHist->GetBinContent( s13BFlaw + 1, s12BFlaw + 1 );
	return Back;
}


////////////////////////////////////////////////////////////////////////
double backgroundF_Sij ( double s12, double s13 ){

	int s12BFlaw = int(( s12 - s12MinForBack )*NBins12ForBack / s12ExtForBack );
	int s13BFlaw = int(( s13 - s13MinForBack )*NBins13ForBack / s13ExtForBack );
	double 	Back = BackHist->GetBinContent( s13BFlaw + 1, s12BFlaw + 1 );

	return Back;

}

////////////////////////////////////////////////////////////////////////
double backgroundF( double s12, double s13 )
{

	return backgroundF_Sij( s12, s13 );

}

////////////////////////////////////////////////////////////////////////
void DrawMypull( TH1F *hdata, TH1F *hfunc, TH1F *htoy ) {

	int    nbins = hdata->GetNbinsX();
	double xmin  = hdata->GetXaxis()->GetBinLowEdge( 1 );
	double xmax  = hdata->GetXaxis()->GetBinLowEdge( hdata->GetNbinsX()) + hdata->GetXaxis()->GetBinWidth( 5 );
	double diff, err; 
	TString pullName = "PullHist";
	TH1F *PullHist = new TH1F( pullName, "", nbins, xmin, xmax );
	for( int i = 0; i < nbins + 1; ++i ) {
		diff =  hfunc->GetBinContent( i + 1 ) - hdata->GetBinContent( i + 1 );

		if (htoy->GetBinContent( i + 1 ) != 0) err =  (hfunc->GetBinContent( i + 1 )/htoy->GetBinContent( i + 1 ))*htoy->GetBinError( i + 1 ) + sqrt(hdata->GetBinContent( i + 1 ));
		else err = 0;

		if (diff != 0 && err != 0) PullHist->SetBinContent(i+1,diff/err);
		else PullHist->SetBinContent(i+1,0);  
	}
	PullHist->SetXTitle( "" );
	PullHist->SetFillColor(kBlue+3);
	PullHist->SetBarWidth(0.75);
	PullHist->SetBarOffset(0.1);
	PullHist->SetStats(0);
	PullHist->SetMinimum(-5);
	PullHist->SetMaximum(5);
	PullHist->SetTitleFont(62);
	PullHist->SetTitleOffset(1.1,"y");
	// PullHist->SetTextSize(0.08);
	PullHist->SetTitleSize(0.3,"y");
	PullHist->SetLabelOffset(0.005);
	PullHist->GetXaxis()->SetLabelSize(0.);
	//PullHist->SetTitleSize(0.06,"z");
	PullHist->GetYaxis()->SetNdivisions(2, kTRUE);
	PullHist->GetYaxis()->SetLabelSize(0.2);
	PullHist->SetYTitle( "#bf{Pull}" );

	PullHist->Draw("b");

}

////////////////////////////////////////////////////////////////////////
TH2D* CalculateChiSq(  TH2D* pdf, TH2D* toy, TH2D* dat, int nparams, int &ndof, double &totalChiSq, double &kolmog ) {

	int nbins=toy->GetNbinsX()*toy->GetNbinsY();
	double xmin = toy->GetXaxis()->GetBinLowEdge(1),
	       xmax = toy->GetXaxis()->GetBinLowEdge(toy->GetNbinsX())+toy->GetXaxis()->GetBinWidth(5),
	       ymin = toy->GetYaxis()->GetBinLowEdge(1),
	       ymax = toy->GetYaxis()->GetBinLowEdge(toy->GetNbinsY())+toy->GetYaxis()->GetBinWidth(5);

	UNUSED(kolmog);
	TH2D* histoChi2 = new TH2D( "histoChi2", "#chi^{2} Histogram", toy->GetNbinsX(), xmin, xmax, toy->GetNbinsY(), ymin, ymax );

	int minimum    = 20,
	    nob        = 0;
	totalChiSq = 0;
	double toyVal, 
	       pdfVal, 
	       dataVal, 
	       diff,
	       errsq,
	       errToy,
	       errPDF;

	for( Int_t i = 1; i <= nbins; ++i ) 
	{	
		Double_t chiSq = 0.;
		pdfVal  = pdf->GetBinContent(i);
		errPDF  = pdf->GetBinError(i);
		toyVal  = toy->GetBinContent(i);
		errToy  = toy->GetBinError(i);
		dataVal = dat->GetBinContent(i);
		diff = pdfVal-dataVal;
		errsq = pow((pdfVal/toyVal)*errToy,2) + dataVal;
		if( pdfVal > 0 && dataVal > minimum ) 
		{

			chiSq = (diff*diff)/errsq;
			++nob;
		}
		totalChiSq += chiSq;
		if (diff > 0)   histoChi2->SetBinContent(i, chiSq);
		else    histoChi2->SetBinContent(i, -chiSq);
	}
	//kolmog = toy->KolmogorovTest(dat);
	ndof = nob-nparams-1;
	cout << "number of bins with dataval > 20: " << nob; 
	cout<<"Total ChiSq/nDof = "<<totalChiSq<<"/"<<ndof<<" = "<<totalChiSq/ndof<<endl;
	//totalChiSq/=ndof;
	//delete toy;
	//delete dat;

	return  histoChi2;

}

////////////////////////////////////////////////////////////////////////
TH2Poly* CalculateChiSq_adaptive( TH2Poly* pdf, TH2Poly* toy, TH2Poly* dat, int nparams, int &ndof, double &totalChiSq, double &kolmog ) {

	//int nbins=toy->GetNbinsX()*toy->GetNbinsY();
	int nbins=toy->GetNumberOfBins();
	TH2Poly  *histoChi2 =  new TH2Poly("Chi2", "Chi2", 0.0, 2.1, 0.9, 3.1);
	double Bins12Max, Bins13Min, Bins13Max, Bins12Min;
	UNUSED(kolmog);
	ifstream infile("bins_2025Dkkpi_TISorTOS.txt",ios::in);

	if (!infile.good()) {
		std::cout << "BinAreasAndInt - ERROR: File not found - infile" << std::endl;
		exit(-1);
	}

	for (int i = 0; i<toy->GetNumberOfBins(); i++) {

		infile >> Bins13Min;
		infile >> Bins13Max;
		infile >> Bins12Min;
		infile >> Bins12Max;

		histoChi2->AddBin(Bins13Min, Bins12Min, Bins13Max, Bins12Max);
	}

	//cout<<" CalculateChiSq: cheke clone = "<<endl;
	histoChi2-> ClearBinContents();
	//cout<<" CalculateChiSq: cheke clear = "<<endl;
	histoChi2->SetXTitle("s_low");
	histoChi2->SetYTitle("s_high");
	//cout<<" CalculateChiSq: cheke titles = "<<endl;
	int minimum    = 20,
	    nob        = 0;
	totalChiSq = 0;
	double toyVal, 
	       pdfVal, 
	       dataVal, 
	       diff,
	       errsq,
	       errToy,
	       errPDF;

	for( Int_t i = 1; i <= nbins; ++i ) 
	{	
		Double_t chiSq = 0.;
		toyVal  = toy->GetBinContent(i);
		pdfVal  = pdf->GetBinContent(i);
		dataVal = dat->GetBinContent(i);
		errPDF  = pdf->GetBinError(i);
		errToy  = toy->GetBinError(i);


		diff = pdfVal-dataVal;
		errsq = pow((pdfVal/toyVal)*errToy,2) + dataVal;
		cout << "toyVal = " << toyVal << "dataVal = " << dataVal << endl;
		if( /*toyVal > minimum &&*/ dataVal > minimum ) 
		{			
			chiSq = (diff*diff)/errsq;
			++nob;
			cout << "Bin " << i << " -  toyVal = " << toyVal << ", errToy = " << errToy << ", pdfVal =  " << pdfVal  << ", errPDF = " << errPDF << ", dataVal =  " << dataVal 
				<< ", errsq =  " << errsq  << ", Chi2Sq = " 
				<< chiSq  << ", totalChiSq = " << totalChiSq << endl; 
			//		cout << "Bin " << i << " -  toyVal = " << toyVal << ", pdfVal =  " << pdfVal  << ", dataVal =  " << dataVal  << ", errsq =  " << errsq  << ", Chi2Sq = " 
			//			<< chiSq  << ", totalChiSq = " << totalChiSq << endl; 
		}
		totalChiSq += chiSq;
		if (diff > 0)	histoChi2->SetBinContent(i, chiSq);
		else	histoChi2->SetBinContent(i, -chiSq);
	}
	//	kolmog = toy->KolmogorovTest(dat);
	ndof = nob-nparams-1;
	//totalChiSq/=ndof;			
	cout<<"Total ChiSq/nDof = "<<totalChiSq<<endl;

	return  histoChi2;

}
////////////////////////////////////////////////////////////////////////
