#include <TFitter.h>
#include <TTree.h>
#include <TFile.h>
#include <TColor.h>
#include <TCanvas.h>
#include <TMultiGraph.h>
#include <TGraphErrors.h>
#include "TStyle.h"
#include "TGaxis.h"
#include "TLegend.h"
#include "TLatex.h"
#include <TStopWatch.h>
#include <iomanip>
#include <fstream>
#include <sstream>
#include "Read_Parameters.h"
#include "SigNorm.h"
#include "TotPdf.h"

void plot_moments (double *m, double *bin_pop, double *t0, double *t1, double *t2, double *t3, double *t4, double *t0_err, double *t1_err, double *t2_err, double *t3_err, double *t4_err, TGraphErrors *&gr0, TGraphErrors *&gr1, TGraphErrors *&gr2, TGraphErrors *&gr3, TGraphErrors *&gr4);

void Draw_Graph_Pulls(TGraphErrors * g1, TGraphErrors * g2);

void set_plot_style_Dalitz()
{
	const Int_t NRGBs = 5;
	const Int_t NCont = 255;

	Double_t stops[NRGBs] = { 0.00, 0.25, 0.5, 0.75, 1.00 };
	Double_t red[NRGBs]   = { 1.0, 0.3, 0.0, 0.5, 0.5 };
	Double_t green[NRGBs] = { 1.0, 0.3, 0.0, 0.3, 0.0 };
	Double_t blue[NRGBs]  = { 1.0, 0.5, .5, 0.3, 0.0 };
	TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
	gStyle->SetNumberContours(NCont);
}

void set_plot_style_chi2()
{
	const Int_t NRGBs = 5;
	const Int_t NCont = 255;

	Double_t stops[NRGBs] = { 0.00, 0.25, 0.5, 0.75, 1.00 };
	Double_t red[NRGBs]   = { 0.3, 0.8, 1.0, 0.5, 0.0 };
	Double_t green[NRGBs] = { 0.00, 0.5, 1.0, 0.8, 0.3 };
	Double_t blue[NRGBs]  = { 0.00, 0.5, 1.0, 0.5, 0.0 };
	TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
	gStyle->SetNumberContours(NCont);
}

void MakeMCPlots(string input_txt_file_name){

	double M = D_Mass, s = M*M, m1=mK, m2=mK, m3=mpi, m1sq = m1*m1, m2sq = m2*m2, m3sq = m3*m3, s12, s23, s13, z, signalweight, WeightSum, SqrWeightSum;
	vector<double> S12;
	vector<double> S13;
	vector<double> S23;
	vector<double> Mass;
	vector<double> KK_bin_limit;
	vector< vector<TComplex> > Events_Signal_Amplitudes;
	vector< vector<TComplex> > Events_Bkg_Amplitudes;
	vector<double> sWeight;
	double sFactor;

	vector<double> Acceptance;

	double bkg_fraction, mass_min, mass_max, bkg_par1, bkg_par2, FCN_value;
	vector<int> resonances, bkg_components;
	vector<TComplex> coefs, bkg_coefs;
	int final_state, number_of_resonances, number_of_bkg_components, number_of_pwa_bins, iteration_number = 0;

	Long64_t nentries;

	vector< vector<TComplex> > sig_normalization;
	vector<double> bkg_normalization, last_iteration_pwa_parameters, last_iteration_masses, last_iteration_widths;
	vector< vector<double> > last_iteration_res_extra_pars;
	bool is_gaussian, real_and_imaginary;

	string input_ntuple_name, output_pwa_txt_file;
	vector<string> included_resonant_channel_string;
	int seed, resonance_number = 0, bkg_component_number = 0, npt, number_of_samples, event_number=0;
	TRandom1 * Random = new TRandom1(time(NULL));
	double best_coef1, error_coef1, best_coef2, error_coef2, best_mass, best_width, error_mass, error_width;
	vector<int> fix_parameter_index;
	vector<TComplex> pwa_coefs, pwa_coefs_prime, sig_amp, bkg_amp;
	bool resonances_calc_norm[AVAILABLE_RESONANCES];
	Resonant_Amplitudes res_amp;
	TStopwatch clock;
	clock.Start();

	int moments_npt = 100, bin=0;
	double m_KK[moments_npt], t0_KK[moments_npt], t1_KK[moments_npt], t2_KK[moments_npt], t3_KK[moments_npt], t4_KK[moments_npt], bin_pop_KK[moments_npt];
	double t0sq_KK[moments_npt], t1sq_KK[moments_npt], t2sq_KK[moments_npt], t3sq_KK[moments_npt], t4sq_KK[moments_npt];
	double t0_err_KK[moments_npt], t1_err_KK[moments_npt], t2_err_KK[moments_npt], t3_err_KK[moments_npt], t4_err_KK[moments_npt];

	double m_Kpi[moments_npt], t0_Kpi[moments_npt], t1_Kpi[moments_npt], t2_Kpi[moments_npt], t3_Kpi[moments_npt], t4_Kpi[moments_npt], bin_pop_Kpi[moments_npt];
	double t0sq_Kpi[moments_npt], t1sq_Kpi[moments_npt], t2sq_Kpi[moments_npt], t3sq_Kpi[moments_npt], t4sq_Kpi[moments_npt];
	double t0_err_Kpi[moments_npt], t1_err_Kpi[moments_npt], t2_err_Kpi[moments_npt], t3_err_Kpi[moments_npt], t4_err_Kpi[moments_npt];

	double s12min, s12max, s13min, s13max, sabmin, sabmax, sacmin, sacmax, ds12, ds13, threshold, Integral_KK = 0, IntegralSq_KK = 0, Integral_Kpi = 0, IntegralSq_Kpi = 0;
	double s12_limits[moments_npt], s13_limits[moments_npt], PScorr_KK[moments_npt], PScorr_Kpi[moments_npt], PL[8];
	double E1, E2, E3, p1, p2, p3;
	double acc, w, sigma_w, av_w, av_w2;

	TGraphErrors *g0_KK_Data, *g1_KK_Data, *g2_KK_Data, *g3_KK_Data, *g4_KK_Data;
	TGraphErrors *g0_Kpi_Data, *g1_Kpi_Data, *g2_Kpi_Data, *g3_Kpi_Data, *g4_Kpi_Data;
	TGraphErrors *g0_KK_MC, *g1_KK_MC, *g2_KK_MC, *g3_KK_MC, *g4_KK_MC;
	TGraphErrors *g0_Kpi_MC, *g1_Kpi_MC, *g2_Kpi_MC, *g3_Kpi_MC, *g4_Kpi_MC;

	TMultiGraph *mg0_Kpi = new TMultiGraph();
	TMultiGraph *mg1_Kpi = new TMultiGraph();
	TMultiGraph *mg2_Kpi = new TMultiGraph();
	TMultiGraph *mg3_Kpi = new TMultiGraph();
	TMultiGraph *mg4_Kpi = new TMultiGraph();
	TMultiGraph *mg0_KK = new TMultiGraph();
	TMultiGraph *mg1_KK = new TMultiGraph();
	TMultiGraph *mg2_KK = new TMultiGraph();
	TMultiGraph *mg3_KK = new TMultiGraph();
	TMultiGraph *mg4_KK = new TMultiGraph();

	TGaxis::SetMaxDigits(2);
	gStyle->SetOptStat(0);
	gStyle->SetPadRightMargin(0.18);
	gStyle->SetPadLeftMargin(0.15);
	gStyle->SetPadBottomMargin(0.20);
	gStyle->SetOptTitle(0);
	gStyle->SetLabelSize(0.05,"x");
	gStyle->SetLabelSize(0.05,"y");
	gStyle->SetLabelSize(0.05,"z");
	gStyle->SetTitleSize(0.06,"x");
	gStyle->SetTitleSize(0.06,"y");
	gStyle->SetTitleSize(0.06,"z");

	vector< vector< TComplex > > SigAmps,
		BkgAmps,
		normalization_sig;      

	vector< double > normalization_bkg,  
		res_masses, 
		res_widths,
		tmp_res_extra_pars(4);

	vector< vector<double> > res_extra_pars;

	TComplex temp_pwa_coef, 
		 tmp;

	Total_PDF TotPdf;

	TFitter * minimizer = new TFitter(200);

	double spd_prod = 0, L_s = 0;

	number_of_resonances = 0;
	number_of_bkg_components = 0;
	resonances.clear();
	fix_parameter_index.clear();
	bkg_components.clear();
	KK_bin_limit.clear();
	pwa_coefs.clear();
	iteration_number = 0;

	Read_Parameters(input_txt_file_name, minimizer, real_and_imaginary, final_state, is_gaussian, bkg_par1, bkg_par2, number_of_samples, resonances, bkg_components, number_of_resonances,
			number_of_bkg_components, bkg_fraction, input_ntuple_name, seed, fix_parameter_index, number_of_pwa_bins, output_pwa_txt_file, pwa_coefs, KK_bin_limit, UsesWeights, UseAcceptance, Acceptance_Ntuple_Name, Acceptance_Histo_Name,UseBackHisto, Back_Ntuple_Name, Back_Histo_Name );

	int  nnh = 65,
	     nevts, 
	     res_number     = 0,
	     pwa_par_number = KK_bin_limit.size(),
	     nparams        = 8 * number_of_resonances - fix_parameter_index.size();

	coefs.resize( number_of_resonances );

	double par[ 8 * number_of_resonances + number_of_bkg_components + 2 * pwa_par_number ];
	vector< vector< TComplex > > coefs_product( number_of_resonances, vector< TComplex > ( number_of_resonances ) );

	for (int i = 0; i < 8 * number_of_resonances + number_of_bkg_components +
			2 * pwa_par_number; i++ ) {

		par[ i ] = minimizer->GetParameter( i );

	}


	for ( int i = 0; i < number_of_resonances; i++ ) {

		if ( real_and_imaginary )

			tmp( par[ 8 * i ], par[8 * i + 1 ] );

		else

			tmp( par[ 8 * i ], par[ 8 * i + 1 ], 1 );
		coefs[ i ] = tmp;
	}

	for ( int res_i = 0; res_i < number_of_resonances; res_i++ ) {

		for ( int res_j = 0; res_j < number_of_resonances; res_j++ ) {

			coefs_product[ res_i ][ res_j ] =
				coefs[ res_i ] * coefs[ res_j ].Conjugate( coefs[ res_j ] );
		}
	}

	for ( int i = 0; i < AVAILABLE_RESONANCES; i++ ) {

		if ( resonances[ i ] == 1 ) {

			res_masses.push_back( par[ 8 * res_number + 2 ] );
			res_widths.push_back( par[ 8 * res_number + 3 ] );
			//cout << "res_mass = " << res_masses[res_number] << ", res_width = " << res_widths[res_number] << endl; 
			for ( int extra_par = 0; extra_par < 4; extra_par++ ) {
				tmp_res_extra_pars[extra_par] = par[ 8 * res_number + 4 + extra_par];
			}
			res_number++;
			resonances_calc_norm[ i ] = 1;
			res_extra_pars.push_back(tmp_res_extra_pars);
		}

	}

	TString colNameLatexX1,
		colNameLatexY1,
		colNameLatexY2;

	//cout << "final_state = " << final_state << endl; 

	switch ( final_state ) {
		case 0:
			colNameLatexX1 = "#bf{#it{S}_{K^{-}#pi^{+}} [GeV^{2}/#it{c}^{4}]}";
			colNameLatexY1 = "#bf{#it{S}_{K^{-}K^{+}} [GeV^{2}/#it{c}^{4}]}";
			colNameLatexY2 = "#bf{#it{S}_{K^{+}#pi^{+}} [GeV^{2}/#it{c}^{4}]}";
			//cout << "passou dentro " << endl;
			break;
		case 2: 
			colNameLatexX1 = "#bf{#it{S}_{#pi^{+}#pi^{-},low} [GeV^{2}/#it{c}^{4}]}";
			colNameLatexY1 = "#bf{#it{S}_{#pi^{+}#pi^{-},high}  [GeV^{2}/#it{c}^{4}]}";
			colNameLatexY2 = "#bf{#it{S}_{#pi^{+}#pi^{-},23}  [GeV^{2}/#it{c}^{4}]}";
			break;
		case 3:
			colNameLatexX1 = "#bf{#it{S}_{K^{+}K^{-},low} [GeV^{2}/#it{c}^{4}]}";
			colNameLatexY1 = "#bf{#it{S}_{K^{+}K^{-},high}  [GeV^{2}/#it{c}^{4}]}";
			break;
		default:
			cout << "Invalid final state" << endl;
			exit( -1 );
	}


	TH2D *TotalPDFHist  = 
		new TH2D( "TotalPDFHist", "TotPdf", nnh, .3, 2.1, nnh, .9, 3.1 ),
		    *TotalToyHist = 
			    new TH2D( "TotalToyHist", "TotToy", nnh, 0.3, 2.1, nnh, 0.9, 3.1 ),
		    *TotalDataHist = 
			    new TH2D( "TotalDataHist", "TotData", nnh, 0.3, 2.1, nnh, 0.9, 3.1 );

	TH1F *S12Hist      = 
		new TH1F( "S12Hist", "S12Hist", nnh, 0.9, 3.1 ),
		    *S12HistData  = 
			    new TH1F( "S12HistData", "S12HistData", nnh, 0.9, 3.1 ),
		    *S12HistToy  = 
			    new TH1F( "S12HistToy", "S12HistToy", nnh, 0.9, 3.1 ),
		    *S13Hist       = 
			    new TH1F( "S13Hist", "S13Hist", nnh, 0., 2.1 ),
		    *S13HistData   = 
			    new TH1F( "S13HistData", "S13HistData", nnh, 0., 2.1 ),
		    *S13HistToy   = 
			    new TH1F( "S13HistToy", "S13HistToy", nnh, 0., 2.1 ),
		    *S23Hist       = 
			    new TH1F( "S23Hist", "S23Hist", nnh, 0., 2.1),
		    *S23HistData   = 
			    new TH1F( "S23HistData", "S23HistData", nnh, 0., 2.1 ),
		    *S23HistToy   = 
			    new TH1F( "S23HistToy", "S23HistToy", nnh, 0., 2.1 ),

		    *histo_sij[ 3 ][ number_of_resonances ],
		    *bhisto_sij[ 3 ][ number_of_bkg_components ];

	std::vector<double> sigPdf( number_of_resonances );

	for ( Int_t iii = 0; iii < number_of_resonances; iii++ ) {

		histo_sij[ 0 ][ iii ] = 
			new TH1F( Form( "histo_s12%d", iii ), "", nnh, 0.9, 3.1 );
		histo_sij[ 1 ][ iii ] = 
			new TH1F( Form( "histo_s13%d", iii ), "", nnh, 0., 2.1 );
		histo_sij[ 2 ][ iii ] = 
			new TH1F( Form( "histo_s23%d", iii ), "", nnh, 0., 2.1 );
		sigPdf[ iii ] = 0;

	}

	for ( Int_t bb = 0; bb < number_of_bkg_components; ++bb ) {

		bhisto_sij[ 0 ][ bb ] = 
			new TH1F( Form( "bhisto_s12%d", bb ), "", nnh, 0.9, 3.1 );
		bhisto_sij[ 1 ][ bb ] = 
			new TH1F( Form( "bhisto_s13%d", bb ), "", nnh, 0., 2.1 );
		bhisto_sij[ 2 ][ bb ] = 
			new TH1F( Form( "bhisto_s23%d", bb ), "", nnh, 0., 2.1 );

	}

	int Adaptative_bin_number = 2025;
	double Bins12Min, Bins13Min, Bins12Max, Bins13Max;
	TH2Poly *TotalPDFHist_Chi2  = new TH2Poly("TotalPDF", "TotalPDF", 0.0, 2.1, 0.9, 3.1);
	TH2Poly *TotalToyHist_Chi2  = new TH2Poly("TotalToy", "TotalToy", 0.0, 2.1, 0.9, 3.1);
	TH2Poly  *TotalDataHist_Chi2 =  new TH2Poly("TotalData", "TotalData", 0.0, 2.1, 0.9, 3.1);

	ifstream infile("bins_2025Dkkpi_TIS.txt",ios::in);

	for (int i = 0; i<Adaptative_bin_number; i++) {

		infile >> Bins13Min;
		infile >> Bins13Max;
		infile >> Bins12Min;
		infile >> Bins12Max;

		TotalPDFHist_Chi2->AddBin(Bins13Min, Bins12Min, Bins13Max, Bins12Max);
		TotalToyHist_Chi2->AddBin(Bins13Min, Bins12Min, Bins13Max, Bins12Max);
		TotalDataHist_Chi2->AddBin(Bins13Min, Bins12Min, Bins13Max, Bins12Max);
		//cout << "Bin " << i << " - Bins12Min = " << Bins12Min << ", Bins12Max = " << Bins12Max << ", Bins13Min = " << Bins13Min 
		//	<< ", Bins13Max = " << Bins13Max << endl;
	}

	for ( unsigned int pwa_par_number = 0; pwa_par_number < KK_bin_limit.size(); 
			pwa_par_number++ ) {

		if ( real_and_imaginary ) {

			temp_pwa_coef( par[ 8 * number_of_resonances + number_of_bkg_components +
					2 * pwa_par_number ],
					par[ 8 * number_of_resonances + number_of_bkg_components +
					2 * pwa_par_number + 1 ] );
		} else {

			temp_pwa_coef( par[ 8 * number_of_resonances + number_of_bkg_components +
					2 * pwa_par_number ],
					par[ 8 * number_of_resonances + number_of_bkg_components +
					2 * pwa_par_number + 1 ], 1 );
		}

		pwa_coefs.push_back( temp_pwa_coef );

	}
	npt = KK_bin_limit.size();

	Complex_Derivative( KK_bin_limit, 
			pwa_coefs, 
			npt, 
			pwa_coefs_prime );

	sig_amp.resize( number_of_resonances );
	bkg_amp.resize( number_of_bkg_components );

	clock.Stop();
	cout << "Finished defining variables " << endl;
	clock.Print("u");

	clock.Start();

	ComplexSigNorm( final_state, 
			is_gaussian, 
			mass_min, 
			mass_max, 
			bkg_par1,
			bkg_par2, 
			number_of_resonances, 
			number_of_bkg_components,
			resonances, 
			bkg_components, 
			res_masses, 
			res_widths,
			res_extra_pars,
			resonances_calc_norm, 
			normalization_sig, 
			normalization_bkg,
			KK_bin_limit, 
			pwa_coefs, 
			pwa_coefs_prime );

	double NL_s = SigNorm( normalization_sig, 
			number_of_resonances, 
			par, 
			real_and_imaginary );

	clock.Stop();
	cout << "Finished calculating normalization" << endl;
	clock.Print("u");

	clock.Start();


	s13min = (m1+m3)*(m1+m3);
	s13max = (M - m2)*(M - m2);
	s12min = (m1+m2)*(m1+m2);
	s12max = (M - m3)*(M - m3);

	ds12 = (s12max - s12min)/moments_npt;
	ds13 = (s13max - s13min)/moments_npt;

	for(int i=0;i<moments_npt;i++)  {

		// initialize moments
		t0_KK[i]=0;
		t1_KK[i]=0;
		t2_KK[i]=0;
		t3_KK[i]=0;
		t4_KK[i]=0;
		bin_pop_KK[i] = 0;

		// initialize  variables for errors on moments     
		t0sq_KK[i] = 0;
		t1sq_KK[i] = 0;
		t2sq_KK[i] = 0;
		t3sq_KK[i] = 0;
		t4sq_KK[i] = 0;

		// define bin limits 
		s12_limits[i] = s12min + (i+1)*ds12 - ds12/2;
		m_KK[i] = sqrt(s12_limits[i]);
		//cout << "m[i] = " << m[i] << endl; 

		// initialize moments
		t0_Kpi[i]=0;
		t1_Kpi[i]=0;
		t2_Kpi[i]=0;
		t3_Kpi[i]=0;
		t4_Kpi[i]=0;
		bin_pop_Kpi[i] = 0;

		// initialize  variables for errors on moments     
		t0sq_Kpi[i] = 0;
		t1sq_Kpi[i] = 0;
		t2sq_Kpi[i] = 0;
		t3sq_Kpi[i] = 0;
		t4sq_Kpi[i] = 0;

		// define bin limits 
		s13_limits[i] = s13min + (i+1)*ds13 - ds13/2;
		m_Kpi[i] = sqrt(s13_limits[i]);


		//compute phase space correction
		s13 = s13_limits[i];
		p1=0.5*sqrt( ( (s13 - m3sq - m1sq)*(s13 - m3sq - m1sq) - 4*m3sq*m1sq )/s13 );
		p2=0.5*sqrt( ( (s - s13 - m2sq)*(s - s13 - m2sq) - 4*s13*m2sq )/s13 );
		E1 = (s13 + m1sq - m3sq) / (2*sqrt(s13));
		E2 = (s - s13 - m2sq)/(2*sqrt(s13));
		sacmax =  m1sq + m2sq + 2*E1*E2 + 2*p1*p2;
		sacmin =   m1sq + m2sq + 2*E1*E2 - 2*p1*p2;
		PScorr_Kpi[i] = 1/(sacmax - sacmin);

		s12 = s12_limits[i];
		p1=0.5*sqrt( ( (s12 - m1sq - m2sq)*(s12 - m1sq - m2sq) - 4*m1sq*m2sq )/s12 );
		p3=0.5*sqrt( ( (s - s12 - m3sq)*(s - s12 - m3sq) - 4*s12*m3sq )/s12 );
		E1 = (s12 + m1sq - m2sq) / (2*sqrt(s12));
		E3 = (s - s12 - m3sq)/(2*sqrt(s12));
		sacmax =  m1sq + m3sq + 2*E1*E3 + 2*p1*p3;
		sacmin =   m1sq + m3sq + 2*E1*E3 - 2*p1*p3;
		PScorr_KK[i] = 1/(sacmax - sacmin);
		//cout << "PScorr[i] = " << PScorr[i] << endl;
		//cout << "PScorr_Kpi[i] = " << PScorr_Kpi[i] << ", PScorr_KK[i] = " << PScorr_KK[i] << endl; 
		//cout << ", s12_limits[i] = " << s12_limits[i]  << ", s13_limits[i] = " << s13_limits[i]  << endl;

	}


	//TFile f0( "/Users/dvieira/trabalho/charm/Data2012_stripping20/KKPfittCMTIS.root" );
	TFile f0( "/Users/dvieira/trabalho/charm/Data2012_stripping20/KKPfitTIS.root" );
	//TFile f0( "/Users/dvieira/trabalho/charm/Data2012_stripping20/KKPfitTISorTOS.root" );
	//TFile f0( "/Users/dvieira/Dropbox/toyMCgenerator/ntuple_kkpi_model4_kappaCleo_0.root" );
	TTree *t0 = (TTree*)f0.Get( "DecayTree" );

	t0->SetBranchAddress( "s12_KK_DTF",    &s12 );
	t0->SetBranchAddress( "s13_Kpi_DTF",    &s13 );
	t0->SetBranchAddress( "s23_Kpi_DTF",    &s23 );
	//t0->SetBranchAddress( "s12",    &s12 );
	//t0->SetBranchAddress( "s13",    &s13 );
	//t0->SetBranchAddress( "s23",    &s23 );
	nentries = t0->GetEntries();
	//nentries = 10000;

	for ( Long64_t jentry = 0; jentry < nentries; ++jentry ) {

		t0->GetEntry( jentry );

		TotalDataHist->Fill( s13, s12 );
		TotalDataHist_Chi2->Fill( s13, s12 );
		S12HistData ->Fill( s12 );
		S13HistData  ->Fill( s13 );
		S23HistData  ->Fill( s23 );


		p1=0.5*sqrt( ( (s12 - m1sq - m2sq)*(s12 - m1sq - m2sq) - 4*m1sq*m2sq )/s12 );
		p3=0.5*sqrt( ( (s - s12 - m3sq)*(s - s12 - m3sq) - 4*s12*m3sq )/s12 );
		E1 = (s12 + m1sq - m2sq) / (2*sqrt(s12));
		E3 = (s - s12 - m3sq)/(2*sqrt(s12));
		z = (m1sq + m3sq + 2*E1*E3 - s13)/(2*p1*p3);

//		if(fabs(z)>1) cout << "z = " << z << ", s12 = " << s12 << ", s13 = " << s13 << endl;

		if (fabs(z)<1){

			// Legendre polynomials for this event
			PL[0] = 1.;
			PL[1] = z;
			PL[2] = (3.*z*z-1.)/2.;
			PL[3] = (5.*pow(z,3)-3.*z)/2.;
			PL[4] = (35.*pow(z,4)-30.*z*z+3.)/8.;

			// check if s23 is in the appropriate mass region
			if((s12 > s12min)&&(s12 < s12max)){

				// acceptance for this event       
				if (s12>1.6) acc = acceptance(s12,s13);

			//	cout << " acc = " << acc << endl;

				// find in which mKK bin this event belongs to;
				if(s12 < s12_limits[0]) bin=0;
				for(int i=1;i<moments_npt;i++) {if(s12>s12_limits[i-1] && s12<s12_limits[i]) bin=i;}

				// this is the event weight
				w = PScorr_KK[bin]/acc;
				// Increment mass bin population and compute the moments for this event 
				bin_pop_KK[bin] = bin_pop_KK[bin] + 1.;
				t0_KK[bin]= t0_KK[bin] + w*PL[0];
				t1_KK[bin]= t1_KK[bin] + w*PL[1];
				t2_KK[bin]= t2_KK[bin] + w*PL[2];
				t3_KK[bin]= t3_KK[bin] + w*PL[3];
				t4_KK[bin]= t4_KK[bin] + w*PL[4];

				//Integral += t0[bin];

				//if (s12>1.6) cout << "w = " << w << ", t0 = " << t0_KK[bin] << ", Integral = " << Integral_KK << ", bin = " << bin  << ", s12_limits[i-1] = " << s12_limits[bin-1]  << ", s12_limits[i] = " << s12_limits[bin]  << ", s12 = " << s12 << endl;  

				// variables for errors on moments
				t0sq_KK[bin] = t0sq_KK[bin] + w*PL[0]*w*PL[0];
				t1sq_KK[bin] = t1sq_KK[bin] + w*PL[1]*w*PL[1];
				t2sq_KK[bin] = t2sq_KK[bin] + w*PL[2]*w*PL[2];
				t3sq_KK[bin] = t3sq_KK[bin] + w*PL[3]*w*PL[3];
				t4sq_KK[bin] = t4sq_KK[bin] + w*PL[4]*w*PL[4];

				//IntegralSq += t0sq[bin];
			}
		}

		p1=0.5*sqrt( ( (s13 - m3sq - m1sq)*(s13 - m3sq - m1sq) - 4*m3sq*m1sq )/s13 );
		p2=0.5*sqrt( ( (s - s13 - m2sq)*(s - s13 - m2sq) - 4*s13*m2sq )/s13 );
		E1 = (s13 + m1sq - m3sq) / (2*sqrt(s13));
		E2 = (s - s13 - m2sq)/(2*sqrt(s13));
		z = (m1sq + m2sq + 2*E1*E2 - s12)/(2*p1*p2);

		if (fabs(z)<1){

			// Legendre polynomials for this event
			PL[0] = 1.;
			PL[1] = z;
			PL[2] = (3.*z*z-1.)/2.;
			PL[3] = (5.*pow(z,3)-3.*z)/2.;
			PL[4] = (35.*pow(z,4)-30.*z*z+3.)/8.;

			// check if s23 is in the appropriate mass region
			if((s13 > s13min)&&(s13 < s13max)){

				// acceptance for this event       
				if (s13>0.9)	acc = acceptance(s12,s13);

				//cout << " acc = " << acc << endl;

				// find in which mKpi bin this event belongs to;
				if(s13 < s13_limits[0]) bin=0;
				for(int i=1;i<moments_npt;i++) {if(s13>s13_limits[i-1] && s13<s13_limits[i]) bin=i;}

				// this is the event weight
				w = PScorr_Kpi[bin]/acc;
				// Increment mass bin population and compute the moments for this event 
				bin_pop_Kpi[bin] = bin_pop_Kpi[bin] + 1.;
				t0_Kpi[bin]= t0_Kpi[bin] + w*PL[0];
				t1_Kpi[bin]= t1_Kpi[bin] + w*PL[1];
				t2_Kpi[bin]= t2_Kpi[bin] + w*PL[2];
				t3_Kpi[bin]= t3_Kpi[bin] + w*PL[3];
				t4_Kpi[bin]= t4_Kpi[bin] + w*PL[4];

				//Integral += t0[bin];

				//if (s13>0.9) cout << "w = " << w << ", t0 = " << t0_Kpi[bin] << ", Integral = " << Integral_Kpi << ", bin = " << bin << ", s13_limits[i-1] = " << s13_limits[bin-1]  << ", s13_limits[i] = " << s13_limits[bin]  << ", s13 = " << s13 << endl;  

				// variables for errors on moments
				t0sq_Kpi[bin] = t0sq_Kpi[bin] + w*PL[0]*w*PL[0];
				t1sq_Kpi[bin] = t1sq_Kpi[bin] + w*PL[1]*w*PL[1];
				t2sq_Kpi[bin] = t2sq_Kpi[bin] + w*PL[2]*w*PL[2];
				t3sq_Kpi[bin] = t3sq_Kpi[bin] + w*PL[3]*w*PL[3];
				t4sq_Kpi[bin] = t4sq_Kpi[bin] + w*PL[4]*w*PL[4];

				//IntegralSq += t0sq[bin];
			}
		}
	}

	for (int i=0; i<moments_npt; i++) {
		Integral_KK += t0_KK[i];
		IntegralSq_KK += t0sq_KK[i];
	}

	for (int i=0; i<moments_npt; i++) {

//		              cout << "t0_KK[" << i <<"] = " << t0_KK[i] << ", t1_KK[" << i <<"] = " << t1_KK[i] << ", t2_KK[" << i <<"] = " << t2_KK[i] << ", t3_KK[" << i <<"] = " << t3_KK[i] << ", t4_KK[" << i <<"] = " << t4_KK[i] << ", Integral = " << Integral_KK << "bin pop = "<< bin_pop_KK[i] << endl; 

		t0_KK[i] = t0_KK[i]/Integral_KK;
		t1_KK[i] = t1_KK[i]/Integral_KK;
		t2_KK[i] = t2_KK[i]/Integral_KK;
		t3_KK[i] = t3_KK[i]/Integral_KK;
		t4_KK[i] = t4_KK[i]/Integral_KK;

		t0sq_KK[i] = t0sq_KK[i]/IntegralSq_KK;
		t1sq_KK[i] = t1sq_KK[i]/IntegralSq_KK;
		t2sq_KK[i] = t2sq_KK[i]/IntegralSq_KK;
		t3sq_KK[i] = t3sq_KK[i]/IntegralSq_KK;
		t4sq_KK[i] = t4sq_KK[i]/IntegralSq_KK;

		av_w = t0_KK[i]/bin_pop_KK[i];
		av_w2 = t0sq_KK[i]/bin_pop_KK[i];
		sigma_w = av_w2 - av_w*av_w;
		t0_err_KK[i] = fabs(t0_KK[i])/sqrt(bin_pop_KK[i])*sqrt(1.+(sigma_w/av_w)*(sigma_w/av_w));

		av_w = t1_KK[i]/bin_pop_KK[i];
		av_w2 = t1sq_KK[i]/bin_pop_KK[i];
		sigma_w = av_w2 - av_w*av_w;
		t1_err_KK[i] = fabs(t1_KK[i])/sqrt(bin_pop_KK[i])*sqrt(1.+(sigma_w/av_w)*(sigma_w/av_w));

		av_w = t2_KK[i]/bin_pop_KK[i];
		av_w2 = t2sq_KK[i]/bin_pop_KK[i];
		sigma_w = av_w2 - av_w*av_w;
		t2_err_KK[i] = fabs(t2_KK[i])/sqrt(bin_pop_KK[i])*sqrt(1.+(sigma_w/av_w)*(sigma_w/av_w));

		av_w = t3_KK[i]/bin_pop_KK[i];
		av_w2 = t3sq_KK[i]/bin_pop_KK[i];
		sigma_w = av_w2 - av_w*av_w;
		t3_err_KK[i] = fabs(t3_KK[i])/sqrt(bin_pop_KK[i])*sqrt(1.+(sigma_w/av_w)*(sigma_w/av_w));

		av_w = t4_KK[i]/bin_pop_KK[i];
		av_w2 = t4sq_KK[i]/bin_pop_KK[i];
		sigma_w = av_w2 - av_w*av_w;
		t4_err_KK[i] = fabs(t4_KK[i])/sqrt(bin_pop_KK[i])*sqrt(1.+(sigma_w/av_w)*(sigma_w/av_w));

//		cout << "t0_err_KK = " << t0_err_KK[i] << "t1_err_KK = " << t1_err_KK[i] << "t2_err_KK = " << t2_err_KK[i] << "t3_err_KK = " << t3_err_KK[i] << "t4_err_KK = " << t4_err_KK[i] << ", bin pop = "<< bin_pop_KK[i] << endl;
	}

	plot_moments(m_KK, bin_pop_KK, t0_KK, t1_KK, t2_KK, t3_KK, t4_KK, t0_err_KK, t1_err_KK, t2_err_KK, t3_err_KK, t4_err_KK, g0_KK_Data, g1_KK_Data, g2_KK_Data, g3_KK_Data, g4_KK_Data);

	for (int i=0; i<moments_npt; i++) {
		Integral_Kpi += t0_Kpi[i];
		IntegralSq_Kpi += t0sq_Kpi[i];
	}

	for (int i=0; i<moments_npt; i++) {

//		              cout << "t0_Kpi[" << i <<"] = " << t0_Kpi[i] << ", t1[" << i <<"] = " << t1_Kpi[i] << ", t2_Kpi[" << i <<"] = " << t2_Kpi[i] << ", t3_Kpi[" << i <<"] = " << t3_Kpi[i] << ", t4_Kpi[" << i <<"] = " << t4_Kpi[i] << ", Integral = " << Integral_Kpi << "bin pop = "<< bin_pop_Kpi[i] << endl; 

		t0_Kpi[i] = t0_Kpi[i]/Integral_Kpi;
		t1_Kpi[i] = t1_Kpi[i]/Integral_Kpi;
		t2_Kpi[i] = t2_Kpi[i]/Integral_Kpi;
		t3_Kpi[i] = t3_Kpi[i]/Integral_Kpi;
		t4_Kpi[i] = t4_Kpi[i]/Integral_Kpi;

		t0sq_Kpi[i] = t0sq_Kpi[i]/IntegralSq_Kpi;
		t1sq_Kpi[i] = t1sq_Kpi[i]/IntegralSq_Kpi;
		t2sq_Kpi[i] = t2sq_Kpi[i]/IntegralSq_Kpi;
		t3sq_Kpi[i] = t3sq_Kpi[i]/IntegralSq_Kpi;
		t4sq_Kpi[i] = t4sq_Kpi[i]/IntegralSq_Kpi;

		av_w = t0_Kpi[i]/bin_pop_Kpi[i];
		av_w2 = t0sq_Kpi[i]/bin_pop_Kpi[i];
		sigma_w = av_w2 - av_w*av_w;
		t0_err_Kpi[i] = fabs(t0_Kpi[i])/sqrt(bin_pop_Kpi[i])*sqrt(1.+(sigma_w/av_w)*(sigma_w/av_w));

		av_w = t1_Kpi[i]/bin_pop_Kpi[i];
		av_w2 = t1sq_Kpi[i]/bin_pop_Kpi[i];
		sigma_w = av_w2 - av_w*av_w;
		t1_err_Kpi[i] = fabs(t1_Kpi[i])/sqrt(bin_pop_Kpi[i])*sqrt(1.+(sigma_w/av_w)*(sigma_w/av_w));

		av_w = t2_Kpi[i]/bin_pop_Kpi[i];
		av_w2 = t2sq_Kpi[i]/bin_pop_Kpi[i];
		sigma_w = av_w2 - av_w*av_w;
		t2_err_Kpi[i] = fabs(t2_Kpi[i])/sqrt(bin_pop_Kpi[i])*sqrt(1.+(sigma_w/av_w)*(sigma_w/av_w));

		av_w = t3_Kpi[i]/bin_pop_Kpi[i];
		av_w2 = t3sq_Kpi[i]/bin_pop_Kpi[i];
		sigma_w = av_w2 - av_w*av_w;
		t3_err_Kpi[i] = fabs(t3_Kpi[i])/sqrt(bin_pop_Kpi[i])*sqrt(1.+(sigma_w/av_w)*(sigma_w/av_w));

		av_w = t4_Kpi[i]/bin_pop_Kpi[i];
		av_w2 = t4sq_Kpi[i]/bin_pop_Kpi[i];
		sigma_w = av_w2 - av_w*av_w;
		t4_err_Kpi[i] = fabs(t4_Kpi[i])/sqrt(bin_pop_Kpi[i])*sqrt(1.+(sigma_w/av_w)*(sigma_w/av_w));

//		cout << "t0_err_Kpi = " << t0_err_Kpi[i] << "t1_err_Kpi = " << t1_err_Kpi[i] << "t2_err_Kpi = " << t2_err_Kpi[i] << "t3_err_Kpi = " << t3_err_Kpi[i] << "t4_err_Kpi = " << t4_err_Kpi[i] << "bin pop = "<< bin_pop_KK[i] << endl;
	}

	plot_moments(m_Kpi, bin_pop_Kpi, t0_Kpi, t1_Kpi, t2_Kpi, t3_Kpi, t4_Kpi, t0_err_Kpi, t1_err_Kpi, t2_err_Kpi, t3_err_Kpi, t4_err_Kpi, g0_Kpi_Data, g1_Kpi_Data, g2_Kpi_Data, g3_Kpi_Data, g4_Kpi_Data);

	for(int i=0;i<moments_npt;i++)  {

		// initialize moments
		t0_KK[i]=0;
		t1_KK[i]=0;
		t2_KK[i]=0;
		t3_KK[i]=0;
		t4_KK[i]=0;
		bin_pop_KK[i] = 0;

		// initialize  variables for errors on moments     
		t0sq_KK[i] = 0;
		t1sq_KK[i] = 0;
		t2sq_KK[i] = 0;
		t3sq_KK[i] = 0;
		t4sq_KK[i] = 0;

		// initialize moments
		t0_Kpi[i]=0;
		t1_Kpi[i]=0;
		t2_Kpi[i]=0;
		t3_Kpi[i]=0;
		t4_Kpi[i]=0;
		bin_pop_Kpi[i] = 0;

		// initialize  variables for errors on moments     
		t0sq_Kpi[i] = 0;
		t1sq_Kpi[i] = 0;
		t2sq_Kpi[i] = 0;
		t3sq_Kpi[i] = 0;
		t4sq_Kpi[i] = 0;

	}

	Integral_KK = 0;
	IntegralSq_KK = 0;
	Integral_Kpi = 0;
	IntegralSq_Kpi = 0;

	TFile f( "/Users/dvieira/Dropbox/toyMCGenerator/ntuple_flat_kkpi_100M_0.root" );
	TTree *t1 = (TTree*)f.Get( "DecayTree" );

	t1->SetBranchAddress( "s12",    &s12 );
	t1->SetBranchAddress( "s13",    &s13 );
	t1->SetBranchAddress( "s23",    &s23 );
	nevts = t1->GetEntries();
	//nevts = 10000;

	double totalPdf = 0;
	reading_ntuple = 1;
	for ( Long64_t jentry = 0; jentry < nevts; ++jentry ) {
		t1->GetEntry( jentry );

		//			Slo.push_back(s_low);
		//			S13.push_back(s_high);

		//cout << " s12 = " << s12 << ", s13 = " << s13 << ", s23 = " << s23 << endl;

		res_amp.Amplitudes( final_state, 
				M , 
				s12, 
				s13, 
				s23, 
				sig_amp,
				bkg_amp, 
				resonances, 
				bkg_components, 
				res_masses, 
				res_widths,
				res_extra_pars,
				KK_bin_limit, 
				pwa_coefs, 
				pwa_coefs_prime );

		//		SigAmps.push_back( sig_amp );
		//		BkgAmps.push_back( bkg_amp );

		L_s = TotPdf.Eval( bkg_fraction, 
				sig_amp, 
				coefs_product, 
				s12,
				s13, 
				NL_s, 
				bkg_amp, 
				par, 
				normalization_bkg,
				number_of_resonances, 
				number_of_bkg_components );

		if ( L_s <= 0 )
			continue;


		double accTerm = ( UseAcceptance ? acceptance( s12, s13 ) : 1 );

		//		cout << "s12 = " << s12 << ", s13 = " << s13 << ", L_s = " << L_s << ", accTerm = " << accTerm << endl;   

		p1=0.5*sqrt( ( (s12 - m1sq - m2sq)*(s12 - m1sq - m2sq) - 4*m1sq*m2sq )/s12 );
		p3=0.5*sqrt( ( (s - s12 - m3sq)*(s - s12 - m3sq) - 4*s12*m3sq )/s12 );
		E1 = (s12 + m1sq - m2sq) / (2*sqrt(s12));
		E3 = (s - s12 - m3sq)/(2*sqrt(s12));
		z = (m1sq + m3sq + 2*E1*E3 - s13)/(2*p1*p3);

		//if(fabs(z)>1 && (s12 < 1.4)) cout << "z = " << z << ", s12 = " << s12 << ", s13 = " << s13 << endl;

		if (fabs(z)<1){

			// Legendre polynomials for this event
			PL[0] = 1.;
			PL[1] = z;
			PL[2] = (3.*z*z-1.)/2.;
			PL[3] = (5.*pow(z,3)-3.*z)/2.;
			PL[4] = (35.*pow(z,4)-30.*z*z+3.)/8.;

			// check if s23 is in the appropriate mass region
			if((s12 > s12min)&&(s12 < s12max)){

				// acceptance for this event       
				acc = acceptance(s12,s13);

				//cout << " acc = " << acc << endl;

				// find in which mKK bin this event belongs to;
				if(s12 < s12_limits[0]) bin=0;
				for(int i=1;i<moments_npt;i++) {if(s12>s12_limits[i-1] && s12<s12_limits[i]) bin=i;}

				// this is the event weight
				w = PScorr_KK[bin]*L_s/acc;
				// Increment mass bin population and compute the moments for this event 
				bin_pop_KK[bin] = bin_pop_KK[bin] + 1.;
				t0_KK[bin]= t0_KK[bin] + w*PL[0];
				t1_KK[bin]= t1_KK[bin] + w*PL[1];
				t2_KK[bin]= t2_KK[bin] + w*PL[2];
				t3_KK[bin]= t3_KK[bin] + w*PL[3];
				t4_KK[bin]= t4_KK[bin] + w*PL[4];

				//Integral += t0[bin];

				//cout << "w = " << w << ", t0 = " << t0[bin] << ", Integral = " << Integral << endl;  

				// variables for errors on moments
				t0sq_KK[bin] = t0sq_KK[bin] + w*PL[0]*w*PL[0];
				t1sq_KK[bin] = t1sq_KK[bin] + w*PL[1]*w*PL[1];
				t2sq_KK[bin] = t2sq_KK[bin] + w*PL[2]*w*PL[2];
				t3sq_KK[bin] = t3sq_KK[bin] + w*PL[3]*w*PL[3];
				t4sq_KK[bin] = t4sq_KK[bin] + w*PL[4]*w*PL[4];

				//IntegralSq += t0sq[bin];
			}
		}

		p1=0.5*sqrt( ( (s13 - m3sq - m1sq)*(s13 - m3sq - m1sq) - 4*m3sq*m1sq )/s13 );
		p2=0.5*sqrt( ( (s - s13 - m2sq)*(s - s13 - m2sq) - 4*s13*m2sq )/s13 );
		E1 = (s13 + m1sq - m3sq) / (2*sqrt(s13));
		E2 = (s - s13 - m2sq)/(2*sqrt(s13));
		z = (m1sq + m2sq + 2*E1*E2 - s12)/(2*p1*p2);

		if (fabs(z)<1){

			// Legendre polynomials for this event
			PL[0] = 1.;
			PL[1] = z;
			PL[2] = (3.*z*z-1.)/2.;
			PL[3] = (5.*pow(z,3)-3.*z)/2.;
			PL[4] = (35.*pow(z,4)-30.*z*z+3.)/8.;

			// check if s23 is in the appropriate mass region
			if((s13 > s13min)&&(s13 < s13max)){

				// acceptance for this event       
				acc = acceptance(s12,s13);

				//cout << " acc = " << acc << endl;

				// find in which mKpi bin this event belongs to;
				if(s13 < s13_limits[0]) bin=0;
				for(int i=1;i<moments_npt;i++) {if(s13>s13_limits[i-1] && s13<s13_limits[i]) bin=i;}

				// this is the event weight
				w = PScorr_Kpi[bin]*L_s/acc;
				// Increment mass bin population and compute the moments for this event 
				bin_pop_Kpi[bin] = bin_pop_Kpi[bin] + 1.;
				t0_Kpi[bin]= t0_Kpi[bin] + w*PL[0];
				t1_Kpi[bin]= t1_Kpi[bin] + w*PL[1];
				t2_Kpi[bin]= t2_Kpi[bin] + w*PL[2];
				t3_Kpi[bin]= t3_Kpi[bin] + w*PL[3];
				t4_Kpi[bin]= t4_Kpi[bin] + w*PL[4];

				//Integral += t0[bin];

				//cout << "w = " << w << ", t0 = " << t0[bin] << ", Integral = " << Integral << endl;  

				// variables for errors on moments
				t0sq_Kpi[bin] = t0sq_Kpi[bin] + w*PL[0]*w*PL[0];
				t1sq_Kpi[bin] = t1sq_Kpi[bin] + w*PL[1]*w*PL[1];
				t2sq_Kpi[bin] = t2sq_Kpi[bin] + w*PL[2]*w*PL[2];
				t3sq_Kpi[bin] = t3sq_Kpi[bin] + w*PL[3]*w*PL[3];
				t4sq_Kpi[bin] = t4sq_Kpi[bin] + w*PL[4]*w*PL[4];

				//IntegralSq += t0sq[bin];
			}
		}

		for ( int xx = 0; xx < number_of_resonances; ++xx ) {

			spd_prod = coefs_product[ xx ][ xx ].Rho() * sig_amp[ xx ].Rho2();

			histo_sij[ 0 ][ xx ]->Fill( s12,
					spd_prod * ( 1 - bkg_fraction ) * accTerm / NL_s );

			histo_sij[ 1 ][ xx ]->Fill( s13,
					spd_prod * ( 1 - bkg_fraction ) * accTerm / NL_s ); 

			histo_sij[ 2 ][ xx ]->Fill( s23,
					spd_prod * ( 1 - bkg_fraction ) * accTerm / NL_s );                            
		}


		for ( int i = 0; i < number_of_bkg_components; ++i ) {

			bhisto_sij[ 0 ][ i ]->Fill( s12, 
					( par[ ( i + 8 * number_of_resonances ) ] *
					  bkg_amp[ i ].Re() /
					  normalization_bkg.at( i ) ) *
					bkg_fraction );

			bhisto_sij[ 1 ][ i ]->Fill( s13, 
					( par[ ( i + 8 * number_of_resonances ) ] *
					  bkg_amp[ i ].Re() /
					  normalization_bkg.at(i)) *
					bkg_fraction );

			bhisto_sij[ 2 ][ i ]->Fill( s23, 
					( par[ ( i + 8 * number_of_resonances ) ] *
					  bkg_amp[ i ].Re() /
					  normalization_bkg.at(i)) *
					bkg_fraction );                                            

		}

		TotalPDFHist->Fill( s13,  s12, L_s );
		TotalToyHist->Fill( s13,  s12);
		TotalPDFHist_Chi2->Fill( s13,  s12, L_s );
		TotalToyHist_Chi2->Fill( s13,  s12);
		S12Hist    ->Fill( s12,  L_s );
		S12HistToy ->Fill( s12 );
		S13Hist     ->Fill( s13, L_s );
		S13HistToy  ->Fill( s13 );
		S23Hist     ->Fill( s23,    L_s );
		S23HistToy  ->Fill( s23 );
		totalPdf += L_s;

	}

	for (int i=0; i<moments_npt; i++) {
		Integral_KK += t0_KK[i];
		IntegralSq_KK += t0sq_KK[i];
	}

	for (int i=0; i<moments_npt; i++) {

		//              cout << "t0[" << i <<"] = " << t0[i] << ", t1[" << i <<"] = " << t1[i] << ", t2[" << i <<"] = " << t2[i] << ", t3[" << i <<"] = " << t3[i] << ", t4[" << i <<"] = " << t4[i] << ", Integral = " << Integral << "bin pop = "<< bin_pop[i] << endl; 

		t0_KK[i] = t0_KK[i]/Integral_KK;
		t1_KK[i] = t1_KK[i]/Integral_KK;
		t2_KK[i] = t2_KK[i]/Integral_KK;
		t3_KK[i] = t3_KK[i]/Integral_KK;
		t4_KK[i] = t4_KK[i]/Integral_KK;

		t0sq_KK[i] = t0sq_KK[i]/IntegralSq_KK;
		t1sq_KK[i] = t1sq_KK[i]/IntegralSq_KK;
		t2sq_KK[i] = t2sq_KK[i]/IntegralSq_KK;
		t3sq_KK[i] = t3sq_KK[i]/IntegralSq_KK;
		t4sq_KK[i] = t4sq_KK[i]/IntegralSq_KK;

		av_w = t0_KK[i]/bin_pop_KK[i];
		av_w2 = t0sq_KK[i]/bin_pop_KK[i];
		sigma_w = av_w2 - av_w*av_w;
		t0_err_KK[i] = fabs(t0_KK[i])/sqrt(bin_pop_KK[i])*sqrt(1.+(sigma_w/av_w)*(sigma_w/av_w));

		av_w = t1_KK[i]/bin_pop_KK[i];
		av_w2 = t1sq_KK[i]/bin_pop_KK[i];
		sigma_w = av_w2 - av_w*av_w;
		t1_err_KK[i] = fabs(t1_KK[i])/sqrt(bin_pop_KK[i])*sqrt(1.+(sigma_w/av_w)*(sigma_w/av_w));

		av_w = t2_KK[i]/bin_pop_KK[i];
		av_w2 = t2sq_KK[i]/bin_pop_KK[i];
		sigma_w = av_w2 - av_w*av_w;
		t2_err_KK[i] = fabs(t2_KK[i])/sqrt(bin_pop_KK[i])*sqrt(1.+(sigma_w/av_w)*(sigma_w/av_w));

		av_w = t3_KK[i]/bin_pop_KK[i];
		av_w2 = t3sq_KK[i]/bin_pop_KK[i];
		sigma_w = av_w2 - av_w*av_w;
		t3_err_KK[i] = fabs(t3_KK[i])/sqrt(bin_pop_KK[i])*sqrt(1.+(sigma_w/av_w)*(sigma_w/av_w));

		av_w = t4_KK[i]/bin_pop_KK[i];
		av_w2 = t4sq_KK[i]/bin_pop_KK[i];
		sigma_w = av_w2 - av_w*av_w;
		t4_err_KK[i] = fabs(t4_KK[i])/sqrt(bin_pop_KK[i])*sqrt(1.+(sigma_w/av_w)*(sigma_w/av_w));

//		cout << "t0_err_KK = " << t0_err_KK[i] << "t1_err_KK = " << t1_err_KK[i] << "t2_err_KK = " << t2_err_KK[i] << "t3_err_KK = " << t3_err_KK[i] << "t4_err_KK = " << t4_err_KK[i] << endl;
	}

	plot_moments(m_KK, bin_pop_KK, t0_KK, t1_KK, t2_KK, t3_KK, t4_KK, t0_err_KK, t1_err_KK, t2_err_KK, t3_err_KK, t4_err_KK, g0_KK_MC, g1_KK_MC, g2_KK_MC, g3_KK_MC, g4_KK_MC);

	for (int i=0; i<moments_npt; i++) {
		Integral_Kpi += t0_Kpi[i];
		IntegralSq_Kpi += t0sq_Kpi[i];
	}

	for (int i=0; i<moments_npt; i++) {

		//              cout << "t0[" << i <<"] = " << t0[i] << ", t1[" << i <<"] = " << t1[i] << ", t2[" << i <<"] = " << t2[i] << ", t3[" << i <<"] = " << t3[i] << ", t4[" << i <<"] = " << t4[i] << ", Integral = " << Integral << "bin pop = "<< bin_pop[i] << endl; 

		t0_Kpi[i] = t0_Kpi[i]/Integral_Kpi;
		t1_Kpi[i] = t1_Kpi[i]/Integral_Kpi;
		t2_Kpi[i] = t2_Kpi[i]/Integral_Kpi;
		t3_Kpi[i] = t3_Kpi[i]/Integral_Kpi;
		t4_Kpi[i] = t4_Kpi[i]/Integral_Kpi;

		t0sq_Kpi[i] = t0sq_Kpi[i]/IntegralSq_Kpi;
		t1sq_Kpi[i] = t1sq_Kpi[i]/IntegralSq_Kpi;
		t2sq_Kpi[i] = t2sq_Kpi[i]/IntegralSq_Kpi;
		t3sq_Kpi[i] = t3sq_Kpi[i]/IntegralSq_Kpi;
		t4sq_Kpi[i] = t4sq_Kpi[i]/IntegralSq_Kpi;

		av_w = t0_Kpi[i]/bin_pop_Kpi[i];
		av_w2 = t0sq_Kpi[i]/bin_pop_Kpi[i];
		sigma_w = av_w2 - av_w*av_w;
		t0_err_Kpi[i] = fabs(t0_Kpi[i])/sqrt(bin_pop_Kpi[i])*sqrt(1.+(sigma_w/av_w)*(sigma_w/av_w));

		av_w = t1_Kpi[i]/bin_pop_Kpi[i];
		av_w2 = t1sq_Kpi[i]/bin_pop_Kpi[i];
		sigma_w = av_w2 - av_w*av_w;
		t1_err_Kpi[i] = fabs(t1_Kpi[i])/sqrt(bin_pop_Kpi[i])*sqrt(1.+(sigma_w/av_w)*(sigma_w/av_w));

		av_w = t2_Kpi[i]/bin_pop_Kpi[i];
		av_w2 = t2sq_Kpi[i]/bin_pop_Kpi[i];
		sigma_w = av_w2 - av_w*av_w;
		t2_err_Kpi[i] = fabs(t2_Kpi[i])/sqrt(bin_pop_Kpi[i])*sqrt(1.+(sigma_w/av_w)*(sigma_w/av_w));

		av_w = t3_Kpi[i]/bin_pop_Kpi[i];
		av_w2 = t3sq_Kpi[i]/bin_pop_Kpi[i];
		sigma_w = av_w2 - av_w*av_w;
		t3_err_Kpi[i] = fabs(t3_Kpi[i])/sqrt(bin_pop_Kpi[i])*sqrt(1.+(sigma_w/av_w)*(sigma_w/av_w));

		av_w = t4_Kpi[i]/bin_pop_Kpi[i];
		av_w2 = t4sq_Kpi[i]/bin_pop_Kpi[i];
		sigma_w = av_w2 - av_w*av_w;
		t4_err_Kpi[i] = fabs(t4_Kpi[i])/sqrt(bin_pop_Kpi[i])*sqrt(1.+(sigma_w/av_w)*(sigma_w/av_w));

		//cout << "t0_err_Kpi = " << t0_err_Kpi[i] << "t1_err_Kpi = " << t1_err_Kpi[i] << "t2_err_Kpi = " << t2_err_Kpi[i] << "t3_err_Kpi = " << t3_err_Kpi[i] << "t4_err_Kpi = " << t4_err_Kpi[i] << endl;
	}

	plot_moments(m_Kpi, bin_pop_Kpi, t0_Kpi, t1_Kpi, t2_Kpi, t3_Kpi, t4_Kpi, t0_err_Kpi, t1_err_Kpi, t2_err_Kpi, t3_err_Kpi, t4_err_Kpi, g0_Kpi_MC, g1_Kpi_MC, g2_Kpi_MC, g3_Kpi_MC, g4_Kpi_MC);

	double scale = double(nentries) / totalPdf;

	TotalPDFHist->Sumw2();
	TotalToyHist->Sumw2();
	S12HistToy->Sumw2();
	S13HistToy->Sumw2();
	S23HistToy->Sumw2();

	S12Hist->Scale( scale );
	TotalPDFHist->Scale( scale );
	TotalPDFHist_Chi2->Scale( scale );
	S13Hist->Scale( scale );
	S23Hist->Scale( scale );

	clock.Stop();
	cout << "Finished filling vectors and histograms" << endl;
	clock.Print("u");

	set_plot_style_Dalitz();

	TCanvas * c0 = new TCanvas(); 
	TotalPDFHist->Draw("colz");
	c0->Print("TotalPDFHist.pdf");

	TCanvas * c1 = new TCanvas(); 
	TotalPDFHist_Chi2->Draw("colz");
	c1->Print("TotalPDFHist_Chi2.pdf"); 

	for ( int xxx = 0; xxx < number_of_resonances; ++xxx ) {
		for ( int yy = 0; yy < 3; ++yy ) {

			histo_sij[ yy ][ xxx ]->SetLineWidth( 2 );
			if ( xxx == 0 )
				histo_sij[ yy ][ xxx ]->SetLineColor( 46 );
			if ( xxx == 1 )
				histo_sij[ yy ][ xxx ]->SetLineColor( 1 );
			if ( xxx == 2 )
				histo_sij[ yy ][ xxx ]->SetLineColor( 2 );


			if ( xxx == 3 )
				histo_sij[ yy ][ xxx ]->SetLineColor( 3 );

			if ( xxx == 4 )
				histo_sij[ yy ][ xxx ]->SetLineColor( kBlue - 3 );
			if ( xxx == 5 )
				histo_sij[ yy ][ xxx ]->SetLineColor( 5 );
			if ( xxx == 6 )
				histo_sij[ yy ][ xxx ]->SetLineColor( 6 );
			if ( xxx == 7 )
				histo_sij[ yy ][ xxx ]->SetLineColor( 7 );
			if ( xxx == 8 )
				histo_sij[ yy ][ xxx ]->SetLineColor( 8 );
			if ( xxx == 9 )
				histo_sij[ yy ][ xxx ]->SetLineColor( 9 );
			if ( xxx == 10 )
				histo_sij[ yy ][ xxx ]->SetLineColor( 30 );
			if ( xxx == 11 )
				histo_sij[ yy ][ xxx ]->SetLineColor( 31 );
			if ( xxx == 12 )
				histo_sij[ yy ][ xxx ]->SetLineColor( 12 );
			if ( xxx == 13 )
				histo_sij[ yy ][ xxx ]->SetLineColor( 28 );
			if ( xxx == 14 )
				histo_sij[ yy ][ xxx ]->SetLineColor( 39 );
			if ( xxx == 15 )
				histo_sij[ yy ][ xxx ]->SetLineColor( 45 );
			if ( xxx == 16 )
				histo_sij[ yy ][ xxx ]->SetLineColor( 26 );


			histo_sij[ yy ][ xxx ]->Scale( scale );

		}

		//cout << " sij["<<0<<"]["<<xxx<<"] = "<< histo_sij[ 0 ][ xxx ]->Integral() << endl;
	}
	for ( int bb = 0; bb < number_of_bkg_components; ++bb ) {
		for ( int yy = 0; yy < 3; ++yy ) {

			bhisto_sij[ yy ][ bb ]->SetLineWidth( 2 );
			bhisto_sij[ yy ][ bb ]->SetLineColor( kRed );
			bhisto_sij[ yy ][ bb ]->SetLineStyle( 3 );

			bhisto_sij[ yy ][ bb ]->Scale(scale);
		}
	}


	string output_histo_file_name = input_txt_file_name.substr(0,input_txt_file_name.size() - 4);
	TString PDFName = "Results_PDF_fit/" + output_histo_file_name + ".pdf"; 
	TString MomentsPDFName = "Results_PDF_fit/" + output_histo_file_name; 
	output_histo_file_name += ".root";


	S12Hist->SetXTitle( colNameLatexY1 );
	S12Hist->SetYTitle( "Entries" );
	S13Hist->SetXTitle ( colNameLatexX1 );
	S13Hist->SetYTitle ( "Entries" );
	S23Hist->SetXTitle ( colNameLatexY2 );
	S23Hist->SetYTitle ( "Entries" );


	TotalDataHist->SetXTitle( colNameLatexX1 );
	TotalDataHist->SetYTitle( colNameLatexY1 );
	TotalDataHist_Chi2->SetXTitle( colNameLatexX1 );
	TotalDataHist_Chi2->SetYTitle( colNameLatexY1 );
	TotalPDFHist->SetXTitle ( colNameLatexX1 );
	TotalPDFHist->SetYTitle ( colNameLatexY1 );
	TotalPDFHist_Chi2->SetXTitle ( colNameLatexX1 );
	TotalPDFHist_Chi2->SetYTitle ( colNameLatexY1 );


	TotalDataHist->GetYaxis()->SetTitleOffset( 1.1 );
	TotalDataHist->GetXaxis()->SetTitleOffset( 0.9 );
	TotalDataHist_Chi2->GetYaxis()->SetTitleOffset( 1.1 );
	TotalDataHist_Chi2->GetXaxis()->SetTitleOffset( 0.9 );
	TotalPDFHist->GetYaxis()->SetTitleOffset( 1.1 );
	TotalPDFHist->GetXaxis()->SetTitleOffset( 0.9 );
	TotalPDFHist_Chi2->GetYaxis()->SetTitleOffset( 1.1 );
	TotalPDFHist_Chi2->GetXaxis()->SetTitleOffset( 0.9 );


	S12Hist->GetYaxis()->SetTitleOffset( 0.9 );
	S12Hist->GetXaxis()->SetTitleOffset( 0.9 );
	S13Hist->GetYaxis()->SetTitleOffset( 0.9);
	S13Hist->GetXaxis()->SetTitleOffset( 0.9 );
	S23Hist->GetYaxis()->SetTitleOffset( 0.9 );
	S23Hist->GetXaxis()->SetTitleOffset( 0.9 );


	S12Hist->SetLineColor( kBlue );
	S12Hist->SetLineWidth( 2 );
	S13Hist->SetLineColor( kBlue );
	S13Hist->SetLineWidth( 2 );
	S23Hist->SetLineColor( kBlue );
	S23Hist->SetLineWidth( 2 );

	S13HistData->SetLineWidth(1);
	S12HistData->SetLineWidth(1);
	S23HistData->SetLineWidth(1);


	TCanvas *AccCanv = new TCanvas( "AccCanv", " ", 0, 0, 1000, 500 );
	AccCanv->SetRightMargin(0.12);
	AccCanv->SetLeftMargin(0.2);
	AccCanv->SetBottomMargin(0.1);
	AccCanv->Divide( 2, 1);
	AccCanv->Print( PDFName + "[" );
	AccCanv->Update();
	AccCanv->cd( 1 );

	TotalPDFHist->Draw( "colz" );
	AccCanv->cd( 2 );
	TotalDataHist->Draw( "colz" );
	AccCanv->Print( PDFName );
	AccCanv->Update();

	AccCanv->Clear();

	AccCanv->Divide( 2, 1 );

	AccCanv->cd( 1 );
	TotalPDFHist_Chi2->Draw( "colz" );
	AccCanv->cd( 2 );
	TotalDataHist_Chi2->Draw( "colz" );
	AccCanv->Print( PDFName );
	AccCanv->Update();

	AccCanv->Clear();



	TPad *pad1 = new TPad("pad1", "",0.03,0.19,0.5,0.92);
	TPad *pad2 = new TPad("pad2", "",0.03,0.02,0.5,0.18);

	TPad *pad3 = new TPad("pad3", "",0.51,0.20,0.97,0.92);
	TPad *pad4 = new TPad("pad4", "",0.51,0.02,0.97,0.18);

	pad1->Draw();
	pad2->SetBorderSize(0);
	pad2->Draw();
	pad3->Draw();
	pad4->Draw();
	pad1->cd();

	AccCanv->cd( 1 );
	gPad->SetLogy();

	S13Hist->Draw( "HIST C " );
	S13HistData->Draw( "E1 SAME" );  

	cout << "Finished setting histograms" << endl;  

	for( int xx = 0; xx < number_of_resonances; ++xx ) 
	{

		histo_sij[ 1 ][ xx ]->Draw( "HIST C SAME" ); 

	}

	for( int bb0 = 0; bb0 < number_of_bkg_components; ++bb0 ) 
	{

		bhisto_sij[ 1 ][ bb0 ]->Draw( "HIST C SAME" ); 

	}

	cout << "Finished Drawing" << endl;  

	TLegend *ls12 = 
		new TLegend( 0.65, 0.55, 0.8, 0.85 ), // ver https://root.cern.ch/root/html/TPave.html#TPave:TPave@1
		    *ls13 = 
			    new TLegend( 0.15, 0.55, 0.3, 0.85),
		    *ls23 = 
			    new TLegend( 0.65, 0.6, 0.8, 0.9);

	ls12->SetBorderSize ( 0 );
	ls12->SetTextSize(0.02);
	ls13->SetBorderSize ( 0 );
	ls13->SetTextSize(0.02);
	ls23->SetTextSize(0.02);
	ls23->SetBorderSize ( 0 );

	for (int xx = 0; xx < AVAILABLE_RESONANCES; ++xx ) {
		if(resonances[xx]){
			//cout << "xx = " << xx <<  ", resonance_number = " << resonance_number << endl;
			//cout << "number_of_resonances = " << number_of_resonances << endl;
			//cout << resonant_channel_string_tex[ xx ].c_str() << endl;

			ls12->AddEntry( histo_sij[ 0 ][ resonance_number ], 
					resonant_channel_string_tex[ xx ].c_str(), 
					"l" );
			ls13->AddEntry( histo_sij[ 1 ][ resonance_number ], 
					resonant_channel_string_tex[ xx ].c_str(), 
					"l" );
			ls23->AddEntry( histo_sij[ 2 ][ resonance_number ], 
					resonant_channel_string_tex[ xx ].c_str(), 
					"l" ); 
			resonance_number++;              
		}
	}
	for (int xx = 0; xx < AVAILABLE_BKG_COMPONENTS; ++xx ) {
		if(bkg_components[xx]){
			//cout << "xx = " << xx <<  ", bkg_component_number = " << bkg_component_number << endl;
			ls12->AddEntry( bhisto_sij[ 0 ][ bkg_component_number ], 
					bkg_component_string[ xx ].c_str(), 
					"l" );
			ls13->AddEntry( bhisto_sij[ 1 ][ bkg_component_number ], 
					bkg_component_string[ xx ].c_str(), 
					"l" );
			ls23->AddEntry( bhisto_sij[ 2 ][ bkg_component_number ], 
					bkg_component_string[ xx ].c_str(), 
					"l" );              
			bkg_component_number++;              
		}
	}


	//ls13->Draw();

	pad2->cd();

	DrawMypull( S13HistData, S13Hist, S13HistToy );

	pad3->cd();

	gPad->SetLogy();
	S12Hist->Draw( "C HIST" );
	S12HistData->Draw( "E same" );

	for( int xx = 0; xx < number_of_resonances; ++xx )

		histo_sij[ 0 ][ xx ]->Draw("csame");

	for( int bb0 = 0; bb0 < number_of_bkg_components; ++bb0 )

		bhisto_sij[ 0 ][ bb0 ]->Draw( "csame" );

	//ls12->Draw();

	pad4->cd();

	DrawMypull( S12HistData, S12Hist, S12HistToy );


	AccCanv->Print( PDFName );
	AccCanv->Update();

	pad1->cd();
	S23Hist->Draw( "C HIST" );
	S23HistData->Draw( "Esame" );

	for( int xx = 0; xx < number_of_resonances; ++xx )

		histo_sij[ 2 ][ xx ]->Draw("csame");

	for( int bb0 = 0; bb0 < number_of_bkg_components; ++bb0 )

		bhisto_sij[ 2 ][ bb0 ]->Draw( "csame" );

	pad2->cd();

	DrawMypull( S23HistData, S23Hist, S23HistToy );

	pad3->cd();
	gPad->SetLogy(0); 

	set_plot_style_chi2();

	double chisq, 
	       kolm,
	       pval;

	int nDof;
	/*
	   TH2D *chi2 = CalculateChiSq( TotalPDFHist, 
	   TotalToyHist, 
	   TotalDataHist, 
	   nparams,
	   nDof, 
	   chisq,
	   kolm );
	 */

	TH2Poly *chi2 = CalculateChiSq_adaptive( TotalPDFHist_Chi2,
			TotalToyHist_Chi2,
			TotalDataHist_Chi2,
			nparams,
			nDof,
			chisq,
			kolm );


	//chi2->GetZaxis()->SetRange(-40,40);
	chi2->SetMaximum(20);
	chi2->SetMinimum(-20);

	pval = TMath::Prob(chisq, nDof);
	chi2->SetXTitle( colNameLatexX1 );
	chi2->SetYTitle( colNameLatexY1 );
	cout<<"Total ChiSq/nDof = "<<chisq<<"/"<<nDof<<" = "<< chisq/nDof <<" Prob = "<<pval<<endl;
	chi2->Draw( "COLZ" );

	char string_chi[ 30 ],
	     string_p[30], 
	     string_kol[ 30 ];

	sprintf( string_chi, "#chi^{2}/ndof =  %0.2f", chisq/nDof );
	sprintf( string_p, "Pval = %0.2f", pval);
	sprintf( string_kol, "Kolm Test =  %0.2f", kolm );

	TLatex latex;
	latex.SetNDC( kTRUE );
	latex.SetTextSize( 0.03 );
	latex.DrawLatex( 0.65, 0.85, string_chi  );
	latex.DrawLatex( 0.64, 0.8, string_p );
	latex.DrawLatex( 0.64, 0.75, string_kol );

	pad4->Clear();
	AccCanv->Print( PDFName );
	AccCanv->Print( PDFName + "]" );

	AccCanv->Close();


	TFile * histo_file = new TFile(output_histo_file_name.c_str(), "RECREATE");
	TotalPDFHist->SetDirectory(histo_file);
	TotalDataHist->SetDirectory(histo_file);
	TotalToyHist->SetDirectory(histo_file);
	TotalPDFHist_Chi2->SetDirectory(histo_file);
	TotalDataHist_Chi2->SetDirectory(histo_file);
	TotalToyHist_Chi2->SetDirectory(histo_file);
	S12Hist     ->SetDirectory(histo_file);
	S12HistData     ->SetDirectory(histo_file);
	S12HistToy  ->SetDirectory(histo_file);
	S13Hist     ->SetDirectory(histo_file);
	S13HistData     ->SetDirectory(histo_file);
	S13HistToy  ->SetDirectory(histo_file);
	S23Hist     ->SetDirectory(histo_file);
	S23HistData     ->SetDirectory(histo_file);
	S23HistToy  ->SetDirectory(histo_file);

	TotalPDFHist->Write();
	TotalDataHist->Write();
	TotalToyHist->Write();
	TotalPDFHist_Chi2->Write();
	TotalDataHist_Chi2->Write();
	TotalToyHist_Chi2->Write();
	S12Hist     ->Write();
	S12HistData     ->Write();
	S12HistToy  ->Write();
	S13Hist     ->Write();
	S13HistData     ->Write();
	S13HistToy  ->Write();
	S23Hist     ->Write();
	S23HistData     ->Write();
	S23HistToy  ->Write();

	for ( int xx = 0; xx < number_of_resonances; ++xx ) {
		histo_sij[ 0 ][ xx ]->Write();
		histo_sij[ 1 ][ xx ]->Write();
		histo_sij[ 2 ][ xx ]->Write();
	} 

	for ( int i = 0; i < number_of_bkg_components; ++i ) {
		bhisto_sij[ 0 ][ i ]->Write();
		bhisto_sij[ 1 ][ i ]->Write();
		bhisto_sij[ 2 ][ i ]->Write();
	}
	histo_file->Close();

	cout << "saved histofile" << ", colName = "  << colNameLatexX1 << ", colName = "  << colNameLatexY1 << ", colName = "  << colNameLatexY2 << endl;  

	gStyle->SetPadRightMargin(0.1);

	//g0_KK_Data->SetLineColor(2);
	g0_KK_Data->SetLineWidth(2);
	g0_KK_Data->SetMarkerColor(1);
	g0_KK_Data->SetMarkerSize(0.5);
	g0_KK_Data->SetMarkerStyle(21);
	g0_KK_Data->GetXaxis()->SetTitle("s_{K#pi} [GeV/c^{2}");
	g0_KK_Data->GetYaxis()->SetTitle("t_{0}^{0} (arbitrary units)");
	g0_KK_Data->GetYaxis()->SetTitleOffset(1.2);
	g0_KK_Data->SetTitle("");

	//g1_KK_Data->SetLineColor(2);
	g1_KK_Data->SetLineWidth(2);
	g1_KK_Data->SetMarkerColor(1);
	g1_KK_Data->SetMarkerSize(0.5);
	g1_KK_Data->SetMarkerStyle(21);
	g1_KK_Data->GetXaxis()->SetTitle("s_{K#pi} [GeV/c^{2}");
	g1_KK_Data->GetYaxis()->SetTitle("t_{1}^{0} (arbitrary units)");
	g1_KK_Data->GetYaxis()->SetTitleOffset(1.2);
	g1_KK_Data->SetTitle("");

	//g2_KK_Data->SetLineColor(2);
	g2_KK_Data->SetLineWidth(2);
	g2_KK_Data->SetMarkerColor(1);
	g2_KK_Data->SetMarkerSize(0.5);
	g2_KK_Data->SetMarkerStyle(21);
	g2_KK_Data->GetXaxis()->SetTitle("s_{K#pi} [GeV/c^{2}");
	g2_KK_Data->GetYaxis()->SetTitle("t_{2}^{0} (arbitrary units)");
	g2_KK_Data->GetYaxis()->SetTitleOffset(1.2);
	g2_KK_Data->SetTitle("");

	//g3_KK_Data->SetLineColor(2);
	g3_KK_Data->SetLineWidth(2);
	g3_KK_Data->SetMarkerColor(1);
	g3_KK_Data->SetMarkerSize(0.5);
	g3_KK_Data->SetMarkerStyle(21);
	g3_KK_Data->GetXaxis()->SetTitle("s_{K#pi} [GeV/c^{2}");
	g3_KK_Data->GetYaxis()->SetTitle("t_{3}^{0} (arbitrary units)");
	g3_KK_Data->GetYaxis()->SetTitleOffset(1.2);
	g3_KK_Data->SetTitle("");

	//g4_KK_Data->SetLineColor(2);
	g4_KK_Data->SetLineWidth(2);
	g4_KK_Data->SetMarkerColor(1);
	g4_KK_Data->SetMarkerSize(0.5);
	g4_KK_Data->SetMarkerStyle(21);
	g4_KK_Data->GetXaxis()->SetTitle("s_{K#pi} [GeV/c^{2}");
	g4_KK_Data->GetYaxis()->SetTitle("t_{4}^{0} (arbitrary units)");
	g4_KK_Data->GetYaxis()->SetTitleOffset(1.2);
	g4_KK_Data->SetTitle("");

	//g0_Kpi_Data->SetLineColor(2);
	g0_Kpi_Data->SetLineWidth(2);
	g0_Kpi_Data->SetMarkerColor(1);
	g0_Kpi_Data->SetMarkerSize(0.5);
	g0_Kpi_Data->SetMarkerStyle(21);
	g0_Kpi_Data->GetXaxis()->SetTitle("s_{K#pi} [GeV/c^{2}");
	g0_Kpi_Data->GetYaxis()->SetTitle("t_{0}^{0} (arbitrary units)");
	g0_Kpi_Data->GetYaxis()->SetTitleOffset(1.2);
	g0_Kpi_Data->SetTitle("");

	//g1_Kpi_Data->SetLineColor(2);
	g1_Kpi_Data->SetLineWidth(2);
	g1_Kpi_Data->SetMarkerColor(1);
	g1_Kpi_Data->SetMarkerSize(0.5);
	g1_Kpi_Data->SetMarkerStyle(21);
	g1_Kpi_Data->GetXaxis()->SetTitle("s_{K#pi} [GeV/c^{2}");
	g1_Kpi_Data->GetYaxis()->SetTitle("t_{1}^{0} (arbitrary units)");
	g1_Kpi_Data->GetYaxis()->SetTitleOffset(1.2);
	g1_Kpi_Data->SetTitle("");

	//g2_Kpi_Data->SetLineColor(2);
	g2_Kpi_Data->SetLineWidth(2);
	g2_Kpi_Data->SetMarkerColor(1);
	g2_Kpi_Data->SetMarkerSize(0.5);
	g2_Kpi_Data->SetMarkerStyle(21);
	g2_Kpi_Data->GetXaxis()->SetTitle("s_{K#pi} [GeV/c^{2}");
	g2_Kpi_Data->GetYaxis()->SetTitle("t_{2}^{0} (arbitrary units)");
	g2_Kpi_Data->GetYaxis()->SetTitleOffset(1.2);
	g2_Kpi_Data->SetTitle("");

	//g3_Kpi_Data->SetLineColor(2);
	g3_Kpi_Data->SetLineWidth(2);
	g3_Kpi_Data->SetMarkerColor(1);
	g3_Kpi_Data->SetMarkerSize(0.5);
	g3_Kpi_Data->SetMarkerStyle(21);
	g3_Kpi_Data->GetXaxis()->SetTitle("s_{K#pi} [GeV/c^{2}");
	g3_Kpi_Data->GetYaxis()->SetTitle("t_{3}^{0} (arbitrary units)");
	g3_Kpi_Data->GetYaxis()->SetTitleOffset(1.2);
	g3_Kpi_Data->SetTitle("");

	//g4_Kpi_Data->SetLineColor(2);
	g4_Kpi_Data->SetLineWidth(2);
	g4_Kpi_Data->SetMarkerColor(1);
	g4_Kpi_Data->SetMarkerSize(0.5);
	g4_Kpi_Data->SetMarkerStyle(21);
	g4_Kpi_Data->GetXaxis()->SetTitle("s_{K#pi} [GeV/c^{2}");
	g4_Kpi_Data->GetYaxis()->SetTitle("t_{4}^{0} (arbitrary units)");
	g4_Kpi_Data->GetYaxis()->SetTitleOffset(1.2);
	g4_Kpi_Data->SetTitle("");

	//g0_KK_MC->SetLineColor(3);
	//g0_KK_MC->SetLineWidth(2);
	//g0_KK_MC->SetMarkerColor(1);
	//g0_KK_MC->SetMarkerSize(1);
	//g0_KK_MC->SetMarkerStyle(2);
	g0_KK_MC->SetFillColor(38);
	g0_KK_MC->GetXaxis()->SetTitle("s_{K#pi} [GeV/c^{2}");
	g0_KK_MC->GetYaxis()->SetTitle("t_{0}^{0} (arbitrary units)");
	g0_KK_MC->GetYaxis()->SetTitleOffset(1.2);
	g0_KK_MC->SetTitle("");

	//g1_KK_MC->SetLineColor(3);
	//g1_KK_MC->SetLineWidth(2);
	//g1_KK_MC->SetMarkerColor(1);
	//g1_KK_MC->SetMarkerSize(1);
	//g1_KK_MC->SetMarkerStyle(2);
	g1_KK_MC->SetFillColor(38);
	g1_KK_MC->GetXaxis()->SetTitle("s_{K#pi} [GeV/c^{2}");
	g1_KK_MC->GetYaxis()->SetTitle("t_{1}^{0} (arbitrary units)");
	g1_KK_MC->GetYaxis()->SetTitleOffset(1.2);
	g1_KK_MC->SetTitle("");

	//g2_KK_MC->SetLineColor(3);
	//g2_KK_MC->SetLineWidth(2);
	//g2_KK_MC->SetMarkerColor(1);
	//g2_KK_MC->SetMarkerSize(1);
	//g2_KK_MC->SetMarkerStyle(2);
	g2_KK_MC->SetFillColor(38);
	g2_KK_MC->GetXaxis()->SetTitle("s_{K#pi} [GeV/c^{2}");
	g2_KK_MC->GetYaxis()->SetTitle("t_{2}^{0} (arbitrary units)");
	g2_KK_MC->GetYaxis()->SetTitleOffset(1.2);
	g2_KK_MC->SetTitle("");

	//g3_KK_MC->SetLineColor(3);
	//g3_KK_MC->SetLineWidth(2);
	//g3_KK_MC->SetMarkerColor(1);
	//g3_KK_MC->SetMarkerSize(1);
	//g3_KK_MC->SetMarkerStyle(2);
	g3_KK_MC->SetFillColor(38);
	g3_KK_MC->GetXaxis()->SetTitle("s_{K#pi} [GeV/c^{2}");
	g3_KK_MC->GetYaxis()->SetTitle("t_{3}^{0} (arbitrary units)");
	g3_KK_MC->GetYaxis()->SetTitleOffset(1.2);
	g3_KK_MC->SetTitle("");

	//g4_KK_MC->SetLineColor(3);
	//g4_KK_MC->SetLineWidth(2);
	//g4_KK_MC->SetMarkerColor(1);
	//g4_KK_MC->SetMarkerSize(1);
	//g4_KK_MC->SetMarkerStyle(2);
	g4_KK_MC->SetFillColor(38);
	g4_KK_MC->GetXaxis()->SetTitle("s_{K#pi} [GeV/c^{2}");
	g4_KK_MC->GetYaxis()->SetTitle("t_{4}^{0} (arbitrary units)");
	g4_KK_MC->GetYaxis()->SetTitleOffset(1.2);
	g4_KK_MC->SetTitle("");

	//g0_Kpi_MC->SetLineColor(3);
	//g0_Kpi_MC->SetLineWidth(2);
	//g0_Kpi_MC->SetMarkerColor(1);
	//g0_Kpi_MC->SetMarkerSize(1);
	//g0_Kpi_MC->SetMarkerStyle(2);
	g0_Kpi_MC->SetFillColor(38);
	g0_Kpi_MC->GetXaxis()->SetTitle("s_{K#pi} [GeV/c^{2}");
	g0_Kpi_MC->GetYaxis()->SetTitle("t_{0}^{0} (arbitrary units)");
	g0_Kpi_MC->GetYaxis()->SetTitleOffset(1.2);
	g0_Kpi_MC->SetTitle("");

	//g1_Kpi_MC->SetLineColor(3);
	//g1_Kpi_MC->SetLineWidth(2);
	//g1_Kpi_MC->SetMarkerColor(1);
	//g1_Kpi_MC->SetMarkerSize(1);
	//g1_Kpi_MC->SetMarkerStyle(2);
	g1_Kpi_MC->SetFillColor(38);
	g1_Kpi_MC->GetXaxis()->SetTitle("s_{K#pi} [GeV/c^{2}");
	g1_Kpi_MC->GetYaxis()->SetTitle("t_{1}^{0} (arbitrary units)");
	g1_Kpi_MC->GetYaxis()->SetTitleOffset(1.2);
	g1_Kpi_MC->SetTitle("");

	//g2_Kpi_MC->SetLineColor(3);
	//g2_Kpi_MC->SetLineWidth(2);
	//g2_Kpi_MC->SetMarkerColor(1);
	//g2_Kpi_MC->SetMarkerSize(1);
	//g2_Kpi_MC->SetMarkerStyle(2);
	g2_Kpi_MC->SetFillColor(38);
	g2_Kpi_MC->GetXaxis()->SetTitle("s_{K#pi} [GeV/c^{2}");
	g2_Kpi_MC->GetYaxis()->SetTitle("t_{2}^{0} (arbitrary units)");
	g2_Kpi_MC->GetYaxis()->SetTitleOffset(1.2);
	g2_Kpi_MC->SetTitle("");

	//g3_Kpi_MC->SetLineColor(3);
	//g3_Kpi_MC->SetLineWidth(2);
	//g3_Kpi_MC->SetMarkerColor(1);
	//g3_Kpi_MC->SetMarkerSize(1);
	//g3_Kpi_MC->SetMarkerStyle(2);
	g3_Kpi_MC->SetFillColor(38);
	g3_Kpi_MC->GetXaxis()->SetTitle("s_{K#pi} [GeV/c^{2}");
	g3_Kpi_MC->GetYaxis()->SetTitle("t_{3}^{0} (arbitrary units)");
	g3_Kpi_MC->GetYaxis()->SetTitleOffset(1.2);
	g3_Kpi_MC->SetTitle("");

	//g4_Kpi_MC->SetLineColor(3);
	//g4_Kpi_MC->SetLineWidth(2);
	//g4_Kpi_MC->SetMarkerColor(1);
	//g4_Kpi_MC->SetMarkerSize(1);
	//g4_Kpi_MC->SetMarkerStyle(2);
	g4_Kpi_MC->SetFillColor(38);
	g4_Kpi_MC->GetXaxis()->SetTitle("s_{K#pi} [GeV/c^{2}");
	g4_Kpi_MC->GetYaxis()->SetTitle("t_{4}^{0} (arbitrary units)");
	g4_Kpi_MC->GetYaxis()->SetTitleOffset(1.2);
	g4_Kpi_MC->SetTitle("");

	mg0_Kpi->Add(g0_Kpi_Data);
	mg0_Kpi->Add(g0_Kpi_MC);

	mg1_Kpi->Add(g1_Kpi_Data);
	mg1_Kpi->Add(g1_Kpi_MC);

	mg2_Kpi->Add(g2_Kpi_Data);
	mg2_Kpi->Add(g2_Kpi_MC);

	mg3_Kpi->Add(g3_Kpi_Data);
	mg3_Kpi->Add(g3_Kpi_MC);

	mg4_Kpi->Add(g4_Kpi_Data);
	mg4_Kpi->Add(g4_Kpi_MC);

	mg0_KK->Add(g0_KK_Data);
	mg0_KK->Add(g0_KK_MC);

	mg1_KK->Add(g1_KK_Data);
	mg1_KK->Add(g1_KK_MC);

	mg2_KK->Add(g2_KK_Data);
	mg2_KK->Add(g2_KK_MC);

	mg3_KK->Add(g3_KK_Data);
	mg3_KK->Add(g3_KK_MC);

	mg4_KK->Add(g4_KK_Data);
	mg4_KK->Add(g4_KK_MC);

	cout << "-----------------------------------------------------------------------------------------------------------------------------------" << endl;
	cout << endl;
	cout << "Drawing" <<  endl;
	cout << endl;

	TCanvas * c0_Kpi = new TCanvas();

	TPad *pad0_Kpi_0 = new TPad("padKpi0_1","top pad",0,0.21,.98,.98);
        pad0_Kpi_0->Draw();
        TPad *pad0_Kpi_1 = new TPad("padKpi0_2","bottom pad",0,0,.98,0.2);
        pad0_Kpi_1->Draw();       

	pad0_Kpi_0->cd();
	g0_Kpi_MC->Draw("APBX");
	g0_Kpi_Data->Draw("P");

	pad0_Kpi_1->cd();
	Draw_Graph_Pulls(g0_Kpi_Data, g0_Kpi_MC);

	c0_Kpi->SaveAs(MomentsPDFName + "t0_Kpi_DataxMC.pdf");

	TCanvas * c1_Kpi = new TCanvas();

	TPad *pad1_Kpi_0 = new TPad("padKpi1_1","top pad",0,0.21,.98,.98);
        pad1_Kpi_0->Draw();
        TPad *pad1_Kpi_1 = new TPad("padKpi1_2","bottom pad",0,0,.98,0.2);
        pad1_Kpi_1->Draw();       

	pad1_Kpi_0->cd();
	g1_Kpi_MC->Draw("APBX");
	g1_Kpi_Data->Draw("P");

	pad1_Kpi_1->cd();
	Draw_Graph_Pulls(g1_Kpi_Data, g1_Kpi_MC);

	c1_Kpi->SaveAs(MomentsPDFName + "t1_Kpi_DataxMC.pdf");

	TCanvas * c2_Kpi = new TCanvas();

	TPad *pad2_Kpi_0 = new TPad("padKpi2_1","top pad",0,0.21,.98,.98);
        pad2_Kpi_0->Draw();
        TPad *pad2_Kpi_1 = new TPad("padKpi2_2","bottom pad",0,0,.98,0.2);
        pad2_Kpi_1->Draw();       

	pad2_Kpi_0->cd();
	g2_Kpi_MC->Draw("APBX");
	g2_Kpi_Data->Draw("P");

	pad2_Kpi_1->cd();
	Draw_Graph_Pulls(g2_Kpi_Data, g2_Kpi_MC);

	c2_Kpi->SaveAs(MomentsPDFName + "t2_Kpi_DataxMC.pdf");

	TCanvas * c3_Kpi = new TCanvas();

	TPad *pad3_Kpi_0 = new TPad("padKpi3_1","top pad",0,0.21,.98,.98);
        pad3_Kpi_0->Draw();
        TPad *pad3_Kpi_1 = new TPad("padKpi3_2","bottom pad",0,0,.98,0.2);
        pad3_Kpi_1->Draw();       

	pad3_Kpi_0->cd();
	g3_Kpi_MC->Draw("APBX");
	g3_Kpi_Data->Draw("P");

	pad3_Kpi_1->cd();
	Draw_Graph_Pulls(g3_Kpi_Data, g3_Kpi_MC);

	c3_Kpi->SaveAs(MomentsPDFName + "t3_Kpi_DataxMC.pdf");

	TCanvas * c4_Kpi = new TCanvas();

	TPad *pad4_Kpi_0 = new TPad("padKpi4_1","top pad",0,0.21,.98,.98);
        pad4_Kpi_0->Draw();
        TPad *pad4_Kpi_1 = new TPad("padKpi4_2","bottom pad",0,0,.98,0.2);
        pad4_Kpi_1->Draw();       

	pad4_Kpi_0->cd();
	g4_Kpi_MC->Draw("APBX");
	g4_Kpi_Data->Draw("P");

	pad4_Kpi_1->cd();
	Draw_Graph_Pulls(g4_Kpi_Data, g4_Kpi_MC);

	c4_Kpi->SaveAs(MomentsPDFName + "t4_Kpi_DataxMC.pdf");

	TCanvas * c0_KK = new TCanvas();

	TPad *pad0_KK_0 = new TPad("padKK0_1","top pad",0,0.21,.98,.98);
        pad0_KK_0->Draw();
        TPad *pad0_KK_1 = new TPad("padKK0_2","bottom pad",0,0,.98,0.2);
        pad0_KK_1->Draw();       

	pad0_KK_0->cd();
	g0_KK_MC->Draw("APBX");
	g0_KK_Data->Draw("P");

	pad0_KK_1->cd();
	Draw_Graph_Pulls(g0_KK_Data, g0_KK_MC);

	c0_KK->SaveAs(MomentsPDFName + "t0_KK_DataxMC.pdf");

	TCanvas * c1_KK = new TCanvas();

	TPad *pad1_KK_0 = new TPad("padKK1_1","top pad",0,0.21,.98,.98);
        pad1_KK_0->Draw();
        TPad *pad1_KK_1 = new TPad("padKK1_2","bottom pad",0,0,.98,0.2);
        pad1_KK_1->Draw();       

	pad1_KK_0->cd();
	g1_KK_MC->Draw("APBX");
	g1_KK_Data->Draw("P");

	pad1_KK_1->cd();
	Draw_Graph_Pulls(g1_KK_Data, g1_KK_MC);

	c1_KK->SaveAs(MomentsPDFName + "t1_KK_DataxMC.pdf");

	TCanvas * c2_KK = new TCanvas();

	TPad *pad2_KK_0 = new TPad("padKK2_1","top pad",0,0.21,.98,.98);
        pad2_KK_0->Draw();
        TPad *pad2_KK_1 = new TPad("padKK2_2","bottom pad",0,0,.98,0.2);
        pad2_KK_1->Draw();       

	pad2_KK_0->cd();
	g2_KK_MC->Draw("APBX");
	g2_KK_Data->Draw("P");

	pad2_KK_1->cd();
	Draw_Graph_Pulls(g2_KK_Data, g2_KK_MC);

	c2_KK->SaveAs(MomentsPDFName + "t2_KK_DataxMC.pdf");

	TCanvas * c3_KK = new TCanvas();

	TPad *pad3_KK_0 = new TPad("padKK3_1","top pad",0,0.21,.98,.98);
        pad3_KK_0->Draw();
        TPad *pad3_KK_1 = new TPad("padKK3_2","bottom pad",0,0,.98,0.2);
        pad3_KK_1->Draw();       

	pad3_KK_0->cd();
	g3_KK_MC->Draw("APBX");
	g3_KK_Data->Draw("P");

	pad3_KK_1->cd();
	Draw_Graph_Pulls(g3_KK_Data, g3_KK_MC);

	c3_KK->SaveAs(MomentsPDFName + "t3_KK_DataxMC.pdf");

	TCanvas * c4_KK = new TCanvas();

	TPad *pad4_KK_0 = new TPad("padKK4_1","top pad",0,0.21,.98,.98);
        pad4_KK_0->Draw();
        TPad *pad4_KK_1 = new TPad("padKK4_2","bottom pad",0,0,.98,0.2);
        pad4_KK_1->Draw();       

	pad4_KK_0->cd();
	g4_KK_MC->Draw("APBX");
	g4_KK_Data->Draw("P");

	pad4_KK_1->cd();
	Draw_Graph_Pulls(g4_KK_Data, g4_KK_MC);

	c4_KK->SaveAs(MomentsPDFName + "t4_KK_DataxMC.pdf");
}

