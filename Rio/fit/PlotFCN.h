#include <iostream>
#include <TTree.h>
#include <TFile.h>
#include <TColor.h>
#include <vector>
#include <TStopwatch.h>
#include <TMultiGraph.h>
#include <TGraphErrors.h>
#include "TGaxis.h"
#include "TColor.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "RioStyle.h"
#include "TH2Poly.h"
#include "Plot_Moments.h"
#include "Draw_Graph_Pulls.h"

using namespace std;

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

void PlotFCN(  TFitter *minimizer, 
		int final_state,
		bool is_gaussian, 
		double mass_min, 
		double mass_max,
		double bkg_par1, 
		double bkg_par2, 
		int number_of_resonances,
		int number_of_bkg_components, 
		int number_of_pwa_bins,
		vector<int> resonances, 
		vector<int> bkg_components,
		vector<double> S12, 
		vector<double> S13, 
		vector<double> S23,
		string rootfile,
	      	string oss2 ) {

	// Do you want to draw the interference terms  ?
	bool Draw_interf = false;
	//////////////////////////////////////////////////////////////////////////////////

	RioStyle();

	double M, 
	       s, 
	       m1sq, 
	       m2sq, 
	       m3sq, 
	       m1, 
	       m2, 
	       m3,
	       s_high, 
	       s_low, 
	       s23,
	       s12,
	       s13,
	       z;

	double spd_prod,
	       interf_prod;

	vector< vector< TComplex > > SigAmps,
		BkgAmps,
		normalization_sig;

	vector< TComplex > coefs, 
		bkg_coefs,
		sig_amp, 
		bkg_amp, 
		pwa_coefs, 
		pwa_coefs_prime;

	vector< double > normalization_bkg,  
		res_masses, 
		res_widths,
		tmp_res_extra_pars;

	vector< vector< double > > res_extra_pars;  

	TComplex temp_pwa_coef, 
		 tmp;


	vector< int > recalc_res_index, 
		recalc_resonances;

	double slow_min;
	double slow_max;

	double shigh_min;
	double shigh_max;

	double s12_min;
	double s12_max;

	double s13_min;
	double s13_max;

	double s23_min;
	double s23_max;

	TString colNameLatex_Slow,
		colNameLatex_Shi,
		colNameLatex_S12,
		colNameLatex_S13,
		colNameLatex_S23,
		colNameLatex_yslo,
		colNameLatex_yshi,
		colNameLatex_ys12,
		colNameLatex_ys13,
		colNameLatex_ys23;

	bool resonances_calc_norm[AVAILABLE_RESONANCES] = {};

	int  npt, 
	     res_number       = 0,
	     resonance_number = 0,
	     bkg_component_number = 0,
	     nparams          = minimizer->GetNumberFreeParameters() ;

	int moments_npt = 100, 
	    bin=0;
	double m_12[moments_npt], 
	       t0_12[moments_npt], 
	       t1_12[moments_npt], 
	       t2_12[moments_npt], 
	       t3_12[moments_npt], 
	       t4_12[moments_npt], 
	       bin_pop_12[moments_npt];
	double t0sq_12[moments_npt], 
	       t1sq_12[moments_npt], 
	       t2sq_12[moments_npt], 
	       t3sq_12[moments_npt], 
	       t4sq_12[moments_npt];
	double t0_err_12[moments_npt], 
	       t1_err_12[moments_npt], 
	       t2_err_12[moments_npt], 
	       t3_err_12[moments_npt], 
	       t4_err_12[moments_npt];

	double m_13[moments_npt], 
	       t0_13[moments_npt], 
	       t1_13[moments_npt], 
	       t2_13[moments_npt], 
	       t3_13[moments_npt], 
	       t4_13[moments_npt], 
	       bin_pop_13[moments_npt];
	double t0sq_13[moments_npt], 
	       t1sq_13[moments_npt], 
	       t2sq_13[moments_npt], 
	       t3sq_13[moments_npt], 
	       t4sq_13[moments_npt];
	double t0_err_13[moments_npt], 
	       t1_err_13[moments_npt], 
	       t2_err_13[moments_npt], 
	       t3_err_13[moments_npt], 
	       t4_err_13[moments_npt];

	double s12min, 
	       s12max, 
	       s13min, 
	       s13max, 
	       sacmin, 
	       sacmax, 
	       ds12, 
	       ds13, 
	       Integral_12 = 0, 
	       IntegralSq_12 = 0, 
	       Integral_13 = 0, 
	       IntegralSq_13 = 0;
	double s12_limits[moments_npt], 
	       s13_limits[moments_npt], 
	       PScorr_12[moments_npt], 
	       PScorr_13[moments_npt], 
	       PL[8];
	double E1, 
	       E2, 
	       E3, 
	       p1, 
	       p2, 
	       p3;
	double acc, 
	       w, 
	       sigma_w, 
	       av_w, 
	       av_w2;

	TGraphErrors *g0_12_Data, 
		     *g1_12_Data, 
		     *g2_12_Data, 
		     *g3_12_Data, 
		     *g4_12_Data;
	TGraphErrors *g0_13_Data, 
		     *g1_13_Data, 
		     *g2_13_Data, 
		     *g3_13_Data, 
		     *g4_13_Data;
	TGraphErrors *g0_12_MC, 
		     *g1_12_MC, 
		     *g2_12_MC, 
		     *g3_12_MC, 
		     *g4_12_MC;
	TGraphErrors *g0_13_MC, 
		     *g1_13_MC, 
		     *g2_13_MC, 
		     *g3_13_MC, 
		     *g4_13_MC;

	TMultiGraph *mg0_13 = new TMultiGraph();
	TMultiGraph *mg1_13 = new TMultiGraph();
	TMultiGraph *mg2_13 = new TMultiGraph();
	TMultiGraph *mg3_13 = new TMultiGraph();
	TMultiGraph *mg4_13 = new TMultiGraph();
	TMultiGraph *mg0_12 = new TMultiGraph();
	TMultiGraph *mg1_12 = new TMultiGraph();
	TMultiGraph *mg2_12 = new TMultiGraph();
	TMultiGraph *mg3_12 = new TMultiGraph();
	TMultiGraph *mg4_12 = new TMultiGraph();

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

	int  nnh = 100,
	     nevts;

	coefs.resize( number_of_resonances );

	double par[ 8 * number_of_resonances + number_of_bkg_components + 2 * number_of_pwa_bins ];
	vector< vector< TComplex > > coefs_product( number_of_resonances, vector< TComplex > ( number_of_resonances ) );

	for (int i = 0; i < 8 * number_of_resonances + number_of_bkg_components +
			2 * number_of_pwa_bins; i++ ) {

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

	switch ( final_state ) {
		case 0:
			M = D_Mass;
			m1 = mK;
			m2 = mK;
			m3 = mpi;
			s12_min = 0.8;
			s12_max = 3.1;
			s13_min = 0.3;
			s13_max = 2.1;
			s23_min = 0.3;
			s23_max = 2.1;
			slow_min = 0; 
			slow_max = 0;
			shigh_min = 0; 
			shigh_max = 0;
			colNameLatex_Slow = "#bf{#it{S}_{K^{+}K^{-},low} [GeV^{2}/#it{c}^{4}]}";
			colNameLatex_Shi = "#bf{#it{S}_{K^{+}K^{-},high}  [GeV^{2}/#it{c}^{4}]}";
			colNameLatex_S12 = "#it{S}_{K^{-}K^{+},12}   [GeV^{2}/#it{c}^{4}]";
			colNameLatex_S13 = "#it{S}_{K^{-}#pi^{+},13}   [GeV^{2}/#it{c}^{4}]";
			colNameLatex_S23 = "#it{S}_{K^{-}#pi^{+},23}   [GeV^{2}/#it{c}^{4}]";
			colNameLatex_yslo   = TString::Format("Entries/(%f GeV^{2}/#it{c}^{4})", (slow_max - slow_min)/nnh);
			colNameLatex_yshi   = TString::Format("Entries/(%f GeV^{2}/#it{c}^{4})", (shigh_max - shigh_min)/nnh);
			colNameLatex_ys12   = TString::Format("Entries/(%f GeV^{2}/#it{c}^{4})", (s12_max - s12_min)/nnh);
			colNameLatex_ys13   = TString::Format("Entries/(%f GeV^{2}/#it{c}^{4})", (s13_max - s13_min)/nnh);
			colNameLatex_ys23   = TString::Format("Entries/(%f GeV^{2}/#it{c}^{4})", (s23_max - s23_min)/nnh);
			break;
		case 1:
			M = D_Mass;
			m1 = mK;
			m2 = mpi;
			m3 = mpi;
			s12_min = 0.3;
			s12_max = 3.1;
			s13_min = 0.3;
			s13_max = 3.1;
			s23_min = 0.;
			s23_max = 2.1;
			slow_min = 0; 
			slow_max = 0;
			shigh_min = 0; 
			shigh_max = 0;
			colNameLatex_Slow = "#bf{#it{S}_{K^{+}K^{-},low} [GeV^{2}/#it{c}^{4}]}";
			colNameLatex_Shi = "#bf{#it{S}_{K^{+}K^{-},high}  [GeV^{2}/#it{c}^{4}]}";
			colNameLatex_S12 = "#it{S}_{K^{-}K^{+},12}   [GeV^{2}/#it{c}^{4}]";
			colNameLatex_S13 = "#it{S}_{K^{-}#pi^{+},13}   [GeV^{2}/#it{c}^{4}]";
			colNameLatex_S23 = "#it{S}_{K^{-}#pi^{+},23}   [GeV^{2}/#it{c}^{4}]";
			colNameLatex_yslo   = TString::Format("Entries/(%f GeV^{2}/#it{c}^{4})", (slow_max - slow_min)/nnh);
			colNameLatex_yshi   = TString::Format("Entries/(%f GeV^{2}/#it{c}^{4})", (shigh_max - shigh_min)/nnh);
			colNameLatex_ys12   = TString::Format("Entries/(%f GeV^{2}/#it{c}^{4})", (s12_max - s12_min)/nnh);
			colNameLatex_ys13   = TString::Format("Entries/(%f GeV^{2}/#it{c}^{4})", (s13_max - s13_min)/nnh);
			colNameLatex_ys23   = TString::Format("Entries/(%f GeV^{2}/#it{c}^{4})", (s23_max - s23_min)/nnh);
			break;
		case 2:
			M = D_Mass;
			m1 = mpi;
			m2 = mpi;
			m3 = mpi;
			s12_min = 0.;
			s12_max = 3.1;
			s13_min = 0.;
			s13_max = 3.1;
			s23_min = 0.;
			s23_max = 3.1;
			slow_min = 0; 
			slow_max = 0;
			shigh_min = 0; 
			shigh_max = 0;
			colNameLatex_Slow = "#bf{#it{S}_{K^{+}K^{-},low} [GeV^{2}/#it{c}^{4}]}";
			colNameLatex_Shi = "#bf{#it{S}_{K^{+}K^{-},high}  [GeV^{2}/#it{c}^{4}]}";
			colNameLatex_S12 = "#it{S}_{K^{-}K^{+},12}   [GeV^{2}/#it{c}^{4}]";
			colNameLatex_S13 = "#it{S}_{K^{-}#pi^{+},13}   [GeV^{2}/#it{c}^{4}]";
			colNameLatex_S23 = "#it{S}_{K^{-}#pi^{+},23}   [GeV^{2}/#it{c}^{4}]";
			colNameLatex_yslo   = TString::Format("Entries/(%f GeV^{2}/#it{c}^{4})", (slow_max - slow_min)/nnh);
			colNameLatex_yshi   = TString::Format("Entries/(%f GeV^{2}/#it{c}^{4})", (shigh_max - shigh_min)/nnh);
			colNameLatex_ys12   = TString::Format("Entries/(%f GeV^{2}/#it{c}^{4})", (s12_max - s12_min)/nnh);
			colNameLatex_ys13   = TString::Format("Entries/(%f GeV^{2}/#it{c}^{4})", (s13_max - s13_min)/nnh);
			colNameLatex_ys23   = TString::Format("Entries/(%f GeV^{2}/#it{c}^{4})", (s23_max - s23_min)/nnh);
			break;
		case 3:
			M = D_Mass;
			m1 = mK;
			m2 = mK;
			m3 = mK;
			s12_min = 0.95;
			s12_max = 1.9;
			s13_min = 0.95;
			s13_max = 1.9;
			s23_min = 0.95;
			s23_max = 1.9;
			slow_min = 0.95; 
			slow_max = 1.65;
			shigh_min = 1.1; 
			shigh_max = 1.9;
			colNameLatex_Slow = "#bf{#it{S}_{K^{+}K^{-},low} [GeV^{2}/#it{c}^{4}]}";
			colNameLatex_Shi = "#bf{#it{S}_{K^{+}K^{-},high}  [GeV^{2}/#it{c}^{4}]}";
			colNameLatex_S12 = "#it{S}_{K^{-}K^{+},12}   [GeV^{2}/#it{c}^{4}]";
			colNameLatex_S13 = "#it{S}_{K^{-}#pi^{+},13}   [GeV^{2}/#it{c}^{4}]";
			colNameLatex_S23 = "#it{S}_{K^{-}#pi^{+},23}   [GeV^{2}/#it{c}^{4}]";
			colNameLatex_yslo   = TString::Format("Entries/(%f GeV^{2}/#it{c}^{4})", (slow_max - slow_min)/nnh);
			colNameLatex_yshi   = TString::Format("Entries/(%f GeV^{2}/#it{c}^{4})", (shigh_max - shigh_min)/nnh);
			colNameLatex_ys12   = TString::Format("Entries/(%f GeV^{2}/#it{c}^{4})", (s12_max - s12_min)/nnh);
			colNameLatex_ys13   = TString::Format("Entries/(%f GeV^{2}/#it{c}^{4})", (s13_max - s13_min)/nnh);
			colNameLatex_ys23   = TString::Format("Entries/(%f GeV^{2}/#it{c}^{4})", (s23_max - s23_min)/nnh);
			break;
		default:
			cout << "Invalid final state" << endl;
			exit( -1 );
	}
	s = M*M;
	m1sq = m1*m1;
	m2sq = m2*m2;
	m3sq = m3*m3;

	TH2D *TotalPDFHist  = 
		new TH2D( "TotalPDFHist", "TotPdf", nnh, s13_min, s13_max, nnh, .9, 3.1 ),
		    *TotalToyHist = 
			    new TH2D( "TotalToyHist", "TotToy", nnh, s13_min, s13_max, nnh, s12_min, s12_max ),
		    *TotalDataHist = 
			    new TH2D( "TotalDataHist", "TotData", nnh, s13_min, s13_max, nnh, s12_min, s12_max );

	TH1F *SlowHist      =
		new TH1F( "SlowHist", "SlowHist", nnh, slow_min, slow_max ),
		    *SlowHistData  =
			    new TH1F( "SlowHistData", "SlowHistData", nnh, slow_min, slow_max ),
		    *SlowHistToy  =
			    new TH1F( "SlowHistToy", "SlowHistToy", nnh, slow_min, slow_max ),

		    *ShiHist       =
			    new TH1F( "ShiHist", "ShiHist", nnh ,shigh_min, shigh_max ),
		    *ShiHistData   =
			    new TH1F( "ShiHistData", "ShiHistData", nnh , shigh_min, shigh_max ),
		    *ShiHistToy   =
			    new TH1F( "ShiHistToy", "ShiHistToy", nnh, shigh_min, shigh_max ),

		    *S12Hist      = 
			    new TH1F( "S12Hist", "S12Hist", nnh, s12_min, s12_max ),
		    *S12HistData  = 
			    new TH1F( "S12HistData", "S12HistData", nnh, s12_min, s12_max ),
		    *S12HistToy  = 
			    new TH1F( "S12HistToy", "S12HistToy", nnh, s12_min, s12_max ),
		    *S13Hist       = 
			    new TH1F( "S13Hist", "S13Hist", nnh, s13_min, s13_max ),
		    *S13HistData   = 
			    new TH1F( "S13HistData", "S13HistData", nnh, s13_min, s13_max ),
		    *S13HistToy   = 
			    new TH1F( "S13HistToy", "S13HistToy", nnh, s13_min, s13_max ),
		    *S23Hist       = 
			    new TH1F( "S23Hist", "S23Hist", nnh, s23_min, s23_max),
		    *S23HistData   = 
			    new TH1F( "S23HistData", "S23HistData", nnh, s23_min, s23_max ),
		    *S23HistToy   = 
			    new TH1F( "S23HistToy", "S23HistToy", nnh, s23_min, s23_max ),

		    *histo_sij[ 5 ][ number_of_resonances ],
		    *histo_interf_sij[ 5 ][ number_of_resonances*(number_of_resonances - 1) ],
		    *bhisto_sij[ 5 ][ number_of_bkg_components ];

	std::vector<double> sigPdf( number_of_resonances );

	for ( Int_t iii = 0; iii < number_of_resonances; iii++ ) {

		histo_sij[ 0 ][ iii ] = 
			new TH1F( Form( "histo_slo%d", iii ), "", nnh, slow_min, slow_max );
		histo_sij[ 1 ][ iii ] = 
			new TH1F( Form( "histo_shi%d", iii ), "", nnh, shigh_min, shigh_max );
		histo_sij[ 2 ][ iii ] = 
			new TH1F( Form( "histo_s12%d", iii ), "", nnh, s12_min, s12_max );
		histo_sij[ 3 ][ iii ] = 
			new TH1F( Form( "histo_s13%d", iii ), "", nnh, s13_min, s13_max );
		histo_sij[ 4 ][ iii ] = 
			new TH1F( Form( "histo_s23%d", iii ), "", nnh, s23_min, s23_max );
		sigPdf[ iii ] = 0;

	}

	int interf_count = 0;
	for ( int interI = 0; interI < number_of_resonances; ++interI )
	{
		for (int interJ=interI; interJ < number_of_resonances; ++interJ )
		{
			if(interJ==interI)continue;
			histo_interf_sij[ 0 ][ interf_count ] = 
				new TH1F( Form( "histo_interf_slo%d%i", interI, interJ ), "", nnh, slow_min, slow_max );
			histo_interf_sij[ 1 ][ interf_count ] = 
				new TH1F( Form( "histo_interf_shi%d%i", interI, interJ ), "", nnh, shigh_min, shigh_max );
			histo_interf_sij[ 2 ][ interf_count ] = 
				new TH1F( Form( "histo_interf_s12%d%i", interI, interJ ), "", nnh, s12_min, s12_max );
			histo_interf_sij[ 3 ][ interf_count ] = 
				new TH1F( Form( "histo_interf_s13%d%i", interI, interJ ), "", nnh, s13_min, s13_max );
			histo_interf_sij[ 4 ][ interf_count ] = 
				new TH1F( Form( "histo_interf_s23%d%i", interI, interJ ), "", nnh, s23_min, s23_max );
			++interf_count;
		}
	}
	cout<<" counter of interf = "<< interf_count<<endl;
	interf_count = 0;

	for ( Int_t bb = 0; bb < number_of_bkg_components; ++bb ) {

		bhisto_sij[ 0 ][ bb ] = 
			new TH1F( Form( "bhisto_slo%d", bb ), "", nnh, slow_min, slow_max );
		bhisto_sij[ 1 ][ bb ] = 
			new TH1F( Form( "bhisto_shi%d", bb ), "", nnh, shigh_min, shigh_max );
		bhisto_sij[ 2 ][ bb ] = 
			new TH1F( Form( "bhisto_s12%d", bb ), "", nnh, s12_min, s12_max );
		bhisto_sij[ 3 ][ bb ] = 
			new TH1F( Form( "bhisto_s13%d", bb ), "", nnh, s13_min, s13_max );
		bhisto_sij[ 4 ][ bb ] = 
			new TH1F( Form( "bhisto_s23%d", bb ), "", nnh, s23_min, s23_max );

	}

	int Adaptative_bin_number = 2025;
	double Bins12Min, Bins13Min, Bins12Max, Bins13Max;
	TH2Poly *TotalPDFHist_Chi2  = new TH2Poly("TotalPDF", "TotalPDF", s13_min, s13_max, s12_min, s12_max);
	TH2Poly *TotalToyHist_Chi2  = new TH2Poly("TotalToy", "TotalToy", s13_min, s13_max, s12_min, s12_max);
	TH2Poly  *TotalDataHist_Chi2 =  new TH2Poly("TotalData", "TotalData", s13_min, s13_max, s12_min, s12_max);

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

	for ( unsigned int number_of_pwa_bins = 0; number_of_pwa_bins < KK_bin_limit.size(); 
			number_of_pwa_bins++ ) {

		if ( real_and_imaginary ) {

			temp_pwa_coef( par[ 8 * number_of_resonances + number_of_bkg_components +
					2 * number_of_pwa_bins ],
					par[ 8 * number_of_resonances + number_of_bkg_components +
					2 * number_of_pwa_bins + 1 ] );
		} else {

			temp_pwa_coef( par[ 8 * number_of_resonances + number_of_bkg_components +
					2 * number_of_pwa_bins ],
					par[ 8 * number_of_resonances + number_of_bkg_components +
					2 * number_of_pwa_bins + 1 ], 1 );
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

	s13min = (m1+m3)*(m1+m3);
	s13max = (M - m2)*(M - m2);
	s12min = (m1+m2)*(m1+m2);
	s12max = (M - m3)*(M - m3);

	ds12 = (s12max - s12min)/moments_npt;
	ds13 = (s13max - s13min)/moments_npt;

	for(int i=0;i<moments_npt;i++)  {

		// initialize moments
		t0_12[i]=0;
		t1_12[i]=0;
		t2_12[i]=0;
		t3_12[i]=0;
		t4_12[i]=0;
		bin_pop_12[i] = 0;

		// initialize  variables for errors on moments     
		t0sq_12[i] = 0;
		t1sq_12[i] = 0;
		t2sq_12[i] = 0;
		t3sq_12[i] = 0;
		t4sq_12[i] = 0;

		// define bin limits 
		s12_limits[i] = s12min + (i+1)*ds12 - ds12/2;
		m_12[i] = sqrt(s12_limits[i]);
		//cout << "m[i] = " << m[i] << endl; 

		// initialize moments
		t0_13[i]=0;
		t1_13[i]=0;
		t2_13[i]=0;
		t3_13[i]=0;
		t4_13[i]=0;
		bin_pop_13[i] = 0;

		// initialize  variables for errors on moments     
		t0sq_13[i] = 0;
		t1sq_13[i] = 0;
		t2sq_13[i] = 0;
		t3sq_13[i] = 0;
		t4sq_13[i] = 0;

		// define bin limits 
		s13_limits[i] = s13min + (i+1)*ds13 - ds13/2;
		m_13[i] = sqrt(s13_limits[i]);


		//compute phase space correction
		s13 = s13_limits[i];
		p1=0.5*sqrt( ( (s13 - m3sq - m1sq)*(s13 - m3sq - m1sq) - 4*m3sq*m1sq )/s13 );
		p2=0.5*sqrt( ( (s - s13 - m2sq)*(s - s13 - m2sq) - 4*s13*m2sq )/s13 );
		E1 = (s13 + m1sq - m3sq) / (2*sqrt(s13));
		E2 = (s - s13 - m2sq)/(2*sqrt(s13));
		sacmax =  m1sq + m2sq + 2*E1*E2 + 2*p1*p2;
		sacmin =   m1sq + m2sq + 2*E1*E2 - 2*p1*p2;
		PScorr_13[i] = 1/(sacmax - sacmin);

		s12 = s12_limits[i];
		p1=0.5*sqrt( ( (s12 - m1sq - m2sq)*(s12 - m1sq - m2sq) - 4*m1sq*m2sq )/s12 );
		p3=0.5*sqrt( ( (s - s12 - m3sq)*(s - s12 - m3sq) - 4*s12*m3sq )/s12 );
		E1 = (s12 + m1sq - m2sq) / (2*sqrt(s12));
		E3 = (s - s12 - m3sq)/(2*sqrt(s12));
		sacmax =  m1sq + m3sq + 2*E1*E3 + 2*p1*p3;
		sacmin =   m1sq + m3sq + 2*E1*E3 - 2*p1*p3;
		PScorr_12[i] = 1/(sacmax - sacmin);
		//cout << "PScorr[i] = " << PScorr[i] << endl;
		//cout << "PScorr_13[i] = " << PScorr_13[i] << ", PScorr_12[i] = " << PScorr_12[i] << endl; 
		//cout << ", s12_limits[i] = " << s12_limits[i]  << ", s13_limits[i] = " << s13_limits[i]  << endl;

	}

	nentries  = S12.size();

	for ( Long64_t jentry = 0; jentry < nentries; ++jentry ) {

		s12 = S12[jentry];
		s13 = S13[jentry];
		s23 = S23[jentry];

		if (s12 > s13){
			s_high = s12;
			s_low = s13;
		} else {
			s_high = s13;
			s_low = s12;
		} 

		TotalDataHist     ->Fill( s13, s12 );
		TotalDataHist_Chi2->Fill( s13, s12 );
		S12HistData       ->Fill( s12 );
		S13HistData       ->Fill( s13 );
		S23HistData       ->Fill( s23 );
		SlowHistData      ->Fill( s_low );
		ShiHistData       ->Fill( s_high );

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

				// find in which m12 bin this event belongs to;
				if(s12 < s12_limits[0]) bin=0;
				for(int i=1;i<moments_npt;i++) {if(s12>s12_limits[i-1] && s12<s12_limits[i]) bin=i;}

				// this is the event weight
				w = PScorr_12[bin]/acc;
				// Increment mass bin population and compute the moments for this event 
				bin_pop_12[bin] = bin_pop_12[bin] + 1.;
				t0_12[bin]= t0_12[bin] + w*PL[0];
				t1_12[bin]= t1_12[bin] + w*PL[1];
				t2_12[bin]= t2_12[bin] + w*PL[2];
				t3_12[bin]= t3_12[bin] + w*PL[3];
				t4_12[bin]= t4_12[bin] + w*PL[4];

				//Integral += t0[bin];

				//if (s12>1.6) cout << "w = " << w << ", t0 = " << t0_12[bin] << ", Integral = " << Integral_12 << ", bin = " << bin  << ", s12_limits[i-1] = " << s12_limits[bin-1]  << ", s12_limits[i] = " << s12_limits[bin]  << ", s12 = " << s12 << endl;  

				// variables for errors on moments
				t0sq_12[bin] = t0sq_12[bin] + w*PL[0]*w*PL[0];
				t1sq_12[bin] = t1sq_12[bin] + w*PL[1]*w*PL[1];
				t2sq_12[bin] = t2sq_12[bin] + w*PL[2]*w*PL[2];
				t3sq_12[bin] = t3sq_12[bin] + w*PL[3]*w*PL[3];
				t4sq_12[bin] = t4sq_12[bin] + w*PL[4]*w*PL[4];

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

				// find in which m13 bin this event belongs to;
				if(s13 < s13_limits[0]) bin=0;
				for(int i=1;i<moments_npt;i++) {if(s13>s13_limits[i-1] && s13<s13_limits[i]) bin=i;}

				// this is the event weight
				w = PScorr_13[bin]/acc;
				// Increment mass bin population and compute the moments for this event 
				bin_pop_13[bin] = bin_pop_13[bin] + 1.;
				t0_13[bin]= t0_13[bin] + w*PL[0];
				t1_13[bin]= t1_13[bin] + w*PL[1];
				t2_13[bin]= t2_13[bin] + w*PL[2];
				t3_13[bin]= t3_13[bin] + w*PL[3];
				t4_13[bin]= t4_13[bin] + w*PL[4];

				//Integral += t0[bin];

				//if (s13>0.9) cout << "w = " << w << ", t0 = " << t0_13[bin] << ", Integral = " << Integral_13 << ", bin = " << bin << ", s13_limits[i-1] = " << s13_limits[bin-1]  << ", s13_limits[i] = " << s13_limits[bin]  << ", s13 = " << s13 << endl;  

				// variables for errors on moments
				t0sq_13[bin] = t0sq_13[bin] + w*PL[0]*w*PL[0];
				t1sq_13[bin] = t1sq_13[bin] + w*PL[1]*w*PL[1];
				t2sq_13[bin] = t2sq_13[bin] + w*PL[2]*w*PL[2];
				t3sq_13[bin] = t3sq_13[bin] + w*PL[3]*w*PL[3];
				t4sq_13[bin] = t4sq_13[bin] + w*PL[4]*w*PL[4];

				//IntegralSq += t0sq[bin];
			}
		}
	}

	for (int i=0; i<moments_npt; i++) {
		Integral_12 += t0_12[i];
		IntegralSq_12 += t0sq_12[i];
	}

	for (int i=0; i<moments_npt; i++) {

		//		              cout << "t0_12[" << i <<"] = " << t0_12[i] << ", t1_12[" << i <<"] = " << t1_12[i] << ", t2_12[" << i <<"] = " << t2_12[i] << ", t3_12[" << i <<"] = " << t3_12[i] << ", t4_12[" << i <<"] = " << t4_12[i] << ", Integral = " << Integral_12 << "bin pop = "<< bin_pop_12[i] << endl; 

		t0_12[i] = t0_12[i]/Integral_12;
		t1_12[i] = t1_12[i]/Integral_12;
		t2_12[i] = t2_12[i]/Integral_12;
		t3_12[i] = t3_12[i]/Integral_12;
		t4_12[i] = t4_12[i]/Integral_12;

		t0sq_12[i] = t0sq_12[i]/IntegralSq_12;
		t1sq_12[i] = t1sq_12[i]/IntegralSq_12;
		t2sq_12[i] = t2sq_12[i]/IntegralSq_12;
		t3sq_12[i] = t3sq_12[i]/IntegralSq_12;
		t4sq_12[i] = t4sq_12[i]/IntegralSq_12;

		av_w = t0_12[i]/bin_pop_12[i];
		av_w2 = t0sq_12[i]/bin_pop_12[i];
		sigma_w = av_w2 - av_w*av_w;
		t0_err_12[i] = fabs(t0_12[i])/sqrt(bin_pop_12[i])*sqrt(1.+(sigma_w/av_w)*(sigma_w/av_w));

		av_w = t1_12[i]/bin_pop_12[i];
		av_w2 = t1sq_12[i]/bin_pop_12[i];
		sigma_w = av_w2 - av_w*av_w;
		t1_err_12[i] = fabs(t1_12[i])/sqrt(bin_pop_12[i])*sqrt(1.+(sigma_w/av_w)*(sigma_w/av_w));

		av_w = t2_12[i]/bin_pop_12[i];
		av_w2 = t2sq_12[i]/bin_pop_12[i];
		sigma_w = av_w2 - av_w*av_w;
		t2_err_12[i] = fabs(t2_12[i])/sqrt(bin_pop_12[i])*sqrt(1.+(sigma_w/av_w)*(sigma_w/av_w));

		av_w = t3_12[i]/bin_pop_12[i];
		av_w2 = t3sq_12[i]/bin_pop_12[i];
		sigma_w = av_w2 - av_w*av_w;
		t3_err_12[i] = fabs(t3_12[i])/sqrt(bin_pop_12[i])*sqrt(1.+(sigma_w/av_w)*(sigma_w/av_w));

		av_w = t4_12[i]/bin_pop_12[i];
		av_w2 = t4sq_12[i]/bin_pop_12[i];
		sigma_w = av_w2 - av_w*av_w;
		t4_err_12[i] = fabs(t4_12[i])/sqrt(bin_pop_12[i])*sqrt(1.+(sigma_w/av_w)*(sigma_w/av_w));

		//		cout << "t0_err_12 = " << t0_err_12[i] << "t1_err_12 = " << t1_err_12[i] << "t2_err_12 = " << t2_err_12[i] << "t3_err_12 = " << t3_err_12[i] << "t4_err_12 = " << t4_err_12[i] << ", bin pop = "<< bin_pop_12[i] << endl;
	}

	plot_moments(m_12, t0_12, t1_12, t2_12, t3_12, t4_12, t0_err_12, t1_err_12, t2_err_12, t3_err_12, t4_err_12, g0_12_Data, g1_12_Data, g2_12_Data, g3_12_Data, g4_12_Data);

	for (int i=0; i<moments_npt; i++) {
		Integral_13 += t0_13[i];
		IntegralSq_13 += t0sq_13[i];
	}

	for (int i=0; i<moments_npt; i++) {

		//		              cout << "t0_13[" << i <<"] = " << t0_13[i] << ", t1[" << i <<"] = " << t1_13[i] << ", t2_13[" << i <<"] = " << t2_13[i] << ", t3_13[" << i <<"] = " << t3_13[i] << ", t4_13[" << i <<"] = " << t4_13[i] << ", Integral = " << Integral_13 << "bin pop = "<< bin_pop_13[i] << endl; 

		t0_13[i] = t0_13[i]/Integral_13;
		t1_13[i] = t1_13[i]/Integral_13;
		t2_13[i] = t2_13[i]/Integral_13;
		t3_13[i] = t3_13[i]/Integral_13;
		t4_13[i] = t4_13[i]/Integral_13;

		t0sq_13[i] = t0sq_13[i]/IntegralSq_13;
		t1sq_13[i] = t1sq_13[i]/IntegralSq_13;
		t2sq_13[i] = t2sq_13[i]/IntegralSq_13;
		t3sq_13[i] = t3sq_13[i]/IntegralSq_13;
		t4sq_13[i] = t4sq_13[i]/IntegralSq_13;

		av_w = t0_13[i]/bin_pop_13[i];
		av_w2 = t0sq_13[i]/bin_pop_13[i];
		sigma_w = av_w2 - av_w*av_w;
		t0_err_13[i] = fabs(t0_13[i])/sqrt(bin_pop_13[i])*sqrt(1.+(sigma_w/av_w)*(sigma_w/av_w));

		av_w = t1_13[i]/bin_pop_13[i];
		av_w2 = t1sq_13[i]/bin_pop_13[i];
		sigma_w = av_w2 - av_w*av_w;
		t1_err_13[i] = fabs(t1_13[i])/sqrt(bin_pop_13[i])*sqrt(1.+(sigma_w/av_w)*(sigma_w/av_w));

		av_w = t2_13[i]/bin_pop_13[i];
		av_w2 = t2sq_13[i]/bin_pop_13[i];
		sigma_w = av_w2 - av_w*av_w;
		t2_err_13[i] = fabs(t2_13[i])/sqrt(bin_pop_13[i])*sqrt(1.+(sigma_w/av_w)*(sigma_w/av_w));

		av_w = t3_13[i]/bin_pop_13[i];
		av_w2 = t3sq_13[i]/bin_pop_13[i];
		sigma_w = av_w2 - av_w*av_w;
		t3_err_13[i] = fabs(t3_13[i])/sqrt(bin_pop_13[i])*sqrt(1.+(sigma_w/av_w)*(sigma_w/av_w));

		av_w = t4_13[i]/bin_pop_13[i];
		av_w2 = t4sq_13[i]/bin_pop_13[i];
		sigma_w = av_w2 - av_w*av_w;
		t4_err_13[i] = fabs(t4_13[i])/sqrt(bin_pop_13[i])*sqrt(1.+(sigma_w/av_w)*(sigma_w/av_w));

		//		cout << "t0_err_13 = " << t0_err_13[i] << "t1_err_13 = " << t1_err_13[i] << "t2_err_13 = " << t2_err_13[i] << "t3_err_13 = " << t3_err_13[i] << "t4_err_13 = " << t4_err_13[i] << "bin pop = "<< bin_pop_12[i] << endl;
	}

	plot_moments(m_13, t0_13, t1_13, t2_13, t3_13, t4_13, t0_err_13, t1_err_13, t2_err_13, t3_err_13, t4_err_13, g0_13_Data, g1_13_Data, g2_13_Data, g3_13_Data, g4_13_Data);

	for(int i=0;i<moments_npt;i++)  {

		// initialize moments
		t0_12[i]=0;
		t1_12[i]=0;
		t2_12[i]=0;
		t3_12[i]=0;
		t4_12[i]=0;
		bin_pop_12[i] = 0;

		// initialize  variables for errors on moments     
		t0sq_12[i] = 0;
		t1sq_12[i] = 0;
		t2sq_12[i] = 0;
		t3sq_12[i] = 0;
		t4sq_12[i] = 0;

		// initialize moments
		t0_13[i]=0;
		t1_13[i]=0;
		t2_13[i]=0;
		t3_13[i]=0;
		t4_13[i]=0;
		bin_pop_13[i] = 0;

		// initialize  variables for errors on moments     
		t0sq_13[i] = 0;
		t1sq_13[i] = 0;
		t2sq_13[i] = 0;
		t3sq_13[i] = 0;
		t4sq_13[i] = 0;

	}

	Integral_12 = 0;
	IntegralSq_12 = 0;
	Integral_13 = 0;
	IntegralSq_13 = 0;

	TFile f( "/Users/dvieira/Dropbox/toyMCGenerator/ntuple_flat_kkpi_100M_0.root" );
	TTree *t1 = (TTree*)f.Get( "DecayTree" );

	t1->SetBranchAddress( "s12",    &s12 );
	t1->SetBranchAddress( "s13",    &s13 );
	t1->SetBranchAddress( "s23",    &s23 );
	t1->SetBranchAddress( "s_low",  &s_low );
	t1->SetBranchAddress( "s_high", &s_high );
	nevts = t1->GetEntries();
	//nevts = 10000;

	double totalPdf = 0;
	double L_s = 0;
	for ( Long64_t jentry = 0; jentry < nevts; ++jentry ) {
		t1->GetEntry( jentry );

		//			Slo.push_back(s_low);
		//			S13.push_back(s_high);

		//cout << " s12 = " << s12 << ", s13 = " << s13 << ", s23 = " << s23 << endl;

		Amplitudes( final_state, 
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

		L_s = TotPdf( bkg_fraction, 
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

				// find in which m12 bin this event belongs to;
				if(s12 < s12_limits[0]) bin=0;
				for(int i=1;i<moments_npt;i++) {if(s12>s12_limits[i-1] && s12<s12_limits[i]) bin=i;}

				// this is the event weight
				w = PScorr_12[bin]*L_s/acc;
				// Increment mass bin population and compute the moments for this event 
				bin_pop_12[bin] = bin_pop_12[bin] + 1.;
				t0_12[bin]= t0_12[bin] + w*PL[0];
				t1_12[bin]= t1_12[bin] + w*PL[1];
				t2_12[bin]= t2_12[bin] + w*PL[2];
				t3_12[bin]= t3_12[bin] + w*PL[3];
				t4_12[bin]= t4_12[bin] + w*PL[4];

				//Integral += t0[bin];

				//cout << "w = " << w << ", t0 = " << t0[bin] << ", Integral = " << Integral << endl;  

				// variables for errors on moments
				t0sq_12[bin] = t0sq_12[bin] + w*PL[0]*w*PL[0];
				t1sq_12[bin] = t1sq_12[bin] + w*PL[1]*w*PL[1];
				t2sq_12[bin] = t2sq_12[bin] + w*PL[2]*w*PL[2];
				t3sq_12[bin] = t3sq_12[bin] + w*PL[3]*w*PL[3];
				t4sq_12[bin] = t4sq_12[bin] + w*PL[4]*w*PL[4];

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

				// find in which m13 bin this event belongs to;
				if(s13 < s13_limits[0]) bin=0;
				for(int i=1;i<moments_npt;i++) {if(s13>s13_limits[i-1] && s13<s13_limits[i]) bin=i;}

				// this is the event weight
				w = PScorr_13[bin]*L_s/acc;
				// Increment mass bin population and compute the moments for this event 
				bin_pop_13[bin] = bin_pop_13[bin] + 1.;
				t0_13[bin]= t0_13[bin] + w*PL[0];
				t1_13[bin]= t1_13[bin] + w*PL[1];
				t2_13[bin]= t2_13[bin] + w*PL[2];
				t3_13[bin]= t3_13[bin] + w*PL[3];
				t4_13[bin]= t4_13[bin] + w*PL[4];

				//Integral += t0[bin];

				//cout << "w = " << w << ", t0 = " << t0[bin] << ", Integral = " << Integral << endl;  

				// variables for errors on moments
				t0sq_13[bin] = t0sq_13[bin] + w*PL[0]*w*PL[0];
				t1sq_13[bin] = t1sq_13[bin] + w*PL[1]*w*PL[1];
				t2sq_13[bin] = t2sq_13[bin] + w*PL[2]*w*PL[2];
				t3sq_13[bin] = t3sq_13[bin] + w*PL[3]*w*PL[3];
				t4sq_13[bin] = t4sq_13[bin] + w*PL[4]*w*PL[4];

				//IntegralSq += t0sq[bin];
			}
		}

		// Fill Interference histograms 
		for (int i = 0; i < number_of_resonances; i++){
			for (int j = i; j < number_of_resonances; j++)
			{
				if(j==i)continue;

				interf_prod = 2*(coefs_product[i][j].Rho() * SigAmps[jentry][i]*SigAmps[jentry][j].Conjugate(SigAmps[jentry][j])).Re();

				histo_interf_sij[ 0 ][ interf_count ]->Fill( s_low,  interf_prod * ( 1 - bkg_fraction ) * accTerm / NL_s );
				histo_interf_sij[ 1 ][ interf_count ]->Fill( s_high, interf_prod * ( 1 - bkg_fraction ) * accTerm / NL_s );
				histo_interf_sij[ 2 ][ interf_count ]->Fill( s12,    interf_prod * ( 1 - bkg_fraction ) * accTerm / NL_s );
				histo_interf_sij[ 3 ][ interf_count ]->Fill( s13,    interf_prod * ( 1 - bkg_fraction ) * accTerm / NL_s );
				histo_interf_sij[ 4 ][ interf_count ]->Fill( s23,    interf_prod * ( 1 - bkg_fraction ) * accTerm / NL_s );
				++interf_count;
			}
		}

		interf_count = 0;


		for ( int xx = 0; xx < number_of_resonances; ++xx ) {

			spd_prod = coefs_product[ xx ][ xx ].Rho() * sig_amp[ xx ].Rho2();

			histo_sij[ 0 ][ xx ]->Fill( s_low,  
					spd_prod * ( 1 - bkg_fraction ) * accTerm / NL_s );

			histo_sij[ 1 ][ xx ]->Fill( s_high, 
					spd_prod * ( 1 - bkg_fraction ) * accTerm / NL_s ); 

			histo_sij[ 2 ][ xx ]->Fill( s12,
					spd_prod * ( 1 - bkg_fraction ) * accTerm / NL_s );

			histo_sij[ 3 ][ xx ]->Fill( s13,
					spd_prod * ( 1 - bkg_fraction ) * accTerm / NL_s ); 

			histo_sij[ 4 ][ xx ]->Fill( s23,
					spd_prod * ( 1 - bkg_fraction ) * accTerm / NL_s );                            
		}


		for ( int i = 0; i < number_of_bkg_components; ++i ) {

			bhisto_sij[ 0 ][ i ]->Fill( s_low,  
					(par[ ( i + 8 * number_of_resonances ) ] *
					 bkg_amp[ i ].Re() /
					 normalization_bkg.at( i ) ) *
					bkg_fraction );

			bhisto_sij[ 1 ][ i ]->Fill( s_high, 
					(par[ ( i + 8 * number_of_resonances ) ] *
					 bkg_amp[ i ].Re() /
					 normalization_bkg.at( i ) ) *
					bkg_fraction );

			bhisto_sij[ 2 ][ i ]->Fill( s12, 
					( par[ ( i + 8 * number_of_resonances ) ] *
					  bkg_amp[ i ].Re() /
					  normalization_bkg.at( i ) ) *
					bkg_fraction );

			bhisto_sij[ 3 ][ i ]->Fill( s13, 
					( par[ ( i + 8 * number_of_resonances ) ] *
					  bkg_amp[ i ].Re() /
					  normalization_bkg.at(i)) *
					bkg_fraction );

			bhisto_sij[ 4 ][ i ]->Fill( s23, 
					( par[ ( i + 8 * number_of_resonances ) ] *
					  bkg_amp[ i ].Re() /
					  normalization_bkg.at(i)) *
					bkg_fraction );                                            

		}

		TotalPDFHist->Fill( s13,  s12, L_s );
		TotalToyHist->Fill( s13,  s12);
		TotalPDFHist_Chi2->Fill( s13,  s12, L_s );
		TotalToyHist_Chi2->Fill( s13,  s12);
		SlowHist    ->Fill( s_low,  L_s );
		SlowHistToy ->Fill( s_low );
		ShiHist     ->Fill( s_high, L_s );
		ShiHistToy  ->Fill( s_high);
		S12Hist    ->Fill( s12,  L_s );
		S12HistToy ->Fill( s12 );
		S13Hist     ->Fill( s13, L_s );
		S13HistToy  ->Fill( s13 );
		S23Hist     ->Fill( s23,    L_s );
		S23HistToy  ->Fill( s23 );
		totalPdf += L_s;

	}

	for (int i=0; i<moments_npt; i++) {
		Integral_12 += t0_12[i];
		IntegralSq_12 += t0sq_12[i];
	}

	for (int i=0; i<moments_npt; i++) {

		//              cout << "t0[" << i <<"] = " << t0[i] << ", t1[" << i <<"] = " << t1[i] << ", t2[" << i <<"] = " << t2[i] << ", t3[" << i <<"] = " << t3[i] << ", t4[" << i <<"] = " << t4[i] << ", Integral = " << Integral << "bin pop = "<< bin_pop[i] << endl; 

		t0_12[i] = t0_12[i]/Integral_12;
		t1_12[i] = t1_12[i]/Integral_12;
		t2_12[i] = t2_12[i]/Integral_12;
		t3_12[i] = t3_12[i]/Integral_12;
		t4_12[i] = t4_12[i]/Integral_12;

		t0sq_12[i] = t0sq_12[i]/IntegralSq_12;
		t1sq_12[i] = t1sq_12[i]/IntegralSq_12;
		t2sq_12[i] = t2sq_12[i]/IntegralSq_12;
		t3sq_12[i] = t3sq_12[i]/IntegralSq_12;
		t4sq_12[i] = t4sq_12[i]/IntegralSq_12;

		av_w = t0_12[i]/bin_pop_12[i];
		av_w2 = t0sq_12[i]/bin_pop_12[i];
		sigma_w = av_w2 - av_w*av_w;
		t0_err_12[i] = fabs(t0_12[i])/sqrt(bin_pop_12[i])*sqrt(1.+(sigma_w/av_w)*(sigma_w/av_w));

		av_w = t1_12[i]/bin_pop_12[i];
		av_w2 = t1sq_12[i]/bin_pop_12[i];
		sigma_w = av_w2 - av_w*av_w;
		t1_err_12[i] = fabs(t1_12[i])/sqrt(bin_pop_12[i])*sqrt(1.+(sigma_w/av_w)*(sigma_w/av_w));

		av_w = t2_12[i]/bin_pop_12[i];
		av_w2 = t2sq_12[i]/bin_pop_12[i];
		sigma_w = av_w2 - av_w*av_w;
		t2_err_12[i] = fabs(t2_12[i])/sqrt(bin_pop_12[i])*sqrt(1.+(sigma_w/av_w)*(sigma_w/av_w));

		av_w = t3_12[i]/bin_pop_12[i];
		av_w2 = t3sq_12[i]/bin_pop_12[i];
		sigma_w = av_w2 - av_w*av_w;
		t3_err_12[i] = fabs(t3_12[i])/sqrt(bin_pop_12[i])*sqrt(1.+(sigma_w/av_w)*(sigma_w/av_w));

		av_w = t4_12[i]/bin_pop_12[i];
		av_w2 = t4sq_12[i]/bin_pop_12[i];
		sigma_w = av_w2 - av_w*av_w;
		t4_err_12[i] = fabs(t4_12[i])/sqrt(bin_pop_12[i])*sqrt(1.+(sigma_w/av_w)*(sigma_w/av_w));

		//		cout << "t0_err_12 = " << t0_err_12[i] << "t1_err_12 = " << t1_err_12[i] << "t2_err_12 = " << t2_err_12[i] << "t3_err_12 = " << t3_err_12[i] << "t4_err_12 = " << t4_err_12[i] << endl;
	}

	plot_moments(m_12, t0_12, t1_12, t2_12, t3_12, t4_12, t0_err_12, t1_err_12, t2_err_12, t3_err_12, t4_err_12, g0_12_MC, g1_12_MC, g2_12_MC, g3_12_MC, g4_12_MC);

	for (int i=0; i<moments_npt; i++) {
		Integral_13 += t0_13[i];
		IntegralSq_13 += t0sq_13[i];
	}

	for (int i=0; i<moments_npt; i++) {

		//              cout << "t0[" << i <<"] = " << t0[i] << ", t1[" << i <<"] = " << t1[i] << ", t2[" << i <<"] = " << t2[i] << ", t3[" << i <<"] = " << t3[i] << ", t4[" << i <<"] = " << t4[i] << ", Integral = " << Integral << "bin pop = "<< bin_pop[i] << endl; 

		t0_13[i] = t0_13[i]/Integral_13;
		t1_13[i] = t1_13[i]/Integral_13;
		t2_13[i] = t2_13[i]/Integral_13;
		t3_13[i] = t3_13[i]/Integral_13;
		t4_13[i] = t4_13[i]/Integral_13;

		t0sq_13[i] = t0sq_13[i]/IntegralSq_13;
		t1sq_13[i] = t1sq_13[i]/IntegralSq_13;
		t2sq_13[i] = t2sq_13[i]/IntegralSq_13;
		t3sq_13[i] = t3sq_13[i]/IntegralSq_13;
		t4sq_13[i] = t4sq_13[i]/IntegralSq_13;

		av_w = t0_13[i]/bin_pop_13[i];
		av_w2 = t0sq_13[i]/bin_pop_13[i];
		sigma_w = av_w2 - av_w*av_w;
		t0_err_13[i] = fabs(t0_13[i])/sqrt(bin_pop_13[i])*sqrt(1.+(sigma_w/av_w)*(sigma_w/av_w));

		av_w = t1_13[i]/bin_pop_13[i];
		av_w2 = t1sq_13[i]/bin_pop_13[i];
		sigma_w = av_w2 - av_w*av_w;
		t1_err_13[i] = fabs(t1_13[i])/sqrt(bin_pop_13[i])*sqrt(1.+(sigma_w/av_w)*(sigma_w/av_w));

		av_w = t2_13[i]/bin_pop_13[i];
		av_w2 = t2sq_13[i]/bin_pop_13[i];
		sigma_w = av_w2 - av_w*av_w;
		t2_err_13[i] = fabs(t2_13[i])/sqrt(bin_pop_13[i])*sqrt(1.+(sigma_w/av_w)*(sigma_w/av_w));

		av_w = t3_13[i]/bin_pop_13[i];
		av_w2 = t3sq_13[i]/bin_pop_13[i];
		sigma_w = av_w2 - av_w*av_w;
		t3_err_13[i] = fabs(t3_13[i])/sqrt(bin_pop_13[i])*sqrt(1.+(sigma_w/av_w)*(sigma_w/av_w));

		av_w = t4_13[i]/bin_pop_13[i];
		av_w2 = t4sq_13[i]/bin_pop_13[i];
		sigma_w = av_w2 - av_w*av_w;
		t4_err_13[i] = fabs(t4_13[i])/sqrt(bin_pop_13[i])*sqrt(1.+(sigma_w/av_w)*(sigma_w/av_w));

		//cout << "t0_err_13 = " << t0_err_13[i] << "t1_err_13 = " << t1_err_13[i] << "t2_err_13 = " << t2_err_13[i] << "t3_err_13 = " << t3_err_13[i] << "t4_err_13 = " << t4_err_13[i] << endl;
	}

	plot_moments(m_13, t0_13, t1_13, t2_13, t3_13, t4_13, t0_err_13, t1_err_13, t2_err_13, t3_err_13, t4_err_13, g0_13_MC, g1_13_MC, g2_13_MC, g3_13_MC, g4_13_MC);

	double scale = double(nentries) / totalPdf;

	TotalPDFHist->Sumw2();
	TotalToyHist->Sumw2();
	SlowHistToy->Sumw2();
	ShiHistToy->Sumw2();
	S12HistToy->Sumw2();
	S13HistToy->Sumw2();
	S23HistToy->Sumw2();

	TotalPDFHist->Scale( scale );
	TotalPDFHist_Chi2->Scale( scale );
	SlowHist->Scale( scale );
	ShiHist->Scale( scale );
	S12Hist->Scale( scale );
	S13Hist->Scale( scale );
	S23Hist->Scale( scale );


	set_plot_style_Dalitz();

	TCanvas * c0 = new TCanvas(); 
	TotalPDFHist->Draw("colz");
	c0->Print("TotalPDFHist.pdf");

	TCanvas * c1 = new TCanvas(); 
	TotalPDFHist_Chi2->Draw("colz");
	c1->Print("TotalPDFHist_Chi2.pdf"); 

	for ( int xxx = 0; xxx < number_of_resonances; ++xxx ) {
		for ( int yy = 0; yy < 5; ++yy ) {

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
		for ( int yy = 0; yy < 5; ++yy ) {

			bhisto_sij[ yy ][ bb ]->SetLineWidth( 2 );
			bhisto_sij[ yy ][ bb ]->SetLineColor( kRed );
			bhisto_sij[ yy ][ bb ]->SetLineStyle( 3 );

			bhisto_sij[ yy ][ bb ]->Scale(scale);
		}
	}


	string output_histo_file_name = rootfile.substr(0,rootfile.size() - 5);
	TString PDFName = "Results_PDF_fit/" + output_histo_file_name + ".pdf"; 
	TString MomentsPDFName = "Results_PDF_fit/" + output_histo_file_name; 
	output_histo_file_name += ".root";


	SlowHist->SetXTitle ( colNameLatex_Slow );
	SlowHist->SetYTitle ( colNameLatex_yslo );
	ShiHist->SetXTitle ( colNameLatex_Shi );
	ShiHist->SetYTitle ( colNameLatex_yshi );
	S12Hist->SetXTitle( colNameLatex_S12 );
	S12Hist->SetYTitle( colNameLatex_ys12);
	S13Hist->SetXTitle ( colNameLatex_S13 );
	S13Hist->SetYTitle ( colNameLatex_ys13 );
	S23Hist->SetXTitle ( colNameLatex_S23 );
	S23Hist->SetYTitle ( colNameLatex_ys23 );


	TotalDataHist->SetXTitle( colNameLatex_S13 );
	TotalDataHist->SetYTitle( colNameLatex_S12 );
	TotalDataHist_Chi2->SetXTitle( colNameLatex_S13 );
	TotalDataHist_Chi2->SetYTitle( colNameLatex_S12 );
	TotalPDFHist->SetXTitle ( colNameLatex_S13 );
	TotalPDFHist->SetYTitle ( colNameLatex_S12 );
	TotalPDFHist_Chi2->SetXTitle ( colNameLatex_S12 );
	TotalPDFHist_Chi2->SetYTitle ( colNameLatex_S12 );


	TotalDataHist->GetYaxis()->SetTitleOffset( 1.1 );
	TotalDataHist->GetXaxis()->SetTitleOffset( 0.9 );
	TotalDataHist_Chi2->GetYaxis()->SetTitleOffset( 1.1 );
	TotalDataHist_Chi2->GetXaxis()->SetTitleOffset( 0.9 );
	TotalPDFHist->GetYaxis()->SetTitleOffset( 1.1 );
	TotalPDFHist->GetXaxis()->SetTitleOffset( 0.9 );
	TotalPDFHist_Chi2->GetYaxis()->SetTitleOffset( 1.1 );
	TotalPDFHist_Chi2->GetXaxis()->SetTitleOffset( 0.9 );


	SlowHist->GetYaxis()->SetTitleOffset( 0.9 );
	SlowHist->GetXaxis()->SetTitleOffset( 0.9 );
	ShiHist->GetYaxis()->SetTitleOffset( 0.9 );
	ShiHist->GetXaxis()->SetTitleOffset( 0.9 );
	S12Hist->GetYaxis()->SetTitleOffset( 0.9 );
	S12Hist->GetXaxis()->SetTitleOffset( 0.9 );
	S13Hist->GetYaxis()->SetTitleOffset( 0.9);
	S13Hist->GetXaxis()->SetTitleOffset( 0.9 );
	S23Hist->GetYaxis()->SetTitleOffset( 0.9 );
	S23Hist->GetXaxis()->SetTitleOffset( 0.9 );


	SlowHist->SetLineColor( kBlue );
	SlowHist->SetLineWidth( 2 );
	ShiHist->SetLineColor( kBlue );
	ShiHist->SetLineWidth( 2 );
	S12Hist->SetLineColor( kBlue );
	S12Hist->SetLineWidth( 2 );
	S13Hist->SetLineColor( kBlue );
	S13Hist->SetLineWidth( 2 );
	S23Hist->SetLineColor( kBlue );
	S23Hist->SetLineWidth( 2 );

	SlowHistData->SetLineWidth(1);
	ShiHistData->SetLineWidth(1);
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

	TLegend     *lslo =
		new TLegend( 0.65, 0.7, 0.8, 0.9 ), // ver https://root.cern.ch/root/html/TPave.html#TPave:TPave@1
		    *lshi =
			    new TLegend( 0.65, 0.6, 0.85, 0.85),
		    *ls12 = 
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

			lslo->AddEntry( histo_sij[ 0 ][ resonance_number ], 
					resonant_channel_string_tex[ xx ].c_str(), 
					"l" );
			lshi->AddEntry( histo_sij[ 1 ][ resonance_number ], 
					resonant_channel_string_tex[ xx ].c_str(), 
					"l" );
			ls12->AddEntry( histo_sij[ 2 ][ resonance_number ], 
					resonant_channel_string_tex[ xx ].c_str(), 
					"l" );
			ls13->AddEntry( histo_sij[ 3 ][ resonance_number ], 
					resonant_channel_string_tex[ xx ].c_str(), 
					"l" );
			ls23->AddEntry( histo_sij[ 4 ][ resonance_number ], 
					resonant_channel_string_tex[ xx ].c_str(), 
					"l" ); 
			resonance_number++;              
		}
	}
	for (int xx = 0; xx < AVAILABLE_BKG_COMPONENTS; ++xx ) {
		if(bkg_components[xx]){
			//cout << "xx = " << xx <<  ", bkg_component_number = " << bkg_component_number << endl;
			lslo->AddEntry( bhisto_sij[ 0 ][ bkg_component_number ], 
					bkg_component_string[ xx ].c_str(), 
					"l" );
			lshi->AddEntry( bhisto_sij[ 1 ][ bkg_component_number ], 
					bkg_component_string[ xx ].c_str(), 
					"l" );
			ls12->AddEntry( bhisto_sij[ 2 ][ bkg_component_number ], 
					bkg_component_string[ xx ].c_str(), 
					"l" );
			ls13->AddEntry( bhisto_sij[ 3 ][ bkg_component_number ], 
					bkg_component_string[ xx ].c_str(), 
					"l" );
			ls23->AddEntry( bhisto_sij[ 4 ][ bkg_component_number ], 
					bkg_component_string[ xx ].c_str(), 
					"l" );              
			bkg_component_number++;              
		}
	}

	////////////////////////////////////////////////// Draw S_low and S_high ///////////////////////////////////////////////////////////////

	// Draw Slo 
	AccCanv->cd( 1 );
	pad1->cd();
	SlowHist->Draw( "C HIST" );
	SlowHistData->Draw( "E1 same" );

	for( int xx = 0; xx < number_of_resonances; ++xx ){
		histo_sij[ 0 ][ xx ]->Draw("csame");
	}

	if(Draw_interf){
		for (int i = 0; i < number_of_resonances; i++){
			for (int j = i; j < number_of_resonances; j++){
				if(j==i)continue;
				histo_interf_sij[ 0 ][ interf_count ]->Draw("csame");
				++interf_count;
			}
		}
		interf_count = 0;          
	}
	for( int bb0 = 0; bb0 < number_of_bkg_components; ++bb0 ){
		bhisto_sij[ 0 ][ bb0 ]->Draw( "csame" );
	}    


	lslo->Draw();
	pad2->cd();
	DrawMypull( SlowHistData, SlowHist, SlowHistToy );


	// Draw Shi 
	AccCanv->cd( 2 );
	pad3->cd();
	ShiHist->Draw( "C HIST" );  
	ShiHistData->Draw( "E1 same" );  

	for( int xx = 0; xx < number_of_resonances; ++xx ){
		histo_sij[ 1 ][ xx ]->Draw( "HIST C SAME" ); 
	}

	if(Draw_interf){     
		for (int i = 0; i < number_of_resonances; i++){
			for (int j = i; j < number_of_resonances; j++){
				if(j==i)continue;
				histo_interf_sij[ 1 ][ interf_count ]->Draw("csame");
				++interf_count;
			}
		}
		interf_count = 0;          
	}

	for( int bb0 = 0; bb0 < number_of_bkg_components; ++bb0 ) {
		bhisto_sij[ 1 ][ bb0 ]->Draw( "HIST C SAME" ); 
	}

	pad4->cd();
	DrawMypull( ShiHistData, ShiHist,ShiHistToy );

	AccCanv->Print( PDFName );
	AccCanv->Update();


	////////////////////////////////////////////////// End of Draw S_low and S_high //////////////////////////////////////////////////////

	///////////////////////////////////////////////////// Draw S13 and S12 ///////////////////////////////////////////////////////////////

	AccCanv->cd( 1 );
	gPad->SetLogy();

	S13Hist->Draw( "HIST C " );
	S13HistData->Draw( "E1 SAME" );  

	cout << "Finished setting histograms" << endl;  

	for( int xx = 0; xx < number_of_resonances; ++xx ) 
	{

		histo_sij[ 3 ][ xx ]->Draw( "HIST C SAME" ); 

	}

	for( int bb0 = 0; bb0 < number_of_bkg_components; ++bb0 ) 
	{

		bhisto_sij[ 3 ][ bb0 ]->Draw( "HIST C SAME" ); 

	}

	cout << "Finished Drawing" << endl;  


	//ls13->Draw();

	pad2->cd();

	DrawMypull( S13HistData, S13Hist, S13HistToy );

	pad3->cd();

	gPad->SetLogy();
	S12Hist->Draw( "C HIST" );
	S12HistData->Draw( "E same" );

	for( int xx = 0; xx < number_of_resonances; ++xx )

		histo_sij[ 2 ][ xx ]->Draw("csame");

	for( int bb0 = 0; bb0 < number_of_bkg_components; ++bb0 )

		bhisto_sij[ 2 ][ bb0 ]->Draw( "csame" );

	//ls12->Draw();

	pad4->cd();

	DrawMypull( S12HistData, S12Hist, S12HistToy );


	AccCanv->Print( PDFName );
	AccCanv->Update();

	////////////////////////////////////////////////// End of Draw S13 and S12 //////////////////////////////////////////////////////

	///////////////////////////////////////////////////// Draw S23 and Chi2 ///////////////////////////////////////////////////////////////
	pad1->cd();
	S23Hist->Draw( "C HIST" );
	S23HistData->Draw( "Esame" );

	for( int xx = 0; xx < number_of_resonances; ++xx )

		histo_sij[ 4 ][ xx ]->Draw("csame");

	for( int bb0 = 0; bb0 < number_of_bkg_components; ++bb0 )

		bhisto_sij[ 4 ][ bb0 ]->Draw( "csame" );

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
	chi2->SetXTitle( colNameLatex_S13 );
	chi2->SetYTitle( colNameLatex_S12 );
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

	double tChi2 = double(chisq/nDof);
	AccCanv->Close();

        ofstream myfile;  //file to output the fit result - ofstream : Stream class to write on files
        myfile.open (oss2.c_str(),ios::app); //ios:app :  append to the end of the file
        myfile << "chi2 " << tChi2 << endl; 
        myfile.close();

	////////////////////////////////////////////////// End of Draw S23 and Chi2 //////////////////////////////////////////////////////

	TFile * histo_file = new TFile(output_histo_file_name.c_str(), "RECREATE");
	TotalPDFHist      ->SetDirectory(histo_file);
	TotalDataHist     ->SetDirectory(histo_file);
	TotalToyHist      ->SetDirectory(histo_file);
	TotalPDFHist_Chi2 ->SetDirectory(histo_file);
	TotalDataHist_Chi2->SetDirectory(histo_file);
	TotalToyHist_Chi2 ->SetDirectory(histo_file);
	S12Hist           ->SetDirectory(histo_file);
	S12HistData       ->SetDirectory(histo_file);
	S12HistToy        ->SetDirectory(histo_file);
	S13Hist           ->SetDirectory(histo_file);
	S13HistData       ->SetDirectory(histo_file);
	S13HistToy        ->SetDirectory(histo_file);
	S23Hist           ->SetDirectory(histo_file);
	S23HistData       ->SetDirectory(histo_file);
	S23HistToy        ->SetDirectory(histo_file);

	TotalPDFHist      ->Write();
	TotalDataHist     ->Write();
	TotalToyHist      ->Write();
	TotalPDFHist_Chi2 ->Write();
	TotalDataHist_Chi2->Write();
	TotalToyHist_Chi2 ->Write();
	S12Hist           ->Write();
	S12HistData       ->Write();
	S12HistToy        ->Write();
	S13Hist           ->Write();
	S13HistData       ->Write();
	S13HistToy        ->Write();
	S23Hist           ->Write();
	S23HistData       ->Write();
	S23HistToy        ->Write();

	for ( int xx = 0; xx < number_of_resonances; ++xx ) {
		histo_sij[ 0 ][ xx ]->Write();
		histo_sij[ 1 ][ xx ]->Write();
		histo_sij[ 2 ][ xx ]->Write();
		histo_sij[ 3 ][ xx ]->Write();
		histo_sij[ 4 ][ xx ]->Write();
	} 

	for ( int i = 0; i < number_of_bkg_components; ++i ) {
		bhisto_sij[ 0 ][ i ]->Write();
		bhisto_sij[ 1 ][ i ]->Write();
		bhisto_sij[ 2 ][ i ]->Write();
		bhisto_sij[ 3 ][ i ]->Write();
		bhisto_sij[ 4 ][ i ]->Write();
	}
	histo_file->Close();

	gStyle->SetPadRightMargin(0.1);

	//g0_12_Data->SetLineColor(2);
	g0_12_Data->SetLineWidth(2);
	g0_12_Data->SetMarkerColor(1);
	g0_12_Data->SetMarkerSize(0.5);
	g0_12_Data->SetMarkerStyle(21);
	g0_12_Data->GetXaxis()->SetTitle("s_{K#pi} [GeV/c^{2}");
	g0_12_Data->GetYaxis()->SetTitle("t_{0}^{0} (arbitrary units)");
	g0_12_Data->GetYaxis()->SetTitleOffset(1.2);
	g0_12_Data->SetTitle("");

	//g1_12_Data->SetLineColor(2);
	g1_12_Data->SetLineWidth(2);
	g1_12_Data->SetMarkerColor(1);
	g1_12_Data->SetMarkerSize(0.5);
	g1_12_Data->SetMarkerStyle(21);
	g1_12_Data->GetXaxis()->SetTitle("s_{K#pi} [GeV/c^{2}");
	g1_12_Data->GetYaxis()->SetTitle("t_{1}^{0} (arbitrary units)");
	g1_12_Data->GetYaxis()->SetTitleOffset(1.2);
	g1_12_Data->SetTitle("");

	//g2_12_Data->SetLineColor(2);
	g2_12_Data->SetLineWidth(2);
	g2_12_Data->SetMarkerColor(1);
	g2_12_Data->SetMarkerSize(0.5);
	g2_12_Data->SetMarkerStyle(21);
	g2_12_Data->GetXaxis()->SetTitle("s_{K#pi} [GeV/c^{2}");
	g2_12_Data->GetYaxis()->SetTitle("t_{2}^{0} (arbitrary units)");
	g2_12_Data->GetYaxis()->SetTitleOffset(1.2);
	g2_12_Data->SetTitle("");

	//g3_12_Data->SetLineColor(2);
	g3_12_Data->SetLineWidth(2);
	g3_12_Data->SetMarkerColor(1);
	g3_12_Data->SetMarkerSize(0.5);
	g3_12_Data->SetMarkerStyle(21);
	g3_12_Data->GetXaxis()->SetTitle("s_{K#pi} [GeV/c^{2}");
	g3_12_Data->GetYaxis()->SetTitle("t_{3}^{0} (arbitrary units)");
	g3_12_Data->GetYaxis()->SetTitleOffset(1.2);
	g3_12_Data->SetTitle("");

	//g4_12_Data->SetLineColor(2);
	g4_12_Data->SetLineWidth(2);
	g4_12_Data->SetMarkerColor(1);
	g4_12_Data->SetMarkerSize(0.5);
	g4_12_Data->SetMarkerStyle(21);
	g4_12_Data->GetXaxis()->SetTitle("s_{K#pi} [GeV/c^{2}");
	g4_12_Data->GetYaxis()->SetTitle("t_{4}^{0} (arbitrary units)");
	g4_12_Data->GetYaxis()->SetTitleOffset(1.2);
	g4_12_Data->SetTitle("");

	//g0_13_Data->SetLineColor(2);
	g0_13_Data->SetLineWidth(2);
	g0_13_Data->SetMarkerColor(1);
	g0_13_Data->SetMarkerSize(0.5);
	g0_13_Data->SetMarkerStyle(21);
	g0_13_Data->GetXaxis()->SetTitle("s_{K#pi} [GeV/c^{2}");
	g0_13_Data->GetYaxis()->SetTitle("t_{0}^{0} (arbitrary units)");
	g0_13_Data->GetYaxis()->SetTitleOffset(1.2);
	g0_13_Data->SetTitle("");

	//g1_13_Data->SetLineColor(2);
	g1_13_Data->SetLineWidth(2);
	g1_13_Data->SetMarkerColor(1);
	g1_13_Data->SetMarkerSize(0.5);
	g1_13_Data->SetMarkerStyle(21);
	g1_13_Data->GetXaxis()->SetTitle("s_{K#pi} [GeV/c^{2}");
	g1_13_Data->GetYaxis()->SetTitle("t_{1}^{0} (arbitrary units)");
	g1_13_Data->GetYaxis()->SetTitleOffset(1.2);
	g1_13_Data->SetTitle("");

	//g2_13_Data->SetLineColor(2);
	g2_13_Data->SetLineWidth(2);
	g2_13_Data->SetMarkerColor(1);
	g2_13_Data->SetMarkerSize(0.5);
	g2_13_Data->SetMarkerStyle(21);
	g2_13_Data->GetXaxis()->SetTitle("s_{K#pi} [GeV/c^{2}");
	g2_13_Data->GetYaxis()->SetTitle("t_{2}^{0} (arbitrary units)");
	g2_13_Data->GetYaxis()->SetTitleOffset(1.2);
	g2_13_Data->SetTitle("");

	//g3_13_Data->SetLineColor(2);
	g3_13_Data->SetLineWidth(2);
	g3_13_Data->SetMarkerColor(1);
	g3_13_Data->SetMarkerSize(0.5);
	g3_13_Data->SetMarkerStyle(21);
	g3_13_Data->GetXaxis()->SetTitle("s_{K#pi} [GeV/c^{2}");
	g3_13_Data->GetYaxis()->SetTitle("t_{3}^{0} (arbitrary units)");
	g3_13_Data->GetYaxis()->SetTitleOffset(1.2);
	g3_13_Data->SetTitle("");

	//g4_13_Data->SetLineColor(2);
	g4_13_Data->SetLineWidth(2);
	g4_13_Data->SetMarkerColor(1);
	g4_13_Data->SetMarkerSize(0.5);
	g4_13_Data->SetMarkerStyle(21);
	g4_13_Data->GetXaxis()->SetTitle("s_{K#pi} [GeV/c^{2}");
	g4_13_Data->GetYaxis()->SetTitle("t_{4}^{0} (arbitrary units)");
	g4_13_Data->GetYaxis()->SetTitleOffset(1.2);
	g4_13_Data->SetTitle("");

	//g0_12_MC->SetLineColor(3);
	//g0_12_MC->SetLineWidth(2);
	//g0_12_MC->SetMarkerColor(1);
	//g0_12_MC->SetMarkerSize(1);
	//g0_12_MC->SetMarkerStyle(2);
	g0_12_MC->SetFillColor(38);
	g0_12_MC->GetXaxis()->SetTitle("s_{K#pi} [GeV/c^{2}");
	g0_12_MC->GetYaxis()->SetTitle("t_{0}^{0} (arbitrary units)");
	g0_12_MC->GetYaxis()->SetTitleOffset(1.2);
	g0_12_MC->SetTitle("");

	//g1_12_MC->SetLineColor(3);
	//g1_12_MC->SetLineWidth(2);
	//g1_12_MC->SetMarkerColor(1);
	//g1_12_MC->SetMarkerSize(1);
	//g1_12_MC->SetMarkerStyle(2);
	g1_12_MC->SetFillColor(38);
	g1_12_MC->GetXaxis()->SetTitle("s_{K#pi} [GeV/c^{2}");
	g1_12_MC->GetYaxis()->SetTitle("t_{1}^{0} (arbitrary units)");
	g1_12_MC->GetYaxis()->SetTitleOffset(1.2);
	g1_12_MC->SetTitle("");

	//g2_12_MC->SetLineColor(3);
	//g2_12_MC->SetLineWidth(2);
	//g2_12_MC->SetMarkerColor(1);
	//g2_12_MC->SetMarkerSize(1);
	//g2_12_MC->SetMarkerStyle(2);
	g2_12_MC->SetFillColor(38);
	g2_12_MC->GetXaxis()->SetTitle("s_{K#pi} [GeV/c^{2}");
	g2_12_MC->GetYaxis()->SetTitle("t_{2}^{0} (arbitrary units)");
	g2_12_MC->GetYaxis()->SetTitleOffset(1.2);
	g2_12_MC->SetTitle("");

	//g3_12_MC->SetLineColor(3);
	//g3_12_MC->SetLineWidth(2);
	//g3_12_MC->SetMarkerColor(1);
	//g3_12_MC->SetMarkerSize(1);
	//g3_12_MC->SetMarkerStyle(2);
	g3_12_MC->SetFillColor(38);
	g3_12_MC->GetXaxis()->SetTitle("s_{K#pi} [GeV/c^{2}");
	g3_12_MC->GetYaxis()->SetTitle("t_{3}^{0} (arbitrary units)");
	g3_12_MC->GetYaxis()->SetTitleOffset(1.2);
	g3_12_MC->SetTitle("");

	//g4_12_MC->SetLineColor(3);
	//g4_12_MC->SetLineWidth(2);
	//g4_12_MC->SetMarkerColor(1);
	//g4_12_MC->SetMarkerSize(1);
	//g4_12_MC->SetMarkerStyle(2);
	g4_12_MC->SetFillColor(38);
	g4_12_MC->GetXaxis()->SetTitle("s_{K#pi} [GeV/c^{2}");
	g4_12_MC->GetYaxis()->SetTitle("t_{4}^{0} (arbitrary units)");
	g4_12_MC->GetYaxis()->SetTitleOffset(1.2);
	g4_12_MC->SetTitle("");

	//g0_13_MC->SetLineColor(3);
	//g0_13_MC->SetLineWidth(2);
	//g0_13_MC->SetMarkerColor(1);
	//g0_13_MC->SetMarkerSize(1);
	//g0_13_MC->SetMarkerStyle(2);
	g0_13_MC->SetFillColor(38);
	g0_13_MC->GetXaxis()->SetTitle("s_{K#pi} [GeV/c^{2}");
	g0_13_MC->GetYaxis()->SetTitle("t_{0}^{0} (arbitrary units)");
	g0_13_MC->GetYaxis()->SetTitleOffset(1.2);
	g0_13_MC->SetTitle("");

	//g1_13_MC->SetLineColor(3);
	//g1_13_MC->SetLineWidth(2);
	//g1_13_MC->SetMarkerColor(1);
	//g1_13_MC->SetMarkerSize(1);
	//g1_13_MC->SetMarkerStyle(2);
	g1_13_MC->SetFillColor(38);
	g1_13_MC->GetXaxis()->SetTitle("s_{K#pi} [GeV/c^{2}");
	g1_13_MC->GetYaxis()->SetTitle("t_{1}^{0} (arbitrary units)");
	g1_13_MC->GetYaxis()->SetTitleOffset(1.2);
	g1_13_MC->SetTitle("");

	//g2_13_MC->SetLineColor(3);
	//g2_13_MC->SetLineWidth(2);
	//g2_13_MC->SetMarkerColor(1);
	//g2_13_MC->SetMarkerSize(1);
	//g2_13_MC->SetMarkerStyle(2);
	g2_13_MC->SetFillColor(38);
	g2_13_MC->GetXaxis()->SetTitle("s_{K#pi} [GeV/c^{2}");
	g2_13_MC->GetYaxis()->SetTitle("t_{2}^{0} (arbitrary units)");
	g2_13_MC->GetYaxis()->SetTitleOffset(1.2);
	g2_13_MC->SetTitle("");

	//g3_13_MC->SetLineColor(3);
	//g3_13_MC->SetLineWidth(2);
	//g3_13_MC->SetMarkerColor(1);
	//g3_13_MC->SetMarkerSize(1);
	//g3_13_MC->SetMarkerStyle(2);
	g3_13_MC->SetFillColor(38);
	g3_13_MC->GetXaxis()->SetTitle("s_{K#pi} [GeV/c^{2}");
	g3_13_MC->GetYaxis()->SetTitle("t_{3}^{0} (arbitrary units)");
	g3_13_MC->GetYaxis()->SetTitleOffset(1.2);
	g3_13_MC->SetTitle("");

	//g4_13_MC->SetLineColor(3);
	//g4_13_MC->SetLineWidth(2);
	//g4_13_MC->SetMarkerColor(1);
	//g4_13_MC->SetMarkerSize(1);
	//g4_13_MC->SetMarkerStyle(2);
	g4_13_MC->SetFillColor(38);
	g4_13_MC->GetXaxis()->SetTitle("s_{K#pi} [GeV/c^{2}");
	g4_13_MC->GetYaxis()->SetTitle("t_{4}^{0} (arbitrary units)");
	g4_13_MC->GetYaxis()->SetTitleOffset(1.2);
	g4_13_MC->SetTitle("");

	mg0_13->Add(g0_13_Data);
	mg0_13->Add(g0_13_MC);

	mg1_13->Add(g1_13_Data);
	mg1_13->Add(g1_13_MC);

	mg2_13->Add(g2_13_Data);
	mg2_13->Add(g2_13_MC);

	mg3_13->Add(g3_13_Data);
	mg3_13->Add(g3_13_MC);

	mg4_13->Add(g4_13_Data);
	mg4_13->Add(g4_13_MC);

	mg0_12->Add(g0_12_Data);
	mg0_12->Add(g0_12_MC);

	mg1_12->Add(g1_12_Data);
	mg1_12->Add(g1_12_MC);

	mg2_12->Add(g2_12_Data);
	mg2_12->Add(g2_12_MC);

	mg3_12->Add(g3_12_Data);
	mg3_12->Add(g3_12_MC);

	mg4_12->Add(g4_12_Data);
	mg4_12->Add(g4_12_MC);

	cout << "-----------------------------------------------------------------------------------------------------------------------------------" << endl;
	cout << endl;
	cout << "Drawing" <<  endl;
	cout << endl;

	TCanvas * c0_13 = new TCanvas();

	TPad *pad0_13_0 = new TPad("pad130_1","top pad",0,0.21,.98,.98);
	pad0_13_0->Draw();
	TPad *pad0_13_1 = new TPad("pad130_2","bottom pad",0,0,.98,0.2);
	pad0_13_1->Draw();       

	pad0_13_0->cd();
	g0_13_MC->Draw("APBX");
	g0_13_Data->Draw("P");

	pad0_13_1->cd();
	Draw_Graph_Pulls(g0_13_Data, g0_13_MC);

	c0_13->SaveAs(MomentsPDFName + "t0_13_DataxMC.pdf");

	TCanvas * c1_13 = new TCanvas();

	TPad *pad1_13_0 = new TPad("pad131_1","top pad",0,0.21,.98,.98);
	pad1_13_0->Draw();
	TPad *pad1_13_1 = new TPad("pad131_2","bottom pad",0,0,.98,0.2);
	pad1_13_1->Draw();       

	pad1_13_0->cd();
	g1_13_MC->Draw("APBX");
	g1_13_Data->Draw("P");

	pad1_13_1->cd();
	Draw_Graph_Pulls(g1_13_Data, g1_13_MC);

	c1_13->SaveAs(MomentsPDFName + "t1_13_DataxMC.pdf");

	TCanvas * c2_13 = new TCanvas();

	TPad *pad2_13_0 = new TPad("pad132_1","top pad",0,0.21,.98,.98);
	pad2_13_0->Draw();
	TPad *pad2_13_1 = new TPad("pad132_2","bottom pad",0,0,.98,0.2);
	pad2_13_1->Draw();       

	pad2_13_0->cd();
	g2_13_MC->Draw("APBX");
	g2_13_Data->Draw("P");

	pad2_13_1->cd();
	Draw_Graph_Pulls(g2_13_Data, g2_13_MC);

	c2_13->SaveAs(MomentsPDFName + "t2_13_DataxMC.pdf");

	TCanvas * c3_13 = new TCanvas();

	TPad *pad3_13_0 = new TPad("pad133_1","top pad",0,0.21,.98,.98);
	pad3_13_0->Draw();
	TPad *pad3_13_1 = new TPad("pad133_2","bottom pad",0,0,.98,0.2);
	pad3_13_1->Draw();       

	pad3_13_0->cd();
	g3_13_MC->Draw("APBX");
	g3_13_Data->Draw("P");

	pad3_13_1->cd();
	Draw_Graph_Pulls(g3_13_Data, g3_13_MC);

	c3_13->SaveAs(MomentsPDFName + "t3_13_DataxMC.pdf");

	TCanvas * c4_13 = new TCanvas();

	TPad *pad4_13_0 = new TPad("pad134_1","top pad",0,0.21,.98,.98);
	pad4_13_0->Draw();
	TPad *pad4_13_1 = new TPad("pad134_2","bottom pad",0,0,.98,0.2);
	pad4_13_1->Draw();       

	pad4_13_0->cd();
	g4_13_MC->Draw("APBX");
	g4_13_Data->Draw("P");

	pad4_13_1->cd();
	Draw_Graph_Pulls(g4_13_Data, g4_13_MC);

	c4_13->SaveAs(MomentsPDFName + "t4_13_DataxMC.pdf");

	TCanvas * c0_12 = new TCanvas();

	TPad *pad0_12_0 = new TPad("pad120_1","top pad",0,0.21,.98,.98);
	pad0_12_0->Draw();
	TPad *pad0_12_1 = new TPad("pad120_2","bottom pad",0,0,.98,0.2);
	pad0_12_1->Draw();       

	pad0_12_0->cd();
	g0_12_MC->Draw("APBX");
	g0_12_Data->Draw("P");

	pad0_12_1->cd();
	Draw_Graph_Pulls(g0_12_Data, g0_12_MC);

	c0_12->SaveAs(MomentsPDFName + "t0_12_DataxMC.pdf");

	TCanvas * c1_12 = new TCanvas();

	TPad *pad1_12_0 = new TPad("pad121_1","top pad",0,0.21,.98,.98);
	pad1_12_0->Draw();
	TPad *pad1_12_1 = new TPad("pad121_2","bottom pad",0,0,.98,0.2);
	pad1_12_1->Draw();       

	pad1_12_0->cd();
	g1_12_MC->Draw("APBX");
	g1_12_Data->Draw("P");

	pad1_12_1->cd();
	Draw_Graph_Pulls(g1_12_Data, g1_12_MC);

	c1_12->SaveAs(MomentsPDFName + "t1_12_DataxMC.pdf");

	TCanvas * c2_12 = new TCanvas();

	TPad *pad2_12_0 = new TPad("pad122_1","top pad",0,0.21,.98,.98);
	pad2_12_0->Draw();
	TPad *pad2_12_1 = new TPad("pad122_2","bottom pad",0,0,.98,0.2);
	pad2_12_1->Draw();       

	pad2_12_0->cd();
	g2_12_MC->Draw("APBX");
	g2_12_Data->Draw("P");

	pad2_12_1->cd();
	Draw_Graph_Pulls(g2_12_Data, g2_12_MC);

	c2_12->SaveAs(MomentsPDFName + "t2_12_DataxMC.pdf");

	TCanvas * c3_12 = new TCanvas();

	TPad *pad3_12_0 = new TPad("pad123_1","top pad",0,0.21,.98,.98);
	pad3_12_0->Draw();
	TPad *pad3_12_1 = new TPad("pad123_2","bottom pad",0,0,.98,0.2);
	pad3_12_1->Draw();       

	pad3_12_0->cd();
	g3_12_MC->Draw("APBX");
	g3_12_Data->Draw("P");

	pad3_12_1->cd();
	Draw_Graph_Pulls(g3_12_Data, g3_12_MC);

	c3_12->SaveAs(MomentsPDFName + "t3_12_DataxMC.pdf");

	TCanvas * c4_12 = new TCanvas();

	TPad *pad4_12_0 = new TPad("pad124_1","top pad",0,0.21,.98,.98);
	pad4_12_0->Draw();
	TPad *pad4_12_1 = new TPad("pad124_2","bottom pad",0,0,.98,0.2);
	pad4_12_1->Draw();       

	pad4_12_0->cd();
	g4_12_MC->Draw("APBX");
	g4_12_Data->Draw("P");

	pad4_12_1->cd();
	Draw_Graph_Pulls(g4_12_Data, g4_12_MC);

	c4_12->SaveAs(MomentsPDFName + "t4_12_DataxMC.pdf");
}


