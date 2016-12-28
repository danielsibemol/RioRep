#include <stdio.h>
#include <string>

using namespace std;

// Trick to avoid warnings for variables that are not used.

#define UNUSED(p)  ((void)(p))

// Defines how many resonant channels and background components can be read from .txt file 

#define AVAILABLE_RESONANCES 36
#define AVAILABLE_BKG_COMPONENTS 4

int const toy = 1;

// Defines identification strings for resonant channels and background components

// s_ stands for simmetric amplitude.

string const resonant_channel_string[AVAILABLE_RESONANCES] = { // k k pi amplitudes
							      "Kstar892+K",                         // 0 
							      "K0star1430+K",                       // 1
							      "phi+pi",                             // 2 
							      "a0+pi",                              // 3 
							      "kappa+K",                            // 4
	                                                      "NR",                                 // 5
							      "KK_S-Wave",                          // 6
							      "f02KK+pi",                           // 7 
                                                              "f0X2KK+pi",                          // 8
							      "K2+K",                               // 9
                                                              "phi1680+pi",			    // 10 
							      "kappaBW+K",                          // 11
							      "f01500+pi",                          // 12 
							      "K1+K",                               // 13 
                                                              "a2+pi",                              // 14
                                                              "f2+pi",                              // 15
							      "Kstar1680+K",                        // 16
							      "f01710+pi",                          // 17 
							      "f2p+pi",                             // 18
							      "Kpi_S-Wave",                         // 19
							      // pi pi pi amplitudes		    //  
							      "s_rho770+pi",                        // 20
							      "s_f0980+pi",                         // 21 
						              "s_f21270+pi",                        // 22 
							      "s_rho1450+pi",                       // 23
							      "s_f0X+pi",                           // 24
							      "s_sigma+pi",                         // 25
							      "s_KK_S-Wave",                        // 26 
							      "s_f01500+pi",                        // 27
							      "s_sigmaBW+pi",                       // 28
							      "s_f0980_BW+pi",			    // 29
							      // k k k  amplitudes                  //   
							      "s_a01450+K",			    // 30 
							      "s_f0980+K",                          // 31
							      "s_f21270+K",                         // 32 
							      "s_f0X+K",                            // 33 
							      "s_phi+K",                            // 34
						              "s_f01500+K"                          // 35
							     };                              

string const resonant_channel_string_tex[AVAILABLE_RESONANCES] = {// k k pi amplitudes
								  "K*(892)+K",			   // 0 
							          "K_{0}*(1430)+K",                // 1
								  "#phi (1020)+#pi",               // 2 
                                                                  "a_{0}(1450)+#pi",               // 3 
                                                                  "#kappa +K",                     // 4
                                                                  "NR",                            // 5
								  "KK S-Wave",                     // 6
								  "f_{0}(980)+#pi",                // 7 
								  "f_{0}(X)+#pi",                  // 8
							          "K_{2}+K",                       // 9
								  "#phi (1680)+#pi",               // 10 
								  "#kappa (BW)+K",                 // 11
 								  "f_{0}(1500)+#pi",               // 12 
								  "K_{1}(1410)+K",                 // 13 
								  "a_{2}(1320)+#pi",               // 14
								  "f_{2}(1270)+#pi",               // 15
								  "K*(1680)+K",                    // 16
								  "f_{0}(1710)+#pi",               // 17 
								  "f_{2}'(1525)+#pi",              // 18
								  "K #pi S-Wave",                  // 19
								  // pi pi pi amplitudes           //  
                                                                  "#rho(770) +#pi",                // 20
                                                                  "f_{0}(980)+#pi",                // 21 
                                                                  "f_{2}(1270)+#pi",               // 22 
								  "#rho (1450)+#pi",               // 23 
                                                                  "f_{0}(X)+#pi",                  // 24
								  "#sigma + #pi",                  // 25
								  "#pi #pi S-Wave",                // 26
 								  "f_{0}(1500)+#pi",               // 27 
								  "#sigma BW + #pi",               // 28
								  "f_{0}(980)(BW)+#pi",            // 29
								  // k k k amplitudes              //
								  "a_{0}(1450)+K",                 // 30  
								  "f_{0}(980)+K",                  // 31 
								  "f_{2}(1270)+K",                 // 32
								  "f_{0}(X)+K",                    // 33 
								  "#phi +K",                       // 34 
								  "f_{0}(1500)+K"};                // 35

string const bkg_component_string[AVAILABLE_BKG_COMPONENTS] = {"Random_phi","Random_Kstar","COMBINATORIAL","HISTOGRAM"};

int const PWA_INDEX = 12;

// Standard values used in Amplitudes and Generator

double const PI = 3.14159265358979323846;

double const D_Mass = 1.86962;

double const D_width = 0.0096;

double const mK = 0.493677;

double const mpi = 0.13956995;

double const Rekappa = 0.71;

double const Imkappa = - 0.31;

double const m0_K892 = 0.89581;

double const m0_K1430 = 1.425;

double const m0_K11410 = 1.414;

double const m0_K21430 = 1.4324;

double const m0_K1680 = 1.717;

double const m0_phi = 1.019461;

double const m0_a0 = 1.474;

double const m0_a2 = 1.3183;

double const w0_K892 = 0.0474;

double const w0_K1430 = 0.270;

double const w0_K11410 = 0.232;

double const w0_K21430 = 0.109;

double const w0_K1680 = 0.322;

double const w0_phi = 0.004266;

double const w0_a0 = 0.265;

double const w0_a2 = 0.107;

double const m0_kappabugg = 3.3;

double const m0_kappa = 0.797;

double const w0_kappa = 0.410;

double const m0_etap = 0.95778;

double const m0_K2 = 1.4324;

double const w0_K2 = 0.109;

double const m0_phi1680 = 1.680;

double const w0_phi1680 = 0.150;

double const rD2 = 25.0;

double const rR2 = 2.25;

double const m0_f0 = 0.965;

double const w0_f0 = 0.07;	

double const gpi = 0.165;

double const gK = 0.695;

double const Resigma = 0.441;

double const Imsigma = -0.272;

double const m0_rho = 0.769;

double const m0_omega = 0.78265;;

double const m0_f2 = 1.2755;

double const m0_f2p = 1.525;

double const m0_f0X = 1.435;

double const m0_rho2 = 1.465;

double const w0_rho = 0.150;

double const w0_omega = 0.00849;

double const w0_f2 = 0.1867;

double const w0_f2p = 0.073;

double const w0_f0X = 0.135;

double const w0_rho2 = 0.310;

double const m0_f01500 = 1.505;

double const w0_f01500 = 0.109;

double const m0_f01710 = 1.722;

double const w0_f01710 = 0.135;

double const mKsq = 0.24371698;

double const mK0sq = 0.247677419;

double const area_Dp_kkpi = 1.68542042;

double const area_Dp_3pi = 2.70431;

// modified kappa constants

double const b1 = 24.49;

double const b2 = 0.0;

double const A = 2.5;

double const sA = 0.234;

// Flatté K*1430 constants

double const g_Kpi = 0.284;
//double const g_Kpi = 0.304;

double const g_Ketap = 0.22362182362;
//double const g_Ketap = 0.;

double const g2Kpi = 2.5275;

//double const g2Ketap = 0.085*g2Kpi;
double const g2Ketap = 0.;

double const alpha  = 0.5661;

double alpha2 = 4.5;

double const sth = (mK+mpi)*(mK+mpi);

double MKstar = 1.33;

double const C = -0.001;

double const Fkpi = 0.102722;

double const_cm =  0.001;

double const_cd = 0.0352;

double const sA2 = (m0_etap*m0_etap - mK*mK)/2;

double const resonant_channel_mass[AVAILABLE_RESONANCES] = { // k k pi
							    m0_K892,
							    m0_K1430, 
							    m0_phi, 
							    m0_a0, 
							    Rekappa, 
						            0, 
							    0, 
							    m0_f0,
							    m0_f0X,
							    m0_K2,
						            m0_phi1680,
							    m0_kappa,
							    m0_f01500,
							    m0_K11410,
						            m0_a2,	
							    m0_f2,
							    m0_K1680,
							    m0_f01710,
							    m0_f2p,
							    MKstar,
							    // pi pi pi 
							    m0_rho,
							    m0_f0,
							    m0_f2,
							    m0_rho2,
					                    m0_f0X,
							    Resigma,
							    0,
							    m0_f01500,
							    Resigma, // ??
							    m0_f0,
							    // k k k 
							    m0_a0,
							    m0_f0,
							    m0_f2,
							    m0_f0X,
							    m0_phi,
							    m0_f01500};
							   

double const resonant_channel_width[AVAILABLE_RESONANCES] = { // k k pi
							     w0_K892,
							     w0_K1430, 
							     w0_phi, 
							     w0_a0, 
							     Imkappa, 
						             100, 
							     100, 
							     w0_f0,
							     w0_f0X,
							     w0_K2,
						             w0_phi1680,
							     w0_kappa,
							     w0_f01500,
							     w0_K11410,
						             w0_a2,	
							     w0_f2,
							     w0_K1680,
							     w0_f01710,
							     w0_f2p,
							     0,
							     // pi pi pi 
							     w0_rho,
							     w0_f0,
							     w0_f2,
							     w0_rho2,
					                     w0_f0X,
							     Resigma,
							     0,
							     w0_f01500,
							     Imsigma, // ??
							     w0_f0,
							     // k k k 
							     w0_a0,
							     w0_f0,
							     w0_f2,
							     w0_f0X,
							     w0_phi,
							     w0_f01500};

