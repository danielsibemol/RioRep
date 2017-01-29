#include <TGraphErrors.h>

void plot_moments (double *m, double *t0, double *t1, double *t2, double *t3, double *t4, double *t0_err, double *t1_err, double *t2_err, double *t3_err, double *t4_err, TGraphErrors *&gr0, TGraphErrors *&gr1, TGraphErrors *&gr2, TGraphErrors *&gr3, TGraphErrors *&gr4){

	int npt = 100;

	gr0 = new TGraphErrors(npt,m,t0,0,t0_err);
	gr1 = new TGraphErrors(npt,m,t1,0,t1_err);
	gr2 = new TGraphErrors(npt,m,t2,0,t2_err);
	gr3 = new TGraphErrors(npt,m,t3,0,t3_err);
	gr4 = new TGraphErrors(npt,m,t4,0,t4_err);

}

