#include <TGraphErrors.h>

void plot_moments (double *m, double *t0, double *t1, double *t2, double *t3, double *t4, double *t0_err, double *t1_err, double *t2_err, double *t3_err, double *t4_err, TGraphErrors *&gr0, TGraphErrors *&gr1, TGraphErrors *&gr2, TGraphErrors *&gr3, TGraphErrors *&gr4){

	int npt = 100;
	double S2[npt], S2_err[npt], P2[npt], P2_err[npt], cosphis[npt], cosphis_err[npt];
	double y0[npt], y1[npt], y2[npt], y3[npt], y4[npt];
	double y0_err[npt], y1_err[npt], y2_err[npt], y3_err[npt], y4_err[npt];
	double S, P, mmax, mmin, dm;
	double y0max, y0min,  y1max, y1min, y2max, y2min, y3max, y3min, y4max, y4min;
	double yg0max, yg0min, yg1max, yg1min, yg2max, yg2min, yg3max, yg3min, yg4max, yg4min;

	// initialize variables that will set the vertical scale of the graphs   
	y0max = 0;
	y0min = 0;
	y1max = 0;
	y1min = 0;
	y2max = 0;
	y2min = 0;
	y3max = 0;
	y3min = 0;
	y4max = 0;
	y4min = 0;

	for (int n=0;n<npt;n++){

		// assume for the moment that very close to threshold there is only S- and P-wave
		P2[n] = 2.5*t2[n];
		S2[n] = t0[n] - P2[n];
		S = sqrt(S2[n]);
		P=1;
		if ((t2[n]>0)) {P = sqrt(P2[n]);}
		cosphis[n] = 1.5*t1[n]/(S*P); 
		S2_err[n]=sqrt(t0_err[n]*t0_err[n] + 2.5*t2_err[n]*2.5*t2_err[n]);
		P2_err[n]=2.5*t2_err[n];
		cosphis_err[n] = 1.5*1.5*t1_err[n];

		// moments and their errors

		y0[n] = t0[n];
		y1[n] = t1[n];
		y2[n] = t2[n];
		y3[n] = t3[n];
		y4[n] = t4[n];

		y0_err[n] = t0_err[n];
		y1_err[n] = t1_err[n];
		y2_err[n] = t2_err[n];
		y3_err[n] = t3_err[n];
		y4_err[n] = t4_err[n];

		// find the largest and smallest value  of the moments, including errors

		if (y0[n]+y0_err[n] > y0max) {y0max = y0[n]+y0_err[n];}
		if (y0[n]-y0_err[n] < y0min) {y0min = y0[n]-y0_err[n];}
		if (y1[n]+y1_err[n] > y1max) {y1max = y1[n]+y1_err[n];}
		if (y1[n]-y1_err[n] < y1min) {y1min = y1[n]-y1_err[n];}
		if (y2[n]+y2_err[n] > y2max) {y2max = y2[n]+y2_err[n];}
		if (y2[n]-y2_err[n] < y2min) {y2min = y2[n]-y2_err[n];}
		if (y3[n]+y3_err[n] > y3max) {y3max = y3[n]+y3_err[n];}
		if (y3[n]-y3_err[n] < y3min) {y3min = y3[n]-y3_err[n];}
		if (y4[n]+y4_err[n] > y4max) {y4max = y4[n]+y4_err[n];}
		if (y4[n]-y4_err[n] < y4min) {y4min = y4[n]-y4_err[n];}

	}

	// set the vertical scale
	yg0max = y0max + (y0max-y0min)*0.1;
	yg0min = y0min - (y0max-y0min)*0.1;   
	yg1max = y1max + (y1max-y1min)*0.1;
	yg1min = y1min - (y1max-y1min)*0.1;
	yg2max = y2max + (y2max-y2min)*0.1;
	yg2min = y2min - (y2max-y2min)*0.1;   
	yg3max = y3max + (y3max-y3min)*0.1;
	yg3min = y3min - (y3max-y3min)*0.1;
	yg4max = y4max + (y4max-y4min)*0.1;
	yg4min = y4min - (y4max-y4min)*0.1;

	// now the horizontal scale   
	dm = (m[npt-1] - m[0]) / npt;
	mmax = m[npt-1] + dm;
	mmin = m[0] - dm;

	gr0 = new TGraphErrors(npt,m,y0,0,y0_err);
	gr1 = new TGraphErrors(npt,m,y1,0,y1_err);
	gr2 = new TGraphErrors(npt,m,y2,0,y2_err);
	gr3 = new TGraphErrors(npt,m,y3,0,y3_err);
	gr4 = new TGraphErrors(npt,m,y4,0,y4_err);

}

