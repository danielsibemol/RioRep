#include <TGraphErrors.h>

void Draw_Graph_Pulls(TGraphErrors * g1, TGraphErrors * g2){
	int    nbins = g1->GetN();
	double diff, err, x1, y1, x2, y2, x[nbins], y[nbins]; 

	for( int i = 0; i < nbins ; ++i ) {
		g1->GetPoint(i, x1, y1);
		g2->GetPoint(i, x2, y2);

		x[i] = x1;

		err = g1->GetErrorY(i); 
		diff = y1 - y2;

		if (diff != 0 && err != 0) y[i] = diff/err;
		else y[i] = 0;

	}

	TGraphErrors * PullGraph = new TGraphErrors(nbins, x,y,0,0);

        PullGraph->SetFillColor(1);
        PullGraph->SetMinimum(-5);
        PullGraph->SetMaximum(5);
	PullGraph->GetYaxis()->SetNdivisions(2, kTRUE);
        PullGraph->GetYaxis()->SetLabelSize(0.2);
        PullGraph->GetYaxis()->SetTitleSize(0.3);
        PullGraph->GetYaxis()->SetTitleOffset(1.1);
        PullGraph->GetYaxis()->SetTitle( "#bf{Pull}" );

	PullGraph->Draw("ABP");

}
