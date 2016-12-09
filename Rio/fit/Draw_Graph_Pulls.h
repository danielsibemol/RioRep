#include <TGraphErrors.h>

void Draw_Graph_Pulls(TGraphErrors * g1, TGraphErrors * g2){
//cout << " teste 0 " << endl;
	int    nbins = g1->GetN();
	double diff, err, x1, y1, x2, y2, x[nbins], y[nbins]; 
//cout << " teste 1 " << endl;
	//TString pullName = "PullHist";
	//TH1F *PullHist = new TH1F( pullName, "", nbins, min, max);
	for( int i = 0; i < nbins ; ++i ) {
		g1->GetPoint(i, x1, y1);
		g2->GetPoint(i, x2, y2);

		x[i] = x1;

		err = g1->GetErrorY(i); 
		diff = y1 - y2;

		if (diff != 0 && err != 0) y[i] = diff/err;
		else y[i] = 0;
		cout << "y[i] = " << y[i] << ", x1 = " << x1 << ", diff = " << diff << ", err = " << err << endl;

		//if (diff != 0 && err != 0) PullHist->SetBinContent(i+1,diff/err);
                //else PullHist->SetBinContent(i+1,0);  
	}

//cout << " teste 2 " << endl;
	TGraphErrors * PullGraph = new TGraphErrors(nbins, x,y,0,0);
//cout << " teste 4 " << endl;
        PullGraph->SetFillColor(1);
        PullGraph->SetMinimum(-5);
        PullGraph->SetMaximum(5);
	PullGraph->GetYaxis()->SetNdivisions(2, kTRUE);
        PullGraph->GetYaxis()->SetLabelSize(0.2);
        PullGraph->GetYaxis()->SetTitleSize(0.3);
        PullGraph->GetYaxis()->SetTitleOffset(1.1);
        PullGraph->GetYaxis()->SetTitle( "#bf{Pull}" );
	PullGraph->Draw("ABP");

/*        PullHist->SetXTitle( "" );
        PullHist->SetFillColor(kBlue+3);
        PullHist->SetBarWidth(0.75);
        PullHist->SetBarOffset(0.1);
        PullHist->SetStats(0);
        PullHist->SetMinimum(-5);
        PullHist->SetMaximum(5);

        //PullHist->GetXaxis()->SetRangeUser(-5., 5.);
        // PullHist->GetYaxis()->SetTickLength( 2. );
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
*/
}
