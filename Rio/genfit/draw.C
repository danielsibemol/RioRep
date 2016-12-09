



void draw(){

gStyle->SetOptStat(1);
gStyle->SetOptFit(1);

TString output = "stabilitytest.pdf";


// model 3 results
	double f1 = 0.18555,
	       e1 = 0.01136,
	       ef1 = f1 - e1, 
	       df1 = f1 + e1; 
	double f2 = 0.29861, 
	       e2 = 0.019131,
	       ef2 = f2 - e2, 
	       df2 = f2 + e2; 
	double f3 = 0.064888, 
	       e3 = 0.0024469,
	       ef3 = f3 - e3, 
	       df3 = f3 + e3; 

TLine *l1  = new TLine(f1,0,f1,0.25);  // 30.2
TLine *el1 = new TLine(ef1,0,ef1,13);  // 30.2
TLine *dl1 = new TLine(df1,0,df1,13);  // 30.2
       l1->SetLineColor(kRed);
       el1->SetLineColor(kRed);
       dl1->SetLineColor(kRed);
TLine *l2  = new TLine(f2,0,f2,0.18);  // 18.4
TLine *el2 = new TLine(ef2,0,ef2,18);  // 18.4
TLine *dl2 = new TLine(df2,0,df2,18);  // 18.4
       l2->SetLineColor(kRed);
       el2->SetLineColor(kRed);
       dl2->SetLineColor(kRed);
TLine *l3  = new TLine(f3,0,f3,0.18);  // 6.65
TLine *el3 = new TLine(ef3,0,ef3,22);  // 6.65
TLine *dl3 = new TLine(df3,0,df3,22);  // 6.65
       l3->SetLineColor(kRed);
       el3->SetLineColor(kRed);
       dl3->SetLineColor(kRed);

// stability results

TCanvas *c = new TCanvas("c","",1500,500);
c->Print(output+"[");
TFile f("tree.root");
//TFile f("check1_withacc/tree.root");
  TTree *tt =(TTree*)f.Get("tree");

  c->Divide(3,1);
  c->cd(1); tt->Draw("res1_frac >> h1(25,0.145,0.22)","chi2 < 2");
	    h1->SetTitle("a_{0}(1450) Fraction");
   	    h1->SetFillColor(kBlue);
   	    h1->SetFillStyle(3005);
	    h1->DrawNormalized("");
	    l1->Draw("");
	    h1->Fit("gaus");
	    //el1->Draw("");
	    //dl1->Draw("");

  c->cd(2); tt->Draw("res2_frac >> h2(25,0.26,0.34)","chi2 < 2");
	    h2->SetTitle("f_{0}(980) Fraction");
   	    h2->SetFillColor(kBlue);
   	    h2->SetFillStyle(3005);
	    h2->DrawNormalized("");
	    l2->Draw("");
	    h2->Fit("gaus");
	    //el2->Draw("");
	    //dl2->Draw("");

  c->cd(3); tt->Draw("res3_frac >> h3(25,0.06,0.07)","chi2 < 2");
	    h3->SetTitle("#phi (1020) Fraction");
   	    h3->SetFillColor(kBlue);
   	    h3->SetFillStyle(3005);
	    h3->DrawNormalized("");
	    l3->Draw("");
	    h3->Fit("gaus");
	    //el3->Draw("");
	    //dl3->Draw("");
  c->Print(output);

  c->cd(1); tt->Draw("res1_frac_err >> h1(25,0.0,0.02)","chi2 < 2");
	    h1->SetTitle("a_{0}(1450) Fraction Error");
   	    h1->SetFillColor(kBlue);
   	    h1->SetFillStyle(3005);
	    h1->DrawNormalized("");
	    l1->Draw("");
	    h1->Fit("gaus");
	    //el1->Draw("");
	    //dl1->Draw("");

  c->cd(2); tt->Draw("res2_frac_err >> h2(25,0.0,0.03)","chi2 < 2");
	    h2->SetTitle("f_{0}(980) Fraction Error ");
   	    h2->SetFillColor(kBlue);
   	    h2->SetFillStyle(3005);
	    h2->DrawNormalized("");
	    l2->Draw("");
	    h2->Fit("gaus");
	    //el2->Draw("");
	    //dl2->Draw("");

  c->cd(3); tt->Draw("res3_frac_err >> h3(25,0.001,0.006)","chi2 < 2");
	    h3->SetTitle("#phi (1020) Fraction Error");
   	    h3->SetFillColor(kBlue);
   	    h3->SetFillStyle(3005);
	    h3->DrawNormalized("");
	    l3->Draw("");
	    h3->Fit("gaus");
	    //el3->Draw("");
	    //dl3->Draw("");
  c->Print(output);


double fa1 = 13.669/4.,
       fp1 = 44.629,
       fa2 = 13.615/4.,
       fp2 = -71.753; 


TLine *la1  = new TLine(fa1,0,fa1,0.3);  // 30.2
       la1->SetLineColor(kRed);
TLine *lp1  = new TLine(fp1,0,fp1,0.18);  // 30.2
       lp1->SetLineColor(kRed);
TLine *la2  = new TLine(fa2,0,fa2,0.3);  // 30.2
       la2->SetLineColor(kRed);
TLine *lp2  = new TLine(fp2,0,fp2,0.16);  // 30.2
       lp2->SetLineColor(kRed);
TCanvas *cc = new TCanvas("cc","",1500,500);
  cc->Divide(2,1);
  cc->cd(1); tt->Draw("res1_amp/4. >> h1(25,3,4)","chi2 < 2");
	    h1->SetTitle("a_{0}(1450) Amplitude");
   	    h1->SetFillColor(kBlue);
   	    h1->SetFillStyle(3005);
	    h1->DrawNormalized("");
	    h1->Fit("gaus");
	    la1->Draw("");

  cc->cd(2); tt->Draw("res1_phs >> h2(25,30,60)","chi2 < 2");
	    h2->SetTitle("a_{0}(1450) Phase");
   	    h2->SetFillColor(kBlue);
   	    h2->SetFillStyle(3005);
	    h2->DrawNormalized("");
	    h2->Fit("gaus");
	    lp1->Draw("");
  cc->Print(output);

  cc->cd(1); tt->Draw("res1_amp_err/4. >> h1(25,0.05,0.08)","chi2 < 2");
	    h1->SetTitle("a_{0}(1450) Amplitude Error");
   	    h1->SetFillColor(kBlue);
   	    h1->SetFillStyle(3005);
	    h1->DrawNormalized("");
	    h1->Fit("gaus");

  cc->cd(2); tt->Draw("res1_phs_err >> h2(25,3,4)","chi2 < 2");
	    h2->SetTitle("a_{0}(1450) Phase Error");
   	    h2->SetFillColor(kBlue);
   	    h2->SetFillStyle(3005);
	    h2->DrawNormalized("");
	    h2->Fit("gaus");
	    lp1->Draw("");
  cc->Print(output);

  cc->cd(1); tt->Draw("res2_amp/4. >> h1(25,3,4)","chi2 < 2");
	    h1->SetTitle("f_{0}(980) Amplitude");
   	    h1->SetFillColor(kBlue);
   	    h1->SetFillStyle(3005);
	    h1->DrawNormalized("");
	    h1->Fit("gaus");
	    la2->Draw("");

  cc->cd(2); tt->Draw("res2_phs >> h2(25,-78,-66)","chi2 < 2");
	    h2->SetTitle("f_{0}(980) Phase");
   	    h2->SetFillColor(kBlue);
   	    h2->SetFillStyle(3005);
	    h2->DrawNormalized("");
	    h2->Fit("gaus");
	    lp2->Draw("");
  cc->Print(output);

  cc->cd(1); tt->Draw("res2_amp_err/4. >> h1(25,0.05,0.07)","chi2 < 2");
	    h1->SetTitle("f_{0}(980) Amplitude");
   	    h1->SetFillColor(kBlue);
   	    h1->SetFillStyle(3005);
	    h1->DrawNormalized("");
	    h1->Fit("gaus");
	    la2->Draw("");

  cc->cd(2); tt->Draw("res2_phs_err >> h2(25,1,2)","chi2 < 2");
	    h2->SetTitle("f_{0}(980) Phase");
   	    h2->SetFillColor(kBlue);
   	    h2->SetFillStyle(3005);
	    h2->DrawNormalized("");
	    h2->Fit("gaus");
	    lp2->Draw("");
  cc->Print(output);

//  cc->cd(1); tt->Draw("res3_amp >> h1(60,-2,2)","chi2 < 2");
//	    h1->SetTitle("#phi(1020) Amplitude");
//   	    h1->SetFillColor(kBlue);
//   	    h1->SetFillStyle(3005);
//	    h1->DrawNormalized("");
//	    //l1->Draw("");
//	    //el1->Draw("");
//	    //dl1->Draw("");
//  cc->cd(2); tt->Draw("res3_phs >> h2(60,-10,10)","chi2 < 2");
//	    h2->SetTitle("#phi(1020) Phase");
//   	    h2->SetFillColor(kBlue);
//   	    h2->SetFillStyle(3005);
//	    h2->DrawNormalized("");
//	    //l2->Draw("");
//	    //el2->Draw("");
//	    //dl2->Draw("");
//  cc->Print(output);

  cc->Print(output+"]"); 




cout << " Generate " << endl;
cout << " a(a0)  = " << fa1 <<  "  delta(a0) = " << fp1 << endl; 
cout << " a(f0)  = " << fa2 <<  "  delta(f0) = " << fp2 << endl; 


}

