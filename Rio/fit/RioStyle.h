#include "TStyle.h"
#include "TROOT.h"
#include "TPaveText.h"
#include "TText.h"
#include "TLatex.h"
#include <iomanip>
#include <iostream>

using namespace std;

TPaveText *lhcbName = new TPaveText(0.65,0.8,0.9,0.9,"BRNDC");
TText *lhcbLabel = new TText();
TLatex *lhcbLatex = new TLatex();

void RioStyle(){

//  gROOT->Reset();

  ////////////////////////////////////////////////////////////////////
  
  cout << "executing RioStyle.h:" << endl;
  cout << "                      " << endl;
  cout << "                      " << endl;
  cout << "                         $$$$$   $   $$$$        $$   " << endl;
  cout << "                         $   $   $   $  $        $$   " << endl;
  cout << "                         $$$$$   $   $  $------$$$$$$ " << endl;
  cout << "                         $ $     $   $  $        $$   " << endl;
  cout << "                         $  $    $   $$$$        $$   " << endl;
  cout << " " << endl;
  cout << "                         Rio Style " << endl;
  cout << " " << endl;
  cout << " " << endl; 
  cout << " " << endl;



 cout<<" Entering RioStyle "<<endl;
  TStyle *RioStyle = new TStyle("RioStyle", "LHCb official plots style");
  Double_t lhcbWidth = 2;
  Int_t RioFont = 60;
  RioStyle->SetPadColor(0);
  RioStyle->SetCanvasColor(0);
  RioStyle->SetStatColor(0);

  RioStyle->SetMarkerStyle(7);
  RioStyle->SetHistMinimumZero();
  RioStyle->SetOptStat(0);
  RioStyle->SetOptFit(1111111);
  RioStyle->SetPalette(1);

///  //RioStyle->SetNumberContours(99);
///
RioStyle->SetFrameBorderMode(0);
RioStyle->SetCanvasBorderMode(0);
RioStyle->SetPadBorderMode(0);

///
RioStyle->SetPadColor(0);
RioStyle->SetCanvasColor(0);
RioStyle->SetStatColor(0);
RioStyle->SetFillColor(0);
///
////*
/////begin comment
///// set the paper & margin sizes
///RioStyle->SetPaperSize(20,26);
///RioStyle->SetPadTopMargin(0.05);
///RioStyle->SetPadRightMargin(0.05); // increase for colz plots!!
///RioStyle->SetPadBottomMargin(0.16);
///RioStyle->SetPadLeftMargin(0.14);
/////end comment
///*/
///
///// use large fonts
///
//RioStyle->SetTextFont(RioFont);
//RioStyle->SetTextSize(0.08);
///
RioStyle->SetLabelFont(RioFont,"x");
RioStyle->SetLabelFont(RioFont,"y");
//RioStyle->SetLabelFont(RioFont,"z");
//
//RioStyle->SetLabelSize(0.05,"x");
//RioStyle->SetLabelSize(0.05,"y");
//RioStyle->SetLabelSize(0.05,"z");
//RioStyle->SetTitleFont(RioFont);
RioStyle->SetTitleSize(0.05,"x");
RioStyle->SetTitleSize(0.05,"y");
//RioStyle->SetTitleSize(0.06,"z");
///
///
///// use bold lines and markers
RioStyle->SetLineWidth(lhcbWidth);
RioStyle->SetFrameLineWidth(lhcbWidth);
RioStyle->SetHistLineWidth(lhcbWidth);
RioStyle->SetFuncWidth(lhcbWidth);
RioStyle->SetGridWidth(lhcbWidth);
//RioStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes
RioStyle->SetMarkerStyle(7);
RioStyle->SetMarkerSize(1.);

///
///// label offsets
/// // RioStyle->SetLabelOffset(0.015);
///
///// by default, do not display histogram decorations:
RioStyle->SetOptStat(0);  
RioStyle->SetOptTitle(0);
RioStyle->SetOptFit(0);
///
///// look of the statistics box:
///RioStyle->SetStatBorderSize(1);
///RioStyle->SetStatFont(RioFont);
///RioStyle->SetStatFontSize(0.05);
///RioStyle->SetStatX(0.9);
///RioStyle->SetStatY(0.9);
///RioStyle->SetStatW(0.25);
///RioStyle->SetStatH(0.15);
///
///// put tick marks on top and RHS of plots
RioStyle->SetPadTickX(1);
RioStyle->SetPadTickY(1);

gROOT->SetStyle("RioStyle");
gROOT->ForceStyle();
}
