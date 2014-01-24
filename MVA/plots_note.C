// Distributions for Summer 2008 approval analysis
// Authors: R. Gonzalez & G. Gomez-Ceballos
#include "TStyle.h"

void plots_note(int nsel = 0, int ReBin = 1) {
   
  gROOT->SetStyle("Plain");
  setTDRStyle();
  gStyle->SetOptStat(0);
  //gStyle->SetOptTitle(0);
  gStyle->SetErrorX(0);

  TFile *_file0;
  
  char plotName[300];
  if(nsel == 1)       sprintf(plotName,"/before_veto_phi_ee");
  else if(nsel == 2)  sprintf(plotName,"/after_veto_phi_ee");
  else if(nsel == 3)  sprintf(plotName,"/before_veto_invmass_ee");
  else if(nsel == 4)  sprintf(plotName,"/after_veto_invmass_ee");
  else if(nsel == 5)  sprintf(plotName,"/before_veto_met_ee");
  else if(nsel == 6)  sprintf(plotName,"/after_veto_met_ee");
  else if(nsel == 7)  sprintf(plotName,"/before_veto_ptmax_ee");
  else if(nsel == 8)  sprintf(plotName,"/after_veto_ptmax_ee");
  else if(nsel == 9)  sprintf(plotName,"/before_veto_ptmin_ee");
  else if(nsel == 10) sprintf(plotName,"/after_veto_ptmin_ee");
  
  else if(nsel == 11)  sprintf(plotName,"/before_veto_phi_mm");
  else if(nsel == 12)  sprintf(plotName,"/after_veto_phi_mm");
  else if(nsel == 13)  sprintf(plotName,"/before_veto_invmass_mm");
  else if(nsel == 14)  sprintf(plotName,"/after_veto_invmass_mm");
  else if(nsel == 15)  sprintf(plotName,"/before_veto_met_mm");
  else if(nsel == 16)  sprintf(plotName,"/after_veto_met_mm");
  else if(nsel == 17)  sprintf(plotName,"/before_veto_ptmax_mm");
  else if(nsel == 18)  sprintf(plotName,"/after_veto_ptmax_mm");
  else if(nsel == 19)  sprintf(plotName,"/before_veto_ptmin_mm");
  else if(nsel == 20)  sprintf(plotName,"/after_veto_ptmin_mm");
  
  else if(nsel == 21)  sprintf(plotName,"/before_veto_phi_em");
  else if(nsel == 22)  sprintf(plotName,"/after_veto_phi_em");
  else if(nsel == 23)  sprintf(plotName,"/before_veto_invmass_em");
  else if(nsel == 24)  sprintf(plotName,"/after_veto_invmass_em");
  else if(nsel == 25)  sprintf(plotName,"/before_veto_met_em");
  else if(nsel == 26)  sprintf(plotName,"/after_veto_met_em");
  else if(nsel == 27)  sprintf(plotName,"/before_veto_ptmax_em");
  else if(nsel == 28)  sprintf(plotName,"/after_veto_ptmax_em");
  else if(nsel == 29)  sprintf(plotName,"/before_veto_ptmin_em");
  else if(nsel == 30)  sprintf(plotName,"/after_veto_ptmin_em");
  
  else if(nsel == 31) sprintf(plotName,"/histo_tmva_130");
  else if(nsel == 32) sprintf(plotName,"/histo_tmva_170");
  
  else if(nsel == 33) sprintf(plotName,"/BDT_WW");
  else if(nsel == 34) sprintf(plotName,"/BDT_tt");
 
  else if(nsel == 35) sprintf(plotName,"/histo_allcutsbut_deltaphi");
  else if(nsel == 36) sprintf(plotName,"/histo_allcutsbut_mass");
  else if(nsel == 37) sprintf(plotName,"/histo_allcutsbut_met");
  else if(nsel == 38) sprintf(plotName,"/histo_allcutsbut_ptmax");
  else if(nsel == 39) sprintf(plotName,"/histo_allcutsbut_ptmin");

  char myRootFile[300];
  sprintf(myRootFile,"rootfiles%s.root",plotName);
  TFile *_file0 = TFile::Open(myRootFile);

  char XTitle[300];
  if(nsel == 1  || nsel == 2 ||  nsel == 11 || nsel == 12 ||
     nsel == 21 || nsel == 22 || nsel == 35) sprintf(XTitle,"#Delta#Phi_{ll} [dg.]");
     
  if(nsel == 3  || nsel == 4 || nsel == 13 || nsel == 14 ||
     nsel == 23 || nsel == 24 || nsel == 36) sprintf(XTitle,"m_{ll} [GeV/c^{2}]");
     
  if(nsel == 5  || nsel == 6 || nsel == 15 || nsel == 16 ||
     nsel == 25 || nsel == 26 || nsel == 37) sprintf(XTitle,"Missing E_{T} [GeV/c]");
  
  if(nsel == 7  || nsel == 8 || nsel == 17 || nsel == 18 ||
     nsel == 27 || nsel == 28 || nsel == 38) sprintf(XTitle,"P_{T}max [GeV/c]");
  
  if(nsel == 9  || nsel == 10 || nsel == 19 || nsel == 20 ||
     nsel == 29 || nsel == 30 || nsel == 39) sprintf(XTitle,"P_{T}min [GeV/c]");
     
  if(nsel == 31  || nsel == 32) sprintf(XTitle,"Neural Network Output");
  if(nsel == 33  || nsel == 34) sprintf(XTitle,"BDT Output");

  char PChannel[300];
  sprintf(PChannel," ");
  if(nsel >= 1 && nsel <= 10)  sprintf(PChannel,"  e^{+}e^{-} Channel",PChannel);
  if(nsel >= 11 && nsel <= 20) sprintf(PChannel,"  #mu^{+}#mu^{-} Channel",PChannel);
  if(nsel >= 21 && nsel <= 30) sprintf(PChannel,"  e^{+/-}#mu^{-\+} Channel",PChannel);

  TCanvas* c1 = new TCanvas();
  c1->SetLogy() ;
  
  if(nsel < 40){   
    TH1F* histo_4 =  histo4->Clone();
    histo3->Add(histo_4);
    TH1F* histo_3 =  histo3->Clone();
    histo2->Add(histo_3);
    TH1F* histo_2 =  histo2->Clone();
    histo1->Add(histo_2);

    histo1->Rebin(ReBin);
    histo1->SetFillColor(kBlue);
    histo1->SetFillStyle(1001);
    histo1->SetLineStyle(0);
    histo1->SetLineWidth(0);

    histo2->Rebin(ReBin);
    histo2->SetFillColor(kMagenta);
    histo2->SetFillStyle(1001);
    histo2->SetLineStyle(0);
    histo2->SetLineWidth(0);

    histo3->Rebin(ReBin);
    histo3->SetFillColor(kGreen);
    histo3->SetFillStyle(1001);
    histo3->SetLineStyle(0);
    histo3->SetLineWidth(0);
    
    histo4->Rebin(ReBin);
    histo4->SetFillColor(kCyan);
    histo4->SetFillStyle(1001);
    histo4->SetLineStyle(0);
    histo4->SetLineWidth(0);

    histo0->Sumw2();
    histo0->Rebin(ReBin);
    histo0->SetMarkerStyle(20);
    histo0->SetMarkerSize(1);
   
    histo0->SetYTitle("normalized events");
    histo1->SetTitleSize(0.05, "Y");
    histo1->GetYaxis()->SetTitleFont(62);

    histo1->SetTitleSize(0.05, "X");
    histo1->GetXaxis()->SetTitleFont(62);
    histo1->GetXaxis()->SetLabelFont(61);
    histo1->GetYaxis()->SetLabelFont(61); 
    histo1->GetYaxis()->SetTitleOffset(1.1);
    histo1->SetLabelSize(0.04, "Y");
    histo1->SetLabelSize(0.04, "X");

    char YTitle[300];
sprintf(YTitle,"events / bin");
   /* if(nsel == 1  || nsel == 2 || nsel == 11 || nsel == 12 ||
       nsel == 21 || nsel == 22 ) sprintf(YTitle," events / 6 dg.");

    if(nsel == 35) sprintf(YTitle," events / 10 dg.");

    if(nsel == 3  || nsel == 4 || nsel == 13 || nsel == 14 ||
       nsel == 23 || nsel == 24) sprintf(YTitle," events / 8 GeV/c^{2}");
    if(nsel == 36) sprintf(YTitle,"  events / 5 GeV/c^{2}");

    if(nsel == 5  || nsel == 6 || nsel == 15 || nsel == 16 ||
       nsel == 25 || nsel == 26) sprintf(YTitle," events / 8 GeV/c");

    if(nsel == 37 || nsel == 38 || nsel == 39) sprintf(YTitle,"  events / 5 GeV/c");


    if(nsel == 7  || nsel == 8 || nsel == 17 || nsel == 18 ||
       nsel == 27 || nsel == 28) sprintf(YTitle," events / 8 GeV/c");

    if(nsel == 9  || nsel == 10 || nsel == 19 || nsel == 20 ||
       nsel == 29 || nsel == 30) sprintf(YTitle,"events / 8 GeV/c");

    if(nsel == 31  || nsel == 32) sprintf(YTitle,"events / bin");
    if(nsel == 33  || nsel == 34) sprintf(YTitle,"events / bin");*/

    histo1->SetYTitle(YTitle);
    histo2->SetYTitle(YTitle);
    histo3->SetYTitle(YTitle);
    histo4->SetYTitle(YTitle);

    histo0->SetXTitle(XTitle);
    histo1->SetXTitle(XTitle);
    histo2->SetXTitle(XTitle);
    histo3->SetXTitle(XTitle);
    histo4->SetXTitle(XTitle);
     
    double max = histo1->GetMaximum();
    histo1->SetMaximum(max*10);
    
  
    //histo1->SetMinimum(0.0);
    histo1->SetMinimum(0.1);
    histo1->Draw("hist,");
    histo2->Draw("hist,same");
    histo3->Draw("hist,same");
    histo4->Draw("hist,same");
    histo0->Draw("E, same");
    cout << histo0->GetSumOfWeights() << " "
         << histo1->GetSumOfWeights() << " "
         << histo2->GetSumOfWeights() << " "
         << histo3->GetSumOfWeights() << " "
         << histo4->GetSumOfWeights() << endl;
         

  } else {
    printf("Wrong option: %d\n",nsel);
    return;
  }

  labelcms  = new TPaveText(0.17,0.90,0.18,0.91,"NDCBR");
  labelcms->SetTextAlign(12);
  labelcms->SetTextSize(0.05);
  labelcms->SetFillColor(kWhite);
  labelcms->AddText("CMS Preliminary");
  labelcms->SetBorderSize(0);
  labelcms->Draw();

  TLegend* leg;
  
  if(nsel == 3  || nsel == 4 || nsel == 13 || nsel == 14 ||
     nsel == 23 || nsel == 24) leg = new TLegend(0.60,0.68,0.95,0.88);

  else if(nsel == 1  || nsel == 2 || nsel == 11 || nsel == 12 ||
          nsel == 21 || nsel == 22) leg = new TLegend(0.18,0.68,0.51,0.87);
     
  else if(nsel == 31 || nsel == 32) leg = new TLegend(0.60,0.70,0.95,0.88);
     
  else  leg = new TLegend(0.60,0.68,0.95,0.92);
  
  leg ->SetFillStyle(0);
  leg ->SetFillColor(kWhite);
  leg ->SetBorderSize(0);
  leg ->SetTextSize(0.035);
  
  if(nsel < 31 || nsel > 32){
    leg ->AddEntry(histo0,"Signal, m_{H}=160 GeV");
    leg ->AddEntry(histo4,"W+Jets, tW","F"); 
    leg ->AddEntry(histo3,"di-boson","F");  
    leg ->AddEntry(histo2,"t#bar{t}","F"); 
    leg ->AddEntry(histo1,"Drell-Yan","F");
  }
  if(nsel == 31){
    leg ->AddEntry(histo0,"Signal, m_{H}=130 GeV");
    leg ->AddEntry(histo4,"W+Jets, tW","F"); 
    leg ->AddEntry(histo3,"di-boson","F");  
    leg ->AddEntry(histo2,"t#bar{t}","F"); 
    leg ->AddEntry(histo1,"Drell-Yan","F"); 
  }
  if(nsel == 32){
    leg ->AddEntry(histo0,"Signal, m_{H}=170 GeV");
    leg ->AddEntry(histo4,"W+Jets, tW","F"); 
    leg ->AddEntry(histo3,"di-boson","F");  
    leg ->AddEntry(histo2,"t#bar{t}","F"); 
    leg ->AddEntry(histo1,"Drell-Yan","F"); 
  }

  leg ->Draw();
  
  if (nsel == 31 || nsel == 32){
    float xpos =0;
    if (nsel == 31) xpos = 0.2;
    else xpos = -0.25;
    l=new TLine(xpos,0.1,xpos,115);
    l->SetLineWidth(4);
    //l->Draw();
  }
  
  bool labelch = false;
  if (nsel < 33 && nsel > 0)labelch = true;
  
  if (labelch){
  if(nsel == 3  || nsel == 4 || nsel == 13 || nsel == 14 ||
     nsel == 23 || nsel == 24)  labelchannel  = new TPaveText(0.58,0.66,0.83,0.62,"NDCBR");
  
  else if(nsel == 1  || nsel == 2 || nsel == 11 || nsel == 12 ||
          nsel == 21 || nsel == 22)  labelchannel  = new TPaveText(0.58,0.88,0.65,0.90,"NDCBR");
	  
  else labelchannel  = new TPaveText(0.58,0.66,0.83,0.62,"NDCBR");
  
  labelchannel->SetTextAlign(12);
  labelchannel->SetTextSize(0.05);
  labelchannel->SetFillColor(0);
  labelchannel->SetFillStyle(0);
  labelchannel->AddText(PChannel);
  labelchannel->SetBorderSize(0);
  labelchannel->Draw();
  
 }

  char myOutputFile[300];

  sprintf(myOutputFile,"plots%s.eps",plotName);
  c1->SaveAs(myOutputFile);
  sprintf(myOutputFile,"plots%s.gif",plotName);
  c1->SaveAs(myOutputFile);
}

void setTDRStyle() {
  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

  tdrStyle->SetCanvasBorderMode(0);
  tdrStyle->SetCanvasColor(kWhite);
  tdrStyle->SetCanvasDefH(600);
  tdrStyle->SetCanvasDefW(600);
  tdrStyle->SetCanvasDefX(0);
  tdrStyle->SetCanvasDefY(0);

  tdrStyle->SetPadBorderMode(0);
  tdrStyle->SetPadColor(kWhite);

  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameBorderSize(1);
  tdrStyle->SetFrameFillColor(0);
  tdrStyle->SetFrameFillStyle(0);
  tdrStyle->SetFrameLineColor(1);
  tdrStyle->SetFrameLineStyle(1);
  tdrStyle->SetFrameLineWidth(1);

  tdrStyle->SetHistLineColor(1);
  tdrStyle->SetHistLineStyle(0);
  tdrStyle->SetHistLineWidth(1);

  tdrStyle->SetEndErrorSize(2);
  tdrStyle->SetErrorX(0.);
  
  tdrStyle->SetMarkerStyle(20);
  tdrStyle->SetMarkerSize(2);

  tdrStyle->SetOptFit(1);
  tdrStyle->SetFitFormat("5.4g");
  tdrStyle->SetFuncColor(2);
  tdrStyle->SetFuncStyle(1);
  tdrStyle->SetFuncWidth(1);

  tdrStyle->SetOptDate(0);

  tdrStyle->SetOptFile(0);
  tdrStyle->SetOptStat(0); 
  tdrStyle->SetStatColor(kWhite);
  tdrStyle->SetStatFont(42);
  tdrStyle->SetStatFontSize(0.025);
  tdrStyle->SetStatTextColor(1);
  tdrStyle->SetStatFormat("6.4g");
  tdrStyle->SetStatBorderSize(1);
  tdrStyle->SetStatH(0.1);
  tdrStyle->SetStatW(0.15);

  tdrStyle->SetPadTopMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.13);
  tdrStyle->SetPadLeftMargin(0.16);
  tdrStyle->SetPadRightMargin(0.027);

  tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.05);

  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(20, "XYZ");
  //tdrStyle->SetTitleSize(0.06, "XYZ");
  tdrStyle->SetTitleXSize(0.02); 
  tdrStyle->SetTitleYSize(0.02);
  tdrStyle->SetTitleXOffset(1.5);
  tdrStyle->SetTitleYOffset(1.7);
  // tdrStyle->SetTitleOffset(1.1, "XYZ"); // Another way to set the Offset

  // For the axis labels:

  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(18, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "XYZ");
  tdrStyle->SetLabelSize(0.04, "XYZ");

  tdrStyle->SetAxisColor(1, "XYZ");
  tdrStyle->SetStripDecimals(kTRUE);
  tdrStyle->SetTickLength(0.03, "XYZ");
  tdrStyle->SetNdivisions(510, "XYZ");
  tdrStyle->SetPadTickX(1);  
  tdrStyle->SetPadTickY(1);

  tdrStyle->SetOptLogx(0);
  tdrStyle->SetOptLogy(0);
  tdrStyle->SetOptLogz(0);

  tdrStyle->SetPaperSize(20.,20.);

  tdrStyle->cd();
}
