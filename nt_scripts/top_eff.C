#include "/home/ceballos/HiggsMVA/NeuralNetworkMaker-3x/factors.h"

void top_eff(int nsel = 0, int option = 1, char myRootFile[300] = "histo_njets.root"){
  Char_t xTitle[]="N_{jets}"; Char_t yTitle[]="efficiency";
  if(nsel == 5) {sprintf(xTitle,"p_{T}^{jet} [GeV/c]"); sprintf(yTitle,"Fraction");}
  if(nsel == 6) {sprintf(xTitle,"|#eta^{jet}|"); sprintf(yTitle,"Fraction");}
  if(nsel == 7) {sprintf(xTitle,"|#eta^{jet}|");}
  if(nsel == 8) {sprintf(xTitle,"|#eta^{jet}|");sprintf(yTitle,"Events/fb^{-1}");}
  TFile *_file0 = TFile::Open(myRootFile);
  _file0->cd();
  char sb[50];
  TH1D* hsww_nj[65];
  TH1D* hsww_ptb[60];
  TH1D* hsww_etab[60];
  for(int nj=0; nj<65; nj++){
    sprintf(sb,"hDnjets_%d",nj);
    hsww_nj[nj]  = (TH1D*) gROOT->FindObject(sb);
  }
  for(int nj=0; nj<60; nj++){
    sprintf(sb,"hDpttag_%d",nj);
    hsww_ptb[nj]  = (TH1D*) gROOT->FindObject(sb);
    sprintf(sb,"hDetatag_%d",nj);
    hsww_etab[nj]  = (TH1D*) gROOT->FindObject(sb);
  }
  atributes(hsww_ptb[20],xTitle,1,yTitle);
  atributes(hsww_ptb[21],xTitle,1,yTitle);
  atributes(hsww_ptb[22],xTitle,2,yTitle);
  atributes(hsww_ptb[23],xTitle,2,yTitle);
  atributes(hsww_ptb[24],xTitle,4,yTitle);  
  atributes(hsww_ptb[25],xTitle,4,yTitle);  
  atributes(hsww_ptb[26],xTitle,6,yTitle);  
  atributes(hsww_ptb[27],xTitle,6,yTitle);  
  atributes(hsww_ptb[28],xTitle,8,yTitle);  
  atributes(hsww_ptb[29],xTitle,8,yTitle);  
  atributes(hsww_etab[20],xTitle,1,yTitle);
  atributes(hsww_etab[21],xTitle,1,yTitle);
  atributes(hsww_etab[22],xTitle,2,yTitle);
  atributes(hsww_etab[23],xTitle,2,yTitle);
  atributes(hsww_etab[24],xTitle,4,yTitle);  
  atributes(hsww_etab[25],xTitle,4,yTitle);  
  atributes(hsww_etab[26],xTitle,6,yTitle);  
  atributes(hsww_etab[27],xTitle,6,yTitle);  
  atributes(hsww_etab[28],xTitle,8,yTitle);  
  atributes(hsww_etab[29],xTitle,8,yTitle);  

  // signal tagging
  hsww_nj[21]->Divide(hsww_nj[20]); // low pt jet tagging
  hsww_nj[22]->Divide(hsww_nj[20]); // soft muon tagging
  hsww_nj[23]->Divide(hsww_nj[20]); // low pt jet and soft muon tagging
  hsww_nj[24]->Divide(hsww_nj[20]); // low pt jet and soft muon tagging, not in jets
  hsww_nj[25]->Divide(hsww_nj[20]); // total jet tagging
  hsww_nj[27]->Divide(hsww_nj[26]); // low pt jet and soft muon tagging, not in jets, with leading jet tagging
  hsww_nj[26]->Divide(hsww_nj[20]); // leading jet tagging
  hsww_nj[29]->Divide(hsww_nj[28]); // leading jet tagging with no low pt tagging
  hsww_nj[31]->Divide(hsww_nj[30]); // low pt jet and soft muon tagging if no high pt tagging
  atributes(hsww_nj[21],xTitle,1,yTitle);
  atributes(hsww_nj[22],xTitle,2,yTitle);
  atributes(hsww_nj[23],xTitle,4,yTitle);
  atributes(hsww_nj[24],xTitle,6,yTitle);
  atributes(hsww_nj[25],xTitle,8,yTitle);
  atributes(hsww_nj[27],xTitle,9,yTitle);
  atributes(hsww_nj[26],xTitle,11,yTitle);
  atributes(hsww_nj[29],xTitle,25,yTitle);
  atributes(hsww_nj[31],xTitle,13,yTitle);
  // data-background tagging
  hsww_nj[41]->Divide(hsww_nj[40]); // low pt jet tagging
  hsww_nj[42]->Divide(hsww_nj[40]); // soft muon tagging
  hsww_nj[43]->Divide(hsww_nj[40]); // low pt jet and soft muon tagging
  hsww_nj[44]->Divide(hsww_nj[40]); // low pt jet and soft muon tagging, not in jets
  hsww_nj[45]->Divide(hsww_nj[40]); // total jet tagging
  hsww_nj[47]->Divide(hsww_nj[46]); // low pt jet and soft muon tagging, not in jets, with leading jet tagging
  hsww_nj[46]->Divide(hsww_nj[40]); // leading jet tagging
  hsww_nj[49]->Divide(hsww_nj[48]); // leading jet tagging with no low pt tagging
  hsww_nj[51]->Divide(hsww_nj[50]); // low pt jet and soft muon tagging if no high pt tagging
  atributes(hsww_nj[41],xTitle,1,yTitle);
  atributes(hsww_nj[42],xTitle,2,yTitle);
  atributes(hsww_nj[43],xTitle,4,yTitle);
  atributes(hsww_nj[44],xTitle,6,yTitle);
  atributes(hsww_nj[45],xTitle,8,yTitle);
  atributes(hsww_nj[47],xTitle,9,yTitle);
  atributes(hsww_nj[46],xTitle,11,yTitle);
  atributes(hsww_nj[49],xTitle,20,yTitle);
  atributes(hsww_nj[51],xTitle,13,yTitle);
  // vbf tagging
  hsww_nj[61]->Divide(hsww_nj[60]); // low pt jet and soft muon tagging
  hsww_nj[63]->Divide(hsww_nj[62]); // low pt jet and soft muon tagging if no high pt tagging
  hsww_nj[64]->Divide(hsww_nj[60]); // total jet tagging
  atributes(hsww_nj[61],xTitle,6,yTitle);
  atributes(hsww_nj[63],xTitle,2,yTitle);
  atributes(hsww_nj[64],xTitle,4,yTitle);

  for(int nj=0; nj<65; nj++){
    hsww_nj[nj]->SetMinimum(0.0);
    hsww_nj[nj]->SetMaximum(1.3);
  }

  if     (nsel == 0){
    hsww_nj[21]->Draw();
    hsww_nj[22]->Draw("same");
    hsww_nj[24]->Draw("same");
    //hsww_nj[27]->Draw("same");
    TLegend* leg = new TLegend(0.2,0.75,0.93,0.9);					      
    leg->SetFillColor(10);								      
    leg->SetTextSize(0.023);
    leg->AddEntry(hsww_nj[21] ,"low p_{T} jet tagging","l");	   
    leg->AddEntry(hsww_nj[22] ,"soft muon tagging","l"); 	   
    leg->AddEntry(hsww_nj[24] ,"low p_{T} jet and soft muon tagging","l");	 
    //leg->AddEntry(hsww_nj[27] ,"low p_{T} jet and soft muon tagging, with tagged highest p_{T} jet","l");
    leg->Draw();	  
  }
  else if(nsel == 1){
    hsww_nj[25]->Draw();
    hsww_nj[26]->Draw("same");
    //hsww_nj[29]->Draw("same");
    TLegend* leg = new TLegend(0.2,0.75,0.85,0.9);					      
    leg->SetFillColor(10);								      
    leg->SetTextSize(0.03);
    leg->SetFillColor(10);								      
    leg->SetTextSize(0.025);
    leg->AddEntry(hsww_nj[25] ,"total event tagging","l");		    
    leg->AddEntry(hsww_nj[26] ,"highest p_{T} jet tagging","l");
    //leg->AddEntry(hsww_nj[29] ,"highest p_{T} jet tagging with no low p_{T} tagged jets","l");
    leg->Draw();	  
  }
  else if(nsel == 2){
    hsww_nj[31]->SetMaximum(1.0);
    hsww_nj[31]->Draw();
    //hsww_nj[23]->Draw("same");
    hsww_nj[63]->Draw("same");
    //hsww_nj[61]->Draw("same");
    TLegend* leg = new TLegend(0.17,0.75,0.94,0.9);					      
    leg->SetFillColor(10);								      
    leg->SetTextSize(0.018);
    leg->AddEntry(hsww_nj[31] ,"low p_{T} jet and soft muon tagging with no high p_{T} tagged jets","l");
    //leg->AddEntry(hsww_nj[23] ,"low p_{T} jet and soft muon tagging","l");		      
    leg->AddEntry(hsww_nj[63] ,"low p_{T} jet and soft muon tagging with no high p_{T} tagged jets for WBF selection","l");
    //leg->AddEntry(hsww_nj[61] ,"low p_{T} jet and soft muon tagging for WBF selection","l");		    
    leg->Draw();	  
  }
  else if(nsel == 3){
    hsww_nj[25]->Draw();
    hsww_nj[64]->Draw("same");
    TLegend* leg = new TLegend(0.2,0.75,0.85,0.9);					      
    leg->SetFillColor(10);								      
    leg->SetTextSize(0.03);
    leg->SetFillColor(10);								      
    leg->SetTextSize(0.025);
    leg->AddEntry(hsww_nj[25] ,"total event tagging","l");		    
    leg->AddEntry(hsww_nj[64] ,"total event tagging with VBF selection","l");		    
    leg->Draw();	  
  }
  else if(nsel == 4){
    char titleName[100] = "";
    if     (option == 1) sprintf(titleName,"low p_{T} jet tagging");
    else if(option == 2) sprintf(titleName,"soft muon tagging");
    else if(option == 3) sprintf(titleName,"low p_{T} jet and soft muon tagging");
    else if(option == 4) sprintf(titleName,"low p_{T} jet and soft muon tagging, not in jets");
    else if(option == 5) sprintf(titleName,"total event tagging");
    else if(option == 6) sprintf(titleName,"leading jet tagging");
    else if(option == 7) sprintf(titleName,"low p_{T} jet and soft muon tagging, with tagged highest p_{T} jet","l");
    else if(option == 9) sprintf(titleName,"highest p_{T} jet tagging with no low p_{T} tagged jets","l");
    else if(option ==11) sprintf(titleName,"low p_{T} jet and soft muon tagging with no high p_{T} tagged jets","l");
    else {printf("Wrong option\n"); return;}
    atributes(hsww_nj[20+option],xTitle,1,yTitle);
    atributes(hsww_nj[40+option],xTitle,4,yTitle);
    hsww_nj[20+option]->Draw("hist,e");
    hsww_nj[40+option]->Draw("same,e");
    TLegend* leg = new TLegend(0.3,0.8,0.5,0.9);					      
    leg->SetFillColor(10);								      
    leg->SetTextSize(0.03);
    leg->AddEntry(hsww_nj[20+option] ,"simulation","l");	
    leg->AddEntry(hsww_nj[40+option] ,"data","l");	
    leg->Draw();	
    labelTitle  = new TPaveText(0.15,0.93,0.5,0.96,"NDCBR");
    labelTitle->SetTextAlign(12);
    labelTitle->SetTextSize(0.05);
    labelTitle->SetFillColor(kWhite);
    labelTitle->AddText(titleName);
    labelTitle->SetBorderSize(0);
    labelTitle->SetTextFont(132);
    labelTitle->SetLineWidth(2);
    labelTitle->Draw();
  }
  else if(nsel == 5){
    atributes(hsww_ptb[36+option],xTitle,1,yTitle);
    atributes(hsww_ptb[38+option],xTitle,4,yTitle);
    if(hsww_ptb[36+option]->GetSumOfWeights() > 0) hsww_ptb[36+option]->Scale(1./hsww_ptb[36+option]->GetSumOfWeights());
    if(hsww_ptb[38+option]->GetSumOfWeights() > 0) hsww_ptb[38+option]->Scale(1./hsww_ptb[38+option]->GetSumOfWeights());
    hsww_ptb[36+option]->Draw("e");
    hsww_ptb[38+option]->Draw("same,e");
    TLegend* leg = new TLegend(0.3,0.8,0.6,0.9);					      
    leg->SetFillColor(10);								      
    leg->SetTextSize(0.03);
    leg->AddEntry( hsww_ptb[36+option] ,"btag one jet","l");	   
    leg->AddEntry( hsww_ptb[38+option] ,"no btag","l");    
    leg->Draw();	
  }
  else if(nsel == 6){
    atributes(hsww_etab[36+option],xTitle,1,yTitle);
    atributes(hsww_etab[38+option],xTitle,4,yTitle);
    if(hsww_etab[36+option]->GetSumOfWeights() > 0) hsww_etab[36+option]->Scale(1./hsww_etab[36+option]->GetSumOfWeights());
    if(hsww_etab[38+option]->GetSumOfWeights() > 0) hsww_etab[38+option]->Scale(1./hsww_etab[38+option]->GetSumOfWeights());
    hsww_etab[38+option]->Draw("e");
    hsww_etab[36+option]->Draw("same,e");
    TLegend* leg = new TLegend(0.3,0.8,0.6,0.9);					      
    leg->SetFillColor(10);								      
    leg->SetTextSize(0.03);
    leg->AddEntry( hsww_etab[36+option] ,"btag one jet","l");
    leg->AddEntry( hsww_etab[38+option] ,"no btag","l");     	   
    leg->Draw();	
  }
  else if(nsel == 7){
    atributes( hsww_etab[24+option],xTitle,2,yTitle);
    atributes( hsww_etab[32+option],xTitle,8,yTitle);
    hsww_etab[24+option]->Divide(hsww_etab[20+option]);
    hsww_etab[32+option]->Divide(hsww_etab[28+option]);
    hsww_etab[24+option]->SetMaximum(1.0);
    hsww_etab[24+option]->Draw("e");
    hsww_etab[32+option]->Draw("same,e");
    TLegend* leg = new TLegend(0.3,0.8,0.8,0.9);					      
    leg->SetFillColor(10);								      
    leg->SetTextSize(0.03);
    leg->AddEntry( hsww_etab[24+option] ,"btag one jet","l");	      
    leg->AddEntry( hsww_etab[32+option] ,"btag one jet, VBF selection","l");		  
    leg->Draw();	
  }
  else if(nsel == 8){
    atributes( hsww_etab[44],xTitle,8,yTitle);
    atributes( hsww_etab[46],xTitle,2,yTitle);
    int theOption = option;
    if(option >= 4) theOption = theOption - 4;
    hsww_etab[44]->Scale(0.0);
    hsww_etab[46]->Scale(0.0);
    hsww_etab[44]->Add(hsww_etab[28+theOption],1.0);
    hsww_etab[44]->Add(hsww_etab[32+theOption],-1.0);
    hsww_etab[24+theOption]->Divide(hsww_etab[20+theOption]);
    for(int i=1; i<=hsww_etab[46]->GetNbinsX(); i++){
      double factor = 1.0;
      if     (option <= 3){
        if(hsww_etab[28+theOption]->GetBinContent(i) > 0 &&
           hsww_etab[32+theOption]->GetBinContent(i) > 0) factor = hsww_etab[32+theOption]->GetBinContent(i) / hsww_etab[28+theOption]->GetBinContent(i);
      }
      else if(option <= 8){
        if(hsww_etab[24+theOption]->GetBinContent(i) > 0 &&
	   hsww_etab[24+theOption]->GetBinContent(i) < 1) factor = hsww_etab[24+theOption]->GetBinContent(i);
      }
      hsww_etab[46]->SetBinContent(i,hsww_etab[32+theOption]->GetBinContent(i)*(1-factor)/factor);
    }
    hsww_etab[44]->Draw("e");
    hsww_etab[46]->Draw("same,e");
    TLegend* leg = new TLegend(0.3,0.8,0.8,0.9);					      
    leg->SetFillColor(10);								      
    leg->SetTextSize(0.03);
    leg->AddEntry( hsww_etab[44] ,"non tagged events","l");         
    leg->AddEntry( hsww_etab[46] ,"non tagged events, prediction","l");              
    leg->Draw();	
  }
}
