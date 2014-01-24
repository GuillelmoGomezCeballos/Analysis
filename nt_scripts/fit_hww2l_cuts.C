void fit_hww2l_cuts(int nsel=0){
 
  const int nBins = 21;
  
  double cutMassLow[nBins]       = { 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12};
  double cutMassHigh[nBins]      = { 40, 40, 45, 45, 50, 50, 50, 60, 80, 90,110,120,130,150,200,250,300,350,400,450,500};
  double cutPtMaxLow[nBins]      = { 20, 20, 25, 25, 27, 30, 34, 36, 38, 40, 44, 48, 52, 55, 70, 80, 90,110,120,130,140};
  double cutPtMinLow[nBins]      = { 10, 10, 10, 15, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25};
  double cutDeltaphilHigh[nBins] = {115,115, 90, 90, 90, 60, 60, 70, 90,100,110,120,130,140,175,175,175,175,175,175,175};
  double cutMTLow[nBins]         = { 70, 70, 75, 80, 80, 90,110,120,120,120,120,120,120,120,120,120,120,120,120,120,120};
  double cutMTHigh[nBins]        = {110,120,125,130,150,160,170,180,190,200,210,220,230,250,300,350,400,450,500,550,600};
  double cutMetLow[nBins]        = { 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35};
  double cutMetLowEM[nBins]      = { 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20};
  double cutsError[nBins];  for(int i=0; i<nBins;i++) cutsError[i] = 0.001;
  double cuts[nBins];
  double mass[nBins]      = {115,120,130,140,150,160,170,180,190,200,210,220,230,250,300,350,400,450,500,550,600};
  double massError[nBins];  for(int i=0; i<nBins;i++) massError[i] = 0.001;

  if     (nsel == 0) for(int i=0; i<nBins;i++) cuts[i] =  cutMassLow[i];    
  else if(nsel == 1) for(int i=0; i<nBins;i++) cuts[i] =  cutMassHigh[i];    
  else if(nsel == 2) for(int i=0; i<nBins;i++) cuts[i] =  cutPtMaxLow[i];    
  else if(nsel == 3) for(int i=0; i<nBins;i++) cuts[i] =  cutPtMinLow[i];
  else if(nsel == 4) for(int i=0; i<nBins;i++) cuts[i] =  cutDeltaphilHigh[i];
  else if(nsel == 5) for(int i=0; i<nBins;i++) cuts[i] =  cutMTLow[i];      
  else if(nsel == 6) for(int i=0; i<nBins;i++) cuts[i] =  cutMTHigh[i];      
  else if(nsel == 7) for(int i=0; i<nBins;i++) cuts[i] =  cutMetLow[i];    
  else if(nsel == 8) for(int i=0; i<nBins;i++) cuts[i] =  cutMetLowEM[i];

  TString title = "";
  if     (nsel == 0) for(int i=0; i<nBins;i++) title = "m_{ll} > X";
  else if(nsel == 1) for(int i=0; i<nBins;i++) title = "m_{ll} < X";
  else if(nsel == 2) for(int i=0; i<nBins;i++) title = "p_{T}^{max} > X";
  else if(nsel == 3) for(int i=0; i<nBins;i++) title = "p_{T}^{min} > X";
  else if(nsel == 4) for(int i=0; i<nBins;i++) title = "#Delta #phi_{ll} < X";
  else if(nsel == 5) for(int i=0; i<nBins;i++) title = "m_{T}^{llE_{T}^{miss}} > X";
  else if(nsel == 6) for(int i=0; i<nBins;i++) title = "m_{T}^{llE_{T}^{miss}} < X";
  else if(nsel == 7) for(int i=0; i<nBins;i++) title = "E_{T}^{miss}-ll > X";
  else if(nsel == 8) for(int i=0; i<nBins;i++) title = "E_{T}^{miss}-e#mu > X";

  c1 = new TCanvas("c1","gerrors2",0,0,500,500);
 
  gr = new TGraphErrors(nBins,mass,cuts,massError,cutsError);
 
  //func = new TF1("func",funcfit,120.,500.,9);
  //func->SetParameter(0,  2.81401e+01);
  //func->SetParameter(1, -1.79109e+01);
  //func->SetParameter(2,  1.70836e+02);
 
  //gr->Fit("func","roe");
 
  gr->Draw("ACP");
  gr->SetTitle("");
  gr->GetXaxis()->SetTitle("Higgs Mass [GeV/c^{2}]");
  gr->GetYaxis()->SetTitle(title);
  gr->GetXaxis()->SetTitleOffset(0.9);
  gr->GetXaxis()->SetTitleSize(0.04);
  gr->GetXaxis()->SetLabelSize(0.04);
  gr->GetYaxis()->SetTitleOffset(1.9);
  gr->GetYaxis()->SetTitleSize(0.04);
  gr->GetYaxis()->SetLabelSize(0.04);
  gr->GetYaxis()->CenterTitle(kTRUE);
  //gr->GetHistogram()->SetMinimum(0.00);
  //gr->GetHistogram()->SetMaximum(0.35);
}
Double_t funcfit(Double_t* x, Double_t* par) {
 
  double val = 0;
  for(int i=0; i<9; i++){
    double aux = par[i];
    for(int j=0; j<i; j++) aux = aux * x[0];
    val = val + aux;
  }
  return val;
}
