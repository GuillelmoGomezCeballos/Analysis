#include "MitNtupleEvent.C"

const int verboseLevel =  0;
const int nDecays      = 30;

//------------------------------------------------------------------------------
// PlotHiggsRes2
//------------------------------------------------------------------------------
void PlotMlpEpochs
(
 int         njets    = 0,
 int         mhiggs   = 140,
 const char* sigFile1 = "data/inputNtuple-train-data-standard-H140_WW_2l.root",
 const char* sigFile2 = "proc/inputNtuple-train-data-standard-H140_WW_2l.epoch-test-140.0jets.tmva.root",
 const char* bgdFile1 = "data/inputNtuple-data-standard-HBCK_WW_2l-train.root",
 const char* bgdFile2 = "proc/inputNtuple-data-standard-HBCK_WW_2l-train.epoch-test-140.0jets.tmva.root"
 )
{
  int myColor[] = {
    1,
    2,
    4,
    8,
    15,
    38,
    50
  };

  char title[200];
  sprintf(title,"pic/mlp_default-vars_H%d-%djets_significance.eps",mhiggs,njets);

  TString cName = title;

  // Recover both the original ntuples and the tmva results for signal
  TTree* sig1 = getTreeFromFile(sigFile1,"all");
  assert(sig1);
  TTree* sig2 = getTreeFromFile(sigFile2,"nt");
  assert(sig2);
  assert(sig1->GetEntries()==sig2->GetEntries());

  // Recover both the original ntuples and the tmva results for background
  TTree* bgd1 = getTreeFromFile(bgdFile1,"all");
  assert(bgd1);
  TTree* bgd2 = getTreeFromFile(bgdFile2,"nt");
  assert(bgd2);
  assert(bgd1->GetEntries()==bgd2->GetEntries());

  // Setup histograms
  TH1D* sigBDT      = new TH1D("sigBDT",   "BDT",   150, -1.0, 0.5);
  TH1D* bgdAllBDT   = new TH1D("bgdAllBDT","BDT",   150, -1.0, 0.5);

  TH1D* sigMLP_0    = new TH1D("sigMLP_0", "MLP_0", 300, -1.5, 1.5);
  TH1D* sigMLP_1    = new TH1D("sigMLP_1", "MLP_1", 300, -1.5, 1.5);
  TH1D* sigMLP_2    = new TH1D("sigMLP_2", "MLP_2", 300, -1.5, 1.5);
  TH1D* sigMLP_3    = new TH1D("sigMLP_3", "MLP_3", 300, -1.5, 1.5);
  TH1D* sigMLP_4    = new TH1D("sigMLP_4", "MLP_4", 300, -1.5, 1.5);
  TH1D* sigMLP_5    = new TH1D("sigMLP_5", "MLP_5", 300, -1.5, 1.5);
  TH1D* sigMLP_6    = new TH1D("sigMLP_6", "MLP_6", 300, -1.5, 1.5);

  TH1D* bgdAllMLP_0 = new TH1D("bgdAllMLP_0","MLP_0", 300, -1.5, 1.5);
  TH1D* bgdAllMLP_1 = new TH1D("bgdAllMLP_1","MLP_1", 300, -1.5, 1.5);
  TH1D* bgdAllMLP_2 = new TH1D("bgdAllMLP_2","MLP_2", 300, -1.5, 1.5);
  TH1D* bgdAllMLP_3 = new TH1D("bgdAllMLP_3","MLP_3", 300, -1.5, 1.5);
  TH1D* bgdAllMLP_4 = new TH1D("bgdAllMLP_4","MLP_4", 300, -1.5, 1.5);
  TH1D* bgdAllMLP_5 = new TH1D("bgdAllMLP_5","MLP_5", 300, -1.5, 1.5);
  TH1D* bgdAllMLP_6 = new TH1D("bgdAllMLP_6","MLP_6", 300, -1.5, 1.5);

  sigBDT     ->Sumw2();
  bgdAllBDT  ->Sumw2();

  sigMLP_0   ->Sumw2();
  sigMLP_1   ->Sumw2();
  sigMLP_2   ->Sumw2();
  sigMLP_3   ->Sumw2();
  sigMLP_4   ->Sumw2();
  sigMLP_5   ->Sumw2();
  sigMLP_6   ->Sumw2();

  bgdAllMLP_0->Sumw2();
  bgdAllMLP_1->Sumw2();
  bgdAllMLP_2->Sumw2();
  bgdAllMLP_3->Sumw2();
  bgdAllMLP_4->Sumw2();
  bgdAllMLP_5->Sumw2();
  bgdAllMLP_6->Sumw2();

  // Start loop over signal events
  //----------------------------------------------------------------------------
  MitNtupleEvent sigEvent(sig1);

  float run,event,bdt,mlp_0,mlp_1,mlp_2,mlp_3,mlp_4,mlp_5,mlp_6;
  sig2->SetBranchAddress("run",  &run  );
  sig2->SetBranchAddress("event",&event);
  sig2->SetBranchAddress("bdt",  &bdt  );
  sig2->SetBranchAddress("mlp_0",&mlp_0);
  sig2->SetBranchAddress("mlp_1",&mlp_1);
  sig2->SetBranchAddress("mlp_2",&mlp_2);
  sig2->SetBranchAddress("mlp_3",&mlp_3);
  sig2->SetBranchAddress("mlp_4",&mlp_4);
  sig2->SetBranchAddress("mlp_5",&mlp_5);
  sig2->SetBranchAddress("mlp_6",&mlp_6);

  double scaleFactorLum = 0.1; // 100 pb-1

  int    nSig        = sig1->GetEntries();
  double nSigAcc     = 0;

  for (int i=0; i<nSig; ++i) {

    if (i%1000 == 0 && verboseLevel > 0)
      printf("--- reading event %5d of %5d\n",i,nSig);

    sigEvent.GetEntry(i);
    sig2->GetEntry(i);
    assert(run == sigEvent.H_run);
  
    if (sigEvent.H_njets != njets) continue;
  
    double weight = scaleFactorLum * sigEvent.H_weight;

    nSigAcc += weight;

    sigBDT  ->Fill(TMath::Max(TMath::Min(bdt,  1.0),-0.99999),weight);
    sigMLP_0->Fill(TMath::Max(TMath::Min(mlp_0,1.5),-1.49999),weight);
    sigMLP_1->Fill(TMath::Max(TMath::Min(mlp_1,1.5),-1.49999),weight);
    sigMLP_2->Fill(TMath::Max(TMath::Min(mlp_2,1.5),-1.49999),weight);
    sigMLP_3->Fill(TMath::Max(TMath::Min(mlp_3,1.5),-1.49999),weight);
    sigMLP_4->Fill(TMath::Max(TMath::Min(mlp_4,1.5),-1.49999),weight);
    sigMLP_5->Fill(TMath::Max(TMath::Min(mlp_5,1.5),-1.49999),weight);
    sigMLP_6->Fill(TMath::Max(TMath::Min(mlp_6,1.5),-1.49999),weight);
  }

  // Start loop over background events
  //----------------------------------------------------------------------------
  MitNtupleEvent bgdEvent(bgd1);

  bgd2->SetBranchAddress("run",  &run  );
  bgd2->SetBranchAddress("event",&event);
  bgd2->SetBranchAddress("bdt",  &bdt);
  bgd2->SetBranchAddress("mlp_0",&mlp_0);
  bgd2->SetBranchAddress("mlp_1",&mlp_1);
  bgd2->SetBranchAddress("mlp_2",&mlp_2);
  bgd2->SetBranchAddress("mlp_3",&mlp_3);
  bgd2->SetBranchAddress("mlp_4",&mlp_4);
  bgd2->SetBranchAddress("mlp_5",&mlp_5);
  bgd2->SetBranchAddress("mlp_6",&mlp_6);

  int nBgd=bgd1->GetEntries();

  for (int i=0; i<nBgd; ++i) {

    if (i%5000 == 0 && verboseLevel > 0)
      printf("--- reading event %5d of %5d\n",i,nBgd);

    bgdEvent.GetEntry(i);
    bgd2->GetEntry(i);
    assert(run == bgdEvent.H_run);
    if (bgdEvent.H_njets != njets) continue;

    double weight = scaleFactorLum * bgdEvent.H_weight;

    //    if (bgdEvent.H_decay == 0 || bgdEvent.H_decay == 1) weight *= 0.01;

    bgdAllBDT  ->Fill(TMath::Max(TMath::Min(bdt,  1.0),-0.99999),weight);
    bgdAllMLP_0->Fill(TMath::Max(TMath::Min(mlp_0,1.5),-1.49999),weight);
    bgdAllMLP_1->Fill(TMath::Max(TMath::Min(mlp_1,1.5),-1.49999),weight);
    bgdAllMLP_2->Fill(TMath::Max(TMath::Min(mlp_2,1.5),-1.49999),weight);
    bgdAllMLP_3->Fill(TMath::Max(TMath::Min(mlp_3,1.5),-1.49999),weight);
    bgdAllMLP_4->Fill(TMath::Max(TMath::Min(mlp_4,1.5),-1.49999),weight);
    bgdAllMLP_5->Fill(TMath::Max(TMath::Min(mlp_5,1.5),-1.49999),weight);
    bgdAllMLP_6->Fill(TMath::Max(TMath::Min(mlp_6,1.5),-1.49999),weight);
  }

  //----------------------------------------------------------------------------
  // Drawing part
  //----------------------------------------------------------------------------
  TCanvas* canvas = new TCanvas("canvas","canvas");
  canvas->SetGrid(1,1);
  TH2F* zone4 = new TH2F("zone4","discriminant vs. significance; xtitle; ytitle",
			 1, -1.5, 0.3, 1, 1.e-3, 0.6);

  axis2F(zone4,"discriminant","significance = #frac{S}{#sqrt{B+0.04B^{2}}}");
  canvas->SetLeftMargin(1.2*canvas->GetLeftMargin());
  zone4->DrawCopy();

  double significance[8];
  for (int i=0; i<8; i++) significance[i] = 0.0;

  printf("double bestSignificance[] = {");

  //  TGraphErrors* gBDT   = makeSignificanceCurve(sigBDT,   bgdAllBDT,   "BDT      ", significance[0]);
  TGraphErrors* gMLP_0 = makeSignificanceCurve(sigMLP_0, bgdAllMLP_0, "N,N      ", significance[1]);
  TGraphErrors* gMLP_1 = makeSignificanceCurve(sigMLP_1, bgdAllMLP_1, "N+1,N    ", significance[2]);
  TGraphErrors* gMLP_2 = makeSignificanceCurve(sigMLP_2, bgdAllMLP_2, "N+2,N    ", significance[3]);
  TGraphErrors* gMLP_3 = makeSignificanceCurve(sigMLP_3, bgdAllMLP_3, "N+2,N+1,N", significance[4]);
  TGraphErrors* gMLP_4 = makeSignificanceCurve(sigMLP_4, bgdAllMLP_4, "N        ", significance[5]);
  TGraphErrors* gMLP_5 = makeSignificanceCurve(sigMLP_5, bgdAllMLP_5, "N+1      ", significance[6]);
  TGraphErrors* gMLP_6 = makeSignificanceCurve(sigMLP_6, bgdAllMLP_6, "N+2      ", significance[7]);

  printf("};\n");

  //  setGraph(gBDT  , 1, 20);
  setGraph(gMLP_0, myColor[0], 21);
  setGraph(gMLP_1, myColor[1], 22);
  setGraph(gMLP_2, myColor[2], 23);
  setGraph(gMLP_3, myColor[3], 24);
  setGraph(gMLP_4, myColor[4], 25);
  setGraph(gMLP_5, myColor[5], 26);
  setGraph(gMLP_6, myColor[6], 27);

  TLegend *leg = SetLegend(0.218,0.523,0.526,0.911);
  leg->SetTextSize(0.033);
  //  sprintf(title," %4.2f #rightarrow BDT      ",significance[0]); leg->AddEntry(gBDT,  title,"L");
  sprintf(title," %4.2f #rightarrow N,N      ",significance[1]); leg->AddEntry(gMLP_0,title,"L");
  sprintf(title," %4.2f #rightarrow N+1,N    ",significance[2]); leg->AddEntry(gMLP_1,title,"L");
  sprintf(title," %4.2f #rightarrow N+2,N    ",significance[3]); leg->AddEntry(gMLP_2,title,"L");
  sprintf(title," %4.2f #rightarrow N+2,N+1,N",significance[4]); leg->AddEntry(gMLP_3,title,"L");
  sprintf(title," %4.2f #rightarrow N        ",significance[5]); leg->AddEntry(gMLP_4,title,"L");
  sprintf(title," %4.2f #rightarrow N+1      ",significance[6]); leg->AddEntry(gMLP_5,title,"L");
  sprintf(title," %4.2f #rightarrow N+2      ",significance[7]); leg->AddEntry(gMLP_6,title,"L");
  leg->Draw("same");

  TPad* pad = (TPad*)canvas->GetPad(0);

  sprintf(title,"%.1f signal events @ %d GeV, %d jets",nSigAcc,mhiggs,njets);

  //  PRLabel(pad,"default variables",0,11,0.05);
  PRLabel(pad,title,1,31,0.05);

  //  canvas->SaveAs(cName.Data());
}

//------------------------------------------------------------------------------
// getTreeFromFile
//------------------------------------------------------------------------------
TTree* getTreeFromFile(const char* infname, const char* tname)
{
  bool verbose = false;

  if (verbose) {
    cout << "--- Open file " << infname << endl;
  }
  
  TFile* inf = new TFile(infname,"read");
  assert(inf);
  TTree* t = (TTree*)inf->Get(tname);
  assert(t);

  if (verbose) {
    cout << "---\tRecovered tree " << t->GetName()
	 << " with "<< t->GetEntries() << " entries" << endl;
  }
  
  return t;
}

//------------------------------------------------------------------------------
// setGraph
//------------------------------------------------------------------------------
void setGraph(TGraphErrors* g, int color, int marker)
{
  g->SetLineColor  ( color);
  g->SetLineWidth  (   2.0);
  g->SetMarkerColor( color);
  g->SetMarkerSize (   0.5);
  g->SetMarkerStyle(marker);

  g->Draw("L");
}

//------------------------------------------------------------------------------
// makeBgdSigCurve
//------------------------------------------------------------------------------
TGraphErrors* makeBgdSigCurve(TH1D* sig, TH1D* bgd, const char* name)
{
  double xs [1000];
  double ys [1000];
  double dys[1000];
  double dxs[1000];

  int n = 0;
	
  double sigI = sig->Integral();
  double bgdI = bgd->Integral();
	
  for (int bin=1; bin<=sig->GetNbinsX(); ++bin) {
    
    double s = sig->Integral(bin, sig->GetNbinsX());
    double b = bgd->Integral(bin, sig->GetNbinsX());
    
    xs[n] = s;
    ys[n] = b;
		
    dxs[n] = 1.e-6;
    dys[n] = 1.e-6;
		
    ++n;
  }
	
  TGraphErrors* g = new TGraphErrors(n, xs, ys, dxs, dys);
  g->SetName(name);
  return g;
}

//------------------------------------------------------------------------------
// makeSignificanceCurve
//------------------------------------------------------------------------------
TGraphErrors* makeSignificanceCurve(TH1D*       sig,
				    TH1D*       bgd,
				    const char* name,
				    double&     significance)
{
  double xs [1000];
  double ys [1000];
  double dys[1000];
  double dxs[1000];

  int n = 0;
	
  double sigI = sig->Integral();
  double bgdI = bgd->Integral();

  double theValue[4] = {0., 0., 0., 0.};

  for (int bin=1; bin<=sig->GetNbinsX(); ++bin) {
		
    double s = sig->Integral(bin, sig->GetNbinsX());
    double b = bgd->Integral(bin, sig->GetNbinsX());
	
    //    if (b != 0) ys[n] = s/sqrt(b+0.04*b*b);
    if (b > 1.0) ys[n] = s/sqrt(b+0.04*b*b);
    else       ys[n] = 0.0;
    xs[n] = sig->GetBinCenter(bin);

    dxs[n] = 1.e-6;
    dys[n] = 1.e-6;
	   
    if (ys[n] > theValue[0] && b > 1.0) {
      theValue[0] = ys[n];
      theValue[1] = s;
      theValue[2] = b;
      theValue[3] = xs[n];	     
    }
    ++n;
  }
  printf("%f,",theValue[0]);

//  printf("%15s -> Sig=%6.3f, S=%6.3f, B=%6.3f, S/B=%6.3f, Bin=%6.3f\n",
//	 name,
//	 theValue[0],theValue[1],
//	 theValue[2],theValue[1]/theValue[2],theValue[3]);

  significance = theValue[0];
	
  TGraphErrors* g = new TGraphErrors(n, xs, ys, dxs, dys);
  g->SetName(name);
  return g;
}

//------------------------------------------------------------------------------
// makeRejEffCurve
//------------------------------------------------------------------------------
TGraphErrors* makeRejEffCurve(TH1D* sig, TH1D* bgd, const char* name)
{
  double xs [1000];
  double ys [1000];
  double dys[1000];
  double dxs[1000];

  int n = 0;
	
  double sigI = sig->Integral();
  double bgdI = bgd->Integral();
	
  for (int bin=1; bin<=sig->GetNbinsX(); ++bin) {
		
    double s = sig->Integral(bin, sig->GetNbinsX());
    double b = bgd->Integral(1, bin);
		
    double eff = s/sigI;
    double rej = b/bgdI;
		
    xs[n] = eff;
    ys[n] = rej;

    dxs[n] = 1.e-6;
    dys[n] = 1.e-6;
		
    ++n;
  }
	
  TGraphErrors* g = new TGraphErrors(n, xs, ys, dxs, dys);
  g->SetName(name);
  return g;
}

//------------------------------------------------------------------------------
// plotHistsInPad
//------------------------------------------------------------------------------
void plotHistsInPad(TH1D* h1, TH1D* h2)
{
  gPad->SetLogy();
  setHist(h1, 4, 20);
  setHist(h2, 2, 24);
  setPair(h1, h2);
  h1->DrawCopy("pe");
  h2->DrawCopy("pe same");

  TLegend* lg = new TLegend(0.65, 0.65, 0.93, 0.90);
  lg->SetFillColor(0);
  lg->AddEntry(h1," signal");
  lg->AddEntry(h2," background");
  lg->Draw("same");
}

//------------------------------------------------------------------------------
// setPair
//------------------------------------------------------------------------------
void setPair(TH1D* h1, TH1D* h2)
{
  double theMax = 10.*TMath::Max(h1->GetMaximum(),h2->GetMaximum());
  h1->SetMaximum(theMax);
  h2->SetMaximum(theMax);
}

//------------------------------------------------------------------------------
// setHist
//------------------------------------------------------------------------------
void setHist(TH1D* h, int color, int style)
{
  h->SetMarkerColor(color        );
  h->SetMarkerStyle(style        );
  h->SetLineColor  (color        );
  h->SetMarkerSize (0.5          );
  h->SetXTitle     (h->GetTitle());
}
