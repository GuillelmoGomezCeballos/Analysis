TAxis *xaxis, *yaxis;

const int nPoints = 8;

void mlp_graph()
{

  TString type[] = {
    "BDT",
    "N,N",
    "N+1,N",
    "N+2,N",
    "N+2,N+1,N",
    "N",
    "N+1",
    "N+2"
  };


  double significance00[nPoints] = {0.048546,0.135713,0.136832,0.131017,0.117603,0.122263,0.117234,0.119644};
  double significance01[nPoints] = {0.132486,0.275378,0.303424,0.272423,0.298035,0.256767,0.277789,0.265437};
  double significance02[nPoints] = {0.226100,0.508630,0.475136,0.480520,0.489636,0.466586,0.477218,0.462684};
  double significance03[nPoints] = {0.416997,0.637742,0.560789,0.605910,0.695875,0.612516,0.622859,0.633298};
  double significance04[nPoints] = {1.995363,1.301982,1.292677,1.308112,1.296675,1.157200,1.243139,1.340399};
  double significance05[nPoints] = {2.368028,1.349330,1.360647,1.477600,1.296656,1.294673,1.409314,1.225513};
  double significance06[nPoints] = {1.045880,0.916360,0.938321,0.952872,1.022146,0.897633,0.975534,0.933211};
  double significance07[nPoints] = {0.582218,0.686911,0.660015,0.611648,0.722242,0.626163,0.640420,0.583328};
  double significance08[nPoints] = {0.332242,0.558672,0.583506,0.527362,0.494861,0.527095,0.513648,0.499829};
  double significance09[nPoints] = {0.179107,0.407321,0.424121,0.393238,0.434051,0.346433,0.368028,0.390279};
  double significance10[nPoints] = {0.074538,0.369159,0.333947,0.361120,0.271573,0.348610,0.325430,0.325058};
  double significance11[nPoints] = {0.057835,0.250245,0.235895,0.269665,0.254727,0.246628,0.256478,0.238247};
  double significance12[nPoints] = {0.127241,0.277091,0.290699,0.270071,0.279029,0.249416,0.233659,0.256198};

  double x[nPoints], xmax[13];

  double max[13];

  for (int i=0; i<13; i++) {max[i] = 0.0;};

  for (int i=0; i<nPoints; i++) {
    x[i]    = i;
    xmax[i] = 0;
  }

  for (int i=0; i<nPoints; i++) {
    if (significance00[i] > max[ 0]) {xmax[ 0] = i; max[ 0] = significance00[i];};
    if (significance01[i] > max[ 1]) {xmax[ 1] = i; max[ 1] = significance01[i];};
    if (significance02[i] > max[ 2]) {xmax[ 2] = i; max[ 2] = significance02[i];};
    if (significance03[i] > max[ 3]) {xmax[ 3] = i; max[ 3] = significance03[i];};
    if (significance04[i] > max[ 4]) {xmax[ 4] = i; max[ 4] = significance04[i];};
    if (significance05[i] > max[ 5]) {xmax[ 5] = i; max[ 5] = significance05[i];};
    if (significance06[i] > max[ 6]) {xmax[ 6] = i; max[ 6] = significance06[i];};
    if (significance07[i] > max[ 7]) {xmax[ 7] = i; max[ 7] = significance07[i];};
    if (significance08[i] > max[ 8]) {xmax[ 8] = i; max[ 8] = significance08[i];};
    if (significance09[i] > max[ 9]) {xmax[ 9] = i; max[ 9] = significance09[i];};
    if (significance10[i] > max[10]) {xmax[10] = i; max[10] = significance10[i];};
    if (significance11[i] > max[11]) {xmax[11] = i; max[11] = significance11[i];};
    if (significance12[i] > max[12]) {xmax[12] = i; max[12] = significance12[i];};
  }

  TGraph *y00Graph = new TGraph(nPoints,x,significance00);
  TGraph *y01Graph = new TGraph(nPoints,x,significance01);
  TGraph *y02Graph = new TGraph(nPoints,x,significance02);
  TGraph *y03Graph = new TGraph(nPoints,x,significance03);
  TGraph *y04Graph = new TGraph(nPoints,x,significance04);
  TGraph *y05Graph = new TGraph(nPoints,x,significance05);
  TGraph *y06Graph = new TGraph(nPoints,x,significance06);
  TGraph *y07Graph = new TGraph(nPoints,x,significance07);
  TGraph *y08Graph = new TGraph(nPoints,x,significance08);
  TGraph *y09Graph = new TGraph(nPoints,x,significance09);
  TGraph *y10Graph = new TGraph(nPoints,x,significance10);
  TGraph *y11Graph = new TGraph(nPoints,x,significance11);
  TGraph *y12Graph = new TGraph(nPoints,x,significance12);

  SetTGraph(y00Graph, 3,  1, 20);
  SetTGraph(y01Graph, 3,  2, 20);
  SetTGraph(y02Graph, 3,  3, 20);
  SetTGraph(y03Graph, 3,  4, 20);
  SetTGraph(y04Graph, 3,  5, 20);
  SetTGraph(y05Graph, 3,  6, 20);
  SetTGraph(y06Graph, 3,  7, 20);
  SetTGraph(y07Graph, 3,  8, 20);
  SetTGraph(y08Graph, 3,  9, 20);
  SetTGraph(y09Graph, 3, 11, 20);
  SetTGraph(y10Graph, 3, 12, 20);
  SetTGraph(y11Graph, 3, 13, 20);
  SetTGraph(y12Graph, 3, 14, 20);

  //----------------------------------------------------------------------------
  // Draw
  //----------------------------------------------------------------------------
  TCanvas *cvs = new TCanvas("cvs","cvs",0,0,750,600);

  cvs->SetGridx();
  cvs->SetGridy();

  TMultiGraph *mgraph = new TMultiGraph();
  mgraph->Add(y00Graph);
  mgraph->Add(y01Graph);
  mgraph->Add(y02Graph);
  mgraph->Add(y03Graph);
  mgraph->Add(y04Graph);
  mgraph->Add(y05Graph);
  mgraph->Add(y06Graph);
  mgraph->Add(y07Graph);
  mgraph->Add(y08Graph);
  mgraph->Add(y09Graph);
  mgraph->Add(y10Graph);
  mgraph->Add(y11Graph);
  mgraph->Add(y12Graph);

  mgraph->Draw("alp");
  SetMultiGraph(mgraph,xaxis,yaxis,"method","significance = #frac{S}{#sqrt{B+0.04B^{2}}}");

  mgraph->GetXaxis()->SetLabelSize(0.0);


 for (int i=0; i<13; i++) {
   TMarker *maxpoint = new TMarker(xmax[i],max[i],4);
   maxpoint->SetMarkerSize(2.0);
   maxpoint->SetMarkerColor(1);
   maxpoint->Draw();
 }

  TLatex *xtype = new TLatex();
  xtype->SetTextFont (  42);  
  xtype->SetTextSize (0.04);
  xtype->SetTextAlign(  12);
  xtype->SetTextAngle( -45);

  for (i=0; i<nPoints; i++) {
    xtype->DrawLatex(x[i],-0.09,type[i].Data());
  }

  xaxis = (TAxis*)mgraph->GetXaxis();
  xaxis->SetTitleOffset(1.15*xaxis->GetTitleOffset());
  
  yaxis = (TAxis*)mgraph->GetYaxis();
  yaxis->SetTitleOffset(0.85*yaxis->GetTitleOffset());

  //----------------------------------------------------------------------------
  // Titles
  //----------------------------------------------------------------------------
  TPad *pad = (TPad*)cvs->GetPad(0);
  pad->Update();

  PRLabel(pad,"default variables, 0 jets",0.0,11,0.05);
  PRLabel(pad,"BDT vs. MLP",              1.0,31,0.05);

  TLegend *leg = SetLegend(0.74,0.65,0.82,0.91);
  //  leg->SetFillStyle(1001);
  leg->SetTextSize(0.03);

  TH1F *kk = new TH1F();
  kk->SetLineColor(10);

  leg->AddEntry(kk,"m_{H}","l");
  leg->AddEntry(y00Graph," 120","l");
  leg->AddEntry(y01Graph," 130","l");
  leg->AddEntry(y02Graph," 140","l");
  leg->AddEntry(y03Graph," 150","l");
  leg->AddEntry(y04Graph," 160","l");
  leg->AddEntry(y05Graph," 170","l");

  leg->Draw();

  TLegend *leg2 = SetLegend(0.83,0.65,0.91,0.91);
  //  leg2->SetFillStyle(1001);
  leg2->SetTextSize(0.03);

  leg2->AddEntry(y06Graph," 180","l");
  leg2->AddEntry(y07Graph," 190","l");
  leg2->AddEntry(y08Graph," 200","l");
  leg2->AddEntry(y09Graph," 210","l");
  leg2->AddEntry(y10Graph," 220","l");
  leg2->AddEntry(y11Graph," 250","l");
  leg2->AddEntry(y12Graph," 300","l");

  leg2->Draw();

  cvs->SaveAs("pic/mlp-summary.eps");
}

void SetTGraph(TGraph *gr, int width, int color, int style)
{
  gr->SetLineWidth  (width);
  gr->SetLineColor  (color);
  gr->SetMarkerColor(color);
  gr->SetMarkerSize ( 1.10);
  gr->SetMarkerStyle(style);
}
