void test_top_eff_mc(){
int reBin = 5;

TFile *_file1d = TFile::Open("rootfiles/histo_nice1d.root");
_file1d->cd();
TH1D* histo1d  = (TH1D*) histo2->Clone("histo1d");
histo1d->SetDirectory(0);
_file1d->Close();

TFile *_file1n = TFile::Open("rootfiles/histo_nice1n.root");
_file1n->cd();
TH1D* histo1n  = (TH1D*) histo2->Clone("histo1n");
TH1D* histoD1n = (TH1D*) histo2->Clone("histoD1n");
histo1n ->SetDirectory(0);
histoD1n->SetDirectory(0);
_file1n->Close();

TFile *_file2d = TFile::Open("rootfiles/histo_nice2d.root");
_file2d->cd();
TH1D* histo2d  = (TH1D*) histo2->Clone("histo2d");
histo2d->SetDirectory(0);
_file2d->Close();

TFile *_file2n = TFile::Open("rootfiles/histo_nice2n.root");
_file2n->cd();
TH1D* histo2n  = (TH1D*) histo2->Clone("histo2n");
TH1D* histoD2n = (TH1D*) histo2->Clone("histoD2n");
histo2n ->SetDirectory(0);
histoD2n->SetDirectory(0);
_file2n->Close();

histo1n ->Rebin(reBin);
histo1d ->Rebin(reBin);
histo2n ->Rebin(reBin);
histo2d ->Rebin(reBin);
histoD1n->Rebin(reBin);
histoD2n->Rebin(reBin);

histoD1n->Divide(histo1d);
histoD2n->Divide(histo2d);

printf("1j: %6.1f / %6.2f = %5.3f\n",histo1n->GetSumOfWeights(),histo1d->GetSumOfWeights(),histo1n->GetSumOfWeights()/histo1d->GetSumOfWeights());
printf("2j: %6.1f / %6.2f = %5.3f\n",histo2n->GetSumOfWeights(),histo2d->GetSumOfWeights(),histo2n->GetSumOfWeights()/histo2d->GetSumOfWeights());

histoD1n->SetLineColor(1);
histoD2n->SetLineColor(2);
histoD1n->SetMarkerColor(1);
histoD2n->SetMarkerColor(2);
histoD1n->SetLineWidth(3);
histoD2n->SetLineWidth(3);
//TFile* outFile = new TFile("test_top_eff_mc.root","recreate");
//outFile->cd();histoD1n->SetMarkerColor(1);
//  histo1n ->Write();
//  histo1d ->Write();
//  histo2n ->Write();
//  histo2d ->Write();
//  histoD1n->Write();
//  histoD2n->Write();
//outFile->Close();
TCanvas *canvas = new TCanvas("c1","c1",500,500);
canvas->cd();
histoD1n->Draw("e");
histoD2n->Draw("e,same");
canvas->Update();
}
