double makeNvtx(){
TFile* outputFile = TFile::Open( "puReweighting.root", "RECREATE" );
outputFile->cd();
TH1D *puWeights  = new TH1D("puWeights","puWeights",19,0.5,19.5);
for(int i=1; i<=19; i++){
  puWeights->SetBinContent(i,PileupReweightFactor(i));
}
puWeights->Write();
outputFile->Close();

}


Double_t PileupReweightFactor(Int_t nvtx) {
  Double_t weight = 1.0;
  if     (nvtx == 0) weight = 0.00000;
  else if(nvtx == 1) weight = 0.23111;
  else if(nvtx == 2) weight = 0.82486;
  else if(nvtx == 3) weight = 1.47045;
  else if(nvtx == 4) weight = 1.83450;
  else if(nvtx == 5) weight = 1.79663;
  else if(nvtx == 6) weight = 1.46496;
  else if(nvtx == 7) weight = 1.05682;
  else if(nvtx == 8) weight = 0.70823;
  else if(nvtx == 9) weight = 0.47386;
  else if(nvtx ==10) weight = 0.32382;
  else if(nvtx ==11) weight = 0.22383;
  else if(nvtx ==12) weight = 0.17413;
  else if(nvtx ==13) weight = 0.10930;
  else if(nvtx ==14) weight = 0.09563;
  else if(nvtx ==15) weight = 0.08367;
  else if(nvtx ==16) weight = 0.05418;
  else if(nvtx ==17) weight = 0.04891;
  else if(nvtx ==18) weight = 0.03515;
  else if(nvtx >=19) weight = 0.01000;
  
  return weight;
}
