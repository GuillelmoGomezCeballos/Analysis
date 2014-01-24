void HwwNoteFigs(TString filename = "../histo_tmva_new-ntuples-1_160_0.root",
                 int mass = 170, int ReBin = 4, int ysel = 0)
{
	TFile* infile = new TFile(filename.Data(),"READ");
	setTDRStyle(0);
	
	//BDTD dists
	TCanvas* c1 = new TCanvas("c1","c1",0,-200,500,500);
	
        TH1F* histoBg =  histo4->Clone();
        histoBg->Add(histo3);
        histoBg->Add(histo2);
        histoBg->Add(histo1);

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

        char YTitle[300];
        sprintf(YTitle,"events / bin");
        char XTitle[300];
        sprintf(XTitle,"BDT Output");
    	histo0->SetYTitle(YTitle);
    	histo1->SetYTitle(YTitle);
    	histo2->SetYTitle(YTitle);
    	histo3->SetYTitle(YTitle);
    	histo4->SetYTitle(YTitle);

    	histo0->SetXTitle(XTitle);
    	histo1->SetXTitle(XTitle);
    	histo2->SetXTitle(XTitle);
    	histo3->SetXTitle(XTitle);
    	histo4->SetXTitle(XTitle);
    	histo1->SetTitleSize(0.05, "X");
    	histo1->GetXaxis()->SetTitleFont(62);
    	histo1->GetXaxis()->SetLabelFont(61);
    	histo1->GetYaxis()->SetLabelFont(61); 
    	histo1->GetYaxis()->SetTitleOffset(1.3);
    	histo1->SetLabelSize(0.04, "Y");
    	histo1->SetLabelSize(0.04, "X");

	int min = histoBg->FindBin(-1.5);
	int max = histoBg->FindBin(1.5);

	histoBg->GetXaxis()->SetRange(min,max);
	
	histoBg->SetMarkerStyle(20);
	histoBg->SetMarkerSize(1.0);
	histoBg->GetYaxis()->SetTitleOffset(1.40);
	
	histo0->SetMarkerStyle(21);
	histo0->SetMarkerSize(1.0);
	
	histoBg->Rebin(ReBin);
	histo0->Rebin(ReBin);

	histoBg->SetLineColor(4);
	histo0->SetLineColor(1);

	scaleHist(histoBg);
	scaleHist(histo0);
	cout << "bg events: " << histoBg->GetSumOfWeights() << endl;
	cout << "si events: " << histo0->GetSumOfWeights() << endl;

	histoBg->SetYTitle("Events");
    	histo1->SetMinimum(0.01);
	if(ysel == 0)  {
	  histo1->Draw("hist");
	}
	else          {
	  histo0->Draw("E");
	  histo1->Draw("hist,same");
    	}
	histo2->Draw("hist,same");
    	histo3->Draw("hist,same");
    	histo4->Draw("hist,same");
    	histo0->Draw("E, same");
	//histoBg->DrawCopy("hist");
	//histo0->DrawCopy("hist,same");

	TLegend* leg = new TLegend(0.63, 0.75, 0.92, 0.92);
	leg->SetFillColor(0);
	char theSLine[100];
	if(mass != 999)	sprintf(theSLine,"Signal, m_{H}=%d GeV",mass);
	else    	sprintf(theSLine,"WW");
        cout << theSLine << endl;
	leg ->AddEntry(histo0,theSLine);
    	leg ->AddEntry(histo4,"W+Jets, W#gamma","F"); 
    	leg ->AddEntry(histo3,"di-boson","F");  
    	leg ->AddEntry(histo2,"t#bar{t}, tW","F"); 
    	leg ->AddEntry(histo1,"Drell-Yan","F"); 
	leg->Draw("same");
	
        TString fileOutput1(filename.Data());
	TString theLine = "_";
	theLine = theLine + "plot.eps";
        fileOutput1.ReplaceAll(".root",theLine.Data());
        fileOutput1.ReplaceAll("../","");
        fileOutput1.ReplaceAll("rootfiles_fastsim/","");
        fileOutput1.ReplaceAll("rootfiles_fullsim/","");
	c1->SaveAs(fileOutput1.Data());
        //return;
	//--------------------------------
	TCanvas* c3 = new TCanvas("c3","c3",550,-200,500,500);
	c3->SetLogx();
	c3->SetLogy();
	gPad->SetGrid(1,1);
	TGraphErrors* gBDTD = makeGraphFromHists(histo0, histoBg);
	gBDTD->Draw("APXl");
	TH1* zBDTD = gBDTD->GetHistogram();
	zBDTD->SetXTitle("Signal Events");
	zBDTD->SetYTitle("Signal/Background");
	//zBDTD->SetYTitle("Background Events");
	//zBDTD->DrawCopy();

	TLegend* leg3 = new TLegend(0.7, 0.7, 0.9, 0.9);
	leg3->SetFillColor(0);
	leg3->AddEntry(gBDTD,theSLine,"lp");
	leg3->Draw("same");
        TString fileOutput2(filename.Data());
	theLine = "_";
	theLine = theLine + "counts.eps";
        fileOutput2.ReplaceAll(".root",theLine.Data());
        fileOutput2.ReplaceAll("../","");
        fileOutput2.ReplaceAll("rootfiles_fastsim/","");
        fileOutput2.ReplaceAll("rootfiles_fullsim/","");
	c3->SaveAs(fileOutput2.Data());
	
}

TGraphErrors* makeGraphFromHists(TH1* hsig, TH1* hbgd)
{
	const int nbins = hsig->GetNbinsX();
	double xs[1000], ys[1000], dxs[1000], dys[1000];

	int i=0;
	for (int bin=1; bin <= nbins; ++bin) {
		double bgds = hbgd->Integral(bin, nbins);
		double sigs = hsig->Integral(bin, nbins);
		if (sigs==0. || bgds==0.) continue;
		xs[i] = sigs;
		dxs[i] = sqrt(sigs);
		ys[i] = sigs/bgds;
                dys[i] = sqrt(bgds);
		//ys[i] = bgds;
		//dys[i] = sqrt(bgds);
		++i;
	}
	TGraphErrors* g = new TGraphErrors(i, xs, ys, dxs, dys);
	return g;
}

void scaleHist(TH1* h)
{
	//h->SetMaximum(500.);
	//h->Scale(1./h->GetBinWidth(1));
}

void setTDRStyle(Int_t ylog) 
{
	
	TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");
	
	// For the canvas:
	tdrStyle->SetCanvasBorderMode(0);
	tdrStyle->SetCanvasColor(kWhite);
	tdrStyle->SetCanvasDefH(600); //Height of canvass
	tdrStyle->SetCanvasDefW(600); //Width of canvas
	tdrStyle->SetCanvasDefX(0);   //POsition on screen
	tdrStyle->SetCanvasDefY(0);
	
	// For the Pad:
	tdrStyle->SetPadBorderMode(0);
	// tdrStyle->SetPadBorderSize(Width_t size = 1);
	tdrStyle->SetPadColor(kWhite);
	tdrStyle->SetPadGridX(false);
	tdrStyle->SetPadGridY(false);
	tdrStyle->SetGridColor(0);
	tdrStyle->SetGridStyle(3);
	tdrStyle->SetGridWidth(1);
	
	// For the frame:
	tdrStyle->SetFrameBorderMode(0);
	tdrStyle->SetFrameBorderSize(1);
	tdrStyle->SetFrameFillColor(0);
	tdrStyle->SetFrameFillStyle(0);
	tdrStyle->SetFrameLineColor(1);
	tdrStyle->SetFrameLineStyle(1);
	tdrStyle->SetFrameLineWidth(1);
	
	// For the histo:
	// tdrStyle->SetHistFillColor(1);
	// tdrStyle->SetHistFillStyle(0);
	tdrStyle->SetHistLineColor(1);
	tdrStyle->SetHistLineStyle(0);
	tdrStyle->SetHistLineWidth(2);
	// tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
	// tdrStyle->SetNumberContours(Int_t number = 20);
	
	tdrStyle->SetEndErrorSize(2);
	//  tdrStyle->SetErrorMarker(20);
	tdrStyle->SetErrorX(0.);
	
	tdrStyle->SetMarkerStyle(20);
	
	//For the fit/function:
	//  tdrStyle->SetOptFit(1);
	tdrStyle->SetOptFit(0);
	tdrStyle->SetFitFormat("5.4g");
	tdrStyle->SetFuncColor(1);
	tdrStyle->SetFuncStyle(1);
	tdrStyle->SetFuncWidth(1);
	
	//For the date:
	tdrStyle->SetOptDate(0);
	// tdrStyle->SetDateX(Float_t x = 0.01);
	// tdrStyle->SetDateY(Float_t y = 0.01);
	
	// For the statistics box:
	tdrStyle->SetOptFile(0);
	tdrStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
	tdrStyle->SetStatColor(kWhite);
	tdrStyle->SetStatFont(42);
	tdrStyle->SetStatFontSize(0.025);
	tdrStyle->SetStatTextColor(1);
	tdrStyle->SetStatFormat("6.4g");
	tdrStyle->SetStatBorderSize(1);
	tdrStyle->SetStatH(0.1);
	tdrStyle->SetStatW(0.15);
	// tdrStyle->SetStatStyle(Style_t style = 1001);
	// tdrStyle->SetStatX(Float_t x = 0);
	// tdrStyle->SetStatY(Float_t y = 0);
	
	// Margins:
	tdrStyle->SetPadTopMargin(0.05);
	tdrStyle->SetPadBottomMargin(0.13);
	tdrStyle->SetPadLeftMargin(0.13);
	tdrStyle->SetPadRightMargin(0.05);
	
	// For the Global title:
	
	tdrStyle->SetOptTitle(0);
	tdrStyle->SetTitleFont(42);
	tdrStyle->SetTitleColor(1);
	tdrStyle->SetTitleTextColor(1);
	tdrStyle->SetTitleFillColor(10);
	tdrStyle->SetTitleFontSize(0.05);
	// tdrStyle->SetTitleH(0); // Set the height of the title box
	// tdrStyle->SetTitleW(0); // Set the width of the title box
	// tdrStyle->SetTitleX(0); // Set the position of the title box
	// tdrStyle->SetTitleY(0.985); // Set the position of the title box
	// tdrStyle->SetTitleStyle(Style_t style = 1001);
	// tdrStyle->SetTitleBorderSize(2);
	
	// For the axis titles:
	
	tdrStyle->SetTitleColor(1, "XYZ");
	tdrStyle->SetTitleFont(42, "XYZ");
	tdrStyle->SetTitleSize(0.06, "XYZ");
	// tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
	// tdrStyle->SetTitleYSize(Float_t size = 0.02);
	tdrStyle->SetTitleXOffset(0.9);
	tdrStyle->SetTitleYOffset(1.05);
	// tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset
	
	// For the axis labels:
	
	tdrStyle->SetLabelColor(1, "XYZ");
	tdrStyle->SetLabelFont(42, "XYZ");
	tdrStyle->SetLabelOffset(0.007, "XYZ");
	tdrStyle->SetLabelSize(0.05, "XYZ");
	
	// For the axis:
	
	tdrStyle->SetAxisColor(1, "XYZ");
	tdrStyle->SetStripDecimals(kTRUE);
	tdrStyle->SetTickLength(0.03, "XYZ");
	tdrStyle->SetNdivisions(510, "XYZ");
	tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
	tdrStyle->SetPadTickY(1);
	
	// Change for log plots:
	tdrStyle->SetOptLogx(0);
	tdrStyle->SetOptLogy(ylog);
	tdrStyle->SetOptLogz(0);
	
	// Postscript options:
	
	//  tdrStyle->SetPaperSize(7.5,7.5);
	
	tdrStyle->SetPaperSize(15.,15.);
	
	//  tdrStyle->SetPaperSize(20.,20.);
	
	// tdrStyle->SetLineScalePS(Float_t scale = 3);
	// tdrStyle->SetLineStyleString(Int_t i, const char* text);
	// tdrStyle->SetHeaderPS(const char* header);
	// tdrStyle->SetTitlePS(const char* pstitle);
	
	// tdrStyle->SetBarOffset(Float_t baroff = 0.5);
	// tdrStyle->SetBarWidth(Float_t barwidth = 0.5);
	// tdrStyle->SetPaintTextFormat(const char* format = "g");
	// tdrStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
	// tdrStyle->SetTimeOffset(Double_t toffset);
	// tdrStyle->SetHistMinimumZero(kTRUE);
	
	tdrStyle->cd();
	
}
