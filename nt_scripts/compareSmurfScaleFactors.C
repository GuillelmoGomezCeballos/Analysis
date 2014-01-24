
void printline(TH2F* h2)
{

    printf("\n");
    printf("%s\n", h2->GetName());
    printf("\t\t & ");
    for (unsigned int y = 1; y <= h2->GetYaxis()->GetNbins(); ++y)
    {
        Float_t min = h2->GetYaxis()->GetBinLowEdge(y);
        Float_t max = min + h2->GetYaxis()->GetBinWidth(y);
        printf("  %4.1f - %4.1f  & ", min, max);
    }
    printf("\n");

    for (unsigned int x = 1; x <= h2->GetXaxis()->GetNbins(); ++x) {

        Float_t min = h2->GetXaxis()->GetBinLowEdge(x);
        Float_t max = min + h2->GetXaxis()->GetBinWidth(x);
        printf("  %4.1f - %4.1f  & ", min, max);

        for (unsigned int y = 1; y <= h2->GetYaxis()->GetNbins(); ++y)
        {
            Float_t eff = h2->GetBinContent(x, y);
            Float_t err = h2->GetBinError(x, y);
            if (y == h2->GetYaxis()->GetNbins())
                printf("\t%4.4f $\\pm$ %4.4f \\\\", eff, err);
            else
                printf("\t%4.4f $\\pm$ %4.4f & ", eff, err);
        }
        printf("\n");
    }
}

void compareSmurfScaleFactors()
{


    TFile f_new("/data/smurf/dlevans/Efficiencies/V00-02-07_trigNameFix_HCP_V1/summary.root");
    TFile f_old("/data/smurf/dlevans/Efficiencies/V00-02-06_V1/summary.root");

    //gROOT->cd();
    //gROOT->ProcessLine(".L ~/tdrStyle.C");
    //gROOT->ProcessLine("setTDRStyle()");
    //gROOT->ForceStyle();
    gStyle->SetPalette(1);

    //
    // muons
    //

    TH2F *ratio_h2_results_muon_selection           = (TH2F*)f_new.Get("h2_results_muon_selection");
    TH2F *ratio_h2_results_muon_single              = (TH2F*)f_new.Get("h2_results_muon_single");
    TH2F *ratio_h2_results_muon_double_leadingleg   = (TH2F*)f_new.Get("h2_results_muon_double_leadingleg");
    TH2F *ratio_h2_results_muon_double_trailingleg  = (TH2F*)f_new.Get("h2_results_muon_double_trailingleg");

    TH2F *old_h2_results_muon_selection           = (TH2F*)f_old.Get("h2_results_muon_selection");
    TH2F *old_h2_results_muon_single              = (TH2F*)f_old.Get("h2_results_muon_single");
    TH2F *old_h2_results_muon_double_leadingleg   = (TH2F*)f_old.Get("h2_results_muon_double_leadingleg");
    TH2F *old_h2_results_muon_double_trailingleg  = (TH2F*)f_old.Get("h2_results_muon_double_trailingleg");

    ratio_h2_results_muon_selection->Divide(old_h2_results_muon_selection);
    ratio_h2_results_muon_single->Divide(old_h2_results_muon_single);
    ratio_h2_results_muon_double_leadingleg->Divide(old_h2_results_muon_double_leadingleg);
    ratio_h2_results_muon_double_trailingleg->Divide(old_h2_results_muon_double_trailingleg);

    printline(ratio_h2_results_muon_selection);
    printline(ratio_h2_results_muon_single);
    printline(ratio_h2_results_muon_double_leadingleg);
    printline(ratio_h2_results_muon_double_trailingleg);


    TCanvas *c1 = new TCanvas();
    c1->cd();
    c1->SetRightMargin(0.15);
    ratio_h2_results_muon_selection->Draw(" COLZ");
    c1->SaveAs("ratio_h2_results_muon_selection.png");
    ratio_h2_results_muon_single->Draw(" COLZ");
    c1->SaveAs("ratio_h2_results_muon_single.png");
    ratio_h2_results_muon_double_leadingleg->Draw(" COLZ");
    c1->SaveAs("ratio_h2_results_muon_leadingleg.png");
    ratio_h2_results_muon_double_trailingleg->Draw(" COLZ");
    c1->SaveAs("ratio_h2_results_muon_trailingleg.png");

    //
    // electrons
    //

    TH2F *ratio_h2_results_electron_selection           = (TH2F*)f_new.Get("h2_results_electron_selection");
    TH2F *ratio_h2_results_electron_single              = (TH2F*)f_new.Get("h2_results_electron_single");
    TH2F *ratio_h2_results_electron_double_leadingleg   = (TH2F*)f_new.Get("h2_results_electron_double_leadingleg");
    TH2F *ratio_h2_results_electron_double_trailingleg  = (TH2F*)f_new.Get("h2_results_electron_double_trailingleg");

    TH2F *old_h2_results_electron_selection           = (TH2F*)f_old.Get("h2_results_electron_selection");
    TH2F *old_h2_results_electron_single              = (TH2F*)f_old.Get("h2_results_electron_single");
    TH2F *old_h2_results_electron_double_leadingleg   = (TH2F*)f_old.Get("h2_results_electron_double_leadingleg");
    TH2F *old_h2_results_electron_double_trailingleg  = (TH2F*)f_old.Get("h2_results_electron_double_trailingleg");

    ratio_h2_results_electron_selection->Divide(old_h2_results_electron_selection);
    ratio_h2_results_electron_single->Divide(old_h2_results_electron_single);
    ratio_h2_results_electron_double_leadingleg->Divide(old_h2_results_electron_double_leadingleg);
    ratio_h2_results_electron_double_trailingleg->Divide(old_h2_results_electron_double_trailingleg);

    printline(ratio_h2_results_electron_selection);
    printline(ratio_h2_results_electron_single);
    printline(ratio_h2_results_electron_double_leadingleg);
    printline(ratio_h2_results_electron_double_trailingleg);

    c1->cd();
    ratio_h2_results_electron_selection->Draw(" COLZ");
    c1->SaveAs("ratio_h2_results_electron_selection.png");
    ratio_h2_results_electron_single->Draw(" COLZ");
    c1->SaveAs("ratio_h2_results_electron_single.png");
    ratio_h2_results_electron_double_leadingleg->Draw(" COLZ");
    c1->SaveAs("ratio_h2_results_electron_leadingleg.png");
    ratio_h2_results_electron_double_trailingleg->Draw(" COLZ");
    c1->SaveAs("ratio_h2_results_electron_trailingleg.png");

}


