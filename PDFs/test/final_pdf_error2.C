#include <iostream>
#include <fstream>
#include "TMath.h"
#include "TString.h"

void mstw_cteq(double *iA, int nPoints, double results[2]) { 
  double lVal = iA[0];
  double lDiffP = 0;
  double lDiffM = 0;
  for(int i0 = 1; i0 < nPoints; i0++) {
    double pDiff = (iA[i0] - lVal);
    if(i0 % 2 == 1)lDiffP += pDiff*pDiff;
    if(i0 % 2 == 0)lDiffM += pDiff*pDiff;
  }
  printf("MSTW/CTEQ syst. => + %8.5f - %8.5f\n",sqrt(lDiffP)*100/lVal,sqrt(lDiffM)*100/lVal);
  results[0] = sqrt(lDiffP)*100/lVal;
  results[1] = sqrt(lDiffM)*100/lVal;
}

void nnpdf(double *iA, int nPoints, double results[2]) { 
  double lVal = iA[0]; 
  double lDiffP = 0; int lNP = 0;
  double lDiffM = 0; int lNM = 0;
  for(int i0 = 1; i0 < nPoints; i0++) {
    double pDiff = (iA[i0] - lVal);
    if(pDiff > 0) {lDiffP += pDiff*pDiff; lNP++;}
    if(pDiff < 0) {lDiffM += pDiff*pDiff; lNM++;}
  }
  double lEMinus = sqrt(lNM/(lNM-1.)*lDiffM/lNM);
  double lEPlus  = sqrt(lNP/(lNP-1.)*lDiffP/lNP);
  printf("NNPDF     syst. => + %8.5f - %8.5f\n",(lEPlus/lVal)*100,(lEMinus/lVal)*100);
  results[0] = (lEPlus/lVal)*100;
  results[1] = (lEMinus/lVal)*100;
}

void combined_syst(double *iA, int nPoints, bool nnpdfFlag, double results[2]) {
  unsigned int npairs = (nPoints-1)/2;
  double events_central = iA[0];
  double wplus = 0.;
  double wminus = 0.;
  unsigned int nplus = 0;
  unsigned int nminus = 0;
  for (unsigned int j=0; j<npairs; ++j) {
    double wa = iA[2*j+1]/events_central-1.;
    double wb = iA[2*j+2]/events_central-1.; 
    if (nnpdfFlag) {
      if (wa>0.) {
    	wplus += wa*wa; 
    	nplus++;
      } else {
    	wminus += wa*wa;
    	nminus++;
      }
      if (wb>0.) {
    	wplus += wb*wb; 
    	nplus++;
      } else {
    	wminus += wb*wb;
    	nminus++;
      }
    } else {
      if (wa>wb) {
        if (wa<0.) wa = 0.;
        if (wb>0.) wb = 0.;
        wplus += wa*wa;
        wminus += wb*wb;
      } else {
    	if (wb<0.) wb = 0.;
    	if (wa>0.) wa = 0.;
    	wplus += wb*wb;
    	wminus += wa*wa;
      }
    }
  }
  if (wplus>0) wplus = sqrt(wplus);
  if (wminus>0) wminus = sqrt(wminus);
  if (nnpdfFlag) {
    if (nplus>0) wplus /= sqrt(nplus);
    if (nminus>0) wminus /= sqrt(nminus);
  }
  printf("PDF syst.     => + %8.5f - %8.5f\n",100.*wplus,100.*wminus);
  results[0] = wplus*100;
  results[1] = wminus*100;
}

void combined_syst_sigma(double *iA, int nPoints, bool nnpdfFlag, double results[2]) {
  unsigned int npairs = (nPoints-1)/2;
  double events_central = iA[0];
  double wplus = 0.;
  double wminus = 0.;
  unsigned int nplus = 0;
  unsigned int nminus = 0;
  for (unsigned int j=0; j<npairs; ++j) {
    double wa = iA[2*j+1]/events_central-1.;
    double wb = iA[2*j+2]/events_central-1.;
    if (nnpdfFlag) {
      if (wa>0.) {
        wplus += wa*wa;
        nplus++;
      } else {
        wminus += wa*wa;
        nminus++;
      }
      if (wb>0.) {
        wplus += wb*wb;
        nplus++;
      } else {
        wminus += wb*wb;
        nminus++;
      }
    } else {
      if (wa>wb) {
        if (wa<0.) wa = 0.;
        if (wb>0.) wb = 0.;
        wplus += wa*wa;
        wminus += wb*wb;
      } else {
        if (wb<0.) wb = 0.;
        if (wa>0.) wa = 0.;
        wplus += wb*wb;
        wminus += wa*wa;
      }
    }
  }
  if (wplus>0) wplus = sqrt(wplus);
  if (wminus>0) wminus = sqrt(wminus);
  if (nnpdfFlag) {
    if (nplus>0) wplus /= sqrt(nplus);
    if (nminus>0) wminus /= sqrt(nminus);
  }
  printf("PDF syst.     => + %8.5f - %8.5f\n",100.*wplus,100.*wminus);
  results[0] = (wplus+wminus)/2*events_central;
  results[1] = 1;
}

void final_pdf_error2(const char* suffix = "") 
{
TString cteq66txt=TString("pdf_cteq66")+TString(suffix)+TString(".txt");
TString nnpdftxt=TString("pdf_nnpdf")+TString(suffix)+TString(".txt");
TString mstwtxt=TString("pdf_mstw")+TString(suffix)+TString(".txt");

TString cteq66astxt=TString("pdf_cteq66_alphas")+TString(suffix)+TString(".txt");
TString nnpdfastxt=TString("pdf_nnpdf_alphas")+TString(suffix)+TString(".txt");
TString mstwastxt=TString("pdf_mstw_alphas")+TString(suffix)+TString(".txt");

const int nPoints0 = 45;
double pdf00[nPoints0],pdf01[nPoints0];

const int nPoints1 = 101;
double pdf10[nPoints1],pdf11[nPoints1];

const int nPoints2 = 41;
double pdf20[nPoints2],pdf21[nPoints2];

const int nPoints3 = 2;
double pdf30[nPoints3],pdf31[nPoints3];

const int nPoints4 = 2;
double pdf40[nPoints4],pdf41[nPoints4];

const int nPoints5 = 2;
double pdf50[nPoints5],pdf51[nPoints5];

int i = 0;
ifstream infile0(cteq66txt.Data());
while (infile0>>pdf00[i]
              >>pdf01[i]){ i++;}

i = 0;
ifstream infile1(nnpdftxt.Data());
while (infile1>>pdf10[i]
              >>pdf11[i]){ i++;}
	     
i = 0;
ifstream infile2(mstwtxt.Data());
while (infile2>>pdf20[i]
              >>pdf21[i]){ i++;}
	     
i = 0;
ifstream infile3(cteq66astxt.Data());
while (infile3>>pdf30[i]
              >>pdf31[i]){ i++;}

i = 0;
ifstream infile4(nnpdfastxt.Data());
while (infile4>>pdf40[i]
              >>pdf41[i]){ i++;}
	     
i = 0;
ifstream infile5(mstwastxt.Data());
while (infile5>>pdf50[i]
              >>pdf51[i]){ i++;}

//combined_syst and mstw_cteq/nnpdf represent distinct methods of determining the uncertainties due to parameter variation for a given PDF set. The difference hinges on floor/ceiling values for specific up/down variations. See code.

double results00[2],results10[2],results20[2];
mstw_cteq(pdf00, nPoints0, results00);
nnpdf    (pdf10, nPoints1, results10);
mstw_cteq(pdf20, nPoints2, results20);
double results01[2],results11[2],results21[2];
mstw_cteq(pdf01, nPoints0, results01);
nnpdf    (pdf11, nPoints1, results11);
mstw_cteq(pdf21, nPoints2, results21);

//combined_syst_sigma puts sigma values on the acceptance in the output arrays resultsN00, resultsN10, resultsN20

double resultsN00[2],resultsN10[2],resultsN20[2];
combined_syst_sigma(pdf00, nPoints0, false, resultsN00);
combined_syst_sigma(pdf10, nPoints1, true , resultsN10);
combined_syst_sigma(pdf20, nPoints2, false, resultsN20);
double resultsN01[2],resultsN11[2],resultsN21[2];
combined_syst_sigma(pdf01, nPoints0, false, resultsN01);
combined_syst_sigma(pdf11, nPoints1, true , resultsN11);
combined_syst_sigma(pdf21, nPoints2, false, resultsN21);

/*
printf("CTEQ  Acc syst. => + %8.5f - %8.5f\n",results00[0],results00[1]);
printf("NNPDF Acc syst. => + %8.5f - %8.5f\n",results10[0],results10[1]);
printf("MSTW  Acc syst. => + %8.5f - %8.5f\n",results20[0],results20[1]);
printf("CTEQ  Rec syst. => + %8.5f - %8.5f\n",results01[0],results01[1]);
printf("NNPDF Rec syst. => + %8.5f - %8.5f\n",results11[0],results11[1]);
printf("MSTW  Rec syst. => + %8.5f - %8.5f\n",results21[0],results21[1]);
*/
/*printf("CTEQ  Acc syst. => + %8.5f - %8.5f\n",resultsN00[0],resultsN00[1]);
printf("NNPDF Acc syst. => + %8.5f - %8.5f\n",resultsN10[0],resultsN10[1]);
printf("MSTW  Acc syst. => + %8.5f - %8.5f\n",resultsN20[0],resultsN20[1]);
printf("CTEQ  Rec syst. => + %8.5f - %8.5f\n",resultsN01[0],resultsN01[1]);
printf("NNPDF Rec syst. => + %8.5f - %8.5f\n",resultsN11[0],resultsN11[1]);
printf("MSTW  Rec syst. => + %8.5f - %8.5f\n",resultsN21[0],resultsN21[1]);
*/
double alpha_s_acc[3] = {100*(pdf30[0]-pdf30[1])/2.0/pdf00[0],100*(pdf40[0]-pdf40[1])/2.0/pdf10[0],100*(pdf50[0]-pdf50[1])/2.0/pdf20[0]};
double alpha_s_rec[3] = {100*(pdf31[0]-pdf31[1])/2.0/pdf01[0],100*(pdf41[0]-pdf41[1])/2.0/pdf11[0],100*(pdf51[0]-pdf51[1])/2.0/pdf21[0]};

printf("CTEQ/NNPDF/MSTW  Acc-alphas syst. => %8.5f  %8.5f  %8.5f\n",alpha_s_acc[0],alpha_s_acc[1],alpha_s_acc[2]);
printf("CTEQ/NNPDF/MSTW  Rec-alphas syst. => %8.5f  %8.5f  %8.5f\n",alpha_s_rec[0],alpha_s_rec[1],alpha_s_rec[2]);

//Values in this array are not normalized to the actual acceptance value
double alpha_s_acc_sig[3] = {(pdf30[0]-pdf30[1])/2.0,(pdf40[0]-pdf40[1])/2.0,(pdf50[0]-pdf50[1])/2.0};

//Values in sigmas_all represent the total combined sigma for variations in alpha_s and all other parameters. Note the addition in quadrature.
double sigmas_all[3]= {TMath::Sqrt(alpha_s_acc_sig[0]*alpha_s_acc_sig[0]+resultsN00[0]*resultsN00[0]),TMath::Sqrt(alpha_s_acc_sig[1]*alpha_s_acc_sig[1]+resultsN10[0]*resultsN10[0]),TMath::Sqrt(alpha_s_acc_sig[2]*alpha_s_acc_sig[2]+resultsN20[0]*resultsN20[0])};

//Determine the maximum acceptance for one sigma variation
double Vmax_acc = pdf00[0]+sigmas_all[0];
if( (pdf10[0]+sigmas_all[1]) > Vmax_acc) Vmax_acc = pdf10[0]+sigmas_all[1];
 if( (pdf20[0]+sigmas_all[2])  > Vmax_acc) Vmax_acc = pdf20[0]+sigmas_all[2];

 //Determine the minimum acceptance for one sigma variation
double Vmin_acc = pdf00[0]-sigmas_all[0];
 if( (pdf10[0]-sigmas_all[1]) < Vmin_acc) Vmin_acc = pdf10[0]-sigmas_all[1];
 if((pdf20[0]-sigmas_all[2])  < Vmin_acc) Vmin_acc = pdf20[0]-sigmas_all[2];

double sigma_acc     = (Vmax_acc-Vmin_acc)/2.;
double x_central_acc = (Vmax_acc+Vmin_acc)/2.;

 std::cout<<"Central Value: " << x_central_acc <<std::endl;
 std::cout<<"Sigma: " <<sigma_acc<<std::endl;
 std::cout<<"Combined Uncertainty (PDF4LHC) : " << 100*sigma_acc/x_central_acc<<std::endl;

}
