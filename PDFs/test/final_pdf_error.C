#include <iostream>
#include <fstream>

void mstw_cteq(double *iA, int nPoints, double results[2]) { 
  double lVal = iA[0];
  double lDiffP = 0;
  double lDiffM = 0;
  for(int i0 = 1; i0 < nPoints; i0++) {
    double pDiff = (iA[i0] - lVal);
    if(i0 % 2 == 1)lDiffP += pDiff*pDiff;
    if(i0 % 2 == 0)lDiffM += pDiff*pDiff;
  }
  printf("MSTW/CTEQ syst. => + %5.3f - %5.3f\n",sqrt(lDiffP)*100/lVal,sqrt(lDiffM)*100/lVal);
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
  printf("NNPDF     syst. => + %5.3f - %5.3f\n",(lEPlus/lVal)*100,(lEMinus/lVal)*100);
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
  printf("PDF syst.     => + %5.3f - %5.3f\n",100.*wplus,100.*wminus);
  results[0] = wplus*100;
  results[1] = wminus*100;
}

void final_pdf_error(){
const int nPoints0 = 53;
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
ifstream infile0("pdf_cteq66.txt");
while (infile0>>pdf00[i]
              >>pdf01[i]){ i++;}

i = 0;
ifstream infile1("pdf_nnpdf.txt");
while (infile1>>pdf10[i]
              >>pdf11[i]){ i++;}
	     
i = 0;
ifstream infile2("pdf_mstw.txt");
while (infile2>>pdf20[i]
              >>pdf21[i]){ i++;}
	     
i = 0;
ifstream infile3("pdf_cteq66_alphas.txt");
while (infile3>>pdf30[i]
              >>pdf31[i]){ i++;}

i = 0;
ifstream infile4("pdf_nnpdf_alphas.txt");
while (infile4>>pdf40[i]
              >>pdf41[i]){ i++;}
	     
i = 0;
ifstream infile5("pdf_mstw_alphas.txt");
while (infile5>>pdf50[i]
              >>pdf51[i]){ i++;}

double results00[2],results10[2],results20[2];
mstw_cteq(pdf00, nPoints0, results00);
nnpdf    (pdf10, nPoints1, results10);
mstw_cteq(pdf20, nPoints2, results20);
double results01[2],results11[2],results21[2];
mstw_cteq(pdf01, nPoints0, results01);
nnpdf    (pdf11, nPoints1, results11);
mstw_cteq(pdf21, nPoints2, results21);

double resultsN00[2],resultsN10[2],resultsN20[2];
combined_syst(pdf00, nPoints0, false, resultsN00);
combined_syst(pdf10, nPoints1, true , resultsN10);
combined_syst(pdf20, nPoints2, false, resultsN20);
double resultsN01[2],resultsN11[2],resultsN21[2];
combined_syst(pdf01, nPoints0, false, resultsN01);
combined_syst(pdf11, nPoints1, true , resultsN11);
combined_syst(pdf21, nPoints2, false, resultsN21);

printf("CTEQ  Acc syst. => + %5.3f - %5.3f\n",results00[0],results00[1]);
printf("NNPDF Acc syst. => + %5.3f - %5.3f\n",results10[0],results10[1]);
printf("MSTW  Acc syst. => + %5.3f - %5.3f\n",results20[0],results20[1]);
printf("CTEQ  Rec syst. => + %5.3f - %5.3f\n",results01[0],results01[1]);
printf("NNPDF Rec syst. => + %5.3f - %5.3f\n",results11[0],results11[1]);
printf("MSTW  Rec syst. => + %5.3f - %5.3f\n",results21[0],results21[1]);

printf("CTEQ  Acc syst. => + %5.3f - %5.3f\n",resultsN00[0],resultsN00[1]);
printf("NNPDF Acc syst. => + %5.3f - %5.3f\n",resultsN10[0],resultsN10[1]);
printf("MSTW  Acc syst. => + %5.3f - %5.3f\n",resultsN20[0],resultsN20[1]);
printf("CTEQ  Rec syst. => + %5.3f - %5.3f\n",resultsN01[0],resultsN01[1]);
printf("NNPDF Rec syst. => + %5.3f - %5.3f\n",resultsN11[0],resultsN11[1]);
printf("MSTW  Rec syst. => + %5.3f - %5.3f\n",resultsN21[0],resultsN21[1]);

double alpha_s_acc[3] = {100*(pdf30[0]-pdf30[1])/2.0/pdf00[0],100*(pdf40[0]-pdf40[1])/2.0/pdf10[0],100*(pdf50[0]-pdf50[1])/2.0/pdf20[0]};
double alpha_s_rec[3] = {100*(pdf31[0]-pdf31[1])/2.0/pdf01[0],100*(pdf41[0]-pdf41[1])/2.0/pdf11[0],100*(pdf51[0]-pdf51[1])/2.0/pdf21[0]};

printf("CTEQ/NNPDF/MSTW  Acc-alphas syst. => %5.3f  %5.3f  %5.3f\n",alpha_s_acc[0],alpha_s_acc[1],alpha_s_acc[2]);
printf("CTEQ/NNPDF/MSTW  Rec-alphas syst. => %5.3f  %5.3f  %5.3f\n",alpha_s_rec[0],alpha_s_rec[1],alpha_s_rec[2]);

double Vmax_acc = (1.0+resultsN00[0]/100.)*pdf00[0];
if((1.0+resultsN10[0]/100.)*pdf10[0] > Vmax_acc) Vmax_acc = (resultsN10[0]/100.+1)*pdf10[0];
if((1.0+resultsN20[0]/100.)*pdf20[0] > Vmax_acc) Vmax_acc = (resultsN20[0]/100.+1)*pdf20[0];
double Vmin_acc = (1.0-resultsN00[1]/100.)*pdf00[0];
if((1.0-resultsN10[0]/100.)*pdf10[0] < Vmin_acc) Vmin_acc = (1.0-resultsN10[0]/100.)*pdf10[0];
if((1.0-resultsN20[0]/100.)*pdf20[0] < Vmin_acc) Vmin_acc = (1.0-resultsN20[0]/100.)*pdf20[0];

double sigma_acc     = (Vmax_acc-Vmin_acc)/2.;
double x_central_acc = (Vmax_acc+Vmin_acc)/2.;

printf("Acc-additional syst. x,sigma: %7.5f +/- %7.5f ==> %7.5f\n",x_central_acc,sigma_acc,100*sigma_acc/x_central_acc);

double Vmax_rec = (1.0+resultsN01[0]/100.)*pdf01[0];
if((1.0+resultsN11[0]/100.)*pdf11[0] > Vmax_rec) Vmax_rec = (resultsN11[0]/100.+1)*pdf11[0];
if((1.0+resultsN21[0]/100.)*pdf21[0] > Vmax_rec) Vmax_rec = (resultsN21[0]/100.+1)*pdf21[0];
double Vmin_rec = (1.0-resultsN01[1]/100.)*pdf01[0];
if((1.0-resultsN11[0]/100.)*pdf11[0] < Vmin_rec) Vmin_rec = (1.0-resultsN11[0]/100.)*pdf11[0];
if((1.0-resultsN21[0]/100.)*pdf21[0] < Vmin_rec) Vmin_rec = (1.0-resultsN21[0]/100.)*pdf21[0];

double sigma_rec     = (Vmax_rec-Vmin_rec)/2.;
double x_central_rec = (Vmax_rec+Vmin_rec)/2.;

printf("Rec-additional syst. x,sigma: %7.5f +/- %7.5f ==> %7.5f\n",x_central_rec,sigma_rec,100*sigma_rec/x_central_rec);

double total0[2] = {sqrt((resultsN00[0]+resultsN00[1])*(resultsN00[0]+resultsN00[1])/4.+alpha_s_acc[0]*alpha_s_acc[0]+(100*sigma_acc/x_central_acc)*(100*sigma_acc/x_central_acc)),
                    sqrt((resultsN01[0]+resultsN01[1])*(resultsN01[0]+resultsN01[1])/4.+alpha_s_rec[0]*alpha_s_rec[0]+(100*sigma_rec/x_central_rec)*(100*sigma_rec/x_central_rec))};
printf("Acc-total syst CTEQ: %5.3f\n",total0[0]);
printf("Rec-total syst CTEQ: %5.3f\n",total0[1]);

double total1[2] = {sqrt((resultsN10[0]+resultsN10[1])*(resultsN10[0]+resultsN10[1])/4.+alpha_s_acc[1]*alpha_s_acc[1]+(100*sigma_acc/x_central_acc)*(100*sigma_acc/x_central_acc)),
                    sqrt((resultsN11[0]+resultsN11[1])*(resultsN11[0]+resultsN11[1])/4.+alpha_s_rec[1]*alpha_s_rec[1]+(100*sigma_rec/x_central_rec)*(100*sigma_rec/x_central_rec))};
printf("Acc-total syst MSTW: %5.3f\n",total1[0]);
printf("Rec-total syst MSTW: %5.3f\n",total1[1]);

double total2[2] = {sqrt((resultsN20[0]+resultsN20[1])*(resultsN20[0]+resultsN20[1])/4.+alpha_s_acc[2]*alpha_s_acc[2]+(100*sigma_acc/x_central_acc)*(100*sigma_acc/x_central_acc)),
                    sqrt((resultsN21[0]+resultsN21[1])*(resultsN21[0]+resultsN21[1])/4.+alpha_s_rec[2]*alpha_s_rec[2]+(100*sigma_rec/x_central_rec)*(100*sigma_rec/x_central_rec))};
printf("Acc-total syst NNPDF: %5.3f\n",total2[0]);
printf("Rec-total syst NNPDF: %5.3f\n",total2[1]);


}
