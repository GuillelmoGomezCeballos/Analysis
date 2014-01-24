#include <TROOT.h>
#include "TMath.h"
#include "TH1D.h"
#include "Math/VectorUtil.h"
#include "Math/LorentzVector.h"
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector; 

double trilepton_info(int nsel, LorentzVector l1,LorentzVector l2,LorentzVector l3,
                      int q1, int q2, int q3, int ld1, int ld2, int ld3,
		      double mt1, double mt2, double mt3,
		      double met = 0.0, double metPhi = 0.0);

double DeltaR(double phi1, double eta1, double phi2, double eta2);

double DeltaR(double phi1, double eta1, double phi2, double eta2)
{
  // Compute DeltaR between two given points in the eta/phi plane.

  double dphi = DeltaPhi(phi1, phi2);
  double deta = eta1-eta2;
  double dR = TMath::Sqrt(dphi*dphi + deta*deta);
  return(dR);
}

double trilepton_info(int nsel, LorentzVector l1,LorentzVector l2,LorentzVector l3,
                      int q1, int q2, int q3, int ld1, int ld2, int ld3,
		      double mt1, double mt2, double mt3,
		      double met, double metPhi){

if     (nsel == 0 || nsel == 3 || nsel == 4 || nsel == 9) {
  double massz = 999.;
  double mass[3]={9999., 9999., 9999.};
  double mtw = 0.0;
  double type = 0.0;
  double which3rdl = -1;
  if(TMath::Abs(ld1) == TMath::Abs(ld2) && q1 != q2) mass[0] = (l1+l2).M();
  if(TMath::Abs(ld1) == TMath::Abs(ld3) && q1 != q3) mass[1] = (l1+l3).M();
  if(TMath::Abs(ld2) == TMath::Abs(ld3) && q2 != q3) mass[2] = (l2+l3).M();
  if(TMath::Abs(mass[0]-91.1876) < massz) {massz = TMath::Abs(mass[0]-91.1876);mtw=mt3;which3rdl=3;}
  if(TMath::Abs(mass[1]-91.1876) < massz) {massz = TMath::Abs(mass[1]-91.1876);mtw=mt2;which3rdl=2;}
  if(TMath::Abs(mass[2]-91.1876) < massz) {massz = TMath::Abs(mass[2]-91.1876);mtw=mt1;which3rdl=1;}
  /*
  if     (TMath::Abs(ld1) == 13 && TMath::Abs(ld2) == 13 && TMath::Abs(ld3) == 13) type = 1;
  else if(TMath::Abs(ld1) == 13 && TMath::Abs(ld2) == 13 && TMath::Abs(ld3) == 11) type = 2;
  else if(TMath::Abs(ld1) == 13 && TMath::Abs(ld2) == 11 && TMath::Abs(ld3) == 13) type = 2;
  else if(TMath::Abs(ld1) == 11 && TMath::Abs(ld2) == 13 && TMath::Abs(ld3) == 13) type = 2;
  else if(TMath::Abs(ld1) == 13 && TMath::Abs(ld2) == 11 && TMath::Abs(ld3) == 11) type = 3;
  else if(TMath::Abs(ld1) == 11 && TMath::Abs(ld2) == 11 && TMath::Abs(ld3) == 13) type = 3;
  else if(TMath::Abs(ld1) == 11 && TMath::Abs(ld2) == 13 && TMath::Abs(ld3) == 11) type = 3;
  else if(TMath::Abs(ld1) == 11 && TMath::Abs(ld2) == 11 && TMath::Abs(ld3) == 11) type = 4;
  */
  if     (TMath::Abs(ld1) == 13 && TMath::Abs(ld2) == 13 && TMath::Abs(ld3) == 13) type = 1;
  else if(TMath::Abs(ld1) == 13 && TMath::Abs(ld2) == 13 && TMath::Abs(ld3) == 11) type = 2;
  else if(TMath::Abs(ld1) == 13 && TMath::Abs(ld2) == 11 && TMath::Abs(ld3) == 13) type = 3;
  else if(TMath::Abs(ld1) == 11 && TMath::Abs(ld2) == 13 && TMath::Abs(ld3) == 13) type = 4;
  else if(TMath::Abs(ld1) == 13 && TMath::Abs(ld2) == 11 && TMath::Abs(ld3) == 11) type = 5;
  else if(TMath::Abs(ld1) == 11 && TMath::Abs(ld2) == 11 && TMath::Abs(ld3) == 13) type = 6;
  else if(TMath::Abs(ld1) == 11 && TMath::Abs(ld2) == 13 && TMath::Abs(ld3) == 11) type = 7;
  else if(TMath::Abs(ld1) == 11 && TMath::Abs(ld2) == 11 && TMath::Abs(ld3) == 11) type = 8;
  if     (nsel == 0) return massz;
  else if(nsel == 3) return mtw;
  else if(nsel == 4) return type;
  else if(nsel == 9) return which3rdl;
}
else if(nsel == 1) {
  double massmin = 999.;
  double mass[3]={9999., 9999., 9999.};
  if(q1 != q2) mass[0] = (l1+l2).M();
  if(q1 != q3) mass[1] = (l1+l3).M();
  if(q2 != q3) mass[2] = (l2+l3).M();
  if(mass[0] < massmin) massmin = mass[0];
  if(mass[1] < massmin) massmin = mass[1];
  if(mass[2] < massmin) massmin = mass[2];
  return massmin;
}
else if(nsel == 2 || nsel == 10) {
  double dr = 999.;
  double mtw = 0.0;
  if(q1 != q2 && DeltaR(l1.Phi(),l1.Eta(),l2.Phi(),l2.Eta()) < dr) {dr =  DeltaR(l1.Phi(),l1.Eta(),l2.Phi(),l2.Eta()); mtw=mt3;}
  if(q1 != q3 && DeltaR(l1.Phi(),l1.Eta(),l3.Phi(),l3.Eta()) < dr) {dr =  DeltaR(l1.Phi(),l1.Eta(),l3.Phi(),l3.Eta()); mtw=mt2;}
  if(q2 != q3 && DeltaR(l2.Phi(),l2.Eta(),l3.Phi(),l3.Eta()) < dr) {dr =  DeltaR(l2.Phi(),l2.Eta(),l3.Phi(),l3.Eta()); mtw=mt1;}
  if(nsel == 10) return mtw;
  return dr;
}
else if(nsel == 5 || nsel == 7 || nsel == 8) {
  double massminSFOS = 999.; double mtw = 0.0; double ptlw = 0.0;
  double mass[3]={9999., 9999., 9999.};
  if(TMath::Abs(ld1) == TMath::Abs(ld2) && q1 != q2) mass[0] = (l1+l2).M();
  if(TMath::Abs(ld1) == TMath::Abs(ld3) && q1 != q3) mass[1] = (l1+l3).M();
  if(TMath::Abs(ld2) == TMath::Abs(ld3) && q2 != q3) mass[2] = (l2+l3).M();
  if(mass[0] < massminSFOS) {massminSFOS = mass[0]; mtw=mt3; ptlw = l3.pt();}
  if(mass[1] < massminSFOS) {massminSFOS = mass[1]; mtw=mt2; ptlw = l2.pt();}
  if(mass[2] < massminSFOS) {massminSFOS = mass[2]; mtw=mt1; ptlw = l1.pt();}
  if     (nsel == 5) return massminSFOS;
  else if(nsel == 7) return mtw;
  else if(nsel == 8) return ptlw;
}
else if(nsel == 6) {
  double massminSFSS = 999.;
  double mass[3]={9999., 9999., 9999.};
  if(TMath::Abs(ld1) == TMath::Abs(ld2) && q1 == q2) mass[0] = (l1+l2).M();
  if(TMath::Abs(ld1) == TMath::Abs(ld3) && q1 == q3) mass[1] = (l1+l3).M();
  if(TMath::Abs(ld2) == TMath::Abs(ld3) && q2 == q3) mass[2] = (l2+l3).M();
  if(mass[0] < massminSFSS) massminSFSS = mass[0];
  if(mass[1] < massminSFSS) massminSFSS = mass[1];
  if(mass[2] < massminSFSS) massminSFSS = mass[2];
  return massminSFSS;
}
else if(nsel > 10) {
  if     ((nsel == 61 || nsel == 62) && TMath::Abs((l1+l2).M() - 91.1876) < 15.0) {
    double metx = met*cos(metPhi)-l3.px();
    double mety = met*sin(metPhi)-l3.py();
    double newMT = sqrt(2.0*(l1+l2).pt()*sqrt(metx*metx+mety*mety)*(1.0-cos(DeltaPhi((l1+l2).phi(),TMath::ATan2(mety,metx)))));
    if     (nsel == 61) return sqrt(metx*metx+mety*mety);
    else if(nsel == 62) return newMT;
  }
  else if((nsel == 61 || nsel == 62) && TMath::Abs((l1+l3).M() - 91.1876) < 15.0) {
    double metx = met*cos(metPhi)-l2.px();
    double mety = met*sin(metPhi)-l2.py();
    double newMT = sqrt(2.0*(l1+l3).pt()*sqrt(metx*metx+mety*mety)*(1.0-cos(DeltaPhi((l1+l3).phi(),TMath::ATan2(mety,metx)))));
    if     (nsel == 61) return sqrt(metx*metx+mety*mety);
    else if(nsel == 62) return newMT;
  }
  else if((nsel == 61 || nsel == 62) && TMath::Abs((l2+l3).M() - 91.1876) < 15.0) {
    double metx = met*cos(metPhi)-l1.px();
    double mety = met*sin(metPhi)-l1.py();
    double newMT = sqrt(2.0*(l2+l3).pt()*sqrt(metx*metx+mety*mety)*(1.0-cos(DeltaPhi((l2+l3).phi(),TMath::ATan2(mety,metx)))));
    if     (nsel == 61) return sqrt(metx*metx+mety*mety);
    else if(nsel == 62) return newMT;
  }
}
return -10.0;
}
