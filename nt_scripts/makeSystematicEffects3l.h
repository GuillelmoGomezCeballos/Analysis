#include <TROOT.h>
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "Math/VectorUtil.h"
#include "Math/LorentzVector.h"
#include "TRandom.h"
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector; 

void makeSystematicEffects3l(int lid1_, int lid2_, int lid3_, LorentzVector lep1_, LorentzVector lep2_, LorentzVector lep3_, LorentzVector dilep_, 
                             double mt_, double met_, double metPhi_, 
                             double trackMet_, double trackMetPhi_, double njets_, LorentzVector jet1_, LorentzVector jet2_,
			     int year, int nsel, double outputVar[16]);

void makeSystematicEffects3l(int lid1_, int lid2_, int lid3_, LorentzVector lep1_, LorentzVector lep2_, LorentzVector lep3_, LorentzVector dilep_, 
                             double mt_, double met_, double metPhi_, 
                             double trackMet_, double trackMetPhi_, double njets_, LorentzVector jet1_, LorentzVector jet2_,
			     int year, int nsel, double outputVar[16]){
  if (year != 2011 && year != 2012) assert(0);
  double rndMon2012[12] = {gRandom->Gaus(0.00,0.010),gRandom->Gaus(0.00,0.017),gRandom->Gaus(0.00,0.015),gRandom->Gaus(0.00,0.030),
  			   gRandom->Gaus(0.00,0.010),gRandom->Gaus(0.00,0.017),gRandom->Gaus(0.00,0.015),gRandom->Gaus(0.00,0.030),
  			   gRandom->Gaus(0.00,0.010),gRandom->Gaus(0.00,0.017),gRandom->Gaus(0.00,0.015),gRandom->Gaus(0.00,0.030)};
  double rndMon2011[12] = {gRandom->Gaus(0.00,0.010),gRandom->Gaus(0.00,0.010),gRandom->Gaus(0.00,0.020),gRandom->Gaus(0.00,0.006),
  			   gRandom->Gaus(0.00,0.010),gRandom->Gaus(0.00,0.010),gRandom->Gaus(0.00,0.020),gRandom->Gaus(0.00,0.006),
  			   gRandom->Gaus(0.00,0.010),gRandom->Gaus(0.00,0.010),gRandom->Gaus(0.00,0.020),gRandom->Gaus(0.00,0.006)};

  double lep1pt,lep2pt,dilmass,dilpt,met,metPhi,trackMet,trackMetPhi,mt,dPhiDiLepMET,dPhiMETTrkMET,pTFrac,mtH,m2l2j,m2j,lep3pt;
  double llPhi = 0.0;
  if(nsel == 0){ // momentum scale +
    double corr[3] = {1.0, 1.0, 1.0};
    if(year == 2012){
      if     (TMath::Abs(lid1_) == 13 && TMath::Abs(lep1_.eta()) <  1.479){
	corr[0] = 1./0.99920 + rndMon2012[0];
      }
      else if(TMath::Abs(lid1_) == 13 && TMath::Abs(lep1_.eta()) >= 1.479){
	corr[0] = 1./0.99934 + rndMon2012[1];
      }
      else if(TMath::Abs(lid1_) == 11 && TMath::Abs(lep1_.eta()) <  1.479){
	corr[0] = 1./0.99807 + rndMon2012[2];
      }
      else if(TMath::Abs(lid1_) == 11 && TMath::Abs(lep1_.eta()) >= 1.479){
	corr[0] = 1./0.99952 + rndMon2012[3];
      }
      if     (TMath::Abs(lid2_) == 13 && TMath::Abs(lep2_.eta()) <  1.479){
	corr[1] = 1./0.99920 + rndMon2012[4];
      }
      else if(TMath::Abs(lid2_) == 13 && TMath::Abs(lep2_.eta()) >= 1.479){
	corr[1] = 1./0.99934 + rndMon2012[5];
      }
      else if(TMath::Abs(lid2_) == 11 && TMath::Abs(lep2_.eta()) <  1.479){
	corr[1] = 1./0.99807 + rndMon2012[6];
      }
      else if(TMath::Abs(lid2_) == 11 && TMath::Abs(lep2_.eta()) >= 1.479){
	corr[1] = 1./0.99952 + rndMon2012[7];
      }
      if     (TMath::Abs(lid3_) == 13 && TMath::Abs(lep3_.eta()) <  1.479){
	corr[2] = 1./0.99920 + rndMon2012[8];
      }
      else if(TMath::Abs(lid3_) == 13 && TMath::Abs(lep3_.eta()) >= 1.479){
	corr[2] = 1./0.99934 + rndMon2012[9];
      }
      else if(TMath::Abs(lid3_) == 11 && TMath::Abs(lep3_.eta()) <  1.479){
	corr[2] = 1./0.99807 + rndMon2012[10];
      }
      else if(TMath::Abs(lid3_) == 11 && TMath::Abs(lep3_.eta()) >= 1.479){
	corr[2] = 1./0.99952 + rndMon2012[11];
      }
    } else {
      if     (TMath::Abs(lid1_) == 13 && TMath::Abs(lep1_.eta()) <  1.479){
	corr[0] = 1.01 + rndMon2011[0];
      }
      else if(TMath::Abs(lid1_) == 13 && TMath::Abs(lep1_.eta()) >= 1.479){
	corr[0] = 1.01 + rndMon2011[1];
      }
      else if(TMath::Abs(lid1_) == 11 && TMath::Abs(lep1_.eta()) <  1.479){
	corr[0] = 1.01 + rndMon2011[2];
      }
      else if(TMath::Abs(lid1_) == 11 && TMath::Abs(lep1_.eta()) >= 1.479){
	corr[0] = 1.06 + rndMon2011[3];
      }
      if     (TMath::Abs(lid2_) == 13 && TMath::Abs(lep2_.eta()) <  1.479){
	corr[1] = 1.01 + rndMon2011[4];
      }
      else if(TMath::Abs(lid2_) == 13 && TMath::Abs(lep2_.eta()) >= 1.479){
	corr[1] = 1.01 + rndMon2011[5];
      }
      else if(TMath::Abs(lid2_) == 11 && TMath::Abs(lep2_.eta()) <  1.479){
	corr[1] = 1.01 + rndMon2011[6];
      }
      else if(TMath::Abs(lid2_) == 11 && TMath::Abs(lep2_.eta()) >= 1.479){
	corr[1] = 1.06 + rndMon2011[7];
      }
      if     (TMath::Abs(lid3_) == 13 && TMath::Abs(lep3_.eta()) <  1.479){
	corr[2] = 1.01 + rndMon2011[8];
      }
      else if(TMath::Abs(lid3_) == 13 && TMath::Abs(lep3_.eta()) >= 1.479){
	corr[2] = 1.01 + rndMon2011[9];
      }
      else if(TMath::Abs(lid3_) == 11 && TMath::Abs(lep3_.eta()) <  1.479){
	corr[2] = 1.01 + rndMon2011[10];
      }
      else if(TMath::Abs(lid3_) == 11 && TMath::Abs(lep3_.eta()) >= 1.479){
	corr[2] = 1.06 + rndMon2011[11];
      }
    }
    lep1pt = lep1_.pt()*corr[0];
    lep2pt = lep2_.pt()*corr[1];
    lep3pt = lep3_.pt()*corr[2];
    double pllx  = lep1_.px()*corr[0]+lep2_.px()*corr[1];
    double plly  = lep1_.py()*corr[0]+lep2_.py()*corr[1];
    double pllz  = lep1_.pz()*corr[0]+lep2_.pz()*corr[1];
    double ell   = lep1_.E()*corr[0] +lep2_.E() *corr[1];
    llPhi = TMath::ATan2(plly,pllx);
    dilmass = ell*ell -pllx*pllx -plly*plly -pllz*pllz;
    if(dilmass >=0) dilmass = sqrt(dilmass); else dilmass = 0.0;
    dilpt = sqrt(pllx*pllx+plly*plly);
    met 	   = met_;
    metPhi	   = metPhi_;
    trackMet	   = trackMet_;
    trackMetPhi    = trackMetPhi_;
    mt  = mt_*sqrt(dilpt/dilep_.pt());
    dPhiDiLepMET = DeltaPhi(llPhi,metPhi_);
    dPhiMETTrkMET = DeltaPhi(trackMetPhi_ ,metPhi_);
    pTFrac = fabs(met_-dilpt)/dilpt;
    double MT3lTotx = met_*cos(metPhi)-lep1_.px()*corr[0]-lep2_.px()*corr[1]-lep3_.px()*corr[2];
    double MT3lToty = met_*sin(metPhi)-lep1_.py()*corr[0]-lep2_.py()*corr[1]-lep3_.py()*corr[2];
    mtH = sqrt(MT3lTotx*MT3lTotx+MT3lToty*MT3lToty);
    double p2l2jx  = lep1_.px()*corr[0]+lep2_.px()*corr[1]+jet1_.px()+jet2_.px();
    double p2l2jy  = lep1_.py()*corr[0]+lep2_.py()*corr[1]+jet1_.py()+jet2_.py();
    double p2l2jz  = lep1_.pz()*corr[0]+lep2_.pz()*corr[1]+jet1_.pz()+jet2_.pz();
    double e2l2j   = lep1_.E() *corr[0]+lep2_.E() *corr[1]+jet1_.E() +jet2_.E() ;
    m2l2j = e2l2j*e2l2j -p2l2jx*p2l2jx -p2l2jy*p2l2jy -p2l2jz*p2l2jz;
    if(m2l2j >=0) m2l2j = sqrt(m2l2j); else m2l2j = 0.0;
    double p2jx  = jet1_.px()+jet2_.px();
    double p2jy  = jet1_.py()+jet2_.py();
    double p2jz  = jet1_.pz()+jet2_.pz();
    double e2j   = jet1_.E() +jet2_.E() ;
    m2j = e2j*e2j -p2jx*p2jx -p2jy*p2jy -p2jz*p2jz;
    if(m2j >=0) m2j = sqrt(m2j); else m2j = 0.0;
  }
  else if(nsel == 1){ // momentum scale -
    double corr[3] = {1.0, 1.0, 1.0};
    if(year == 2012){
      if   (TMath::Abs(lid1_) == 13 && TMath::Abs(lep1_.eta()) <  1.479){
	corr[0] = 0.99920 - rndMon2012[0];
      }
      else if(TMath::Abs(lid1_) == 13 && TMath::Abs(lep1_.eta()) >= 1.479){
	corr[0] = 0.99934 - rndMon2012[1];
      }
      else if(TMath::Abs(lid1_) == 11 && TMath::Abs(lep1_.eta()) <  1.479){
	corr[0] = 0.99807 - rndMon2012[2];
      }
      else if(TMath::Abs(lid1_) == 11 && TMath::Abs(lep1_.eta()) >= 1.479){
	corr[0] = 0.99952 - rndMon2012[3];
      }
      if   (TMath::Abs(lid2_) == 13 && TMath::Abs(lep2_.eta()) <  1.479){
	corr[1] = 0.99920 - rndMon2012[4];
      }
      else if(TMath::Abs(lid2_) == 13 && TMath::Abs(lep2_.eta()) >= 1.479){
	corr[1] = 0.99934 - rndMon2012[5];
      }
      else if(TMath::Abs(lid2_) == 11 && TMath::Abs(lep2_.eta()) <  1.479){
	corr[1] = 0.99807 - rndMon2012[6];
      }
      else if(TMath::Abs(lid2_) == 11 && TMath::Abs(lep2_.eta()) >= 1.479){
	corr[1] = 0.99952 - rndMon2012[7];
      }
    } else {
      if     (TMath::Abs(lid1_) == 13 && TMath::Abs(lep1_.eta()) <  1.479){
	corr[0] = 1/1.01 - rndMon2011[0];
      }
      else if(TMath::Abs(lid1_) == 13 && TMath::Abs(lep1_.eta()) >= 1.479){
	corr[0] = 1/1.01 - rndMon2011[1];
      }
      else if(TMath::Abs(lid1_) == 11 && TMath::Abs(lep1_.eta()) <  1.479){
	corr[0] = 1/1.01 - rndMon2011[2];
      }
      else if(TMath::Abs(lid1_) == 11 && TMath::Abs(lep1_.eta()) >= 1.479){
	corr[0] = 1/1.06 - rndMon2011[3];
      }
      if     (TMath::Abs(lid2_) == 13 && TMath::Abs(lep2_.eta()) <  1.479){
	corr[1] = 1/1.01 - rndMon2011[4];
      }
      else if(TMath::Abs(lid2_) == 13 && TMath::Abs(lep2_.eta()) >= 1.479){
	corr[1] = 1/1.01 - rndMon2011[5];
      }
      else if(TMath::Abs(lid2_) == 11 && TMath::Abs(lep2_.eta()) <  1.479){
	corr[1] = 1/1.01 - rndMon2011[6];
      }
      else if(TMath::Abs(lid2_) == 11 && TMath::Abs(lep2_.eta()) >= 1.479){
	corr[1] = 1/1.06 - rndMon2011[7];
      }
    }
    lep1pt = lep1_.pt()*corr[0];
    lep2pt = lep2_.pt()*corr[1];
    lep3pt = lep3_.pt()*corr[2];
    double pllx  = lep1_.px()*corr[0]+lep2_.px()*corr[1];
    double plly  = lep1_.py()*corr[0]+lep2_.py()*corr[1];
    double pllz  = lep1_.pz()*corr[0]+lep2_.pz()*corr[1];
    double ell   = lep1_.E() *corr[0]+lep2_.E() *corr[1];
    llPhi = TMath::ATan2(plly,pllx);
    dilmass = ell*ell -pllx*pllx -plly*plly -pllz*pllz;
    if(dilmass >=0) dilmass = sqrt(dilmass); else dilmass = 0.0;
    dilpt = sqrt(pllx*pllx+plly*plly);
    met 	   = met_;
    metPhi	   = metPhi_;
    trackMet	   = trackMet_;
    trackMetPhi    = trackMetPhi_;
    mt  = mt_*sqrt(dilpt/dilep_.pt());
    dPhiDiLepMET = DeltaPhi(llPhi,metPhi_);
    dPhiMETTrkMET = DeltaPhi(trackMetPhi_ ,metPhi_);
    pTFrac = fabs(met_-dilpt)/dilpt;
    double MT3lTotx = met_*cos(metPhi)-lep1_.px()*corr[0]-lep2_.px()*corr[1]-lep3_.px()*corr[2];
    double MT3lToty = met_*sin(metPhi)-lep1_.py()*corr[0]-lep2_.py()*corr[1]-lep3_.py()*corr[2];
    mtH = sqrt(MT3lTotx*MT3lTotx+MT3lToty*MT3lToty);
    double p2l2jx  = lep1_.px()*corr[0]+lep2_.px()*corr[1]+jet1_.px()+jet2_.px();
    double p2l2jy  = lep1_.py()*corr[0]+lep2_.py()*corr[1]+jet1_.py()+jet2_.py();
    double p2l2jz  = lep1_.pz()*corr[0]+lep2_.pz()*corr[1]+jet1_.pz()+jet2_.pz();
    double e2l2j   = lep1_.E() *corr[0]+lep2_.E() *corr[1]+jet1_.E() +jet2_.E() ;
    m2l2j = e2l2j*e2l2j -p2l2jx*p2l2jx -p2l2jy*p2l2jy -p2l2jz*p2l2jz;
    if(m2l2j >=0) m2l2j = sqrt(m2l2j); else m2l2j = 0.0;
    double p2jx  = jet1_.px()+jet2_.px();
    double p2jy  = jet1_.py()+jet2_.py();
    double p2jz  = jet1_.pz()+jet2_.pz();
    double e2j   = jet1_.E() +jet2_.E() ;
    m2j = e2j*e2j -p2jx*p2jx -p2jy*p2jy -p2jz*p2jz;
    if(m2j >=0) m2j = sqrt(m2j); else m2j = 0.0;
  }
  else if(nsel == 2){ // met
    double metx=0.0;double mety=0.0;double trkmetx=0.0;double trkmety=0.0;
    if(year == 2012){
      if    (njets_ == 0){
	metx	= met_*cos(metPhi_)	 -0.45189+gRandom->Gaus(0.0,3.2);
	mety	= met_*sin(metPhi_)	 -0.20148+gRandom->Gaus(0.0,3.2);
	trkmetx = trackMet_*cos(trackMetPhi_)+0.12580+gRandom->Gaus(0.0,1.0);
	trkmety = trackMet_*sin(trackMetPhi_)+0.02615+gRandom->Gaus(0.0,1.0);
      }
      else if(njets_ == 1){
	metx	= met_*cos(metPhi_)	 -0.39040+gRandom->Gaus(0.0,3.6);
	mety	= met_*sin(metPhi_)	 -0.20427+gRandom->Gaus(0.0,3.6);
	trkmetx = trackMet_*cos(trackMetPhi_)+0.07639+gRandom->Gaus(0.0,4.5);
	trkmety = trackMet_*sin(trackMetPhi_)+0.01167+gRandom->Gaus(0.0,4.5);
      }
      else if(njets_ >= 2){
	metx	= met_*cos(metPhi_)	 -0.27127+gRandom->Gaus(0.0,4.3);
	mety	= met_*sin(metPhi_)	 -0.18935+gRandom->Gaus(0.0,4.3);
	trkmetx = trackMet_*cos(trackMetPhi_)+0.13328+gRandom->Gaus(0.0,6.0);
	trkmety = trackMet_*sin(trackMetPhi_)-0.01351+gRandom->Gaus(0.0,6.0);
      }
    } else {
      if    (njets_ == 0){
	metx	= met_*cos(metPhi_)+gRandom->Gaus(0.0,3.2);
	mety	= met_*sin(metPhi_)+gRandom->Gaus(0.0,3.2);
	trkmetx = trackMet_*cos(trackMetPhi_)+gRandom->Gaus(0.0,2.1);
	trkmety = trackMet_*sin(trackMetPhi_)+gRandom->Gaus(0.0,2.1);
      }
      else if(njets_ == 1){
	metx	= met_*cos(metPhi_)+gRandom->Gaus(0.0,3.6);
	mety	= met_*sin(metPhi_)+gRandom->Gaus(0.0,3.6);
	trkmetx = trackMet_*cos(trackMetPhi_)+gRandom->Gaus(0.0,7.6);
	trkmety = trackMet_*sin(trackMetPhi_)+gRandom->Gaus(0.0,7.6);
      }
      else if(njets_ >= 2){
	metx	= met_*cos(metPhi_)+gRandom->Gaus(0.0,4.3);
	mety	= met_*sin(metPhi_)+gRandom->Gaus(0.0,4.3);
	trkmetx = trackMet_*cos(trackMetPhi_)+gRandom->Gaus(0.0,12.4);
	trkmety = trackMet_*sin(trackMetPhi_)+gRandom->Gaus(0.0,12.4);
      }
    }
    double newMet      = sqrt(metx*metx+mety*mety);
    double newTrackMet = sqrt(trkmetx*trkmetx+trkmety*trkmety);
    double deltaPhiA[3] = {TMath::Abs(lep1_.Phi()-TMath::ATan2(mety,metx)),TMath::Abs(lep2_.Phi()-TMath::ATan2(mety,metx)),0.0};
    while(deltaPhiA[0]>TMath::Pi()) deltaPhiA[0] = TMath::Abs(deltaPhiA[0] - 2*TMath::Pi());
    while(deltaPhiA[1]>TMath::Pi()) deltaPhiA[1] = TMath::Abs(deltaPhiA[1] - 2*TMath::Pi());
    deltaPhiA[2] = TMath::Min(deltaPhiA[0],deltaPhiA[1]);
    double pmetA = newMet;
    if(deltaPhiA[2]<TMath::Pi()/2) pmetA = pmetA * sin(deltaPhiA[2]);

    double deltaPhiB[3] = {TMath::Abs(lep1_.Phi()-TMath::ATan2(trkmety,trkmetx)),TMath::Abs(lep2_.Phi()-TMath::ATan2(trkmety,trkmetx)),0.0};
    while(deltaPhiB[0]>TMath::Pi()) deltaPhiB[0] = TMath::Abs(deltaPhiB[0] - 2*TMath::Pi());
    while(deltaPhiB[1]>TMath::Pi()) deltaPhiB[1] = TMath::Abs(deltaPhiB[1] - 2*TMath::Pi());
    deltaPhiB[2] = TMath::Min(deltaPhiB[0],deltaPhiB[1]);
    double pmetB = newTrackMet;
    if(deltaPhiB[2]<TMath::Pi()/2) pmetB = pmetB * sin(deltaPhiB[2]);

    lep1pt  = lep1_.pt();
    lep2pt  = lep2_.pt();
    lep3pt  = lep3_.pt();
    dilmass = dilep_.mass();
    dilpt   = dilep_.pt();
    met 	   = newMet;
    metPhi	   = TMath::ATan2(mety,metx);
    trackMet	   = newTrackMet;
    trackMetPhi    = TMath::ATan2(trkmety,trkmetx);
    mt  = mt_*sqrt(newMet/met_); // 8
    dPhiDiLepMET = DeltaPhi(dilep_.phi(),TMath::ATan2(mety,metx));
    dPhiMETTrkMET = DeltaPhi(TMath::ATan2(trkmety,trkmetx) ,TMath::ATan2(mety,metx));
    pTFrac = fabs(newMet-dilpt)/dilpt;
    double MT3lTotx = met_*cos(metPhi)-lep1_.px()-lep2_.px()-lep3_.px();
    double MT3lToty = met_*sin(metPhi)-lep1_.py()-lep2_.py()-lep3_.py();
    mtH = sqrt(MT3lTotx*MT3lTotx+MT3lToty*MT3lToty);
    double p2l2jx  = lep1_.px()+lep2_.px()+jet1_.px()+jet2_.px();
    double p2l2jy  = lep1_.py()+lep2_.py()+jet1_.py()+jet2_.py();
    double p2l2jz  = lep1_.pz()+lep2_.pz()+jet1_.pz()+jet2_.pz();
    double e2l2j   = lep1_.E() +lep2_.E() +jet1_.E() +jet2_.E() ;
    m2l2j = e2l2j*e2l2j -p2l2jx*p2l2jx -p2l2jy*p2l2jy -p2l2jz*p2l2jz;
    if(m2l2j >=0) m2l2j = sqrt(m2l2j); else m2l2j = 0.0;
    double p2jx  = jet1_.px()+jet2_.px();
    double p2jy  = jet1_.py()+jet2_.py();
    double p2jz  = jet1_.pz()+jet2_.pz();
    double e2j   = jet1_.E() +jet2_.E() ;
    m2j = e2j*e2j -p2jx*p2jx -p2jy*p2jy -p2jz*p2jz;
    if(m2j >=0) m2j = sqrt(m2j); else m2j = 0.0;
  }
  else if(nsel == 3){ // nothing
    lep1pt = lep1_.pt();
    lep2pt = lep2_.pt();
    lep3pt = lep3_.pt();
    double pllx  = lep1_.px()+lep2_.px();
    double plly  = lep1_.py()+lep2_.py();
    double pllz  = lep1_.pz()+lep2_.pz();
    double ell   = lep1_.E() +lep2_.E() ;
    llPhi = TMath::ATan2(plly,pllx);
    dilmass = ell*ell -pllx*pllx -plly*plly -pllz*pllz;
    if(dilmass >=0) dilmass = sqrt(dilmass); else dilmass = 0.0;
    dilpt = sqrt(pllx*pllx+plly*plly);
    met 	   = met_;
    metPhi	   = metPhi_;
    trackMet	   = trackMet_;
    trackMetPhi    = trackMetPhi_;
    mt  = mt_*sqrt(dilpt/dilep_.pt());
    dPhiDiLepMET = DeltaPhi(llPhi,metPhi_);
    dPhiMETTrkMET = DeltaPhi(trackMetPhi_ ,metPhi_);
    pTFrac = fabs(met_-dilpt)/dilpt;
    double MT3lTotx = met_*cos(metPhi)-lep1_.px()-lep2_.px()-lep3_.px();
    double MT3lToty = met_*sin(metPhi)-lep1_.py()-lep2_.py()-lep3_.py();
    mtH = sqrt(MT3lTotx*MT3lTotx+MT3lToty*MT3lToty);
    double p2l2jx  = lep1_.px()+lep2_.px()+jet1_.px()+jet2_.px();
    double p2l2jy  = lep1_.py()+lep2_.py()+jet1_.py()+jet2_.py();
    double p2l2jz  = lep1_.pz()+lep2_.pz()+jet1_.pz()+jet2_.pz();
    double e2l2j   = lep1_.E() +lep2_.E() +jet1_.E() +jet2_.E() ;
    m2l2j = e2l2j*e2l2j -p2l2jx*p2l2jx -p2l2jy*p2l2jy -p2l2jz*p2l2jz;
    if(m2l2j >=0) m2l2j = sqrt(m2l2j); else m2l2j = 0.0;
    double p2jx  = jet1_.px()+jet2_.px();
    double p2jy  = jet1_.py()+jet2_.py();
    double p2jz  = jet1_.pz()+jet2_.pz();
    double e2j   = jet1_.E() +jet2_.E() ;
    m2j = e2j*e2j -p2jx*p2jx -p2jy*p2jy -p2jz*p2jz;
    if(m2j >=0) m2j = sqrt(m2j); else m2j = 0.0;
  }
  else if(nsel == 4){ // jet+
    lep1pt = lep1_.pt();
    lep2pt = lep2_.pt();
    lep3pt = lep3_.pt();
    double pllx  = lep1_.px()+lep2_.px();
    double plly  = lep1_.py()+lep2_.py();
    double pllz  = lep1_.pz()+lep2_.pz();
    double ell   = lep1_.E() +lep2_.E() ;
    llPhi = TMath::ATan2(plly,pllx);
    dilmass = ell*ell -pllx*pllx -plly*plly -pllz*pllz;
    if(dilmass >=0) dilmass = sqrt(dilmass); else dilmass = 0.0;
    dilpt = sqrt(pllx*pllx+plly*plly);
    met 	   = met_;
    metPhi	   = metPhi_;
    trackMet	   = trackMet_;
    trackMetPhi    = trackMetPhi_;
    mt  = mt_*sqrt(dilpt/dilep_.pt());
    dPhiDiLepMET = DeltaPhi(llPhi,metPhi_);
    dPhiMETTrkMET = DeltaPhi(trackMetPhi_ ,metPhi_);
    pTFrac = fabs(met_-dilpt)/dilpt;
    double MT3lTotx = met_*cos(metPhi)-lep1_.px()-lep2_.px()-lep3_.px();
    double MT3lToty = met_*sin(metPhi)-lep1_.py()-lep2_.py()-lep3_.py();
    mtH = sqrt(MT3lTotx*MT3lTotx+MT3lToty*MT3lToty);
    double p2l2jx  = lep1_.px()+lep2_.px()+jet1_.px()*1.04+jet2_.px()*1.04;
    double p2l2jy  = lep1_.py()+lep2_.py()+jet1_.py()*1.04+jet2_.py()*1.04;
    double p2l2jz  = lep1_.pz()+lep2_.pz()+jet1_.pz()*1.04+jet2_.pz()*1.04;
    double e2l2j   = lep1_.E() +lep2_.E() +jet1_.E() *1.04+jet2_.E() *1.04;
    m2l2j = e2l2j*e2l2j -p2l2jx*p2l2jx -p2l2jy*p2l2jy -p2l2jz*p2l2jz;
    if(m2l2j >=0) m2l2j = sqrt(m2l2j); else m2l2j = 0.0;
    double p2jx  = jet1_.px()*1.04+jet2_.px()*1.04;
    double p2jy  = jet1_.py()*1.04+jet2_.py()*1.04;
    double p2jz  = jet1_.pz()*1.04+jet2_.pz()*1.04;
    double e2j   = jet1_.E() *1.04+jet2_.E() *1.04;
    m2j = e2j*e2j -p2jx*p2jx -p2jy*p2jy -p2jz*p2jz;
    if(m2j >=0) m2j = sqrt(m2j); else m2j = 0.0;
  }
  else if(nsel == 5){ // jet-
    lep1pt = lep1_.pt();
    lep2pt = lep2_.pt();
    lep3pt = lep3_.pt();
    double pllx  = lep1_.px()+lep2_.px();
    double plly  = lep1_.py()+lep2_.py();
    double pllz  = lep1_.pz()+lep2_.pz();
    double ell   = lep1_.E() +lep2_.E() ;
    llPhi = TMath::ATan2(plly,pllx);
    dilmass = ell*ell -pllx*pllx -plly*plly -pllz*pllz;
    if(dilmass >=0) dilmass = sqrt(dilmass); else dilmass = 0.0;
    dilpt = sqrt(pllx*pllx+plly*plly);
    met 	   = met_;
    metPhi	   = metPhi_;
    trackMet	   = trackMet_;
    trackMetPhi    = trackMetPhi_;
    mt  = mt_*sqrt(dilpt/dilep_.pt());
    dPhiDiLepMET = DeltaPhi(llPhi,metPhi_);
    dPhiMETTrkMET = DeltaPhi(trackMetPhi_ ,metPhi_);
    pTFrac = fabs(met_-dilpt)/dilpt;
    double MT3lTotx = met_*cos(metPhi)-lep1_.px()-lep2_.px()-lep3_.px();
    double MT3lToty = met_*sin(metPhi)-lep1_.py()-lep2_.py()-lep3_.py();
    mtH = sqrt(MT3lTotx*MT3lTotx+MT3lToty*MT3lToty);
    double p2l2jx  = lep1_.px()+lep2_.px()+jet1_.px()*0.96+jet2_.px()*0.96;
    double p2l2jy  = lep1_.py()+lep2_.py()+jet1_.py()*0.96+jet2_.py()*0.96;
    double p2l2jz  = lep1_.pz()+lep2_.pz()+jet1_.pz()*0.96+jet2_.pz()*0.96;
    double e2l2j   = lep1_.E() +lep2_.E() +jet1_.E() *0.96+jet2_.E() *0.96;
    m2l2j = e2l2j*e2l2j -p2l2jx*p2l2jx -p2l2jy*p2l2jy -p2l2jz*p2l2jz;
    if(m2l2j >=0) m2l2j = sqrt(m2l2j); else m2l2j = 0.0;
    double p2jx  = jet1_.px()*0.96+jet2_.px()*0.96;
    double p2jy  = jet1_.py()*0.96+jet2_.py()*0.96;
    double p2jz  = jet1_.pz()*0.96+jet2_.pz()*0.96;
    double e2j   = jet1_.E() *0.96+jet2_.E() *0.96;
    m2j = e2j*e2j -p2jx*p2jx -p2jy*p2jy -p2jz*p2jz;
    if(m2j >=0) m2j = sqrt(m2j); else m2j = 0.0;
  }
  else { assert(0);}

  outputVar[ 0] = lep1pt;
  outputVar[ 1] = lep2pt;
  outputVar[ 2] = dilmass;
  outputVar[ 3] = dilpt;
  outputVar[ 4] = met;
  outputVar[ 5] = metPhi;
  outputVar[ 6] = trackMet;
  outputVar[ 7] = trackMetPhi;
  outputVar[ 8] = mt;
  outputVar[ 9] = dPhiDiLepMET;
  outputVar[10] = dPhiMETTrkMET;
  outputVar[11] = pTFrac;
  outputVar[12] = mtH;
  outputVar[13] = m2l2j;
  outputVar[14] = m2j;
  outputVar[15] = lep3pt;
}
