	  int type = -1;
	  if     ((TMath::Abs(bgdEvent.lid1_) == 11 && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && bgdEvent.lep2_.Pt() < 20 && TMath::Abs(bgdEvent.lep2_.Eta())< 1.479 && ((bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4) && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) ||
	          (TMath::Abs(bgdEvent.lid2_) == 11 && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection && bgdEvent.lep1_.Pt() < 20 && TMath::Abs(bgdEvent.lep1_.Eta())< 1.479 && ((bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4) && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection)){
	    type = 0;
	  }
	  else if((TMath::Abs(bgdEvent.lid1_) == 11 && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && bgdEvent.lep2_.Pt() >=20 && TMath::Abs(bgdEvent.lep2_.Eta())< 1.479 && ((bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4) && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) ||
	          (TMath::Abs(bgdEvent.lid2_) == 11 && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection && bgdEvent.lep1_.Pt() >=20 && TMath::Abs(bgdEvent.lep1_.Eta())< 1.479 && ((bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4) && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection)){
	    type = 1;
	  }
	  else if((TMath::Abs(bgdEvent.lid1_) == 11 && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && bgdEvent.lep2_.Pt() < 20 && TMath::Abs(bgdEvent.lep2_.Eta())>=1.479 && ((bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4) && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) ||
	          (TMath::Abs(bgdEvent.lid2_) == 11 && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection && bgdEvent.lep1_.Pt() < 20 && TMath::Abs(bgdEvent.lep1_.Eta())>=1.479 && ((bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4) && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection)){
	    type = 2;
	  }
	  else if((TMath::Abs(bgdEvent.lid1_) == 11 && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && bgdEvent.lep2_.Pt() >=20 && TMath::Abs(bgdEvent.lep2_.Eta())>=1.479 && ((bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4) && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) ||
	          (TMath::Abs(bgdEvent.lid2_) == 11 && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection && bgdEvent.lep1_.Pt() >=20 && TMath::Abs(bgdEvent.lep1_.Eta())>=1.479 && ((bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4) && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection)){
	    type = 3;
	  }
	  else if((TMath::Abs(bgdEvent.lid1_) == 13 && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && bgdEvent.lep2_.Pt() < 20 && TMath::Abs(bgdEvent.lep2_.Eta())< 1.479 && ((bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4) && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) ||
	          (TMath::Abs(bgdEvent.lid2_) == 13 && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection && bgdEvent.lep1_.Pt() < 20 && TMath::Abs(bgdEvent.lep1_.Eta())< 1.479 && ((bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4) && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection)){
	    type = 4;
	  }
	  else if((TMath::Abs(bgdEvent.lid1_) == 13 && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && bgdEvent.lep2_.Pt() >=20 && TMath::Abs(bgdEvent.lep2_.Eta())< 1.479 && ((bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4) && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) ||
	          (TMath::Abs(bgdEvent.lid2_) == 13 && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection && bgdEvent.lep1_.Pt() >=20 && TMath::Abs(bgdEvent.lep1_.Eta())< 1.479 && ((bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4) && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection)){
	    type = 5;
	  }
	  else if((TMath::Abs(bgdEvent.lid1_) == 13 && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && bgdEvent.lep2_.Pt() < 20 && TMath::Abs(bgdEvent.lep2_.Eta())>=1.479 && ((bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4) && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) ||
	          (TMath::Abs(bgdEvent.lid2_) == 13 && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection && bgdEvent.lep1_.Pt() < 20 && TMath::Abs(bgdEvent.lep1_.Eta())>=1.479 && ((bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4) && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection)){
	    type = 6;
	  }
	  else if((TMath::Abs(bgdEvent.lid1_) == 13 && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && bgdEvent.lep2_.Pt() >=20 && TMath::Abs(bgdEvent.lep2_.Eta())>=1.479 && ((bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4) && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) ||
	          (TMath::Abs(bgdEvent.lid2_) == 13 && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection && bgdEvent.lep1_.Pt() >=20 && TMath::Abs(bgdEvent.lep1_.Eta())>=1.479 && ((bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4) && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection)){
	    type = 7;
	  }
	  else if((TMath::Abs(bgdEvent.lid1_) == 11 && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && bgdEvent.lep2_.Pt() < 20 && TMath::Abs(bgdEvent.lep2_.Eta())< 1.479 && ((bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV2) == SmurfTree::Lep2LooseMuV2) && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) ||
	          (TMath::Abs(bgdEvent.lid2_) == 11 && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection && bgdEvent.lep1_.Pt() < 20 && TMath::Abs(bgdEvent.lep1_.Eta())< 1.479 && ((bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV2) == SmurfTree::Lep1LooseMuV2) && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection)){
	    type = 8;
	  }
	  else if((TMath::Abs(bgdEvent.lid1_) == 11 && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && bgdEvent.lep2_.Pt() >=20 && TMath::Abs(bgdEvent.lep2_.Eta())< 1.479 && ((bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV2) == SmurfTree::Lep2LooseMuV2) && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) ||
	          (TMath::Abs(bgdEvent.lid2_) == 11 && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection && bgdEvent.lep1_.Pt() >=20 && TMath::Abs(bgdEvent.lep1_.Eta())< 1.479 && ((bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV2) == SmurfTree::Lep1LooseMuV2) && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection)){
	    type = 9;
	  }
	  else if((TMath::Abs(bgdEvent.lid1_) == 11 && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && bgdEvent.lep2_.Pt() < 20 && TMath::Abs(bgdEvent.lep2_.Eta())>=1.479 && ((bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV2) == SmurfTree::Lep2LooseMuV2) && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) ||
	          (TMath::Abs(bgdEvent.lid2_) == 11 && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection && bgdEvent.lep1_.Pt() < 20 && TMath::Abs(bgdEvent.lep1_.Eta())>=1.479 && ((bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV2) == SmurfTree::Lep1LooseMuV2) && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection)){
	    type = 10;
	  }
	  else if((TMath::Abs(bgdEvent.lid1_) == 11 && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && bgdEvent.lep2_.Pt() >=20 && TMath::Abs(bgdEvent.lep2_.Eta())>=1.479 && ((bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV2) == SmurfTree::Lep2LooseMuV2) && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) ||
	          (TMath::Abs(bgdEvent.lid2_) == 11 && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection && bgdEvent.lep1_.Pt() >=20 && TMath::Abs(bgdEvent.lep1_.Eta())>=1.479 && ((bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV2) == SmurfTree::Lep1LooseMuV2) && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection)){
	    type = 11;
	  }
	  else if((TMath::Abs(bgdEvent.lid1_) == 13 && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && bgdEvent.lep2_.Pt() < 20 && TMath::Abs(bgdEvent.lep2_.Eta())< 1.479 && ((bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV2) == SmurfTree::Lep2LooseMuV2) && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) ||
	          (TMath::Abs(bgdEvent.lid2_) == 13 && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection && bgdEvent.lep1_.Pt() < 20 && TMath::Abs(bgdEvent.lep1_.Eta())< 1.479 && ((bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV2) == SmurfTree::Lep1LooseMuV2) && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection)){
	    type = 12;
	  }
	  else if((TMath::Abs(bgdEvent.lid1_) == 13 && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && bgdEvent.lep2_.Pt() >=20 && TMath::Abs(bgdEvent.lep2_.Eta())< 1.479 && ((bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV2) == SmurfTree::Lep2LooseMuV2) && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) ||
	          (TMath::Abs(bgdEvent.lid2_) == 13 && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection && bgdEvent.lep1_.Pt() >=20 && TMath::Abs(bgdEvent.lep1_.Eta())< 1.479 && ((bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV2) == SmurfTree::Lep1LooseMuV2) && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection)){
	    type = 13;
	  }
	  else if((TMath::Abs(bgdEvent.lid1_) == 13 && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && bgdEvent.lep2_.Pt() < 20 && TMath::Abs(bgdEvent.lep2_.Eta())>=1.479 && ((bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV2) == SmurfTree::Lep2LooseMuV2) && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) ||
	          (TMath::Abs(bgdEvent.lid2_) == 13 && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection && bgdEvent.lep1_.Pt() < 20 && TMath::Abs(bgdEvent.lep1_.Eta())>=1.479 && ((bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV2) == SmurfTree::Lep1LooseMuV2) && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection)){
	    type = 14;
	  }
	  else if((TMath::Abs(bgdEvent.lid1_) == 13 && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && bgdEvent.lep2_.Pt() >=20 && TMath::Abs(bgdEvent.lep2_.Eta())>=1.479 && ((bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV2) == SmurfTree::Lep2LooseMuV2) && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) ||
	          (TMath::Abs(bgdEvent.lid2_) == 13 && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection && bgdEvent.lep1_.Pt() >=20 && TMath::Abs(bgdEvent.lep1_.Eta())>=1.479 && ((bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV2) == SmurfTree::Lep1LooseMuV2) && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection)){
	    type = 15;
	  }
	  else {
	    cout << "NOOOOOOOOOOOOOO" << endl;
	    cout << ((bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV2) == SmurfTree::Lep1LooseMuV2) << " "
	         << ((bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV2) == SmurfTree::Lep2LooseMuV2) << " "
		 << ((bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4) << " "
		 << ((bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4) << endl;
	    assert(0);
	  }
	  bgdDecayFake[type]  += theWeight;
	  bgdDecayFakeE[type] += theWeight*theWeight;
