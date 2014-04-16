#include "../include/TauLoader.hh"
#include <iostream> 

using namespace baconhep;

TauLoader::TauLoader(TTree *iTree) { 
  fTaus  = new TClonesArray("baconhep::TTau");
  iTree->SetBranchAddress("Tau",       &fTaus);
  fTauBr  = iTree->GetBranch("Tau");
}
TauLoader::~TauLoader() { 
  delete fTaus;
  delete fTauBr;
}
void TauLoader::reset() { 
  fNTaus = 0; 
  fPtr1  = 0; 
  fPtr2  = 0; 

  fPt1   = 0; 
  fEta1  = 0; 
  fPhi1  = 0; 
  fM1    = 0; 
  
  fPt2   = 0; 
  fEta2  = 0; 
  fPhi2  = 0; 
  fM2    = 0; 
}
void TauLoader::setupTree(TTree *iTree) { 
  reset();
  fTree = iTree;
  fTree->Branch("ntaus",&fNTaus,"fNTaus/I");
  //fTree->Branch("tau1" , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &fPtr1);
  //fTree->Branch("tau2" , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &fPtr2);
  fTree->Branch("tau1" , "TLorentzVector", &fPtr1);
  fTree->Branch("tau2" , "TLorentzVector", &fPtr2);
  /*
  fTree->Branch("taupt_1"  ,&fPt1  ,"fPt1/F");
  fTree->Branch("taueta_1" ,&fEta1 ,"fEta1/F");
  fTree->Branch("tauphi_1" ,&fPhi1 ,"fPhi1/F");
  fTree->Branch("taum_1"   ,&fM1   ,"fM1/F");
  
  fTree->Branch("taupt_2"  ,&fPt2  ,"fPt2/F");
  fTree->Branch("taueta_2" ,&fEta2 ,"fEta2/F");
  fTree->Branch("tauphi_2" ,&fPhi2 ,"fPhi2/F");
  fTree->Branch("taum_2"   ,&fM2   ,"fM2/F");
  */
}
void TauLoader::load(int iEvent) { 
  fTaus   ->Clear();
  fTauBr ->GetEntry(iEvent);
}
void TauLoader::fillVetoes(std::vector<TLorentzVector> &iVec) { 
  TLorentzVector pVec1,pVec2;
  if(fPt1 > 0) pVec1.SetPtEtaPhiM(fPt1,fEta1,fPhi1,fM1);
  if(fPt2 > 0) pVec2.SetPtEtaPhiM(fPt2,fEta2,fPhi2,fM2);
  if(fPt1 > 0) iVec.push_back(pVec1);
  if(fPt2 > 0) iVec.push_back(pVec2);
}
bool TauLoader::selectTaus(std::vector<TLorentzVector> iVetoes) {
  reset(); 
  int lCount = 0; 
  TTau *lTau1 = 0; 
  TTau *lTau2 = 0; 
  for  (int i0 = 0; i0 < fTaus->GetEntriesFast(); i0++) { 
    TTau *pTau = (TTau*)((*fTaus)[i0]);
    if(pTau->pt < 15.)            continue;
    if(!passLoose(pTau))          continue;
    bool pMatch = false;
    for(unsigned int i1 = 0; i1 < iVetoes.size(); i1++) { 
      double pDEta = pTau->eta      - iVetoes[i1].Eta();
      double pDPhi = fabs(pTau->phi - iVetoes[i1].Phi());
      if(fabs(pDPhi) > 2.*TMath::Pi()-fabs(pDPhi)) pDPhi =  2.*TMath::Pi()-fabs(pDPhi);
      if(sqrt(pDPhi*pDPhi+pDEta*pDEta) > 0.5) continue;
      pMatch = true;
    }
    if(pMatch) continue;
    if(lTau1 == 0) lTau1 = pTau;
    if(pTau->pt > lTau1->pt) {lTau2 = pTau; lTau1 = pTau; continue;}
    if(lTau2 == 0)            lTau2 = pTau;
    if(pTau->pt > lTau2->pt) {lTau2 = pTau; }
    lCount++;
  }
  fNTaus = lCount;
  if(lTau1 == 0) return false;
  fPt1  = lTau1->pt;
  fEta1 = lTau1->eta;
  fPhi1 = lTau1->phi;
  fM1   = lTau1->m;
  delete fPtr1;
  fPtr1 = new TLorentzVector();
  fPtr1->SetPtEtaPhiM(fPt1,fEta1,fPhi1,fM1);
  //TLorentzVector lVec1; lVec1.SetPtEtaPhiM(fPt1,fEta1,fPhi1,fM1);
  //fPtr1.SetPxPyPzE(lVec.Px(),lVec.Py(),lVec.Pz(),lVec.E());

  fPt2  = lTau2->pt;
  fEta2 = lTau2->eta;
  fPhi2 = lTau2->phi;
  fM2   = lTau2->m;
  delete fPtr2;
  fPtr2 = new TLorentzVector();
  fPtr2->SetPtEtaPhiM(fPt2,fEta2,fPhi2,fM2);
  //TLorentzVector lVec2; lVec2.SetPtEtaPhiM(fPt2,fEta2,fPhi2,fM2);
  //fPtr2.SetPxPyPzE(lVec2.Px(),lVec2.Py(),lVec2.Pz(),lVec2.E());

  return true;
}
bool TauLoader::vetoTau() {
  for  (int i0 = 0; i0 < fTaus->GetEntriesFast(); i0++) { 
    TTau *pTau = (TTau*)((*fTaus)[i0]);
    if(passVeto(pTau)) return true;
  }
  return false;
}
//Loose Tau Id
bool TauLoader::passLoose(TTau *iTau) { 
  if(!(iTau->hpsDisc & kByMVA3LooseElectronRejection &&
       iTau->hpsDisc & kByLooseMuonRejection  &&
       iTau->hpsDisc & kByDecayModeFinding))   return false;  
  //if(iTau->hpsDisc & kByMVA3LooseElectronRejection) std::cout << "===>1 " << (iTau->hpsDisc & kByMVA3LooseElectronRejection)  << std::endl;
  //if(iTau->hpsDisc & kByLooseMuonRejection        ) std::cout << "===>2 " << (iTau->hpsDisc & kByLooseMuonRejection)          << std::endl;
  //if(iTau->hpsDisc & kByDecayModeFinding          ) std::cout << "===>3 " << (iTau->hpsDisc & kByDecayModeFinding)            << std::endl;
  if(iTau->rawIso3Hits  > 3.0) return false;
  return true;
}
bool TauLoader::passTight(TTau *iTau) { 
  if(!(iTau->hpsDisc & kByMVA3MediumElectronRejection &&
       iTau->hpsDisc & kByLooseMuonRejection  &&
       iTau->hpsDisc & kByDecayModeFinding))   return false;  
  //if(!passAntiEMVA3(iTau->antiEleMVA3Cat,iTau-> antiEleMVA3,"Loose")) return false;
  if(iTau->rawIso3Hits  > 1.5) return false;
  return true;
}
bool TauLoader::passVeto(TTau *iTau) { 
  if(!(//iTau->hpsDiscriminators & TTau::kByMVA3LooseElectronRejection &&
       //iTau->hpsDiscriminators & TTau::kByLooseMuonRejection  && 
       iTau->hpsDisc & kByDecayModeFinding))   return false;  
  if(iTau->rawIso3Hits  > 5.0) return false;
  return true;
}
bool TauLoader::passAntiEMVA3(int iCat, float raw, TString WP) {
  if(iCat<0) return false;
  if(iCat>15) return true;
  float cutsLoose[16]={0.835,0.831,0.849,0.859,0.873,0.823,0.85,0.855,0.816,0.861,0.862,0.847,0.893,0.82,0.845,0.851};
  float cutsMedium[16]={0.933,0.921,0.944,0.945,0.918,0.941,0.981,0.943,0.956,0.947,0.951,0.95,0.897,0.958,0.955,0.942};
  float cutsTight[16]={ 0.96,0.968,0.971,0.972,0.969,0.959,0.981,0.965,0.975,0.972,0.974,0.971,0.897,0.971,0.961,0.97};
  float cutsVeryTight[16]={0.978,0.98,0.982,0.985,0.977,0.974,0.989,0.977,0.986,0.983,0.984,0.983,0.971,0.987,0.977,0.981};
  float cut=0;
  if(WP=="Loose")  cut = cutsLoose[iCat];
  if(WP=="Medium")   cut = cutsMedium[iCat];
  if(WP=="Tight") cut = cutsTight[iCat];
  if(WP=="VeryTight") cut = cutsVeryTight[iCat];
  return (raw>cut);
}
