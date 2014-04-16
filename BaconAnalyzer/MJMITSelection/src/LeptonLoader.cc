#include "../include/LeptonLoader.hh"
#include "TMath.h"

using namespace baconhep;

LeptonLoader::LeptonLoader(TTree *iTree) {
  fPtr1 = new TLorentzVector(0.,0.,0.,0.);
  fPtr2 = new TLorentzVector(0.,0.,0.,0.);
  fPtr3 = new TLorentzVector(0.,0.,0.,0.);
}
LeptonLoader::~LeptonLoader() { 
  delete fPtr1;
  delete fPtr2;
  delete fPtr3;
}
void LeptonLoader::reset() { 
  fNLep = 0; 

  fPtr1->SetPtEtaPhiM(0,0,0,0);
  fPtr2->SetPtEtaPhiM(0,0,0,0);
  fPtr3->SetPtEtaPhiM(0,0,0,0);
  
  fId1 = 0; 
  fId2 = 0; 
  fId3 = 0; 

  fLep1IsTightMuon = 0; 
  fLep2IsTightMuon = 0; 
  fLep3IsTightMuon = 0; 
}
void LeptonLoader::setupTree(TTree *iTree) { 
  reset();
  fTree = iTree;
  fTree->Branch("nlep" , &fNLep, "fNLep/i");

  fTree->Branch("lep1" , "TLorentzVector", &fPtr1);
  fTree->Branch("lep2" , "TLorentzVector", &fPtr2);
  fTree->Branch("lep3" , "TLorentzVector", &fPtr3);

  fTree->Branch("lep1IsTightMuon" , &fLep1IsTightMuon , "fLep1IsTightMuon/I");
  fTree->Branch("lep2IsTightMuon" , &fLep2IsTightMuon , "fLep2IsTightMuon/I");
  fTree->Branch("lep3IsTightMuon" , &fLep3IsTightMuon , "fLep3IsTightMuon/I");

  fTree->Branch("lid1" , &fId1 , "fId1/I");
  fTree->Branch("lid2" , &fId2 , "fId2/I");
  fTree->Branch("lid3" , &fId3 , "fId3/I");
}
void LeptonLoader::fillLeptons(std::vector<TMuon*> iMuons,std::vector<TElectron*> iElectrons) {
  reset(); 
  fNLep = (iMuons.size()+iElectrons.size());
  std::vector<TLorentzVector> lVec; 
  std::vector<int>            lId; 
  for(unsigned int i0 = 0; i0 < iMuons.size(); i0++) { 
    TLorentzVector pVec;
    pVec.SetPtEtaPhiM(iMuons[i0]->pt,iMuons[i0]->eta,iMuons[i0]->phi,0.105); 
    int pVId = 13*iMuons[i0]->q; if(passTight(iMuons[i0])) pVId+=20*iMuons[i0]->q;
    lVec.push_back(pVec);
    lId.push_back(pVId);
  }
  for(unsigned int i0 = 0; i0 < iElectrons.size(); i0++) { 
    TLorentzVector pVec;
    pVec.SetPtEtaPhiM(iElectrons[i0]->pt,iElectrons[i0]->eta,iElectrons[i0]->phi,0.000511); 
    int pVId = 11*iElectrons[i0]->q; 
    bool lFill = false;
    std::vector<int>::iterator pIdIter = lId.begin();
    for(std::vector<TLorentzVector>::iterator pMuIter = lVec.begin(); pMuIter != lVec.end(); pMuIter++) { 
      if(pMuIter->Pt() > pVec.Pt()) {pIdIter++; continue;}
      lVec.insert(pMuIter,pVec);
      lId.insert(pIdIter,pVId);
      lFill = true;
      break;
    }
    if(!lFill) lVec.push_back(pVec);
    if(!lFill) lId .push_back(pVId);
  }
  if(lVec.size() > 0) fPtr1->SetPtEtaPhiM(lVec[0].Pt(),lVec[0].Eta(),lVec[0].Phi(),lVec[0].M());
  if(lVec.size() > 1) fPtr2->SetPtEtaPhiM(lVec[1].Pt(),lVec[1].Eta(),lVec[1].Phi(),lVec[1].M());
  if(lVec.size() > 2) fPtr3->SetPtEtaPhiM(lVec[2].Pt(),lVec[2].Eta(),lVec[2].Phi(),lVec[2].M());
  if(lId .size() > 0) fId1 = lId[0] % 20;
  if(lId .size() > 1) fId2 = lId[1] % 20;
  if(lId .size() > 2) fId3 = lId[2] % 20;
  if(lId .size() > 0) fLep1IsTightMuon = (lId[0] > 20); 
  if(lId .size() > 1) fLep2IsTightMuon = (lId[0] > 20); 
  if(lId .size() > 2) fLep3IsTightMuon = (lId[0] > 20); 
}
bool LeptonLoader::passTight(TMuon *iMuon) { 
  if(!(iMuon->typeBits & kGlobal))  return false;
  if(!(iMuon->typeBits & kTracker)) return false;
  if(iMuon->muNchi2        > 10)    return false;
  if(iMuon->nValidHits     < 1)     return false;
  if(iMuon->nMatchStn      < 2)     return false;
  if(iMuon->nPixHits       < 1)     return false;
  if(!(iMuon->typeBits & kPFMuon))  return false;
  return true;
}
