#include "../include/MuonLoader.hh"
#include "TMath.h"
#include <iostream> 

using namespace baconhep;

MuonLoader::MuonLoader(TTree *iTree) { 
  fMuons  = new TClonesArray("baconhep::TMuon");
  iTree->SetBranchAddress("Muon",       &fMuons);
  fMuonBr  = iTree->GetBranch("Muon");
  fDiMuon  = new TLorentzVector(0.,0.,0.,0.);
  fMassMin = 115;
  fMassMax = 130;
}
MuonLoader::~MuonLoader() { 
  delete fMuons;
  delete fMuonBr;
  delete fDiMuon;
}
void MuonLoader::resetDiMu() { 
  fDiMuon->SetPtEtaPhiM(1e-9,0,0,0);
}
void MuonLoader::reset() { 
  fSelMuons.clear();
}
void MuonLoader::setupTree(TTree *iTree) { 
  reset();
  fTree = iTree;
  fTree->Branch("nmuons" ,&fNMuons     ,"fNMuons/I");
  fTree->Branch("dimu"   ,"TLorentzVector", &fDiMuon);
}
void MuonLoader::load(int iEvent) { 
  fMuons   ->Clear();
  fMuonBr ->GetEntry(iEvent);
}
void MuonLoader::fillVetoes(std::vector<TLorentzVector> &iVec) { 
  for(unsigned int i0 = 0; i0 < fSelMuons.size(); i0++) { 
    TLorentzVector pVec;
    pVec.SetPtEtaPhiM(fSelMuons[i0]->pt,fSelMuons[i0]->eta,fSelMuons[i0]->phi,0.105);
    iVec.push_back(pVec);
  }
}
bool MuonLoader::selectMuons(std::vector<TLorentzVector>& iVetoes) {
  reset(); 
  int  lNCount = 0; 
  for  (int i0 = 0; i0 < fMuons->GetEntriesFast(); i0++) { 
    TMuon *pMuon = (TMuon*)((*fMuons)[i0]);
    if(pMuon->pt > 10 && fabs(pMuon->eta) < 2.4) lNCount++;
    if(!passWW(pMuon)) continue;
    bool pMatch = false;
    for(unsigned int i1 = 0; i1 < iVetoes.size(); i1++) { 
      double pDEta = pMuon->eta      - iVetoes[i1].Eta();
      double pDPhi = fabs(pMuon->phi - iVetoes[i1].Phi());
      if(fabs(pDPhi) > 2.*TMath::Pi()-fabs(pDPhi)) pDPhi =  2.*TMath::Pi()-fabs(pDPhi);
      if(sqrt(pDPhi*pDPhi+pDEta*pDEta) > 0.3) continue;
      pMatch = true;
    }
    if(pMatch) continue;
    bool lFill = false;
    for( std::vector<TMuon*>::iterator pMuonIter = fSelMuons.begin(); pMuonIter != fSelMuons.end(); pMuonIter++) { 
      if((*pMuonIter)->pt > pMuon->pt) continue;
      fSelMuons.insert(pMuonIter,pMuon);
      lFill = true;
      break;
    } 
    if(!lFill)  fSelMuons.push_back(pMuon);
  }
  fNMuons = lNCount;
  if(fSelMuons.size() == 0) return false;
  return true;
}
void MuonLoader::selectDiMuon(std::vector<TLorentzVector> &iVetoes) {
  resetDiMu();
  for  (int i0 = 0; i0 < fMuons->GetEntriesFast(); i0++) {
    TMuon *pMuon1 = (TMuon*)((*fMuons)[i0]);
    if(!passWW(pMuon1)) continue;
    for(int i1 = 0; i1 < fMuons->GetEntriesFast(); i1++) {
      if(i0 == i1) continue;
      TMuon *pMuon2 = (TMuon*)((*fMuons)[i1]);
      if(!passWW(pMuon2)) continue;
      TLorentzVector pVec1;      pVec1.SetPtEtaPhiM(pMuon1->pt,pMuon1->eta,pMuon1->phi,0.105);
      TLorentzVector pVec2;      pVec2.SetPtEtaPhiM(pMuon2->pt,pMuon2->eta,pMuon2->phi,0.105);      
      if(fMassMin > (pVec1+pVec2).M() || fMassMax < (pVec1+pVec2).M()) continue;
      iVetoes.push_back(pVec1);
      iVetoes.push_back(pVec2);
      break;
    }
    if(iVetoes.size() > 1) break;
  }
  TLorentzVector lDi;
  if(iVetoes.size() > 1) lDi =  (iVetoes[0] + iVetoes[1]);
  if(iVetoes.size() > 1) fDiMuon->SetPtEtaPhiM( lDi.Pt(),lDi.Eta(),lDi.Phi(),lDi.M());
}
bool MuonLoader::vetoMu() {
  for  (int i0 = 0; i0 < fMuons->GetEntriesFast(); i0++) { 
    TMuon *pMuon = (TMuon*)((*fMuons)[i0]);
    if(passLoose(pMuon)) return true;
  }
  return false;
}
//H=>ZZ Mu Id
bool MuonLoader::passLoose(TMuon *iMuon) { 
  if(!(iMuon->typeBits     & kGlobal || iMuon->typeBits & kTracker))  return false;
  if(!(iMuon->selectorBits & kAllArbitrated))                        return false;
  if(!(iMuon->typeBits     & kPFMuon))                               return false;
  if(fabs(iMuon->dz)> 1.0)                                           return false;
  double chargedIso = iMuon->chHadIso04;
  double neutralIso = TMath::Max(iMuon->gammaIso04 + iMuon->neuHadIso04 - 0.5 * iMuon->puIso04, 0.0);
  double totalIso   = chargedIso+neutralIso;
  if(totalIso/iMuon->pt > 0.4) return false;
  return true;
}
bool MuonLoader::passWW(TMuon *iMuon) { 
  //if(iMuon->pt > 10) std::cout << "===> " << iMuon->pt << " -- " << iMuon->eta << " -- " << iMuon->nTkLayers << " -- " << iMuon->nPixHits << " -- " << iMuon->trkKink << " -- " << iMuon->dz << " -- " << iMuon->ptErr/iMuon->pt << " -- " << (iMuon->typeBits & kPFMuon) << std::endl;
  if(!(fabs(iMuon->eta)     < 2.4     ))                               return false;       
  if(!(iMuon->pt            > 10      ))                               return false;       
  if(!(iMuon->typeBits     & kPFMuon  ))                               return false;
  if(!(iMuon->nTkLayers      > 0      ))                               return false;
  if(!(iMuon->nPixHits       > 0      ))                               return false;
  if(!(iMuon->ptErr/iMuon->pt < 0.1   ))                               return false;
  if(!(iMuon->trkKink        < 20.    ))                               return false;
  if(fabs(iMuon->dz)> 1.0)                                             return false;
  //replace with very loose delta beta iso
  double chargedIso = iMuon->chHadIso04;
  double neutralIso = TMath::Max(iMuon->gammaIso04 + iMuon->neuHadIso04 - 0.5 * iMuon->puIso04, 0.0);
  double totalIso   = chargedIso+neutralIso;
  if(totalIso/iMuon->pt > 0.4) return false;
  return true;
}
bool MuonLoader::passTight(TMuon *iMuon) { 
  if(!(iMuon->typeBits & kGlobal))  return false;
  if(fabs(iMuon->dz)  > 0.2)        return false;
  if(fabs(iMuon->d0)  > 0.045)      return false;
  if(iMuon->muNchi2        > 10)    return false;
  if(iMuon->nValidHits     < 1)     return false;
  if(iMuon->nMatchStn      < 2)     return false;
  if(iMuon->nPixHits       < 1)     return false;
  if(iMuon->nTkLayers      < 6)     return false;
  if(!(iMuon->typeBits & kPFMuon))  return false;


  double chargedIso = iMuon->chHadIso04;
  double neutralIso = TMath::Max(iMuon->gammaIso04 + iMuon->neuHadIso04 - 0.5 * iMuon->puIso04, 0.0);
  double totalIso   = chargedIso+neutralIso;
  if(totalIso/iMuon->pt > 0.15) return false;
  return true;
}
