#include "../include/PhotonLoader.hh"
#include "TMath.h"

using namespace baconhep;

PhotonLoader::PhotonLoader(TTree *iTree) { 
  fPhotons  = new TClonesArray("baconhep::TPhoton");
  iTree->SetBranchAddress("Photon",       &fPhotons);
  fPhotonBr  = iTree->GetBranch("Photon");
}
PhotonLoader::~PhotonLoader() { 
  delete fPhotons;
  delete fPhotonBr;
}
void PhotonLoader::reset() { 
  fNPhotons = 0; 
  delete fPtr1;
  fPtr1     = 0; 
}
void PhotonLoader::setupTree(TTree *iTree) { 
  reset();
  fTree = iTree;
  fTree->Branch("nphotons" ,&fNPhotons,"fNPhotons/I");
  fTree->Branch("pho1"     ,"TLorentzVector", &fPtr1);
}
void PhotonLoader::load(int iEvent) { 
  fPhotons   ->Clear();
  fPhotonBr ->GetEntry(iEvent);
}
void PhotonLoader::fillVetoes(std::vector<TLorentzVector> &iVec) { 
  TLorentzVector pVec;
  if(fPtr1 != 0 && fPtr1->Pt() > 0) pVec.SetPtEtaPhiM(fPtr1->Pt(),fPtr1->Eta(),fPtr1->Phi(),0);
  iVec.push_back(pVec);
}
bool PhotonLoader::selectPhotons(std::vector<TLorentzVector> iElectrons,float iRho) {
  reset(); 
  int lCount = 0; 
  TPhoton *lPhoton = 0; 
  for  (int i0 = 0; i0 < fPhotons->GetEntriesFast(); i0++) { 
    TPhoton *pPhoton = (TPhoton*)((*fPhotons)[i0]);
    bool pPass = true;
    for(unsigned int i1 = 0; i1 < iElectrons.size(); i1++) { 
      if(fabs(iElectrons[i1].Eta() - pPhoton->scEta) < 0.01) {pPass = false; break;}
    }
    if(!pPass) continue;
    if(!passCiCPFIso(pPhoton,iRho)) continue;
    if(lPhoton !=0 && lPhoton->pt > pPhoton->pt) continue;
    lPhoton = pPhoton;
    lCount++;
  }
  fNPhotons = lCount;
  if(lPhoton == 0) return false;
  delete fPtr1;
  fPtr1 = new TLorentzVector();
  fPtr1->SetPtEtaPhiM(lPhoton->pt,lPhoton->eta,lPhoton->phi,0);
  return true;
}
bool PhotonLoader::vetoPhoton() {
  for  (int i0 = 0; i0 < fPhotons->GetEntriesFast(); i0++) { 
    TPhoton *pPhoton = (TPhoton*)((*fPhotons)[i0]);
    if(passLoose(pPhoton)) return true;
  }
  return false;
}
//Basic POG Id except for Iso, which is just a random guess
bool PhotonLoader::passLoose(TPhoton *photon) { 
  if(photon->pt     < 15)   return false;
  if(photon->hovere > 0.05) return false;
  if(photon->sieie  > 0.01) return false;
  double chargedIso = photon->chHadIso03;
  double neutralIso = photon->gammaIso03 + photon->neuHadIso03;
  double totalIso   = chargedIso+neutralIso;
  if(totalIso/photon->pt > 0.4) return false;
  return true;
}
bool PhotonLoader::passCiCPhotonPre(TPhoton *p,float iRho) { 
  bool lIsBarrel = true; 
  if(fabs(p->scEta)<1.566 && fabs(p->scEta)>1.4443) return false;
  if(fabs(p->scEta)>2.5)                            return false;
  if(fabs(p->scEta)>1.566 && fabs(p->scEta)<2.5) lIsBarrel = false;
  if(p->r9 < 0.9) {
    if(p->hovere > 0.075)              return false;
    if(lIsBarrel  && p->sieie > 0.014) return false; 
    if(!lIsBarrel && p->sieie > 0.034) return false; 
    //Note these iso were not PF 
    if(p->gammaIso03 -0.012*p->pt > 4)  return false; 
    if(p->neuHadIso03-0.005*p->pt > 4)  return false; 
  }
  if(p->r9 > 0.9) { 
    if( lIsBarrel && p->hovere > 0.082) return false;
    if(!lIsBarrel && p->hovere > 0.075) return false;
    if( lIsBarrel && p->sieie  > 0.014) return false;
    if(!lIsBarrel && p->sieie  > 0.034) return false;
    if(p->gammaIso03 -0.012*p->pt > 50)  return false; 
    if(p->neuHadIso03-0.005*p->pt > 50)  return false; 
  }
  if(p->chHadIso03             < 2.8) return false;  //This is not really correct, but close
  if((p->neuHadIso03 + p->gammaIso03)-0.17*iRho < 3  ) return false;
  return true;
}
bool PhotonLoader::passCiCPFIso(TPhoton *p,float iRho) { 
  bool lIsBarrel = true; 
  if(fabs(p->scEta)<1.566 && fabs(p->scEta)>1.4443) return false;
  if(fabs(p->scEta)>2.5)                            return false;
  if(fabs(p->scEta)>1.566 && fabs(p->scEta)<2.5) lIsBarrel = false;

  int lCat = 1;
  if ( !lIsBarrel   ) lCat = 3;
  if ( p->r9 < 0.94 ) lCat++;
  float cic4_allcuts_temp_sublead[] = {
    6.0, 4.7, 5.6, 3.6,
    10.0, 6.5, 5.6, 4.4,
    3.8, 2.5, 3.1, 2.2,
    0.0108, 0.0102, 0.028, 0.028,
    0.124, 0.092, 0.142, 0.063,
    0.94, 0.28, 0.94, 0.24 }; 
  bool lPass = false;
  if((p->gammaIso03+p->hcalIso04+p->chHadIso03 - 0.17*iRho)*50./p->pt < cic4_allcuts_temp_sublead[lCat-1+0*4] && 
     (p->ecalIso04 +p->hcalIso04+p->chHadIso03 - 0.52*iRho)*50./p->pt < cic4_allcuts_temp_sublead[lCat-1+1*4] && //Note the track Iso here is 04 not 03 in reality + worst Vtx is used (!==?)
     p->chHadIso03/p->pt                                              < cic4_allcuts_temp_sublead[lCat-1+2*4] && 
     p->sieie                                                         < cic4_allcuts_temp_sublead[lCat-1+3*4] && 
     p->hovere                                                        < cic4_allcuts_temp_sublead[lCat-1+4*4] && 
     p->r9							      > cic4_allcuts_temp_sublead[lCat-1+5*4])  lPass = true;
  return lPass;
}
