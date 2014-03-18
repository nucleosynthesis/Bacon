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
  fPt   = 0; 
  fEta  = 0; 
  fPhi  = 0; 
}
void PhotonLoader::setupTree(TTree *iTree) { 
  reset();
  fTree = iTree;
  fTree->Branch("pt_1"  ,&fPt ,"fPt/F");
  fTree->Branch("eta_1" ,&fEta,"fEta/F");
  fTree->Branch("phi_1" ,&fPhi,"fPhi/F");
}
void PhotonLoader::load(int iEvent) { 
  fPhotons   ->Clear();
  fPhotonBr ->GetEntry(iEvent);
}
bool PhotonLoader::selectSinglePhoton() {
  reset(); 
  int lCount = 0; 
  TPhoton *lPhoton = 0; 
  for  (int i0 = 0; i0 < fPhotons->GetEntriesFast(); i0++) { 
    TPhoton *pPhoton = (TPhoton*)((*fPhotons)[i0]);
    if(!passLoose(pPhoton)) continue;
    lPhoton = pPhoton;
    lCount++;
    if(lCount > 1) return false;
  }
  if(lPhoton == 0) return false;
  fPt  = lPhoton->pt;
  fEta = lPhoton->eta;
  fPhi = lPhoton->phi;
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
