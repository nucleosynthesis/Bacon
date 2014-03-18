#include "../include/JetLoader.hh"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "TMath.h"

#include <cmath>

using namespace baconhep;

JetLoader::JetLoader(TTree *iTree) { 
  fJets  = new TClonesArray("baconhep::TJet");
  iTree->SetBranchAddress("Jet05",       &fJets);
  fJetBr  = iTree->GetBranch("Jet05");
}
JetLoader::~JetLoader() { 
  delete fJets;
  delete fJetBr;
}
void JetLoader::reset() { 
  fPt   = 0; 
  fEta  = 0; 
  fPhi  = 0; 
}
void JetLoader::setupTree(TTree *iTree) { 
  reset();
  fTree = iTree;
  fTree->Branch("pt_1"  ,&fPt ,"fPt/F");
  fTree->Branch("eta_1" ,&fEta,"fEta/F");
  fTree->Branch("phi_1" ,&fPhi,"fPhi/F");
  fTree->Branch("m_1"   ,&fM  ,"fM/F");
}
void JetLoader::load(int iEvent) { 
  fJets   ->Clear();
  fJetBr ->GetEntry(iEvent);
}
bool JetLoader::selectSingleJet() {
  reset(); 
  TJet *lJet = 0; 
  int lCount = 0;
  for  (int i0 = 0; i0 < fJets->GetEntriesFast(); i0++) { 
    TJet *pJet = (TJet*)((*fJets)[i0]);
    if(!passLoose(pJet)) continue;
    lJet = pJet;
    lCount++;
    if(lCount > 1) return false;
  }
  if(lJet == 0) return false;
  fPt   = lJet->pt;
  fEta  = lJet->eta;
  fPhi  = lJet->phi;
  fM    = lJet->mass;
  return true;
}
bool JetLoader::vetoAdditionalJets(int lead_jet_i){
  int ncounts = 0 ;
  TJet *lJet = (TJet*)((*fJets)[lead_jet_i]);
  TJet *sJet = 0;
  for  (int i0 = 0; i0 < fJets->GetEntriesFast(); i0++) { 
    if (i0 == lead_jet_i) continue;
    TJet *pJet = (TJet*)((*fJets)[i0]);
    if(passVeto(pJet)) {
    	if (!(pJet->pt > 60.)) continue;
    	if (!(fabs(pJet->eta) < 4.5)) continue;
	sJet = pJet;
	ncounts++;
    }
  }
  if (ncounts == 0) return false; // no additional jets
  if (ncounts > 1) return true;  // too many jets
  if (sJet){
     if(!(fabs(reco::deltaPhi(lJet->eta,sJet->eta) < 2.5)))  return true; 
  }
  return false;
  
}
bool JetLoader::selectJetsWithAdditionalJetsVeto() {
  reset(); 
  TJet *lJet = 0; 
  int lCount = 0;
  int selected_jet_index = 0;
  for  (int i0 = 0; i0 < fJets->GetEntriesFast(); i0++) { 
    TJet *pJet = (TJet*)((*fJets)[i0]);
    if(!passLoose(pJet)) continue;
    lJet = pJet; selected_jet_index = i0;
    lCount++;
    //if(lCount > 1) return false;
  }
  if (lJet == 0) return false;
  if (vetoAdditionalJets(selected_jet_index)) return false;
  fPt   = lJet->pt;
  fEta  = lJet->eta;
  fPhi  = lJet->phi;
  fM    = lJet->mass;
  return true;
}
bool JetLoader::vetoJet() {
  for  (int i0 = 0; i0 < fJets->GetEntriesFast(); i0++) { 
    TJet *pJet = (TJet*)((*fJets)[i0]);
    if(passVeto(pJet)) return true;
  }
  return false;
}
//H=>ZZ Mu Id
bool JetLoader::passLoose(TJet *iJet) { 
  if(iJet->neuEmFrac        >  0.99)                         return false;
  if(iJet->neuHadFrac       >  0.99)                         return false;
  if(iJet->nParticles       <  2)                            return false;
  if(iJet->chHadFrac        <= 0     && fabs(iJet->eta) < 2.4 )    return false;
  if(iJet->chEmFrac         >  0.99  && fabs(iJet->eta) < 2.4 )    return false;
  if(iJet->nCharged         < 1      && fabs(iJet->eta) < 2.4 )    return false;
  return true;
}
bool JetLoader::passTight(TJet *iJet) { 
  if(iJet->neuEmFrac        >  0.9)                          return false;
  if(iJet->neuHadFrac       >  0.9)                          return false;
  if(iJet->nParticles       <  2)                            return false;
  if(iJet->chHadFrac        <= 0     && fabs(iJet->eta) < 2.4 )    return false;
  if(iJet->chEmFrac         >  0.9   && fabs(iJet->eta) < 2.4 )    return false;
  if(iJet->nCharged         < 1      && fabs(iJet->eta) < 2.4 )    return false;
  return true;
}
bool JetLoader::passVeto(TJet *iJet) { 
  return passLoose(iJet);
}
TLorentzVector JetLoader::jet() { 
  TLorentzVector lJet; 
  lJet.SetPtEtaPhiM(fPt,fEta,fPhi,fM);
  return lJet;
}
