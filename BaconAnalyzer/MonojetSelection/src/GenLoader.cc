#include <iostream>
#include <assert.h>  
#include "../include/GenLoader.hh"

using namespace baconhep;

GenLoader::GenLoader(TTree *iTree) { 
  fGens  = new TClonesArray("baconhep::TGenParticle");
  iTree->SetBranchAddress("GenParticle",       &fGens);
  fGenBr  = iTree->GetBranch("GenParticle");
}
GenLoader::~GenLoader() { 
  delete fGens;
  delete fGenBr;
}
void GenLoader::reset() { 
  fVPt   = 0; 
  fVEta  = 0; 
  fVPhi  = 0; 
  fVM    = 0; 
  fVId   = 0; 
  
  fPt1   = 0; 
  fEta1  = 0; 
  fPhi1  = 0; 
  fM1    = 0; 
  fId1   = 0; 

  fPt2   = 0; 
  fEta2  = 0; 
  fPhi2  = 0; 
  fM2    = 0; 
  fId2   = 0; 
}
void GenLoader::setupTree(TTree *iTree) { 
  reset();
  fTree = iTree;
  fTree->Branch("genvpt"   ,&fVPt  ,"fVPt/F");
  fTree->Branch("genveta"  ,&fVEta ,"fVEta/F");
  fTree->Branch("genvphi"  ,&fVPhi ,"fVPhi/F");
  fTree->Branch("genvm"    ,&fVM   ,"fVM/F");
  fTree->Branch("genvid"   ,&fVId  ,"fVId/I");

  fTree->Branch("genpt_1"  ,&fPt1 ,"fPt1/F");
  fTree->Branch("geneta_1" ,&fEta1,"fEta1/F");
  fTree->Branch("genphi_1" ,&fPhi1,"fPhi1/F");
  fTree->Branch("genm_1"   ,&fM1  ,"fM1/F");
  fTree->Branch("genid_1"  ,&fId1 ,"fId1/I");

  fTree->Branch("genpt_2"  ,&fPt2 ,"fPt2/F");
  fTree->Branch("geneta_2" ,&fEta2,"fEta2/F");
  fTree->Branch("genphi_2" ,&fPhi2,"fPhi2/F");
  fTree->Branch("genm_2"   ,&fM2  ,"fM2/F");
  fTree->Branch("genid_2"  ,&fId2 ,"fId2/I");
}
void GenLoader::resetRecoil() { 
  fGenWLepPhi = 0; 
}
void GenLoader::setupRecoilTree() { 
  resetRecoil();
  fTree->Branch("genwlepphi"  ,&fGenWLepPhi,"fGenWLepPhi/F");
}
void GenLoader::load(int iEvent) { 
  fGens   ->Clear();
  fGenBr ->GetEntry(iEvent);
}
void GenLoader::selectBoson() {
  reset(); 
  TGenParticle *lBoson   = 0; 
  TGenParticle *lLep1    = 0; 
  TGenParticle *lLep2    = 0; 
  int   lBosonId = -10; 
  for  (int i0 = 0; i0 < fGens->GetEntriesFast(); i0++) { 
    TGenParticle *pGen = (TGenParticle*)((*fGens)[i0]);
    if(fabs(pGen->pdgId) == 23 ||   // Select Z or
       fabs(pGen->pdgId) == 24 ||   // Select W or
       fabs(pGen->pdgId) == 25) {   // Select Higgs
      lBoson   = pGen;
      lBosonId = i0;
    }
    if(pGen->parent == lBosonId) { //All of these guys have two daughters
      //!!!Note this will not work for taus
      if(pGen->status == 1 && lLep1 != 0) lLep2 = pGen;
      if(pGen->status == 1 && lLep1 == 0) lLep1 = pGen;
      if(pGen->status != 1 && lLep1 != 0) lLep2 = getStatus1(i0);  //Obtain the simulation level if not already
      if(pGen->status != 1 && lLep1 == 0) lLep1 = getStatus1(i0); 
    } 
  }
  fVPt  = lBoson->pt;
  fVEta = lBoson->eta;
  fVPhi = lBoson->phi;
  fVM   = lBoson->mass;
  fVId  = lBoson->pdgId;
  if(lLep1 == 0 || lLep2 == 0) return;
  assert(lLep1);
  assert(lLep2);
  if(lLep2->pt > lLep1->pt || isNeutrino(lLep1)) {  
    TGenParticle *lLep = 0; 
    lLep = lLep1; 
    //Swaps
    lLep1 = lLep2;
    lLep2 = lLep;
  }
  fPt1  = lLep1->pt;
  fEta1 = lLep1->eta;
  fPhi1 = lLep1->phi;
  fM1   = lLep1->mass;
  fId1  = lLep1->pdgId;
  
  fPt2  = lLep2->pt;
  fEta2 = lLep2->eta;
  fPhi2 = lLep2->phi;
  fM2   = lLep2->mass;
  fId2  = lLep2->pdgId;
}
void GenLoader::fillRecoil(TLorentzVector &iLep) { 
  resetRecoil();
  TLorentzVector lBoson = boson();
  fGenWLepPhi = lBoson.DeltaPhi(iLep);
}
TLorentzVector GenLoader::boson() { 
  TLorentzVector lBoson; 
  lBoson.SetPtEtaPhiM(fVPt,fVEta,fVPhi,fVM);
  return lBoson;
}
//H=>ZZ Mu Id
TGenParticle* GenLoader::getStatus1(int iId) { 
  int lId = iId;
  TGenParticle *lGen = 0; 
  for  (int i0 = 0; i0 < fGens->GetEntriesFast(); i0++) { 
    TGenParticle *pGen = (TGenParticle*)((*fGens)[i0]);
    if(pGen->parent == lId) { 
      lGen = pGen;
      if(pGen->status == 1) break; 
      lId = i0;  //Keep searching down the chain for status 1 !!! Assumes gen particle list is ordered
    }
  }    
  assert(lGen); 
  //assert(lGen->status == 1); ===> commented out to fix issues with taus 
  return lGen;
}
bool GenLoader::isNeutrino(TGenParticle *iPart) { 
  if(fabs(iPart->pdgId) == 12) return true; 
  if(fabs(iPart->pdgId) == 14) return true; 
  if(fabs(iPart->pdgId) == 16) return true; 
  return false;
}
