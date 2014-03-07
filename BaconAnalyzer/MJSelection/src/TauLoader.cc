#include "../include/TauLoader.hh"

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
  fPt   = 0; 
  fEta  = 0; 
  fPhi  = 0; 
}
void TauLoader::setupTree(TTree *iTree) { 
  reset();
  fTree = iTree;
  fTree->Branch("pt_1"  ,&fPt ,"fPt/F");
  fTree->Branch("eta_1" ,&fEta,"fEta/F");
  fTree->Branch("phi_1" ,&fPhi,"fPhi/F");
}
void TauLoader::load(int iEvent) { 
  fTaus   ->Clear();
  fTauBr ->GetEntry(iEvent);
}
bool TauLoader::selectSingleTau() {
  reset(); 
  int lCount = 0; 
  TTau *lTau = 0; 
  for  (int i0 = 0; i0 < fTaus->GetEntriesFast(); i0++) { 
    TTau *pTau = (TTau*)((*fTaus)[i0]);
    if(!passLoose(pTau)) continue;
    lTau = pTau;
    lCount++;
    if(lCount > 1) return false;
  }
  if(lTau == 0) return false;
  fPt  = lTau->pt;
  fEta = lTau->eta;
  fPhi = lTau->phi;
  return true;
}
bool TauLoader::vetoTau() {
  for  (int i0 = 0; i0 < fTaus->GetEntriesFast(); i0++) { 
    TTau *pTau = (TTau*)((*fTaus)[i0]);
    if(passVeto(pTau)) return true;
  }
  return false;
}
//H=>ZZ Mu Id
bool TauLoader::passLoose(TTau *iTau) { 
  if(!(iTau->hpsDisc & kByMVA3LooseElectronRejection &&
       iTau->hpsDisc & kByLooseMuonRejection  &&
       iTau->hpsDisc & kByDecayModeFinding))   return false;  
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
