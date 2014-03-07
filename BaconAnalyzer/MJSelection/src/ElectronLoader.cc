#include "../include/ElectronLoader.hh"
#include <cmath>


using namespace baconhep;

#define ELE_REFERENCE_IDMVA_CUT_BIN0  0.470   // pT<10, |eta|<0.8                                                                                                                               
#define ELE_REFERENCE_IDMVA_CUT_BIN1  0.004   // pT<10, 0.8<|eta|<1.479                                                                                                                             
#define ELE_REFERENCE_IDMVA_CUT_BIN2  0.295   // pT<10, |eta|>1.479                                                                                                                                  
#define ELE_REFERENCE_IDMVA_CUT_BIN3 -0.340   // pT>10, |eta|<0.8                                                                                                                                    
#define ELE_REFERENCE_IDMVA_CUT_BIN4 -0.650   // pT>10, 0.8<|eta|<1.479                                                                                                                               
#define ELE_REFERENCE_IDMVA_CUT_BIN5  0.600   // pT>10, |eta|>1.479                                                                                                                                   

ElectronLoader::ElectronLoader(TTree *iTree) { 
  fElectrons  = new TClonesArray("baconhep::TElectron");
  iTree->SetBranchAddress("Electron",       &fElectrons);
  fElectronBr  = iTree->GetBranch("Electron");
}
ElectronLoader::~ElectronLoader() { 
  delete fElectrons;
  delete fElectronBr;
}
void ElectronLoader::reset() { 
  fPt   = 0; 
  fEta  = 0; 
  fPhi  = 0; 
}
void ElectronLoader::setupTree(TTree *iTree) { 
  reset();
  fTree = iTree;
  fTree->Branch("pt_1"  ,&fPt ,"fPt/F");
  fTree->Branch("eta_1" ,&fEta,"fEta/F");
  fTree->Branch("phi_1" ,&fPhi,"fPhi/F");
}
void ElectronLoader::load(int iEvent) { 
  fElectrons   ->Clear();
  fElectronBr ->GetEntry(iEvent);
}
bool ElectronLoader::selectSingleEle(float iRho) {
  reset(); 
  TElectron *lElectron = 0; 
  int lCount = 0; 
  for  (int i0 = 0; i0 < fElectrons->GetEntriesFast(); i0++) { 
    TElectron *pElectron = (TElectron*)((*fElectrons)[i0]);
    if(!passLoose(pElectron,iRho)) continue;
    lElectron = pElectron;
    lCount++;
    if(lCount > 1) return false;
  }
  if(lElectron == 0) return false;
  fPt  = lElectron->pt;
  fEta = lElectron->eta;
  fPhi = lElectron->phi;
  return true;
}
bool ElectronLoader::vetoEle(float iRho) {
  for  (int i0 = 0; i0 < fElectrons->GetEntriesFast(); i0++) { 
    TElectron *pElectron = (TElectron*)((*fElectrons)[i0]);
    if(passVeto(pElectron,iRho)) return true;
  }
  return false;
}
//POG based veto id
bool  ElectronLoader::passVeto(const TElectron *electron, const float iRho) {
  const double ECAL_GAP_LOW  = 1.4442;
  const double ECAL_GAP_HIGH = 1.566;
  
  if((fabs(electron->scEta)>ECAL_GAP_LOW) && (fabs(electron->scEta)<ECAL_GAP_HIGH)) return false;
  if(!(electron->typeBits & kEcalDriven)) return false;
  
  if(fabs(electron->d0) > 0.02) return false;
  if(fabs(electron->dz) > 0.1)  return false;
  
  // conversion rejection
  if(electron->nMissingHits > 1) return false;
  if(electron->isConv)            return false;
     
  double ea = getEffArea(electron->scEta);
  // barrel/endcap dependent requirments      
  if(fabs(electron->scEta)<=ECAL_GAP_LOW) {
    // barrel
    double iso =  electron->chHadIso03 + TMath::Max(electron->neuHadIso03 + electron->gammaIso03 - iRho*ea, 0.);
    if(iso > 0.15*(electron->pt)) return false;
     
    if(electron->sieie  > 0.01)                       return false;
    if(fabs(electron->dPhiIn) < 0.06)                   return false;
    if(fabs(electron->dEtaIn) < 0.004)                  return false;
    if(electron->hovere   > 0.12)                           return false;
    if(fabs(1.0-electron->eoverp) > 0.05*(electron->ecalEnergy)) return false;
  
  } else {
    // endcap
    Double_t iso = electron->chHadIso03 + TMath::Max(electron->neuHadIso03 + electron->gammaIso03 - iRho*ea, 0.);
    if(iso > 0.15*(electron->pt))                           return false;
    if(electron->sieie  > 0.03)                       return false;
    if(fabs(electron->dPhiIn) < 0.03)                   return false;
    if(fabs(electron->dEtaIn) < 0.007)                  return false;
    if(electron->hovere   > 0.10)                           return false;
    if(fabs(1.0-electron->eoverp) > 0.05*(electron->ecalEnergy)) return false;
  }

  return kTRUE;
}
//H=>ZZ Ele Id
bool ElectronLoader::passLoose(const TElectron *ele,float iRho) { 
  // missing hits cut for conversion rejection                                                                                                                                                     
  if(ele->nMissingHits > 1) return false;
  
  // impact parameters cuts
  if(fabs(ele->sip3d) >= 100) return false;
  if(fabs(ele->d0)    >= 0.5) return false;
  if(fabs(ele->dz)    >= 1.0) return false;
  
  int ptBin = (ele->ptHZZ4l > 10) ? 1 : 0;
  int etaBin = -1;
  if     (fabs(ele->scEta) < 0.8)   etaBin = 0;
  else if(fabs(ele->scEta) < 1.479) etaBin = 1;
  else                              etaBin = 2;
  if(ptBin == 0 && etaBin == 0) return (ele->mva > ELE_REFERENCE_IDMVA_CUT_BIN0);
  if(ptBin == 0 && etaBin == 1) return (ele->mva > ELE_REFERENCE_IDMVA_CUT_BIN1);
  if(ptBin == 0 && etaBin == 2) return (ele->mva > ELE_REFERENCE_IDMVA_CUT_BIN2);
  if(ptBin == 1 && etaBin == 0) return (ele->mva > ELE_REFERENCE_IDMVA_CUT_BIN3);
  if(ptBin == 1 && etaBin == 1) return (ele->mva > ELE_REFERENCE_IDMVA_CUT_BIN4);
  if(ptBin == 1 && etaBin == 2) return (ele->mva > ELE_REFERENCE_IDMVA_CUT_BIN5);
  double lIso = ele->chHadIso04 + TMath::Max(ele->gammaIso04 + ele->neuHadIso04 - iRho*getEffArea(ele->scEta), 0.);
  if(lIso/ele->pt > 0.4) return kFALSE;
  return true;
}
double ElectronLoader::getEffArea(double eta) {
  if     (fabs(eta) < 1.0)   return 0.19;
  else if(fabs(eta) < 1.479) return 0.25;
  else if(fabs(eta) < 2.0)   return 0.12;
  else if(fabs(eta) < 2.2)   return 0.21;
  else if(fabs(eta) < 2.3)   return 0.27;
  else if(fabs(eta) < 2.4)   return 0.44;
  else                       return 0.52;
  return 0.52;
}

//bool ElectronLoader::passTight(TElectron *iElectron) { 
//  return kTRUE;
//}
