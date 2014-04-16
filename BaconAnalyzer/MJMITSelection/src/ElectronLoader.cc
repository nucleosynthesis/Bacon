#include "../include/ElectronLoader.hh"
#include <cmath>
#include <iostream> 

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
  fNElectrons = 0; 
  fSelElectrons.clear();
}
void ElectronLoader::setupTree(TTree *iTree) { 
  reset();
  fTree = iTree;
  fTree->Branch("nelectrons",&fNElectrons,"fNElectrons/I");
}
void ElectronLoader::load(int iEvent) { 
  fElectrons   ->Clear();
  fElectronBr ->GetEntry(iEvent);
}
void ElectronLoader::fillVetoes(std::vector<TLorentzVector> &iVec) { 
  for(unsigned int i0 = 0; i0 < fSelElectrons.size(); i0++) { 
    TLorentzVector pVec;
    pVec.SetPtEtaPhiM(fSelElectrons[i0]->pt,fSelElectrons[i0]->eta,fSelElectrons[i0]->phi,0.000511);
    iVec.push_back(pVec);
  }
}

bool ElectronLoader::selectElectrons(float iRho,std::vector<TLorentzVector> &iVetoes) {
  reset(); 
  int lCount = 0; 
  for  (int i0 = 0; i0 < fElectrons->GetEntriesFast(); i0++) { 
    TElectron *pElectron = (TElectron*)((*fElectrons)[i0]);
    if(pElectron->pt < 10)     continue;
    if(!passVeto(pElectron,iRho)) continue;
    bool pMatch = false;
    for(unsigned int i1 = 0; i1 < iVetoes.size(); i1++) { 
      double pDEta = pElectron->eta      - iVetoes[i1].Eta();
      double pDPhi = fabs(pElectron->phi - iVetoes[i1].Phi());
      if(fabs(pDPhi) > 2.*TMath::Pi()-fabs(pDPhi)) pDPhi =  2.*TMath::Pi()-fabs(pDPhi);
      if(sqrt(pDPhi*pDPhi+pDEta*pDEta) > 0.3) continue;
      pMatch = true;
    }
    if(pMatch) continue;
    bool lFill = false;
    for( std::vector<TElectron*>::iterator pElectronIter = fSelElectrons.begin(); pElectronIter != fSelElectrons.end(); pElectronIter++) { 
      if((*pElectronIter)->pt > pElectron->pt) continue;
      fSelElectrons.insert(pElectronIter,pElectron);
      lFill = true;
      break;
    } 
    if(!lFill)  fSelElectrons.push_back(pElectron);
    //Limit this to the top 4 Jets
    //if(fSelElectrons.size() > 3) fSelElectrons.pop_back();
    lCount++;
  }
  fNElectrons = lCount;
  if(fSelElectrons.size() == 0) return false;
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
//VBTF95? really?
bool ElectronLoader::passVBTF95(const TElectron *iEle) { 
  const double ECAL_GAP_LOW  = 1.4442;
  const double ECAL_GAP_HIGH = 1.566;
  
  if((fabs(iEle->scEta)>ECAL_GAP_LOW) && (fabs(iEle->scEta)<ECAL_GAP_HIGH)) return false;
  if(!(iEle->typeBits & kEcalDriven)) return false;

  if(fabs(iEle->d0) > 0.02)  return false;
  if(fabs(iEle->dz) > 0.1)   return false;
  
  // conversion rejection
  if(iEle->nMissingHits > 1) return false;
  if(iEle->isConv)           return false;
  double lIso = iEle->chHadIso03 + iEle->neuHadIso03 + iEle->gammaIso03;
  //std::cout << " here : " << iEle->pt << " -- " << lIso/iEle->pt << " - " << std::endl;
  
  if(fabs(iEle->scEta) < ECAL_GAP_LOW && lIso/iEle->pt > 0.13) return false;
  if(fabs(iEle->scEta) > ECAL_GAP_LOW && lIso/iEle->pt > 0.09) return false;  
  if(fabs(iEle->scEta) < ECAL_GAP_LOW && iEle->hovere > 0.15)  return false;
  if(fabs(iEle->scEta) > ECAL_GAP_LOW && iEle->hovere > 0.07)  return false;
  if(fabs(iEle->scEta) < ECAL_GAP_LOW && iEle->sieie  > 0.01)  return false;
  if(fabs(iEle->scEta) > ECAL_GAP_LOW && iEle->sieie  > 0.03)  return false;
  if(fabs(iEle->scEta) < ECAL_GAP_LOW && fabs(iEle->dPhiIn) > 0.8)   return false;
  if(fabs(iEle->scEta) > ECAL_GAP_LOW && fabs(iEle->dPhiIn) > 0.7)   return false;
  if(fabs(iEle->scEta) < ECAL_GAP_LOW && fabs(iEle->dEtaIn) > 0.007)   return false;
  if(fabs(iEle->scEta) > ECAL_GAP_LOW && fabs(iEle->dEtaIn) > 0.009)   return false;
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
std::vector<TLorentzVector> ElectronLoader::conversions() { 
  std::vector<TLorentzVector> lConv;
  for  (int i0 = 0; i0 < fElectrons->GetEntriesFast(); i0++) { 
    TElectron *pElectron = (TElectron*)((*fElectrons)[i0]);
    if(pElectron->nMissingHits == 0) continue;
    if(!pElectron->isConv)           continue;
    TLorentzVector pVec; 
    pVec.SetPtEtaPhiM(pElectron->scEt,pElectron->scEta,pElectron->scPhi,0.000511);
    lConv.push_back(pVec);
  }
  return lConv;
}
std::vector<TLorentzVector> ElectronLoader::nonConversions() { 
  std::vector<TLorentzVector> lConv;
  for  (int i0 = 0; i0 < fElectrons->GetEntriesFast(); i0++) { 
    TElectron *pElectron = (TElectron*)((*fElectrons)[i0]);
    if(pElectron->nMissingHits != 0 || pElectron->isConv) continue;
    TLorentzVector pVec; 
    pVec.SetPtEtaPhiM(pElectron->scEt,pElectron->scEta,pElectron->scPhi,0.000511);
    lConv.push_back(pVec);
  }
  return lConv;
}

//bool ElectronLoader::passTight(TElectron *iElectron) { 
//  return kTRUE;
//}
