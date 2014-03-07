#include "../include/EvtLoader.hh"

using namespace baconhep;

EvtLoader::EvtLoader(TTree *iTree,std::string iHLTFile) { 
  fEvt  = new TEventInfo();
  iTree->SetBranchAddress("Info",       &fEvt);
  fEvtBr  = iTree->GetBranch("Info");
  fTrigger = new TTrigger(iHLTFile);
}
EvtLoader::~EvtLoader() { 
  delete  fEvt;
  delete  fEvtBr;
}
void EvtLoader::reset() { 
  fRun       = 0;
  fEvtV      = 0; 
  fLumi      = 0; 
  fMet       = 0; 
  fMetPhi    = 0; 
  fTKMet     = 0; 
  fTKMetPhi  = 0; 
  fMVAMet    = 0; 
  fMVAMetPhi = 0; 
  fRho       = 0; 
}
void EvtLoader::setupTree(TTree *iTree) { 
  reset();
  fTree = iTree;
  fTree->Branch("run"  ,&fRun   ,"fRun/i");
  fTree->Branch("lumi" ,&fLumi  ,"fLumi/i");
  fTree->Branch("evt"  ,&fEvtV  ,"fEvtV/i");
  fTree->Branch("met"  ,&fMet   ,"fMet/F");
  fTree->Branch("metphi"     ,&fMetPhi   ,"fMetPhi/F");
  fTree->Branch("tkmet"      ,&fTKMet    ,"fTKMet/F");
  fTree->Branch("tkmetphi"   ,&fTKMetPhi ,"fTKMetPhi/F");
  fTree->Branch("mvamet"     ,&fMVAMet   ,"fMVAMet/F");
  fTree->Branch("mvametphi"  ,&fMVAMetPhi,"fMVAMetPhi/F");
  fTree->Branch("rho"        ,&fRho      ,"fRho/F");
}
void EvtLoader::resetRecoil() { 
  fU1     = 0;
  fU2     = 0; 
  fTKU1   = 0; 
  fTKU2   = 0; 
  fMVAMt  = 0; 
  fMt     = 0; 
  fTKMt   = 0; 
  fTKUDPhi = 0; 
  fUDPhi   = 0; 
}
void EvtLoader::setupRecoilTree() { 
  resetRecoil();
  fTree->Branch("u1"       ,&fU1      ,"fU1/F");
  fTree->Branch("u2"       ,&fU2      ,"fU2/F");
  fTree->Branch("tku1"     ,&fTKU1    ,"fEvt/F");
  fTree->Branch("tku2"     ,&fTKU2    ,"fMet/F");
  fTree->Branch("mvamt"    ,&fMVAMt   ,"fMVAMt/F");
  fTree->Branch("mt"       ,&fMt      ,"fTKMet/F");
  fTree->Branch("tkmt"     ,&fTKMt    ,"fTKMetPhi/F");
  fTree->Branch("tkudphi"  ,&fTKUDPhi ,"fTKUDPhi/F");
  fTree->Branch("udphi"    ,&fUDPhi   ,"fUDPhi/F");
}
void EvtLoader::load(int iEvent) { 
  fEvtBr ->GetEntry(iEvent);
}
void EvtLoader::addTrigger(std::string iName) { 
  fTrigString.push_back(iName);
}
bool EvtLoader::passTrigger() {
  bool lPass = false;
  for(unsigned int i0 = 0; i0 < fTrigString.size(); i0++) { 
    if(fTrigger->pass(fTrigString[i0],fEvt->triggerBits)) lPass = true;
  }
  return lPass;
}
bool EvtLoader::passTrigger(std::string iTrigger) { 
  return fTrigger->pass(iTrigger,fEvt->triggerBits);
}
void EvtLoader::fillEvent() {
  reset(); 
  fRun     = fEvt->runNum;
  fLumi    = fEvt->lumiSec;
  fEvtV    = fEvt->evtNum;
  fMet     = fEvt->pfMET;
  fMetPhi  = fEvt->pfMETphi;
  fMVAMet  = fEvt->mvaMETU;
  fMVAMetPhi  = fEvt->mvaMETUphi;
  fTKMet      = fEvt->trkMET;
  fTKMetPhi   = fEvt->trkMETphi;
  fRho        = fEvt->rhoJet;
  return;
}
void EvtLoader::fillRecoil(TLorentzVector &iLep) {
  resetRecoil(); 
  TLorentzVector lMet    = Met(0);
  TLorentzVector lTKMet  = Met(1);
  TLorentzVector lMVAMet = Met(3);
  TLorentzVector lU      = recoil(iLep,lMet);
  TLorentzVector lTKU    = recoil(iLep,lTKMet);
  fU1   = lU  .Px();
  fU2   = lU  .Py();
  fTKU1 = lTKU.Px();
  fTKU2 = lTKU.Py();
  fMt      = (lMet    + iLep).M();
  fMVAMt   = (lMVAMet + iLep).M();
  fTKMt    = (lTKMet  + iLep).M();
  //Un rotate the system
  lU  .RotateZ(iLep.Phi());
  lTKU.RotateZ(iLep.Phi());
  fUDPhi   = lU  .DeltaPhi(iLep);
  fTKUDPhi = lTKU.DeltaPhi(iLep);
  return;
}
TLorentzVector EvtLoader::Met(int iOption) { 
  TLorentzVector lVec;
  if(iOption == 0) {lVec.SetPtEtaPhiM(fEvt->pfMET  ,0.,fEvt->pfMETphi  ,0.);}
  if(iOption == 1) {lVec.SetPtEtaPhiM(fEvt->trkMET ,0.,fEvt->trkMETphi ,0.);}
  if(iOption == 2) {lVec.SetPtEtaPhiM(fEvt->mvaMET ,0.,fEvt->mvaMETphi ,0.);}
  if(iOption == 3) {lVec.SetPtEtaPhiM(fEvt->mvaMETU,0.,fEvt->mvaMETUphi,0.);}
  return lVec;
}
TLorentzVector EvtLoader::recoil(TLorentzVector &iLep,TLorentzVector &iMet) { 
  //Subtract Lepton to make recoil (note lepton is negavtive vector sum)
  TLorentzVector lU; lU.SetPtEtaPhiM(iLep.Pt(),0.,iLep.Eta(),0.);
  lU+=iMet;
  lU.RotateZ(TMath::Pi());
  //Now Project it on the Lepton
  lU.RotateZ(-iLep.Phi());
  return lU;
}
