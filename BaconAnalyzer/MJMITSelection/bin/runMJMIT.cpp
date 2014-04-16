//Driver to run Bacon 

#include "../include/GenLoader.hh"
#include "../include/EvtLoader.hh"
#include "../include/ElectronLoader.hh"
#include "../include/MuonLoader.hh"
#include "../include/LeptonLoader.hh"
#include "../include/PhotonLoader.hh"
#include "../include/TauLoader.hh"
#include "../include/JetLoader.hh"
#include "../include/RunLumiRangeMap.h"

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include <string>
#include <iostream>

//Object Processors
GenLoader       *fGen      = 0; 
EvtLoader       *fEvt      = 0; 
MuonLoader      *fMuon     = 0; 
ElectronLoader  *fElectron = 0; 
LeptonLoader    *fLepton   = 0; 
TauLoader       *fTau      = 0; 
PhotonLoader    *fPhoton   = 0; 
JetLoader       *fJet      = 0; 
RunLumiRangeMap *fRangeMap = 0; 

TH1F *fHist                = 0; 

TTree* load(std::string iName) { 
  TFile *lFile = TFile::Open(iName.c_str());
  TTree *lTree = (TTree*) lFile->FindObjectAny("Events");
  fHist        = (TH1F* ) lFile->FindObjectAny("TotalEvents");
  return lTree;
}
bool passEvent(unsigned int iRun,unsigned int iLumi) { 
  RunLumiRangeMap::RunLumiPairType lRunLumi(iRun,iLumi);
  return fRangeMap->HasRunLumi(lRunLumi);
}
int main( int argc, char **argv ) {
  gROOT->ProcessLine("#include <vector>");          
  int maxEvents     = atoi(argv[1]);
  std::string lName = argv[2];
  bool        lGen  = atoi(argv[3]);
  std::string lJSON = argv[4];
  int         lDMu  = atoi(argv[5]);
  double      lXS   = atof(argv[6]);
  fRangeMap = new RunLumiRangeMap();
  if(lJSON.size() > 0) fRangeMap->AddJSONFile(lJSON.c_str());

  TTree *lTree = load(lName); 
  if(lTree->GetEntries() < maxEvents || maxEvents == -1) maxEvents = lTree->GetEntries(); 
  //Declare Readers
  fEvt      = new EvtLoader     (lTree);
  fMuon     = new MuonLoader    (lTree);
  fElectron = new ElectronLoader(lTree);
  fLepton   = new LeptonLoader  (lTree);
  fTau      = new TauLoader     (lTree); 
  fPhoton   = new PhotonLoader  (lTree); 
  fJet      = new JetLoader     (lTree);
  if(lGen) fGen      = new GenLoader     (lTree);
  if(lDMu == 2) {fMuon->fMassMin = 60; fMuon->fMassMax = 120;}

  TFile *lFile = new TFile("Output.root","RECREATE");
  TTree *lOut  = new TTree("Tree","Tree");
  //Setup Tree
  float lWeight = float(lXS)/float(fHist->Integral())/1000.; if(!lGen) lWeight = 1.;
  fEvt     ->setupTree      (lOut,lWeight);
  fJet     ->setupTree      (lOut); 
  fMuon    ->setupTree      (lOut); 
  fElectron->setupTree      (lOut); 
  fLepton  ->setupTree      (lOut); 
  fTau     ->setupTree      (lOut); 
  fPhoton  ->setupTree      (lOut); 
  if(lGen) fGen ->setupTree      (lOut);
  //Add the triggers we want
  fEvt ->addTrigger("HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v*");
  fEvt ->addTrigger("HLT_DiPFJet40_PFMETnoMu65_MJJ800VBF_AllJets_v*");

  fEvt ->addTrigger("HLT_MET80_Parked_v*");
  fEvt ->addTrigger("HLT_MET80_Parked_v*");
  fEvt ->addTrigger("HLT_MET100_HBHENoiseCleaned_v*");
  fEvt ->addTrigger("HLT_MET120_HBHENoiseCleaned_v*");

  fJet ->addTrigger("HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v*");
  fJet ->addTrigger("HLT_DiPFJet40_PFMETnoMu65_MJJ800VBF_AllJets_v*");

  for(int i0 = 0; i0 < maxEvents; i0++) { 
    if(i0 % 1000 == 0) std::cout << "===> Processed " << i0 << " - Done : " << (float(i0)/float(maxEvents)) << std::endl;
    //Load event and require trigger
    std::vector<TLorentzVector> lVetoes; 
    if(lDMu > 0) fMuon->load(i0);
    if(lDMu > 0) fMuon->selectDiMuon(lVetoes);
    if(lDMu > 0 && lVetoes.size() == 0) continue;
    fEvt     ->load(i0);
    fEvt     ->fillEvent(lVetoes);
    if(lDMu == 0) if(!fEvt->passSkim()) continue;
    if(!lGen && !passEvent(fEvt->fRun,fEvt->fLumi)) continue;

    if(lDMu == 0) fMuon    ->load(i0);    
    fMuon    ->selectMuons(lVetoes);
    fMuon    ->fillVetoes(lVetoes);
    
    fElectron->load(i0);
    fElectron->selectElectrons(fEvt->fRho,lVetoes); //Add a muon veto?
    fElectron->fillVetoes(lVetoes);

    //std::cout << "Size ===> " << fMuon->fSelMuons.size() << " -- " << fElectron->fSelElectrons.size() << std::endl;
    fLepton  ->fillLeptons(fMuon->fSelMuons,fElectron->fSelElectrons); 

    fTau     ->load(i0);
    fTau     ->selectTaus(lVetoes);
    fTau     ->fillVetoes(lVetoes);

    fPhoton  ->load(i0);
    fPhoton  ->selectPhotons(lVetoes,fEvt->fRho);
    
    fJet->load(i0); 
    fJet->selectJets(lVetoes);

    if(lGen) fGen->load(i0);
    if(lGen) fGen->selectBoson();
    if(lGen) fGen->fillGenEvent();
    lOut->Fill();
  }
  lFile->cd();
  lOut->Write();
  lFile->Close();
}
