//Driver to run Bacon 

#include "../include/GenLoader.hh"
#include "../include/EvtLoader.hh"
#include "../include/ElectronLoader.hh"
#include "../include/MuonLoader.hh"
#include "../include/PhotonLoader.hh"
#include "../include/TauLoader.hh"
#include "../include/JetLoader.hh"

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
TauLoader       *fTau      = 0; 
PhotonLoader    *fPhoton   = 0; 
JetLoader       *fJet      = 0; 

TTree* load(std::string iName) { 
  TFile *lFile = TFile::Open(iName.c_str());
  TTree *lTree = (TTree*) lFile->FindObjectAny("Events");
  return lTree;
}
int main( int argc, char **argv ) {
  gROOT->ProcessLine("#include <vector>");          
  int maxEvents     = atoi(argv[1]);
  std::string lName = argv[2];
  bool        lGen  = atoi(argv[3]);
  //void runBacon(int iNEvents=10,std::string lName="test.root",bool lGen=false) {   
  TTree *lTree = load(lName); 
  if(lTree->GetEntries() < maxEvents || maxEvents == -1) maxEvents = lTree->GetEntries(); 
  //Declare Readers
  fEvt      = new EvtLoader     (lTree);
  fMuon     = new MuonLoader    (lTree);
  fElectron = new ElectronLoader(lTree); 
  fTau      = new TauLoader     (lTree); 
  fPhoton   = new PhotonLoader  (lTree); 
  fJet      = new JetLoader     (lTree);
  if(lGen) fGen      = new GenLoader     (lTree);

  TFile *lFile = new TFile("Output.root","RECREATE");
  TTree *lOut  = new TTree("Tree","Tree");
  //Setup Tree
  fEvt ->setupTree      (lOut);
  fJet ->setupTree      (lOut); 
  if(lGen) fGen ->setupTree      (lOut);
  //Add the triggers we want
  fEvt ->addTrigger("HLT_MET80_Parked_v*");
  fEvt ->addTrigger("HLT_MET80_Parked_v*");
  fEvt ->addTrigger("HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v*");
  fEvt ->addTrigger("HLT_MET100_HBHENoiseCleaned_v*");
  fEvt ->addTrigger("HLT_MET120_HBHENoiseCleaned_v*");
  for(int i0 = 0; i0 < maxEvents; i0++) { 
    if(i0 % 1000 == 0) std::cout << "===> Processed " << i0 << " - Done : " << (float(i0)/float(maxEvents)) << std::endl;
    //Load event and require trigger
    fEvt    ->load(i0);
    if(!fEvt->passTrigger()) continue;
    //Apply all Lepton  Vetoes => do we want a pt cut?
    fTau     ->load(i0);
    if(fTau->vetoTau()) continue;
    fMuon    ->load(i0);
    if(fMuon->vetoMu()) continue;
    fElectron ->load(i0);
    if(fElectron->vetoEle(fEvt->fRho)) continue;
    //Apply a photon veto? (to be checked)
    //fPhoton ->load(i0);
    //if(!fPhoton->vetoPhoton()) continue;

    //Select a Jet (to be fixed by Dan)
    fJet->load(i0); 
    if(!fJet->selectSingleJet()) continue;
    //Make the flat ntuple
    fEvt    ->fillEvent();
    //This needs to be fixed
    if(lGen) fGen->load(i0);
    if(lGen) fGen->selectBoson();
    lOut->Fill();
  }
  lFile->cd();
  lOut->Write();
  lFile->Close();
}
