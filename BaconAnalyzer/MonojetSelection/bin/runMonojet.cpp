//Driver to run Bacon for Monojet SUS-13-009

#include "../include/GenLoader.hh"
#include "../include/EvtLoader.hh"
#include "../include/MuonLoader.hh"
#include "../include/TauLoader.hh"
#include "../include/ElectronLoader.hh"
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
TauLoader       *fTau      = 0; 
ElectronLoader  *fElectron = 0; 
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
  fJet	    = new JetLoader	(lTree); 
  fTau	    = new TauLoader	(lTree); 
  if(lGen) fGen = new GenLoader (lTree);

  TFile *lFile = new TFile("Output.root","RECREATE");
  TTree *lOut  = new TTree("Tree","Tree");

  //Setup Tree
  fEvt ->setupTree      (lOut);  // Event Output variables 
  fJet ->setupTree      (lOut);  // Jet Output variables 
  
  //fMuon->setupTree      (lOut); 
  if(lGen) fGen ->setupTree      (lOut);
  for(int i0 = 0; i0 < maxEvents; i0++) { 
    if(i0 % 1000 == 0) std::cout << "===> Processed " << i0 << " - Done : " << (float(i0)/float(maxEvents)) << std::endl;

    // Load Ntuple 
    fEvt    ->load(i0);
    fEvt    ->fillEvent();
    float rho = fEvt->getEventRho();

    // Select events with 1 Jet PT> 110 |eta|<2.4
    fJet->load(i0);

    // Will veto on these objects
    fElectron->load(i0);
    fMuon->load(i0);
    fTau->load(i0);
    if (fTau->vetoTau())      continue;
    if (fMuon->vetoMu())      continue;
    if (fElectron->vetoEle(rho)) continue;

    if (!fJet->selectJetsWithAdditionalJetsVeto()) continue;
    // Finally apply cuts to Jet (can leave out MET cut for now
    TLorentzVector lJet = fJet->jet();
    if (!(lJet.Pt() > 110)) continue;
    if (!(fabs(lJet.Eta() < 2.4 ))) continue;
    //Load the rest and fill
    lOut->Fill();
  }
  lFile->cd();
  lOut->Write();
  lFile->Close();
}
