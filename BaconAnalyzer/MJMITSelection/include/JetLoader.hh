#include "TTree.h"
#include "TLorentzVector.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "BaconAna/DataFormats/interface/TJet.hh"
#include "BaconAna/DataFormats/interface/TAddJet.hh"
#include "BaconAna/DataFormats/interface/TTrigger.hh"

using namespace baconhep;

class JetLoader { 
public:
  JetLoader(TTree *iTree,std::string iHLTFile="/afs/cern.ch/user/p/pharris/pharris/public/bacon/CMSSW_5_3_13/src/BaconAna/DataFormats/data/HLTFile_v0");
  ~JetLoader();
  void reset();
  void reset(TJet &iJet,TAddJet &iAddJet);

  void setupTree(TTree *iTree);
  void load (int iEvent);
  
  bool selectJets(std::vector<TLorentzVector> &iVetoes);
  void fillVars(TJet *iJet,TLorentzVector *iPtr,TJet &iSaveJet,TAddJet &iASaveJet);
  //Selectors
  bool vetoJet();
  bool passLoose      (TJet *iJet);
  bool passTight      (TJet *iJet);
  bool passVeto       (TJet *iJet);
  bool passPUId       (TJet *iJet);
  TAddJet *addJet     (TJet *iJet);
  //Trigger Stuff
  void addTrigger (std::string iName);
  bool passTrigObj(TJet *iJet,int iId);

protected: 
  TClonesArray *fJets;
  TBranch      *fJetBr;
  TClonesArray *fAddJets;
  TBranch      *fAddJetBr;
  
  TTree        *fTree;
  std::vector<std::string> fTrigString;
  TTrigger     *fTrigger;

  unsigned int fNJets;
  unsigned int fNBTags;
  unsigned int fNBTags10;
  unsigned int fNQTags;
  float        fJDPhi;
  float        fJDEta;

  TLorentzVector *fPtr1;
  TLorentzVector *fPtr2;
  TLorentzVector *fPtr3;
  TLorentzVector *fPtr4;
  TLorentzVector *fDiJet;

  TJet fJet1;
  TJet fJet2;
  TJet fJet3;
  TJet fJet4;

  TAddJet fAJet1;
  TAddJet fAJet2;
  TAddJet fAJet3;
  TAddJet fAJet4;

  unsigned int fHLTMatch;
  unsigned int fNoiseClean;
};
