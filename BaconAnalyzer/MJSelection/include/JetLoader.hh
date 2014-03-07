#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "BaconAna/DataFormats/interface/TJet.hh"

using namespace baconhep;

class JetLoader { 
public:
  JetLoader(TTree *iTree);
  ~JetLoader();
  void reset();
  void setupTree(TTree *iTree);
  void load (int iEvent);
  bool selectSingleJet();
  bool vetoJet();
  bool passLoose      (TJet *iJet);
  bool passTight      (TJet *iJet);
  bool passVeto       (TJet *iJet);
protected: 
  TClonesArray *fJets;
  TBranch      *fJetBr;
  TTree        *fTree;
  float         fPt;
  float         fEta;
  float         fPhi;
  float         fM;
};
