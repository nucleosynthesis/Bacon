#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "BaconAna/DataFormats/interface/TJet.hh"
#include "TLorentzVector.h"

using namespace baconhep;

class JetLoader { 
public:
  JetLoader(TTree *iTree);
  ~JetLoader();
  void reset();
  void setupTree(TTree *iTree);
  void load (int iEvent);
  bool selectSingleJet();
  bool vetoAdditionalJets(int);
  bool selectJetsWithAdditionalJetsVeto();
  bool vetoJet();
  bool passLoose      (TJet *iJet);
  bool passTight      (TJet *iJet);
  bool passVeto       (TJet *iJet);
  TLorentzVector      jet();
protected: 
  TClonesArray *fJets;
  TBranch      *fJetBr;
  TTree        *fTree;
  float         fPt;
  float         fEta;
  float         fPhi;
  float         fM;
};
