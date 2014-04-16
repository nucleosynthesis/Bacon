#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "BaconAna/DataFormats/interface/TMuon.hh"
#include "BaconAna/DataFormats/interface/TElectron.hh"

using namespace baconhep;

class LeptonLoader { 
public:

  LeptonLoader(TTree *iTree);
  ~LeptonLoader();
  void reset();
  void setupTree(TTree *iTree);
  void load (int iEvent);
  void fillLeptons(std::vector<TMuon*> iMuons,std::vector<TElectron*> iElectrons);
  bool passTight(TMuon *iMuon);

protected: 
  TTree        *fTree;
  int fNLep;
  TLorentzVector *fPtr1;
  TLorentzVector *fPtr2;
  TLorentzVector *fPtr3;
  int fLep1IsTightMuon;
  int fLep2IsTightMuon;
  int fLep3IsTightMuon;

  int   fId1;
  int   fId2;
  int   fId3;
};
