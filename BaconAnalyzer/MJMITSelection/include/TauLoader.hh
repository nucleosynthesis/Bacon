#include "TTree.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "BaconAna/DataFormats/interface/TTau.hh"
#include "BaconAna/DataFormats/interface/TElectron.hh"
#include "BaconAna/DataFormats/interface/TPhoton.hh"

using namespace baconhep;

class TauLoader { 
public:
  TauLoader(TTree *iTree);
  ~TauLoader();
  void reset();
  void setupTree(TTree *iTree);
  void load (int iEvent);
  bool selectTaus(std::vector<TLorentzVector> lVetoes);
  bool vetoTau  ();
  void fillVetoes(std::vector<TLorentzVector> &iVec);  
  bool passLoose(TTau *iTau);
  bool passTight(TTau *iTau);
  bool passVeto (TTau *iTau);
  bool passAntiEMVA3(int iCat, float raw, TString WP);
protected: 
  TClonesArray *fTaus;
  TBranch      *fTauBr;
  TTree        *fTree;
  int   fNTaus;
  TLorentzVector *fPtr1;
  TLorentzVector *fPtr2;
  
  float fPt1;
  float fEta1;
  float fPhi1;
  float fM1;

  float fPt2;
  float fEta2;
  float fPhi2;
  float fM2;
};
