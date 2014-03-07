#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "BaconAna/DataFormats/interface/TTau.hh"

using namespace baconhep;

class TauLoader { 
public:
  TauLoader(TTree *iTree);
  ~TauLoader();
  void reset();
  void setupTree(TTree *iTree);
  void load (int iEvent);
  bool selectSingleTau();
  bool vetoTau  ();
  bool passLoose(TTau *iTau);
  bool passTight(TTau *iTau);
  bool passVeto (TTau *iTau);
  bool passAntiEMVA3(int iCat, float raw, TString WP);
protected: 
  TClonesArray *fTaus;
  TBranch      *fTauBr;
  TTree        *fTree;
  float fPt;
  float fEta;
  float fPhi;
};
