#include "BaconAna/DataFormats/interface/TElectron.hh"
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"

using namespace baconhep;

class ElectronLoader { 
public:
  ElectronLoader(TTree *iTree);
  ~ElectronLoader();
  void reset();
  void setupTree(TTree *iTree);
  void load(int iEvent);
  bool selectSingleEle(float iRho);
  bool vetoEle        (float iRho);
  bool passVeto       (const TElectron *iElectron,float iRho);
  bool passLoose      (const TElectron *iElectron,float iRho);
  double getEffArea(double eta);
  TLorentzVector electron();

protected: 
  TClonesArray *fElectrons;
  TBranch      *fElectronBr;
  TTree        *fTree;
  float fPt;
  float fEta;
  float fPhi;
};
