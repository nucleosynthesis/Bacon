#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "BaconAna/DataFormats/interface/TElectron.hh"

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

protected: 
  TClonesArray *fElectrons;
  TBranch      *fElectronBr;
  TTree        *fTree;
  float fPt;
  float fEta;
  float fPhi;
};
