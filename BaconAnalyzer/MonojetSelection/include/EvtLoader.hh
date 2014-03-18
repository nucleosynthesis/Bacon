#include "BaconAna/DataFormats/interface/TEventInfo.hh"
#include "BaconAna/DataFormats/interface/TTrigger.hh"

#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"

#include <string>
#include <vector>

using namespace baconhep;

class EvtLoader { 
public:
  EvtLoader(TTree *iTree,std::string iHLTFile="/afs/cern.ch/user/p/pharris/pharris/public/bacon/CMSSW_5_3_13/src/BaconAna/DataFormats/data/HLTFile_v0");
  ~EvtLoader(); 
  void reset();
  void resetRecoil(); 
  void setupTree  (TTree *iTree);
  void setupRecoilTree();
  void load (int iEvent);
  void fillEvent();
  void fillRecoil(TLorentzVector &iLep);
  TLorentzVector Met(int iOption);
  TLorentzVector recoil(TLorentzVector &iLep,TLorentzVector &iMet);
  void addTrigger(std::string iName);
  bool passTrigger();
  bool passTrigger(std::string iTrigger);
  float fRho;

  float getEventRho(){return fEvt->rhoIso; } ;

protected: 
  TEventInfo   *fEvt;
  TBranch      *fEvtBr;
  TTree        *fTree;
  TTrigger     *fTrigger;
  
  std::vector<std::string>   fTrigString;
  unsigned int fRun;
  unsigned int fEvtV;
  unsigned int fLumi;
  float fMet;
  float fMetPhi;
  float fTKMet;
  float fTKMetPhi;
  float fMVAMet;
  float fMVAMetPhi;
  float fU1;
  float fU2;
  float fTKU1;
  float fTKU2;
  float fMVAMt;
  float fMt;
  float fTKMt;
  float fTKUDPhi;
  float fUDPhi;
};
