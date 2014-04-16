#include "TTree.h"
#include "TLorentzVector.h"
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
  bool selectElectrons(float iRho,std::vector<TLorentzVector> &iVetoes);
  //Veto
  bool vetoEle        (float iRho);
  void fillVetoes(std::vector<TLorentzVector> &iVec);
  //Electron Ids
  bool   passVeto  (const TElectron *iElectron,float iRho);
  bool   passLoose (const TElectron *iElectron,float iRho);
  bool   passVBTF95(const TElectron *iEle);
  double getEffArea(double eta);
  //Get Electrons with certain properties
  std::vector<TLorentzVector> conversions();
  std::vector<TLorentzVector> nonConversions();
  //
  std::vector<TElectron*> fSelElectrons;

protected: 
  TClonesArray *fElectrons;
  TBranch      *fElectronBr;
  TTree        *fTree;
  int           fNElectrons;
};
