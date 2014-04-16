#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "BaconAna/DataFormats/interface/TMuon.hh"

using namespace baconhep;

class MuonLoader { 
public:
  MuonLoader(TTree *iTree);
  ~MuonLoader();
  void resetDiMu();
  void reset();
  void setupTree(TTree *iTree);
  void load (int iEvent);
  //Selecters
  bool selectMuons(std::vector<TLorentzVector>& iVetoes);
  void selectDiMuon(std::vector<TLorentzVector> &iVetoes);
  //Build vetoes
  bool vetoMu();
  void fillVetoes(std::vector<TLorentzVector> &iVec);  
  //Ids
  bool passLoose(TMuon *iMuon);
  bool passTight(TMuon *iMuon);
  bool passWW   (TMuon *iMuon);
  std::vector<TMuon *> fSelMuons;
  float            fMassMin;
  float            fMassMax;

protected: 
  TClonesArray    *fMuons;
  TBranch         *fMuonBr;
  TTree           *fTree;
  int              fNMuons;
  TLorentzVector  *fDiMuon;
};
