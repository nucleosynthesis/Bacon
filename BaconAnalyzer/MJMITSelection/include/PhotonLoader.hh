#include "TTree.h"
#include "TLorentzVector.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "BaconAna/DataFormats/interface/TPhoton.hh"

using namespace baconhep;

class PhotonLoader { 
public:
  PhotonLoader(TTree *iTree);
  ~PhotonLoader();
  void reset();
  void setupTree(TTree *iTree);
  void load (int iEvent);
  bool selectPhotons(std::vector<TLorentzVector> iElectrons,float iRho);
  void fillVetoes(std::vector<TLorentzVector> &iVec);  
  //Vetos
  bool vetoPhoton();
  //Ids
  bool passLoose       (TPhoton *iPhoton);
  bool passCiCPhotonPre(TPhoton *p,float iRho);  
  bool passCiCPFIso    (TPhoton *p,float iRho);
protected: 
  TClonesArray *fPhotons;
  TBranch      *fPhotonBr;
  TTree        *fTree;
  int             fNPhotons;
  TLorentzVector *fPtr1;
};
