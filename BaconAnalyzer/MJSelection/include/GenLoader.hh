#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"

using namespace baconhep;

class GenLoader { 
public:
  GenLoader(TTree *iTree);
  ~GenLoader();
  void reset();
  void setupTree(TTree *iTree);
  void load (int iEvent);
  void selectBoson();
  TGenParticle* getStatus1(int iId);
  bool isNeutrino(TGenParticle *iPart);

protected: 
  TClonesArray *fGens;
  TBranch      *fGenBr;
  TTree        *fTree;
  float fVPt;
  float fVEta;
  float fVPhi;
  float fVM;
  int   fVId;
  
  float fPt1;
  float fEta1;
  float fPhi1;
  float fM1;
  int   fId1;

  float fPt2;
  float fEta2;
  float fPhi2;
  float fM2;
  int   fId2;
};
