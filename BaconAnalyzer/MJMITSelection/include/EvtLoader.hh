#include "BaconAna/DataFormats/interface/TEventInfo.hh"
#include "BaconAna/DataFormats/interface/TVertex.hh"
#include "BaconAna/DataFormats/interface/TTrigger.hh"

#include "TH1F.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"

#include <string>
#include <vector>

using namespace baconhep;

class EvtLoader { 
public:
  EvtLoader(TTree *iTree,std::string iHLTFile="/afs/cern.ch/user/p/pharris/pharris/public/bacon/CMSSW_5_3_13/src/BaconAna/DataFormats/data/HLTFile_v0",std::string iPUWeight="/afs/cern.ch/user/p/pharris/pharris/public/bacon/CMSSW_5_3_13/src/BaconAnalyzer/MJMITSelection/data/PUWeights_2012.root");
  ~EvtLoader(); 
  void reset();
  void setupTree  (TTree *iTree,float iWeight);
  void load (int iEvent);
  //Fillers
  void fillEvent(std::vector<TLorentzVector> &iVecCorr);
  bool passSkim();
  TLorentzVector Met(int iOption);
  //Trigger Related Stuff
  void addTrigger(std::string iName);
  bool passFilter();
  bool passTrigger();
  bool passTrigger(std::string iTrigger);
  unsigned int triggerBit();
  //PU relates stuff
  float        puWeight(double iPU);
  unsigned int nVtx();
  void  correctMet(float &iMet,float &iMetPhi,TLorentzVector &iCorr);
  //Met Stuff
  float        metSig(float iMet,float iMetPhi,float iCov00,float iCov01,float iCov10,float iCov11);
  unsigned int metFilter(unsigned int iMetFilter);
  float fRho;
  unsigned int fRun;
  unsigned int fEvtV;
  unsigned int fLumi;

protected: 
  TEventInfo   *fEvt;
  TBranch      *fEvtBr;

  TClonesArray *fVertices;
  TBranch      *fVertexBr;

  TTree        *fTree;
  TTrigger     *fTrigger;
  
  std::vector<std::string>   fTrigString;
  unsigned int fITrigger;
  unsigned int fHLTMatch;
  unsigned int fMetFilters;
  unsigned int fNVtx;
  unsigned int fNPU;
  unsigned int fNPUP;
  unsigned int fNPUM;
  float        fPUWeight;
  float        fScale;

  float fMet;
  float fMetPhi;
  float fMetCov00;
  float fMetCov01;
  float fMetCov10;
  float fMetCov11;
  float fMetSig;
  
  float fTKMet;
  float fTKMetPhi;

  float fRawSumEt;
  float fMetRaw;
  float fMetRawPhi;
  float fMVAMet;
  float fMVAMetPhi;
  float fMVAMetUnity;
  float fMVAMetUnityPhi;
  
  float fMVAMetCov00;
  float fMVAMetCov01;
  float fMVAMetCov10;
  float fMVAMetCov11;
  float fMVAMetSig;
  float fMVAMetUnityCov00;
  float fMVAMetUnityCov01;
  float fMVAMetUnityCov10;
  float fMVAMetUnityCov11;
  float fMVAMetUnitySig;

  TH1F *fPUWeightHist;
};
