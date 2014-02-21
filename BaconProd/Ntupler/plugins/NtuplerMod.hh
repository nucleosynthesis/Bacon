#include "FWCore/Framework/interface/MakerMacros.h"    // definitions for declaring plug-in modules
#include "FWCore/Framework/interface/Frameworkfwd.h"   // declaration of EDM types
#include "FWCore/Framework/interface/EDAnalyzer.h"     // EDAnalyzer class
#include "FWCore/ParameterSet/interface/ParameterSet.h"// Parameters
#include "FWCore/ParameterSet/interface/ParameterSet.h"// Parameters
#include "FWCore/Utilities/interface/InputTag.h"
#include <string>                                      // string class

// forward class declarations
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
class TFile;
class TH1F;
class TTree;
class TClonesArray;
namespace edm {
  class TriggerResults;
  class TriggerNames;
}
namespace baconhep {
  class TEventInfo;
  class TGenEventInfo;
  class TTrigger;
  class FillerEventInfo;
  class FillerGenInfo;
  class FillerVertex;
  class FillerElectron;
  class FillerMuon;
  class FillerPhoton;
  class FillerTau;
  class FillerJet;
}

//
class NtuplerMod : public edm::EDAnalyzer {
  public:
    explicit NtuplerMod(const edm::ParameterSet &iConfig);
    ~NtuplerMod();

    static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

  private:
    virtual void beginJob();
    virtual void endJob();
    virtual void analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup);

    virtual void beginRun(const edm::Run &iRun, const edm::EventSetup &iSetup);
    virtual void endRun  (const edm::Run &iRun, const edm::EventSetup &iSetup);
    virtual void beginLuminosityBlock(const edm::LuminosityBlock &iLumi, const edm::EventSetup &iSetup);
    virtual void endLuminosityBlock  (const edm::LuminosityBlock &iLumi, const edm::EventSetup &iSetup);

    // specify trigger paths of interest
    void setTriggers();
    
    // initialization from HLT menu; needs to be called on every change in HLT menu
    void initHLT(const edm::TriggerResults&, const edm::TriggerNames&);
    
    //
    void separatePileUp(const edm::Event &iEvent, const reco::Vertex &pv);


    //--------------------------------------------------------------------------------------------------
    //  data members
    //==================================================================================================   
    bool fSkipOnHLTFail;
    bool fIsData;
    bool fUseGen;
    
    // variables to handle triggers
    edm::ParameterSetID fTriggerNamesID;
    edm::InputTag	fHLTTag;
    edm::InputTag       fHLTObjTag;
    std::string         fHLTFile;
    
    std::vector<const reco::PFCandidate*> fPFNoPU;
    std::vector<const reco::PFCandidate*> fPFPU;
   
    float fEleMinPt;
    float fMuonMinPt;
    float fTauMinPt;
    float fPhotonMinPt;
    float fJetMinPt;

    // AOD collection names
    std::string fGenEvtInfoName;
    std::string fGenParName;
    bool        fFillAllGen;
    std::string fPFCandName;
    std::string fPVName;
    std::string fBSName;
    std::string fPUInfoName;
    std::string fPFMETName;
    std::string fMVAMETName;
    std::string fMVAMETUName;
    std::string fRhoIsoName;
    std::string fRhoJetName;
    std::string fEleName;
    std::string fMuonName;
    bool        fApplyMuscle;
    std::string fPhotonName;
    std::string fTauName;
    std::string fJetName;
    std::string fGenJetName;
    std::string fJetFlavorName;
    std::string fJetFlavorPhysName;
    std::string fPruneJetName;
    edm::InputTag fSubJetName;
    std::string fRhoName;
    std::string fCSVbtagName;
    std::string fCSVbtagNameSubJets;
    std::string fJettinessName;
    std::string fQGLikelihood;
    std::string fQGLikelihoodSubJets;
    bool        fComputeFullJetInfo;
    std::string fTrackName;
    std::string fConvName;
    std::string fEBSCName;
    std::string fEESCName;
    std::string fEBRecHitName;
    std::string fEERecHitName;

    // bacon fillers
    baconhep::FillerEventInfo *fFillerEvtInfo;
    baconhep::FillerGenInfo   *fFillerGenInfo;
    baconhep::FillerVertex    *fFillerPV;
    baconhep::FillerElectron  *fFillerEle;
    baconhep::FillerMuon      *fFillerMuon;
    baconhep::FillerPhoton    *fFillerPhoton;
    baconhep::FillerTau       *fFillerTau;
    baconhep::FillerJet       *fFillerJet;
    
    baconhep::TTrigger        *fTrigger;
//    bool fIsActiveEvtInfo;
//    bool fIsActiveGenInfo;
//    bool fIsActivePV;
//    bool fIsActiveEle;
//    bool fIsActiveMuon;
//    bool fIsActivePhoton;
//    bool fIsActiveTau;
//    bool fIsActiveJet;
    
    // Objects and arrays for output file
    std::string              fOutputName;
    TFile                   *fOutputFile;
    TH1F                    *fTotalEvents;
    TTree                   *fEventTree;
    baconhep::TEventInfo    *fEvtInfo;
    baconhep::TGenEventInfo *fGenEvtInfo;
    TClonesArray            *fGenParArr;
    TClonesArray	    *fEleArr;
    TClonesArray	    *fMuonArr;
    TClonesArray	    *fTauArr;
    TClonesArray	    *fJetArr;
    TClonesArray	    *fPhotonArr;
    TClonesArray	    *fPVArr;
    TClonesArray	    *fAddJetArr;
};
