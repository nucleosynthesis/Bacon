#include "FWCore/Framework/interface/MakerMacros.h"    // definitions for declaring plug-in modules
#include "FWCore/Framework/interface/Frameworkfwd.h"   // declaration of EDM types
#include "FWCore/Framework/interface/EDAnalyzer.h"     // EDAnalyzer class
#include <string>                                      // string class

// forward class declarations
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
class TFile;
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
  class FillerPF;
}

//
class ExpertMod : public edm::EDAnalyzer {
  public:
    explicit ExpertMod(const edm::ParameterSet &iConfig);
    ~ExpertMod();

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

    // AOD collection names
    std::string fGenEvtInfoName;
    std::string fGenParName;
    std::string fPFCandName;
    std::string fPVName;
    std::string fBSName;
    std::string fPUInfoName;
    std::string fRhoIsoName;
    std::string fRhoJetName;
    
    // bacon fillers
    baconhep::FillerEventInfo *fFillerEvtInfo;
    baconhep::FillerGenInfo   *fFillerGenInfo;
    baconhep::FillerVertex    *fFillerPV;
    baconhep::FillerPF        *fFillerPF;
 
    baconhep::TTrigger        *fTrigger;
//    bool fIsActiveEvtInfo;
//    bool fIsActiveGenInfo;
//    bool fIsActivePV;
    
    // Objects and arrays for output file
    std::string              fOutputName;
    TFile                   *fOutputFile;
    TTree                   *fEventTree;
    baconhep::TEventInfo    *fEvtInfo;
    baconhep::TGenEventInfo *fGenEvtInfo;
    TClonesArray            *fGenParArr;
    TClonesArray            *fPFParArr;
    TClonesArray	    *fPVArr;
};
