#include "ExpertMod.hh"

// bacon classes and constants
#include "BaconAna/DataFormats/interface/BaconAnaDefs.hh"
#include "BaconAna/DataFormats/interface/TTrigger.hh"
#include "BaconAna/DataFormats/interface/TEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"
#include "BaconAna/DataFormats/interface/TPFPart.hh"

#include "BaconProd/Ntupler/interface/FillerEventInfo.hh"
#include "BaconProd/Ntupler/interface/FillerGenInfo.hh"
#include "BaconProd/Ntupler/interface/FillerVertex.hh"
#include "BaconProd/Ntupler/interface/FillerPF.hh"

// tools to parse HLT name patterns
#include <boost/foreach.hpp>
#include "FWCore/Utilities/interface/RegexMatch.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

// data format classes
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

// ROOT classes
#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TMath.h>


//--------------------------------------------------------------------------------------------------
ExpertMod::ExpertMod(const edm::ParameterSet &iConfig):
  fSkipOnHLTFail  (iConfig.getUntrackedParameter<Bool_t>("skipOnHLTFail",kFALSE)),
  fIsData         (iConfig.getUntrackedParameter<Bool_t>("isData",kFALSE)),
  fUseGen         (iConfig.getUntrackedParameter<Bool_t>("useGen",kFALSE)),
  fHLTTag         ("TriggerResults","","HLT"),
  fHLTObjTag      ("hltTriggerSummaryAOD","","HLT"),
  fHLTFile        (iConfig.getUntrackedParameter<std::string>("TriggerFile","HLT")),
  fGenEvtInfoName (iConfig.getUntrackedParameter<std::string>("genEventInfoName", "generator")),
  fGenParName     (iConfig.getUntrackedParameter<std::string>("genParticlesName", "genParticles")),
  fPFCandName     (iConfig.getUntrackedParameter<std::string>("pflowCandidatesName", "particleFlow")),
  fPVName         (iConfig.getUntrackedParameter<std::string>("primaryVerticesName", "offlinePrimaryVertices")),
  fBSName         (iConfig.getUntrackedParameter<std::string>("beamspotName", "offlineBeamSpot")),
  fPUInfoName     (iConfig.getUntrackedParameter<std::string>("pileupInfoName", "addPileupInfo")),
  fFillerEvtInfo  (0),
  fFillerGenInfo  (0),
  fFillerPV       (0),
  fOutputName     (iConfig.getUntrackedParameter<std::string>("outputName", "ntuple.root")),
  fOutputFile     (0),
  fEventTree      (0),
  fEvtInfo        (0),
  fGenEvtInfo     (0),
  fGenParArr      (0),
  fPFParArr       (0),
  fPVArr          (0)
{
  // Don't write TObject part of the objects
  baconhep::TEventInfo::Class()->IgnoreTObjectStreamer();
  baconhep::TGenEventInfo::Class()->IgnoreTObjectStreamer();
  baconhep::TGenParticle::Class()->IgnoreTObjectStreamer();
  baconhep::TPFPart::Class()->IgnoreTObjectStreamer();
}

//--------------------------------------------------------------------------------------------------
ExpertMod::~ExpertMod()
{}

//--------------------------------------------------------------------------------------------------
void ExpertMod::beginJob()
{
  fEvtInfo    = new baconhep::TEventInfo();
  fGenEvtInfo = new baconhep::TGenEventInfo();
  
  //
  // Set up arrays
  //
  fGenParArr = new TClonesArray("baconhep::TGenParticle",5000); assert(fGenParArr);
  fPFParArr  = new TClonesArray("baconhep::TPFPart"     ,5000); assert(fPFParArr);
  fPVArr     = new TClonesArray("baconhep::TVertex"     ,5000); assert(fPVArr);

  //
  // Create output file, trees, and structs
  //
  fOutputFile = new TFile(fOutputName.c_str(), "RECREATE");
  fEventTree = new TTree("Events","Events");
  
  fEventTree->Branch("Info",fEvtInfo);
  if(fUseGen) {
    fEventTree->Branch("GenEvtInfo",fGenEvtInfo);
    fEventTree->Branch("GenParticle",&fGenParArr);
  }
  fEventTree->Branch("PFPart",   &fPFParArr);
  fEventTree->Branch("PV",       &fPVArr);

  //
  // Triggers
  //
  setTriggers();

  //
  // Fillers
  //
  fFillerEvtInfo = new baconhep::FillerEventInfo();
  fFillerEvtInfo->fPFCandName  = fPFCandName;
  fFillerEvtInfo->fPUInfoName  = fPUInfoName;
  fFillerEvtInfo->fBSName      = fBSName;
  fFillerEvtInfo->fFillMET     = false;
  fFillerEvtInfo->fPFMETName   = "";
  fFillerEvtInfo->fMVAMETName  = "";
  fFillerEvtInfo->fMVAMETUName = "";
  //fFillerEvtInfo->fRhoIsoName  = fRhoIsoName;
  //fFillerEvtInfo->fRhoJetName  = fRhoJetName;
  

  fFillerGenInfo = new baconhep::FillerGenInfo();
  fFillerGenInfo->fGenEvtInfoName = fGenEvtInfoName;
  fFillerGenInfo->fGenParName     = fGenParName;


  fFillerPV = new baconhep::FillerVertex();
  fFillerPV->fPVName	    = fPVName;
  fFillerPV->fMinNTracksFit = 0;
  fFillerPV->fMinNdof	    = 4;
  fFillerPV->fMaxAbsZ	    = 24;
  fFillerPV->fMaxRho	    = 2;

  fFillerPF  = new baconhep::FillerPF();
}

//--------------------------------------------------------------------------------------------------
void ExpertMod::endJob() 
{
  //
  // Save to ROOT file
  //
  fEventTree->Print();
  fOutputFile->Write();
  fOutputFile->Close();
  
  delete fFillerEvtInfo;
  delete fFillerGenInfo;
  delete fFillerPV;
  delete fFillerPF;

  delete fEvtInfo;
  delete fGenEvtInfo;
  delete fGenParArr;
  delete fPFParArr;
  delete fPVArr;
}

//--------------------------------------------------------------------------------------------------
void ExpertMod::setTriggers()
{
  fTrigger = new baconhep::TTrigger(fHLTFile);
}

//--------------------------------------------------------------------------------------------------
void ExpertMod::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<edm::TriggerResults> hTrgRes;
  iEvent.getByLabel(fHLTTag,hTrgRes);
  assert(hTrgRes.isValid());  
  const edm::TriggerNames &triggerNames = iEvent.triggerNames(*hTrgRes);
  Bool_t config_changed = false;
  if(fTriggerNamesID != triggerNames.parameterSetID()) {
    fTriggerNamesID = triggerNames.parameterSetID();
    config_changed  = true;
  }
  if(config_changed) {
    initHLT(*hTrgRes, triggerNames);
  }
  TriggerBits triggerBits;
  for(unsigned int irec=0; irec<fTrigger->fRecords.size(); irec++) {
    if(fTrigger->fRecords[irec].hltPathIndex == (unsigned int)-1) continue;
    if(hTrgRes->accept(fTrigger->fRecords[irec].hltPathIndex)) {
      triggerBits [fTrigger->fRecords[irec].baconTrigBit] = 1;
    }
  }
  if(fSkipOnHLTFail && triggerBits == 0) return;  

  if(fUseGen) {
    fGenParArr->Clear();
    fFillerGenInfo->fill(fGenEvtInfo, fGenParArr, iEvent);
  }
  fPVArr->Clear();
  int nvertices = 0;
  const reco::Vertex *pv = fFillerPV->fill(fPVArr, nvertices, iEvent);
  assert(pv);
  
  separatePileUp(iEvent, *pv);
  fFillerEvtInfo->fill(fEvtInfo, iEvent, *pv, (nvertices>0), triggerBits);
  fPFParArr->Clear();
  fFillerPF->fill(fPFParArr,fPVArr,iEvent);
  fEventTree->Fill();
}

//--------------------------------------------------------------------------------------------------
void ExpertMod::initHLT(const edm::TriggerResults& result, const edm::TriggerNames& triggerNames)
{
  for(unsigned int irec=0; irec<fTrigger->fRecords.size(); irec++) {
    fTrigger->fRecords[irec].hltPathName  = "";
    fTrigger->fRecords[irec].hltPathIndex = (unsigned int)-1;
    const std::string pattern = fTrigger->fRecords[irec].hltPattern;
    if(edm::is_glob(pattern)) {  // handle pattern with wildcards (*,?)
      std::vector<std::vector<std::string>::const_iterator> matches = edm::regexMatch(triggerNames.triggerNames(), pattern);
      if(matches.empty()) {
        std::cout << "requested pattern [" << pattern << "] does not match any HLT paths" << std::endl;
      } else {
        BOOST_FOREACH(std::vector<std::string>::const_iterator match, matches) {
          fTrigger->fRecords[irec].hltPathName = *match;
       }
      }
    } else {  // take full HLT path name given
      fTrigger->fRecords[irec].hltPathName = pattern;
    }
    // Retrieve index in trigger menu corresponding to HLT path
    unsigned int index = triggerNames.triggerIndex(fTrigger->fRecords[irec].hltPathName);
    if(index < result.size()) {  // check for valid index
      fTrigger->fRecords[irec].hltPathIndex = index;
    }
  }
}

//--------------------------------------------------------------------------------------------------
void ExpertMod::separatePileUp(const edm::Event &iEvent, const reco::Vertex &pv)
{
  // recipe from Matthew Chan

  // Get PF-candidates collection
  edm::Handle<reco::PFCandidateCollection> hPFCandProduct;
  iEvent.getByLabel(fPFCandName,hPFCandProduct);
  assert(hPFCandProduct.isValid());
  const reco::PFCandidateCollection *pfCandCol = hPFCandProduct.product();  
  
  // Get vertex collection
  edm::Handle<reco::VertexCollection> hVertexProduct;
  iEvent.getByLabel(fPVName,hVertexProduct);
  assert(hVertexProduct.isValid());
  const reco::VertexCollection *pvCol = hVertexProduct.product();
  
  fPFNoPU.clear();
  fPFPU.clear();
  
  for(reco::PFCandidateCollection::const_iterator iP = pfCandCol->begin(); iP!=pfCandCol->end(); ++iP) {
    if(iP->particleId() == reco::PFCandidate::h) {  // charged hadrons
      if(iP->trackRef().isNonnull() && pv.trackWeight(iP->trackRef())>0) {
        // charged hadrons with track used to compute PV
	fPFNoPU.push_back(&(*iP)); 
      
      } else {
        // Find closest vertex to charged hadron's vertex source
	bool vertexFound = false;
	const reco::Vertex *closestVtx = 0;
	double dzmin = 10000;
	
	for(reco::VertexCollection::const_iterator iV = pvCol->begin(); iV!=pvCol->end(); ++iV) {
	  if(iP->trackRef().isNonnull() && iV->trackWeight(iP->trackRef())>0) {
	    vertexFound = true;
	    closestVtx  = &(*iV);
	    break;
	  }
	  
	  double dz = fabs(iP->vertex().z() - iV->z());
	  if(dz < dzmin) {
	    closestVtx = &(*iV);
	    dzmin      = dz;
	  }
	}
	
	if(vertexFound || closestVtx != &pv) {
	  fPFPU.push_back(&(*iP));
	} else {
	  fPFNoPU.push_back(&(*iP));  // Note: when no associated vertex found, assume to come from PV
	}
      }
      
    } else {  // all non-charged-hadron PFCandidates are considered to be from PV
      fPFNoPU.push_back(&(*iP));
    }
  }
}

//--------------------------------------------------------------------------------------------------
void ExpertMod::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
void ExpertMod::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup){}
void ExpertMod::endRun  (const edm::Run& iRun, const edm::EventSetup& iSetup){}
void ExpertMod::beginLuminosityBlock(const edm::LuminosityBlock& iLumi, const edm::EventSetup& iSetup){}
void ExpertMod::endLuminosityBlock  (const edm::LuminosityBlock& iLumi, const edm::EventSetup& iSetup){}


//define this as a plug-in
DEFINE_FWK_MODULE(ExpertMod);
