#include "NtuplerMod.hh"

// bacon classes and constants
#include "BaconAna/DataFormats/interface/BaconAnaDefs.hh"
#include "BaconAna/DataFormats/interface/TEventInfo.hh"
#include "BaconAna/DataFormats/interface/TTrigger.hh"
#include "BaconAna/DataFormats/interface/TGenEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"
#include "BaconAna/DataFormats/interface/TElectron.hh"
#include "BaconAna/DataFormats/interface/TMuon.hh"
#include "BaconAna/DataFormats/interface/TTau.hh"
#include "BaconAna/DataFormats/interface/TJet.hh"
#include "BaconAna/DataFormats/interface/TPhoton.hh"
#include "BaconAna/DataFormats/interface/TVertex.hh"
#include "BaconAna/DataFormats/interface/TAddJet.hh"

#include "BaconProd/Ntupler/interface/FillerEventInfo.hh"
#include "BaconProd/Ntupler/interface/FillerGenInfo.hh"
#include "BaconProd/Ntupler/interface/FillerVertex.hh"
#include "BaconProd/Ntupler/interface/FillerMuon.hh"
#include "BaconProd/Ntupler/interface/FillerElectron.hh"
#include "BaconProd/Ntupler/interface/FillerPhoton.hh"
#include "BaconProd/Ntupler/interface/FillerTau.hh"
#include "BaconProd/Ntupler/interface/FillerJet.hh"

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
#include <TH1F.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TMath.h>


//--------------------------------------------------------------------------------------------------
NtuplerMod::NtuplerMod(const edm::ParameterSet &iConfig):
  fSkipOnHLTFail  (iConfig.getUntrackedParameter<Bool_t>("skipOnHLTFail",kFALSE)),
  fIsData         (iConfig.getUntrackedParameter<Bool_t>("isData",kFALSE)),
  fUseGen         (iConfig.getUntrackedParameter<Bool_t>("useGen",kFALSE)),
  fHLTTag         ("TriggerResults","","HLT"),
  fHLTObjTag      ("hltTriggerSummaryAOD","","HLT"),
  fHLTFile        (iConfig.getUntrackedParameter<std::string>("TriggerFile","HLT")),
  fEleMinPt       (iConfig.getUntrackedParameter<double>("electronMinPt",5)),
  fMuonMinPt      (iConfig.getUntrackedParameter<double>("muonMinPt",0)),
  fTauMinPt       (iConfig.getUntrackedParameter<double>("tauMinPt",20)),
  fPhotonMinPt    (iConfig.getUntrackedParameter<double>("photonMinPt",2)),
  fJetMinPt       (iConfig.getUntrackedParameter<double>("jetMinPt",20)),
  fGenEvtInfoName (iConfig.getUntrackedParameter<std::string>("genEventInfoName", "generator")),
  fGenParName     (iConfig.getUntrackedParameter<std::string>("genParticlesName", "genParticles")),
  fFillAllGen     (iConfig.getUntrackedParameter<bool>("fillAllGen", false)),
  fPFCandName     (iConfig.getUntrackedParameter<std::string>("pflowCandidatesName", "particleFlow")),
  fPVName         (iConfig.getUntrackedParameter<std::string>("primaryVerticesName", "offlinePrimaryVertices")),
  fBSName         (iConfig.getUntrackedParameter<std::string>("beamspotName", "offlineBeamSpot")),
  fPUInfoName     (iConfig.getUntrackedParameter<std::string>("pileupInfoName", "addPileupInfo")),
  fPFMETName      (iConfig.getUntrackedParameter<std::string>("pflowMETName", "pfMet")),
  fMVAMETName     (iConfig.getUntrackedParameter<std::string>("mvaMETName", "pfMEtMVA")),
  fMVAMETUName    (iConfig.getUntrackedParameter<std::string>("mvaMETUnityName", "pfMEtMVAUnity")),
  fRhoIsoName     (iConfig.getUntrackedParameter<std::string>("rhoForIsolationName", "kt6PFJets")),
  fRhoJetName     (iConfig.getUntrackedParameter<std::string>("rhoForJetsName", "kt6PFJets")),
  fEleName        (iConfig.getUntrackedParameter<std::string>("electronName", "gsfElectrons")),
  fMuonName       (iConfig.getUntrackedParameter<std::string>("muonName", "muons")),
  fPhotonName     (iConfig.getUntrackedParameter<std::string>("photonName", "photons")),
  fTauName        (iConfig.getUntrackedParameter<std::string>("tauName", "hpsPFTauProducer")),
  fJetName        (iConfig.getUntrackedParameter<std::string>("jetName", "ak5PFJets")),
  fGenJetName     (iConfig.getUntrackedParameter<std::string>("genJetName"   , "ak5GenJets")),
  fJetFlavorName  (iConfig.getUntrackedParameter<std::string>("jetFlavorName", "jetCombinedSecondaryVertexBJetTagsSJ")),
  fJetFlavorPhysName(iConfig.getUntrackedParameter<std::string>("jetFlavorPhysName", "jetCombinedSecondaryVertexBJetTagsSJ")),
  fPruneJetName   (iConfig.getUntrackedParameter<std::string>("pruneJetName" , "ca5PFJetsPruned")),
  fSubJetName     (iConfig.getParameter<edm::InputTag>("subJetName")),
  fRhoName        (iConfig.getUntrackedParameter<std::string>("rhoName"      , "kt6PFJets")),
  fCSVbtagName    (iConfig.getUntrackedParameter<std::string>("csvBTagName"  , "combinedSecondaryVertexBJetTags")),
  fCSVbtagNameSubJets(iConfig.getUntrackedParameter<std::string>("csvBTagSubJetName", "comibnedSecondaryVertexBJetTagsSJ")),
  fJettinessName  (iConfig.getUntrackedParameter<std::string>("jettiness"            , "NJettiness")),
  fQGLikelihood   (iConfig.getUntrackedParameter<std::string>("QGLikelihood"         , "QGTagger")),
  fQGLikelihoodSubJets(iConfig.getUntrackedParameter<std::string>("QGLikelihoodSubjet", "QGTaggerSubJets")),
  fComputeFullJetInfo(iConfig.getUntrackedParameter<bool>("computeFullJetInfo", false)),
  fTrackName      (iConfig.getUntrackedParameter<std::string>("trackName", "generalTracks")),
  fConvName       (iConfig.getUntrackedParameter<std::string>("conversionName", "allConversions")),
  fEBSCName       (iConfig.getUntrackedParameter<std::string>("ecalBarrelSuperclusterName", "correctedHybridSuperClusters")),
  fEESCName       (iConfig.getUntrackedParameter<std::string>("ecalEndcapSuperclusterName", "correctedMulti5x5SuperClustersWithPreshower")),
  fEBRecHitName   (iConfig.getUntrackedParameter<std::string>("ecalBarrelRecHitName", "reducedEcalRecHitsEB")),
  fEERecHitName   (iConfig.getUntrackedParameter<std::string>("ecalEndcapRecHitName", "reducedEcalRecHitsEE")),
  fFillerEvtInfo  (0),
  fFillerGenInfo  (0),
  fFillerPV       (0),
  fFillerEle      (0),
  fFillerMuon     (0),
  fFillerPhoton   (0),
  fFillerTau      (0),
  fFillerJet      (0),
//  fIsActiveEvtInfo(iConfig.getUntrackedParameter<bool>("isActiveEventInfo", true)),
//  fIsActiveGenInfo(iConfig.getUntrackedParameter<bool>("isActiveGenInfo", true)),
//  fIsActivePV     (iConfig.getUntrackedParameter<bool>("isActivePV", true)),
//  fIsActiveEle    (iConfig.getUntrackedParameter<bool>("isActiveElectron", true)),
//  fIsActiveMuon   (iConfig.getUntrackedParameter<bool>("isActiveMuon", true)),
//  fIsActivePhoton (iConfig.getUntrackedParameter<bool>("isActivePhoton", true)),
//  fIsActiveTau    (iConfig.getUntrackedParameter<bool>("isActiveTau", true)),
//  fIsActiveJet    (iConfig.getUntrackedParameter<bool>("isActiveJet", true)),
  fOutputName     (iConfig.getUntrackedParameter<std::string>("outputName", "ntuple.root")),
  fOutputFile     (0),
  fEventTree      (0),
  fEvtInfo        (0),
  fGenEvtInfo     (0),
  fGenParArr      (0),
  fEleArr         (0),
  fMuonArr        (0),
  fTauArr         (0),
  fJetArr         (0),
  fPhotonArr      (0),
  fPVArr          (0),
  fAddJetArr      (0)
{
  // Don't write TObject part of the objects
  baconhep::TEventInfo::Class()->IgnoreTObjectStreamer();
  baconhep::TGenEventInfo::Class()->IgnoreTObjectStreamer();
  baconhep::TGenParticle::Class()->IgnoreTObjectStreamer();
  baconhep::TMuon::Class()->IgnoreTObjectStreamer();
  baconhep::TElectron::Class()->IgnoreTObjectStreamer();
  baconhep::TTau::Class()->IgnoreTObjectStreamer();
  baconhep::TJet::Class()->IgnoreTObjectStreamer();
  baconhep::TPhoton::Class()->IgnoreTObjectStreamer();
  baconhep::TVertex::Class()->IgnoreTObjectStreamer();
  baconhep::TAddJet::Class()->IgnoreTObjectStreamer();
}

//--------------------------------------------------------------------------------------------------
NtuplerMod::~NtuplerMod()
{}

//--------------------------------------------------------------------------------------------------
void NtuplerMod::beginJob()
{
  fEvtInfo    = new baconhep::TEventInfo();
  fGenEvtInfo = new baconhep::TGenEventInfo();
  
  //
  // Set up arrays
  //
  fGenParArr = new TClonesArray("baconhep::TGenParticle",5000); assert(fGenParArr);
  fEleArr    = new TClonesArray("baconhep::TElectron");         assert(fEleArr);
  fMuonArr   = new TClonesArray("baconhep::TMuon");	        assert(fMuonArr);
  fTauArr    = new TClonesArray("baconhep::TTau");              assert(fTauArr);
  fJetArr    = new TClonesArray("baconhep::TJet");	        assert(fJetArr);
  fPhotonArr = new TClonesArray("baconhep::TPhoton");	        assert(fPhotonArr);
  fPVArr     = new TClonesArray("baconhep::TVertex");           assert(fPVArr);
  fAddJetArr = new TClonesArray("baconhep::TAddJet");	        assert(fAddJetArr);

  //
  // Create output file, trees, and structs
  //
  fOutputFile  = new TFile(fOutputName.c_str(), "RECREATE");
  fTotalEvents = new TH1F("TotalEvents","TotalEvents",1,-10,10);
  fEventTree   = new TTree("Events","Events");
  
  fEventTree->Branch("Info",fEvtInfo);
  if(fUseGen) {
    fEventTree->Branch("GenEvtInfo",fGenEvtInfo);
    fEventTree->Branch("GenParticle",&fGenParArr);
  }
  fEventTree->Branch("Electron", &fEleArr);
  fEventTree->Branch("Muon",     &fMuonArr);
  fEventTree->Branch("Tau",      &fTauArr);
  fEventTree->Branch("Jet",      &fJetArr);
  fEventTree->Branch("Photon",   &fPhotonArr);
  fEventTree->Branch("PV",       &fPVArr);
  if(fComputeFullJetInfo) fEventTree->Branch("JetAdd",   &fAddJetArr);

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
  fFillerEvtInfo->fPFMETName   = fPFMETName;
  fFillerEvtInfo->fMVAMETName  = fMVAMETName;
  fFillerEvtInfo->fMVAMETUName = fMVAMETUName;
  fFillerEvtInfo->fRhoIsoName  = fRhoIsoName;
  fFillerEvtInfo->fRhoJetName  = fRhoJetName;
  

  fFillerGenInfo = new baconhep::FillerGenInfo();
  fFillerGenInfo->fGenEvtInfoName = fGenEvtInfoName;
  fFillerGenInfo->fGenParName     = fGenParName;
  fFillerGenInfo->fFillAll        = fFillAllGen;

  fFillerPV = new baconhep::FillerVertex();
  fFillerPV->fPVName	    = fPVName;
  fFillerPV->fMinNTracksFit = 0;
  fFillerPV->fMinNdof	    = 4;
  fFillerPV->fMaxAbsZ	    = 24;
  fFillerPV->fMaxRho	    = 2;
  
  fFillerEle = new baconhep::FillerElectron();
  fFillerEle->fMinPt	    = fEleMinPt;
  fFillerEle->fEleName      = fEleName;
  fFillerEle->fPFCandName   = fPFCandName;
  fFillerEle->fTrackName    = fTrackName;
  fFillerEle->fConvName     = fConvName;
  fFillerEle->fRhoName      = fRhoIsoName;
  fFillerEle->fEBSCName     = fEBSCName;
  fFillerEle->fEESCName     = fEESCName;
  fFillerEle->fEBRecHitName = fEBRecHitName;
  fFillerEle->fEERecHitName = fEERecHitName;
  
  std::string cmsenv = getenv("CMSSW_BASE"); cmsenv+="/src";
  std::vector<std::string> eleIDFilenames;
  eleIDFilenames.push_back(cmsenv+"/BaconProd/Utils/data/ElectronID_BDTG_EGamma2012NonTrigV0_Cat1.weights.xml");
  eleIDFilenames.push_back(cmsenv+"/BaconProd/Utils/data/ElectronID_BDTG_EGamma2012NonTrigV0_Cat2.weights.xml");
  eleIDFilenames.push_back(cmsenv+"/BaconProd/Utils/data/ElectronID_BDTG_EGamma2012NonTrigV0_Cat3.weights.xml");
  eleIDFilenames.push_back(cmsenv+"/BaconProd/Utils/data/ElectronID_BDTG_EGamma2012NonTrigV0_Cat4.weights.xml");
  eleIDFilenames.push_back(cmsenv+"/BaconProd/Utils/data/ElectronID_BDTG_EGamma2012NonTrigV0_Cat5.weights.xml");
  eleIDFilenames.push_back(cmsenv+"/BaconProd/Utils/data/ElectronID_BDTG_EGamma2012NonTrigV0_Cat6.weights.xml");
  fFillerEle->fEleIDMVA.initialize("BDT", EGammaMvaEleEstimator::kNonTrig, true, eleIDFilenames);
  fFillerEle->fEleCorr.initialize((cmsenv+"/BaconProd/Utils/data/eleEnergyRegWeights_WithSubClusters_VApr15.root").c_str(),
  				  baconhep::ElectronEnergyRegression::kWithSubCluVar,
        			  baconhep::ElectronEnergySmearingScaling::kSummer12_LegacyPaper,
        			  2,
        			  (cmsenv+"/BaconProd/Utils/data/scalesCorr.csv").c_str(),
				  (cmsenv+"/BaconProd/Utils/data/smearsCorrType1.csv").c_str(),
        			  (cmsenv+"/BaconProd/Utils/data/smearsCorrType2.csv").c_str(),
				  (cmsenv+"/BaconProd/Utils/data/smearsCorrType3.csv").c_str(),
				  (cmsenv+"/BaconProd/Utils/data/linearityNewReg-May2013.csv").c_str(),
        			  false);

  fFillerMuon = new baconhep::FillerMuon();
  fFillerMuon->fMinPt	   = fMuonMinPt;
  fFillerMuon->fMuonName   = fMuonName;
  fFillerMuon->fPFCandName = fPFCandName;
  fFillerMuon->fTrackName  = fTrackName;
  fFillerMuon->fSaveTracks = false;
  fFillerMuon->fTrackMinPt = 20;
  fFillerMuon->fMuCorr.initialize(baconhep::MuonMomentumCorrector::kMuScleSummer12_DR53X_smearReReco,
  				  (cmsenv+"/MuScleFit/Calibration/data").c_str(),
        			  false);
  

  fFillerPhoton = new baconhep::FillerPhoton();
  fFillerPhoton->fMinPt        = fPhotonMinPt;
  fFillerPhoton->fPhotonName   = fPhotonName;
  fFillerPhoton->fPFCandName   = fPFCandName;
  fFillerPhoton->fEleName      = fEleName;
  fFillerPhoton->fConvName     = fConvName;
  fFillerPhoton->fEBSCName     = fEBSCName;
  fFillerPhoton->fEESCName     = fEESCName;
  fFillerPhoton->fEBRecHitName = fEBRecHitName;
  fFillerPhoton->fEERecHitName = fEERecHitName;
  
  fFillerTau = new baconhep::FillerTau();
  fFillerTau->fMinPt   = fTauMinPt;
  fFillerTau->fTauName = fTauName;
  
//  fFillerTau->fRingIso.initialize();
//  fFillerTau->fRingIso2.initialize();  
  
  fFillerJet = new baconhep::FillerJet();
  fFillerJet->fMinPt	         = fJetMinPt;
  fFillerJet->fUseGen            = fUseGen;
  fFillerJet->fJetName           = fJetName;
  fFillerJet->fGenJetName        = fGenJetName;
  fFillerJet->fJetFlavorName     = fJetFlavorName;
  fFillerJet->fJetFlavorPhysName = fJetFlavorPhysName;
  fFillerJet->fPruneJetName      = fPruneJetName;
  fFillerJet->fSubJetName        = fSubJetName;//(fSubJetName.label()+"_"+fSubJetName.instance());
  fFillerJet->fPVName	         = fPVName;
  fFillerJet->fRhoName           = fRhoJetName;
  fFillerJet->fCSVbtagName       = fCSVbtagName;
  fFillerJet->fCSVbtagSubJetName = fCSVbtagNameSubJets;
  fFillerJet->fJettinessName     = fJettinessName;
  fFillerJet->fQGLikelihood      = fQGLikelihood;
  fFillerJet->fQGLikelihoodSubJets  = fQGLikelihoodSubJets;
  fFillerJet->fComputeFullInfo   = fComputeFullJetInfo;

  std::vector<std::string> jecFiles;
  std::vector<std::string> jecUncFiles;
  jecFiles.push_back(cmsenv+"/BaconProd/Utils/data/Summer13_V1_MC_L1FastJet_AK5PF.txt");
  jecFiles.push_back(cmsenv+"/BaconProd/Utils/data/Summer13_V1_MC_L2Relative_AK5PF.txt");
  jecFiles.push_back(cmsenv+"/BaconProd/Utils/data/Summer13_V1_MC_L3Absolute_AK5PF.txt");
  if(fIsData)   jecFiles.push_back(cmsenv+"/BaconProd/Utils/data/Summer13_V1_DATA_L2L3Residual_AK5PF.txt");
  if(!fIsData)  jecUncFiles.push_back(cmsenv+"/BaconProd/Utils/data/Summer13_V1_MC_Uncertainty_AK5PF.txt");
  if(fIsData)   jecUncFiles.push_back(cmsenv+"/BaconProd/Utils/data/Summer13_V1_DATA_Uncertainty_AK5PF.txt");
  std::vector<std::string> jecFilesForID;
  jecFilesForID.push_back(cmsenv+"/BaconProd/Utils/data/START53_V15_L1FastJet_AK5PF.txt");
  jecFilesForID.push_back(cmsenv+"/BaconProd/Utils/data/START53_V15_L2Relative_AK5PF.txt");
  jecFilesForID.push_back(cmsenv+"/BaconProd/Utils/data/START53_V15_L3Absolute_AK5PF.txt");
  if(fIsData) jecFilesForID.push_back(cmsenv+"/BaconProd/Utils/data/Summer13_V1_DATA_L2L3Residual_AK5PF.txt");
  fFillerJet->initJetCorr(jecFiles,jecUncFiles, jecFilesForID);
  fFillerJet->fJetPUIDMVACalc.initialize(baconhep::JetPUIDMVACalculator::k53,
  					 "BDT","",
  					 "BDT",cmsenv+"/BaconProd/Utils/data/TMVAClassificationCategory_JetID_53X_Dec2012.weights.xml");
  fFillerJet->fQGLLCalc.initialize("BDT",cmsenv+"/BaconProd/Utils/data/QG.weights.xml"); 
  
}

//--------------------------------------------------------------------------------------------------
void NtuplerMod::endJob() 
{
  //
  // Save to ROOT file
  //
  //fEventTree->Print();
  fOutputFile->cd();
  fTotalEvents->Write();
  fOutputFile->Write();
  fOutputFile->Close();
  
  delete fFillerEvtInfo;
  delete fFillerGenInfo;
  delete fFillerPV;
  delete fFillerEle;
  delete fFillerMuon;
  delete fFillerPhoton;
  delete fFillerTau;
  delete fFillerJet;
  
  delete fEvtInfo;
  delete fGenEvtInfo;
  delete fGenParArr;
  delete fEleArr;
  delete fMuonArr;
  delete fTauArr;
  delete fJetArr;
  delete fPhotonArr;
  delete fPVArr;
  delete fAddJetArr;
}

//--------------------------------------------------------------------------------------------------
void NtuplerMod::setTriggers()
{
  fTrigger = new baconhep::TTrigger(fHLTFile);
}

//--------------------------------------------------------------------------------------------------
void NtuplerMod::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  fTotalEvents->Fill(1);
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
  
  edm::Handle<trigger::TriggerEvent> hTrgEvt;
  iEvent.getByLabel(fHLTObjTag,hTrgEvt);
  
  fEleArr->Clear();
  fFillerEle->fill(fEleArr, iEvent, iSetup, *pv, nvertices, fPFNoPU, fTrigger->fRecords, *hTrgEvt);

  fMuonArr->Clear();  
  fFillerMuon->fill(fMuonArr, iEvent, iSetup, *pv, fPFNoPU, fPFPU, fTrigger->fRecords, *hTrgEvt);

  fPhotonArr->Clear();  
  fFillerPhoton->fill(fPhotonArr, iEvent, iSetup, *pv, fPFNoPU, fPFPU, fTrigger->fRecords, *hTrgEvt);

  fTauArr->Clear();
  fFillerTau->fill(fTauArr, iEvent, iSetup, *pv, fTrigger->fRecords, *hTrgEvt);

  fJetArr->Clear();
  fAddJetArr->Clear();
  fFillerJet->fill(fJetArr,fAddJetArr, iEvent, iSetup, *pv, fTrigger->fRecords, *hTrgEvt);
  fEventTree->Fill();
}

//--------------------------------------------------------------------------------------------------
void NtuplerMod::initHLT(const edm::TriggerResults& result, const edm::TriggerNames& triggerNames)
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
void NtuplerMod::separatePileUp(const edm::Event &iEvent, const reco::Vertex &pv)
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
void NtuplerMod::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
void NtuplerMod::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup){}
void NtuplerMod::endRun  (const edm::Run& iRun, const edm::EventSetup& iSetup){}
void NtuplerMod::beginLuminosityBlock(const edm::LuminosityBlock& iLumi, const edm::EventSetup& iSetup){}
void NtuplerMod::endLuminosityBlock  (const edm::LuminosityBlock& iLumi, const edm::EventSetup& iSetup){}


//define this as a plug-in
DEFINE_FWK_MODULE(NtuplerMod);
