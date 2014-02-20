import FWCore.ParameterSet.Config as cms

process = cms.Process('BaconMaking')

# import of standard configurations
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/EventContent/EventContent_cff')
process.load('TrackingTools/TransientTrack/TransientTrackBuilder_cfi')

process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.GlobalTag.globaltag = 'FT_53_V21_AN3::All'

process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")

# import custom configurations
#process.load('BaconProd/Ntupler/mybtagging_cff')     # produce b-tag discriminators
#process.load('BaconProd/Ntupler/myGenJets_cff')
process.load('BaconProd/Ntupler/myJetExtras_cff')    #include gen jets and b-tagging
process.load('BaconProd/Ntupler/myMETFilters_cff')   # apply MET filters set to tagging mode
process.load('BaconProd/Ntupler/myMVAMetData_cff')


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )
process.source = cms.Source("PoolSource",
#  firstRun   = cms.untracked.uint32(0),
#  firstLumi  = cms.untracked.uint32(0),
#  firstEvent = cms.untracked.uint32(0),
  fileNames  = cms.untracked.vstring('/store/data/Run2012D/MET/AOD/PromptReco-v1/000/203/994/163AA34A-570E-E211-ABD5-001D09F242EF.root')
)
process.source.inputCommands = cms.untracked.vstring("keep *",
                                                     "drop *_MEtoEDMConverter_*_*")

process.options = cms.untracked.PSet(
  wantSummary = cms.untracked.bool(True),
  Rethrow     = cms.untracked.vstring('ProductNotFound'),
  fileMode    = cms.untracked.string('NOMERGE')
)

import os
cmssw_base = os.environ['CMSSW_BASE']

process.ntupler = cms.EDAnalyzer('NtuplerMod',
  skipOnHLTFail = cms.untracked.bool(True),
  outputName    = cms.untracked.string('ntuple.root'),
  TriggerFile   = cms.untracked.string(cmssw_base+"/src/BaconAna/DataFormats/data/HLTFile_v0"),                                  

  useGen = cms.untracked.bool(False),
  isData = cms.untracked.bool(True),


  pflowCandidatesName        = cms.untracked.string('particleFlow'),
  primaryVerticesName        = cms.untracked.string('offlinePrimaryVertices'),
  beamspotName               = cms.untracked.string('offlineBeamSpot'),
  pileupInfoName             = cms.untracked.string('addPileupInfo'),
  pflowMETName               = cms.untracked.string('pfMet'),
  rhoForIsolationName        = cms.untracked.string('kt6PFJets'),
  rhoForJetsName             = cms.untracked.string('kt6PFJets'),
  electronName               = cms.untracked.string('gsfElectrons'),
  muonName                   = cms.untracked.string('muons'),
  photonName                 = cms.untracked.string('photons'),
  tauName                    = cms.untracked.string('hpsPFTauProducer'),
  jetName                    = cms.untracked.string('ak5PFJets'),
  genJetName                 = cms.untracked.string('ak5GenJets'),
  jetFlavorName              = cms.untracked.string('AK5byValAlgo'),
  jetFlavorPhysName          = cms.untracked.string('AK5byValPhys'),
  pruneJetName               = cms.untracked.string('ca5PFJetsPruned'),
  subJetName                 = cms.InputTag('ca5PFJetsPruned','SubJets'),
  rhoName                    = cms.untracked.string('kt6PFJets'),
  csvBTagName                = cms.untracked.string('jetCombinedSecondaryVertexMVABJetTags'),
  csvBTagSubJetName          = cms.untracked.string('jetCombinedSecondaryVertexMVABJetTagsSJ'),
  jettiness                  = cms.untracked.string('Njettiness'),
  QGLikelihood               = cms.untracked.string('QGTagger'),
  QGLikelihoodSubjet         = cms.untracked.string('QGTaggerSubJets'),
  computeFullJetInfo         = cms.untracked.bool(True),
  trackName                  = cms.untracked.string('generalTracks'),
  conversionName             = cms.untracked.string('allConversions'),
  ecalBarrelSuperclusterName = cms.untracked.string('correctedHybridSuperClusters'),
  ecalEndcapSuperclusterName = cms.untracked.string('correctedMulti5x5SuperClustersWithPreshower'),
  ecalBarrelRecHitName       = cms.untracked.string('reducedEcalRecHitsEB'),
  ecalEndcapRecHitName       = cms.untracked.string('reducedEcalRecHitsEE')
)

process.baconSequence = cms.Sequence(
                                     process.jetsequence*
                                     process.recoTauClassicHPSSequence*   ### must come after antiktGenJets otherwise conflict on RecoJets/JetProducers/plugins
				     process.MVAMetSeq*
				     process.metFilters*
                                     process.ntupler)

process.p = cms.Path(process.baconSequence)

process.output = cms.OutputModule("PoolOutputModule",                                                                                                                                                     
                                  outputCommands = cms.untracked.vstring('keep *'),                                                                                                                      
                                  fileName       = cms.untracked.string ("test.root")                                                                                                                    
)

# schedule definition                                                                                                       
#process.outpath  = cms.EndPath(process.output)                                                                                                                                                
