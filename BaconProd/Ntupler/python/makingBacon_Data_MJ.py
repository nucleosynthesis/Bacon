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
process.load('BaconProd/Ntupler/myJetExtras04_cff')    # include gen jets and b-tagging
process.load('BaconProd/Ntupler/myJetExtras05_cff')    # include gen jets and b-tagging
process.load('BaconProd/Ntupler/myJetExtras06_cff')    # include gen jets and b-tagging
process.load('BaconProd/Ntupler/myJetExtras07_cff')    # include gen jets and b-tagging
process.load('BaconProd/Ntupler/myJetExtras08_cff')    # include gen jets and b-tagging
process.load('BaconProd/Ntupler/myJetExtras09_cff')    # include gen jets and b-tagging
process.load('BaconProd/Ntupler/myMETFilters_cff')     # apply MET filters set to tagging mode
process.load('BaconProd/Ntupler/myMVAMetData_cff')


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.source = cms.Source("PoolSource",
  fileNames  = cms.untracked.vstring('AAA')
)
process.source.inputCommands = cms.untracked.vstring("keep *",
                                                     "drop *_MEtoEDMConverter_*_*")

process.options = cms.untracked.PSet(
#  wantSummary = cms.untracked.bool(True),
  Rethrow     = cms.untracked.vstring('ProductNotFound'),
  fileMode    = cms.untracked.string('NOMERGE')
)

import os
cmssw_base = os.environ['CMSSW_BASE']

process.ntupler = cms.EDAnalyzer('NtuplerMod',
  skipOnHLTFail = cms.untracked.bool(True),
  outputName    = cms.untracked.string('BBB'),
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
  applyMuscle                = cms.untracked.bool(False),
  photonName                 = cms.untracked.string('photons'),
  tauName                    = cms.untracked.string('hpsPFTauProducer'),
  NumCones                   = cms.untracked.int32(6),    
  MinCone                    = cms.untracked.double(0.4), 
  ConeIter                   = cms.untracked.double(0.1),
  jetName                    = cms.untracked.string('PFJets'),
  genJetName                 = cms.untracked.string('GenJets'),
  jetFlavorName              = cms.untracked.string('byValAlgo'),
  jetFlavorPhysName          = cms.untracked.string('byValPhys'),
  pruneJetName               = cms.untracked.string('caPFJetsPruned'),
  subJetName                 = cms.untracked.string('caPFJetsPruned'),
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

process.baconSequence = cms.Sequence(process.metFilters*
                                     process.recojetsequence*
				     process.AK4jetsequence*
                                     process.AK5jetsequence*
                                     process.AK6jetsequence*
                                     process.AK7jetsequence*
                                     process.AK8jetsequence*
                                     process.AK9jetsequence*                                     
                                     process.recoTauClassicHPSSequence*   ### must come after antiktGenJets otherwise conflict on RecoJets/JetProducers/plugins
				     process.MVAMetSeq*
				     process.ntupler)

process.p = cms.Path(process.baconSequence)

#process.output = cms.OutputModule("PoolOutputModule",                                                                                                                                                     
#                                  outputCommands = cms.untracked.vstring('keep *'),                                                                                                                      
#                                  fileName       = cms.untracked.string ("test.root")                                                                                                                    
#)

# schedule definition                                                                                                       
#process.outpath  = cms.EndPath(process.output)                                                                                                                                                
