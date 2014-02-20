import FWCore.ParameterSet.Config as cms

process = cms.Process('MakingJetMETExpert')

# import of standard configurations
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/EventContent/EventContent_cff')
process.load('TrackingTools/TransientTrack/TransientTrackBuilder_cfi')

process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.GlobalTag.globaltag = 'START53_V7G::All'


process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
process.load("BaconProd.Ntupler.myRho_cff")


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
process.source = cms.Source("PoolSource",
#  firstRun   = cms.untracked.uint32(0),
#  firstLumi  = cms.untracked.uint32(0),
#  firstEvent = cms.untracked.uint32(0),
#  fileNames  = cms.untracked.vstring('file:/afs/cern.ch/work/k/ksung/private/HZZ4lAna/temp/GluGluToHToZZTo4L_M-125_8TeV-powheg-pythia6_PU_S10_START53_V7A.root')
  fileNames  = cms.untracked.vstring('/store/cmst3/group/cmgtools/CMG/ggH125_HTobb_8TeV_madgraph/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/V5_B/PFAOD_99.root')
#file:/afs/cern.ch/work/k/ksung/private/HZZ4lAna/temp/VBF_HToZZTo4L_M-125_8TeV-powheg-pythia6_PU_S10_START53_V7A_4228BABE-70FA-E111-941B-001A92971B26.root')
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

process.ntupler = cms.EDAnalyzer('ExpertMod',
  skipOnHLTFail = cms.untracked.bool(False),
  outputName    = cms.untracked.string('ntuple.root'),
  TriggerFile   = cms.untracked.string(cmssw_base+"/src/BaconAna/DataFormats/data/HLTFile_v0"),                                  
                                 
  useGen = cms.untracked.bool(True),
  genEventInfoName = cms.untracked.string('generator'),
  genParticlesName = cms.untracked.string('genParticles'),
  
  pflowCandidatesName        = cms.untracked.string('particleFlow'),
  primaryVerticesName        = cms.untracked.string('offlinePrimaryVertices'),
  beamspotName               = cms.untracked.string('offlineBeamSpot'),
  pileupInfoName             = cms.untracked.string('addPileupInfo'),
  pflowMETName               = cms.untracked.string('pfMet'),
  rhoForIsolationName        = cms.untracked.string('kt6PFJets'),
  rhoForJetsName             = cms.untracked.string('kt6PFJets'),
  trackName                  = cms.untracked.string('generalTracks'),
)

process.baconSequence = cms.Sequence(process.kt6PFJets*
                                     process.ntupler)
				     
process.p = cms.Path(process.baconSequence)

process.output = cms.OutputModule("PoolOutputModule",                                                                                                                                                      
                                  outputCommands = cms.untracked.vstring('keep *'),                                                                                                                        
                                  fileName       = cms.untracked.string ("test.root")                                                                                                                      
)                                                                                                                                                                                                          

# schedule definition                                                                                                                                                                                      
#process.schedule = cms.Schedule(process.bambu_step)
#process.outpath  = cms.EndPath(process.output)                                                                                                                                                             

