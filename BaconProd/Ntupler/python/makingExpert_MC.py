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
process.load('RecoParticleFlow/PFClusterProducer/particleFlowCluster_cff')
process.load('RecoLocalCalo.EcalRecAlgos.EcalSeverityLevelESProducer_cfi')
process.load('RecoLocalCalo.HcalRecAlgos.hcalRecAlgoESProd_cfi')

process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.GlobalTag.globaltag = 'START53_V7G::All'


process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
process.load("BaconProd.Ntupler.myRho_cff")


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.source = cms.Source("PoolSource",
#  firstRun   = cms.untracked.uint32(0),
#  firstLumi  = cms.untracked.uint32(0),
#  firstEvent = cms.untracked.uint32(0),
#  fileNames  = cms.untracked.vstring('file:/afs/cern.ch/work/k/ksung/private/HZZ4lAna/temp/GluGluToHToZZTo4L_M-125_8TeV-powheg-pythia6_PU_S10_START53_V7A.root')
  fileNames  = cms.untracked.vstring('/store/relval/CMSSW_5_3_6/RelValTTbar/GEN-SIM-RECO/PU_START53_V14-v1/0003/3E3EDF4A-E92C-E211-A1BF-003048D2BD66.root')
#/store/relval/CMSSW_5_3_14/RelValTTbar/GEN-SIM-RECO/START53_LV6_Jan31-v1/00000/348BFC23-468B-E311-AAC5-0025905A6118.root')#/store/cmst3/group/cmgtools/CMG/ggH125_HTobb_8TeV_madgraph/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/V5_B/PFAOD_99.root')
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
  outputName    = cms.untracked.string('expert.root'),
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

process.baconSequence = cms.Sequence(process.particleFlowCluster*
                                     process.kt6PFJets*
                                     process.ntupler)
				     
process.p = cms.Path(process.baconSequence)

process.output = cms.OutputModule("PoolOutputModule",                                                                                                                                                      
                                  outputCommands = cms.untracked.vstring('keep *'),                                                                                                                        
                                  fileName       = cms.untracked.string ("test.root")                                                                                                                      
)                                                                                                                                                                                                          

# schedule definition                                                                                                                                                                                      
#process.schedule = cms.Schedule(process.bambu_step)
#process.outpath  = cms.EndPath(process.output)                                                                                                                                                             

