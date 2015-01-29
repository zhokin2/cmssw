import FWCore.ParameterSet.Config as cms

process = cms.Process('DIJETANALYSIS')

process.load('FWCore.MessageService.MessageLogger_cfi')

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.autoCond import autoCond
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.GlobalTag.globaltag=autoCond['startup']

#load the response corrections calculator
process.load('Calibration.HcalCalibAlgos.diJetAnalyzer_cfi')
process.load('JetMETCorrections.Configuration.JetCorrectionProducers_cff')

# run over files

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(
        'file:/eos/uscms/store/user/dgsheffi/QCD_Pt-120To170_13TeV_0002E63C-FC89-E411-B8D6-003048FFCBA8.root'
        ))

#print readFiles

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery=cms.untracked.int32(100)
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# Load pfNoPileUP

#process.load("CommonTools.ParticleFlow.pfNoPileUp_cff")
#process.load("CommonTools.ParticleFlow.PF2PAT_cff")
#from RecoJets.JetProducers.ak5PFJets_cfi import *
#process.ak5PFJetsCHS = ak5PFJets.clone(
#    src = cms.InputTag("pfNoPileUp")
#    )
#process.load('HcalClosureTest.Analyzers.calcrespcorr_CHSJECs_cff')

# timing
#process.Timing = cms.Service('Timing')

#process.p = cms.Path(process.pfNoPileUpSequence+process.PF2PAT+process.ak5PFJetsCHS+process.calcrespcorrdijets)
process.p = cms.Path(process.diJetAnalyzer)
