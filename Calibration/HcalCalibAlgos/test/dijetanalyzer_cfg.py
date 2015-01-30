import FWCore.ParameterSet.Config as cms

process = cms.Process('DIJETANALYSIS')

process.load('FWCore.MessageService.MessageLogger_cfi')

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.autoCond import autoCond
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.GlobalTag.globaltag=autoCond['startup']

#load the response corrections calculator
#process.load('Calibration.HcalCalibAlgos.diJetAnalyzer_cfi')
process.load('JetMETCorrections.Configuration.JetCorrectionProducers_cff')

# run over files

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(
        'file:/eos/uscms/store/user/dgsheffi/QCD_Pt-120To170_13TeV_0002E63C-FC89-E411-B8D6-003048FFCBA8.root'
        ))

#print readFiles

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery=cms.untracked.int32(100)
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.diJetAnalyzer = cms.EDAnalyzer(
    'DiJetAnalyzer',
    pfJetCollName       = cms.string('ak4PFJetsCHS'),
    pfJetCorrName       = cms.string('ak4PFCHSL2L3'),
    hbheRecHitName      = cms.string('hbhereco'),
    hfRecHitName        = cms.string('hfreco'),
    hoRecHitName        = cms.string('horeco'),
    pvCollName          = cms.string('offlinePrimaryVertices'),
    rootHistFilename    = cms.string('dijettree.root'),
    maxDeltaEta         = cms.double(1.5),
    minTagJetEta        = cms.double(0.0),
    maxTagJetEta        = cms.double(5.0),
    minSumJetEt         = cms.double(20.),
    minJetEt            = cms.double(10.),
    maxThirdJetEt       = cms.double(100.),
    maxJetEMF           = cms.double(0.9),
    debug               = cms.untracked.bool(False)
    )

process.p = cms.Path(process.diJetAnalyzer)
