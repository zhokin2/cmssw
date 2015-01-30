import FWCore.ParameterSet.Config as cms

process = cms.Process("DIJETPRODUCER")

#process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
#######process.load('Calibration.HcalAlCaRecoProducers.alcagammajet_cfi')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.autoCond import autoCond
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.GlobalTag.globaltag=autoCond['startup']

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)
process.source = cms.Source("PoolSource",
    fileNames = 
cms.untracked.vstring('file:/eos/uscms/store/user/dgsheffi/QCD_Pt-120To170_13TeV_0002E63C-FC89-E411-B8D6-003048FFCBA8.root')
)

#process.ak5PFJetsCHS = ak5PFJets.clone(
#   src = cms.InputTag("pfNoPileUp")
#)

process.DiJetsProd1 = cms.EDProducer("AlCaDiJetsProducer",
                                      #PhoInput = cms.InputTag("photons"),
                                      PFjetInput = cms.InputTag("ak5PFJets"),
                                      HBHEInput = cms.InputTag("hbhereco"),
                                      HFInput = cms.InputTag("hfreco"),
                                      HOInput = cms.InputTag("horeco"),
                                      #METInput = cms.InputTag("pfMet"),
                                      #Type1METInput = cms.InputTag("pfType1CorrectedMet"),
                                      #gsfeleInput = cms.InputTag("gedGsfElectrons"),
                                      particleFlowInput = cms.InputTag("particleFlow"),
                                      VertexInput = cms.InputTag("offlinePrimaryVertices"),
                                      #ConversionsInput = cms.InputTag("allConversions"),
                                      #rhoInput = cms.InputTag("fixedGridRhoFastjetAll"),
                                      #BeamSpotInput = cms.InputTag("offlineBeamSpot"),
                                      MinPtJet = cms.double(10.0)
                                      #MinPtPhoton = cms.double(10.0)
                                      )

process.DiJetsRecos = cms.OutputModule("PoolOutputModule",
#    outputCommands = cms.untracked.vstring('drop *', 
#        'keep *_GammaJetProd_*_*'),
    fileName = cms.untracked.string('dijet.root')
)

#process.p = cms.Path(process.ak5PFJetsCHS*process.GammaJetProd)
process.p = cms.Path(process.DiJetsProd1)                                                                                       
process.e = cms.EndPath(process.DiJetsRecos)
