import FWCore.ParameterSet.Config as cms

process = cms.Process("OWNPARTICLES")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'root://eoscms////eos/cms/store/mc/Phys14DR/TToLeptons_t-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0260CBE1-9F6A-E411-88C8-E0CB4E29C514.root'
    )
)


# -- PF-Weighted
from CommonTools.ParticleFlow.ParticleSelectors.pfAllChargedHadrons_cfi import *
from CommonTools.ParticleFlow.ParticleSelectors.pfAllNeutralHadrons_cfi import *
from CommonTools.ParticleFlow.ParticleSelectors.pfAllPhotons_cfi import *
#process.load('CommonTools.ParticleFlow.deltaBetaWeights_cff')

## PUPPI weights ###
from CommonTools.PileupAlgos.Puppi_cff import puppi
process.puppi = puppi.clone()
process.puppi.candName=cms.InputTag("packedPFCandidates")
process.puppi.vertexName=cms.InputTag("offlineSlimmedPrimaryVertices")
process.puppiSequence = cms.Sequence(process.puppi)

process.myProducerLabel = cms.EDProducer('PUPPILeptonIsoProducer',
    electrons = cms.InputTag("slimmedElectrons"),
    muons = cms.InputTag("slimmedMuons"),
    jets = cms.InputTag("slimmedJets"),
    pfCands = cms.InputTag("packedPFCandidates"),
    puppi = cms.InputTag("puppi", "PuppiWeights")                                     
)


process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('myOutputFile.root')
)

  
process.p = cms.Path(process.puppiSequence*process.myProducerLabel)

process.e = cms.EndPath(process.out)
