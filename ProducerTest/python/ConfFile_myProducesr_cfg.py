import FWCore.ParameterSet.Config as cms

process = cms.Process("OWNPARTICLES")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        #'root://xrootd.unl.edu//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU4bx50_PHYS14_25_V1-v1/00000/003B199E-0F81-E411-8E76-0025905A60B0.root'
        'root://eoscms////eos/cms/store/mc/Phys14DR/TToLeptons_t-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0260CBE1-9F6A-E411-88C8-E0CB4E29C514.root'
    )
)


#PF-Weighted Candidates
## PF ChargedParticles
process.pfAllChargedParticles = cms.EDFilter("CandPtrSelector",
    src = cms.InputTag("packedPFCandidates"),
    cut = cms.string(
        '''
        charge!=0 && fromPV
        '''
    )
)
## PF Pileup ChargedParticles
process.pfPileUpAllChargedParticles = cms.EDFilter("CandPtrSelector",
    src = cms.InputTag("packedPFCandidates"),
    cut = cms.string(
        '''
        charge!=0 && !fromPV
        '''
    )
)


## PF Photons
process.pfAllPhotons = cms.EDFilter("CandPtrSelector", 
                                    src = cms.InputTag("slimmedPhotons"), 
                                    cut = cms.string("pt>0.5 && pdgId==22"),
                                    filter = cms.bool(True)
                                    )

## PF NeutralHadrons
process.pfAllNeutralHadrons = cms.EDFilter("CandPtrSelector",
    src = cms.InputTag("packedPFCandidates"),
    cut = cms.string(
        ''' 
        pt>0.5 && charge==0 && !pdgId==22
        '''
    )
)

process.PFSequence = cms.Sequence(process.pfAllChargedParticles+process.pfPileUpAllChargedParticles+process.pfAllPhotons+process.pfAllNeutralHadrons)

## PF weights
from CommonTools.ParticleFlow.deltaBetaWeights_cfi import *
process.pfWeightedPhotons = pfWeightedPhotons.clone()
process.pfWeightedPhotons.src  =  cms.InputTag('pfAllPhotons')
process.pfWeightedPhotons.chargedFromPV  = cms.InputTag('pfAllChargedParticles')
process.pfWeightedPhotons.chargedFromPU  = cms.InputTag("pfPileUpAllChargedParticles")

process.pfWeightedNeutralHadrons = pfWeightedNeutralHadrons.clone()
process.pfWeightedNeutralHadrons.src  = cms.InputTag("pfAllNeutralHadrons")
process.pfWeightedNeutralHadrons.chargedFromPV  = cms.InputTag("pfAllChargedParticles")
process.pfWeightedNeutralHadrons.chargedFromPU  = cms.InputTag("pfPileUpAllChargedParticles")

process.pfDeltaBetaWeightingSequence = cms.Sequence(process.pfWeightedPhotons*process.pfWeightedNeutralHadrons)

# PUPPI weights
from CommonTools.PileupAlgos.Puppi_cff import puppi
process.puppi = puppi.clone()
process.puppi.candName=cms.InputTag("packedPFCandidates")
process.puppi.vertexName=cms.InputTag("offlineSlimmedPrimaryVertices")
process.puppiSequence = cms.Sequence(process.puppi)


process.myProducerLabel = cms.EDProducer('PUPPILeptonIsoProducer',
    electrons = cms.InputTag("slimmedElectrons"),
    muons = cms.InputTag("slimmedMuons"),
    pfCands = cms.InputTag("packedPFCandidates"),
    puppi = cms.InputTag("puppi", "PuppiWeights"),
    dRConeSize = cms.untracked.double(0.4),
    writeCandidateSums = cms.untracked.bool(True),
    includeLeptoninIso = cms.untracked.bool(False)
)



process.myProducerLabel1 = cms.EDProducer('PFWeightedLeptonIsoProducer',
    electrons = cms.InputTag("slimmedElectrons"),
    muons = cms.InputTag("slimmedMuons"),
    pfCands = cms.InputTag("packedPFCandidates"),
    pfWeightedHadrons = cms.InputTag("pfWeightedNeutralHadrons"),
    pfWeightedPhotons =cms.InputTag("pfWeightedPhotons"),
    dRConeSize = cms.untracked.double(0.4),
)


process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('myOutputFile.root')
)

  
process.p = cms.Path(process.PFSequence*
                     process.pfDeltaBetaWeightingSequence*
                     process.puppiSequence*
                     process.myProducerLabel*
                     process.myProducerLabel1
                     )

process.e = cms.EndPath(process.out)
