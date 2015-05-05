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


## PUPPI weights ###
from CommonTools.PileupAlgos.Puppi_cff import puppi
process.puppi = puppi.clone()
process.puppi.candName=cms.InputTag("packedPFCandidates")
process.puppi.vertexName=cms.InputTag("offlineSlimmedPrimaryVertices")

process.puppiSequence = cms.Sequence(process.puppi)

### load default PAT sequence
process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
process.load("PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff")
process.patseq = cms.Sequence(process.patCandidates * process.selectedPatCandidates)
process.p = cms.Path(process.patseq)


# remove unnecessary PAT modules
process.p.remove(process.makePatElectrons)
process.p.remove(process.makePatPhotons)
process.p.remove(process.makePatJets)
process.p.remove(process.makePatTaus)
process.p.remove(process.makePatMETs)
process.p.remove(process.patCandidateSummary)
process.p.remove(process.selectedPatElectrons)
process.p.remove(process.selectedPatPhotons)
process.p.remove(process.selectedPatJets)
process.p.remove(process.selectedPatTaus)
process.p.remove(process.selectedPatCandidateSummary)

### muon selection
process.selectedPatMuons.cut = 'pt>10 && abs(eta)<2.4'

# load user-defined muon PF-isolation values
muon_src, cone_size = 'selectedPatMuons', 0.4

from MuonPFIsolationSequence_cff import *

load_muonPFiso_sequence(process, 'MuonPFIsoSequencePUPPI', algo = 'R04PUPPI',
  coneR = cone_size,
  src = muon_src,
  src_charged_hadron = 'pfPUPPIChargedHadrons',
  src_neutral_hadron = 'pfPUPPINeutralHadrons',
  src_photon         = 'pfPUPPIPhotons'
)

process.MuonPFIsoSequences = cms.Sequence(
  process.MuonPFIsoSequencePUPPI
)


process.p.replace(
  process.selectedPatMuons,
  process.selectedPatMuons *
  process.MuonPFIsoSequences
)


process.p = cms.Path(process.puppiSequence)


# --- Output configuration ------------------------------------------------------------
process.out = cms.OutputModule('PoolOutputModule',
  fileName = cms.untracked.string('first.root'),
)


process.e = cms.EndPath(process.out)
