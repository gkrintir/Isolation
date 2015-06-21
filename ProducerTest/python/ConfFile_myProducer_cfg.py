import FWCore.ParameterSet.Config as cms

process = cms.Process("NewIsolation")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(2500) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
          'root://eoscms////eos/cms/store/mc/Phys14DR/TToLeptons_t-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0260CBE1-9F6A-E411-88C8-E0CB4E29C514.root'
      )
)

#  PF-Weighting 
## PF ChargedParticles
process.pfAllChargedParticles = cms.EDFilter("CandPtrSelector",
    src = cms.InputTag("packedPFCandidates"),
    cut = cms.string("abs(charge)>0 && fromPV>=2")
)
## PF Pileup ChargedParticles
process.pfPileUpAllChargedParticles = cms.EDFilter("CandPtrSelector",
    src = cms.InputTag("packedPFCandidates"),
    cut = cms.string("abs(charge)>0 && fromPV<=1")
)

## PF Photons
process.pfAllPhotons = cms.EDFilter("CandPtrSelector", 
    src = cms.InputTag("packedPFCandidates"), 
    cut = cms.string("pt>0.5 && pdgId==22"),
    filter = cms.bool(False)
)

## PF NeutralHadrons
process.pfAllNeutralHadrons = cms.EDFilter("CandPtrSelector",
    src = cms.InputTag("packedPFCandidates"),
    cut = cms.string("pt>0.5 && charge==0 && !pdgId==22")
)

process.PFSequence = cms.Sequence(
    process.pfAllChargedParticles+
    process.pfPileUpAllChargedParticles+
    process.pfAllPhotons+
    process.pfAllNeutralHadrons
)

## test muon inlusion(for the PF-Weighting mitigating Isolation)
process.pfAllChargedParticlesMuonOut = cms.EDFilter("CandPtrSelector",
    src = cms.InputTag("packedPFCandidates"),
    cut = cms.string("abs(charge)>0 && fromPV>=2 && abs(pdgId)!=13")
)

from CommonTools.ParticleFlow.deltaBetaWeights_cfi import *
## with muon included
process.pfWeightedPhotonsMuonIn = pfWeightedPhotons.clone()
process.pfWeightedPhotonsMuonIn.src  =  cms.InputTag('pfAllPhotons')
process.pfWeightedPhotonsMuonIn.chargedFromPV  = cms.InputTag('pfAllChargedParticles')
process.pfWeightedPhotonsMuonIn.chargedFromPU  = cms.InputTag("pfPileUpAllChargedParticles")

process.pfWeightedNeutralHadronsMuonIn = pfWeightedNeutralHadrons.clone()
process.pfWeightedNeutralHadronsMuonIn.src  = cms.InputTag("pfAllNeutralHadrons")
process.pfWeightedNeutralHadronsMuonIn.chargedFromPV  = cms.InputTag("pfAllChargedParticles")
process.pfWeightedNeutralHadronsMuonIn.chargedFromPU  = cms.InputTag("pfPileUpAllChargedParticles")

## with muon excluded
process.pfWeightedPhotonsMuonOut = pfWeightedPhotons.clone()
process.pfWeightedPhotonsMuonOut.src  =  cms.InputTag('pfAllPhotons')
process.pfWeightedPhotonsMuonOut.chargedFromPV  = cms.InputTag('pfAllChargedParticlesMuonOut')
process.pfWeightedPhotonsMuonOut.chargedFromPU  = cms.InputTag("pfPileUpAllChargedParticles")

process.pfWeightedNeutralHadronsMuonOut = pfWeightedNeutralHadrons.clone()
process.pfWeightedNeutralHadronsMuonOut.src  = cms.InputTag("pfAllNeutralHadrons")
process.pfWeightedNeutralHadronsMuonOut.chargedFromPV  = cms.InputTag("pfAllChargedParticlesMuonOut")
process.pfWeightedNeutralHadronsMuonOut.chargedFromPU  = cms.InputTag("pfPileUpAllChargedParticles")


process.pfDeltaBetaWeightingSequence = cms.Sequence(
    process.pfWeightedPhotonsMuonIn * 
    process.pfWeightedNeutralHadronsMuonIn + 
    process.pfAllChargedParticlesMuonOut * 
    process.pfWeightedPhotonsMuonOut * 
    process.pfWeightedNeutralHadronsMuonOut
)

## produce values for the reweighted isolation deposits (if enabled) and (relative) isolation
pfWeightedMuonIsoInfo_list = []
CandidateSum_list = ["sumChargedCandidatePt", "sumNeutralHadronEt", "sumPhotonEt"]
writeCandidateSums_flag = True
testLeptonsInCandidatePt_flag = True 
testMuonInIsoWeights_cfg = ["MuonIn", "MuonOut"]

if testLeptonsInCandidatePt_flag:
    testLeptonsInCandidatePt_cfg = ["wLepton", "woLepton"] 
else:
    testLeptonsInCandidatePt_cfg = ["woLepton"]

for testMuonInIsoWeights in testMuonInIsoWeights_cfg:
    for testLeptonsInCandidatePt in testLeptonsInCandidatePt_cfg:
        if 'wo' in testLeptonsInCandidatePt:
            includeLeptonsInCandidatePt_flag = False
        else:
            includeLeptonsInCandidatePt_flag = True
        for coneSize in [0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5]:
            pfWeightedLeptonIso = cms.EDProducer('PFWeightedMuonIsoProducer',
            muons = cms.InputTag("slimmedMuons"),
            electrons = cms.InputTag("slimmedElectrons"),
            pfCands = cms.InputTag("packedPFCandidates"),
            pfWeightedHadrons = cms.InputTag("pfWeightedNeutralHadrons%s" % testMuonInIsoWeights),
            pfWeightedPhotons = cms.InputTag("pfWeightedPhotons%s" % testMuonInIsoWeights),
            writeCandidateSums = cms.untracked.bool(writeCandidateSums_flag),
            includeLeptoninIso = cms.untracked.bool(includeLeptonsInCandidatePt_flag),
            dRConeSize = cms.untracked.double(coneSize)                                    
            )
            producerName =  "pfWeightedIso%sR%02i%s" % (testMuonInIsoWeights, int(coneSize*100), testLeptonsInCandidatePt)
            setattr(process,producerName,pfWeightedLeptonIso)
            process.pfDeltaBetaWeightingSequence+=getattr(process,producerName)
            pfWeightedMuonIsoInfo_list.append(producerName)
            if (writeCandidateSums_flag):
                for CandidateSum in CandidateSum_list:
                    if ("Charged" in CandidateSum):
                        pfWeightedLeptonIso = cms.EDProducer('PFWeightedMuonIsoProducer',
                        muons = cms.InputTag("slimmedMuons"),
                        electrons = cms.InputTag("slimmedElectrons"),
                        pfCands = cms.InputTag("packedPFCandidates"),
                        pfWeightedHadrons = cms.InputTag("pfWeightedNeutralHadrons%s" % testMuonInIsoWeights),
                        pfWeightedPhotons =cms.InputTag("pfWeightedPhotons%s" % testMuonInIsoWeights),
                        writeCandidateSums = cms.untracked.bool(writeCandidateSums_flag),
                        includeLeptoninIso = cms.untracked.bool(includeLeptonsInCandidatePt_flag),
                        dRConeSize = cms.untracked.double(coneSize)                                                  
                        )
                        if (includeLeptonsInCandidatePt_flag):
                            producerName = "pfWeightedIso%sR%02isumChargedParticlePt" % (testMuonInIsoWeights, int(coneSize*100))
                        else:
                            producerName = "pfWeightedIso%sR%02isumChargedHadronPt" % (testMuonInIsoWeights, int(coneSize*100))
                        setattr(process,producerName,pfWeightedLeptonIso)
                        process.pfDeltaBetaWeightingSequence+=getattr(process,producerName)
                        pfWeightedMuonIsoInfo_list.append(producerName)
                    elif ("Neutral" not in pfWeightedMuonIsoInfo_list or "Photon" not in pfWeightedMuonIsoInfo_list):
                        pfWeightedLeptonIso = cms.EDProducer('PFWeightedMuonIsoProducer',
                        muons = cms.InputTag("slimmedMuons"),
                        electrons = cms.InputTag("slimmedElectrons"),
                        pfCands = cms.InputTag("packedPFCandidates"),
                        pfWeightedHadrons = cms.InputTag("pfWeightedNeutralHadrons%s" % testMuonInIsoWeights),
                        pfWeightedPhotons = cms.InputTag("pfWeightedPhotons%s" % testMuonInIsoWeights),
                        writeCandidateSums = cms.untracked.bool(writeCandidateSums_flag),
                        includeLeptoninIso = cms.untracked.bool(includeLeptonsInCandidatePt_flag),
                        dRConeSize = cms.untracked.double(coneSize)
                        )
                        producerName = "pfWeightedIso%sR%02i%s" % (testMuonInIsoWeights, int(coneSize*100), CandidateSum)
                        setattr(process,producerName,pfWeightedLeptonIso)
                        process.pfDeltaBetaWeightingSequence+=getattr(process,producerName)
                        pfWeightedMuonIsoInfo_list.append(producerName)


# iPUPPI 
## with muon included
from CommonTools.PileupAlgos.Puppi_cff import puppi
process.puppiMuonIn = puppi.clone()
process.puppiMuonIn.candName=cms.InputTag("packedPFCandidates")
process.puppiMuonIn.vertexName=cms.InputTag("offlineSlimmedPrimaryVertices")

## with muon excluded
process.packedPFCandidatesMuonOut  = cms.EDFilter("CandPtrSelector", 
    src = cms.InputTag("packedPFCandidates"), 
    cut = cms.string("fromPV>=2 && abs(pdgId)!=13" ) 
)

process.puppiMuonOut = puppi.clone()
process.puppiMuonOut.candName=cms.InputTag("packedPFCandidatesMuonOut") 
process.puppiMuonOut.vertexName=cms.InputTag("offlineSlimmedPrimaryVertices")

process.puppiSequence = cms.Sequence(
    process.puppiMuonIn+ 
    process.packedPFCandidatesMuonOut+
    process.puppiMuonOut
)

## produce values for the reweighted isolation deposits (if enabled) and (relative) isolation
iPUPPIMuonIsoInfo_list = []
CandidateSum_list = ["sumChargedCandidatePt", "sumNeutralHadronEt", "sumPhotonEt"]
writeCandidateSums_flag = True
testLeptonsInCandidatePt_flag = True
testMuonInIsoWeights_cfg = ["MuonIn", "MuonOut"]

if testLeptonsInCandidatePt_flag:
    testLeptonsInCandidatePt_cfg = ["wLepton", "woLepton"]
else:
    testLeptonsInCandidatePt_cfg = ["woLepton"]

for testMuonInIsoWeights in testMuonInIsoWeights_cfg:
    for testLeptonsInCandidatePt in testLeptonsInCandidatePt_cfg:
        if 'wo' in testLeptonsInCandidatePt:
            includeLeptonsInCandidatePt_flag = False
        else:
            includeLeptonsInCandidatePt_flag = True
        for coneSize in [0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5]:
            puppiIsoLepton = cms.EDProducer('PUPPIMuonIsoProducer',
            muons = cms.InputTag("slimmedMuons"),
            electrons = cms.InputTag("slimmedElectrons"),
            pfCands = cms.InputTag("packedPFCandidates"),
            puppi = cms.InputTag("puppi%s" % testMuonInIsoWeights),
            dRConeSize = cms.untracked.double(coneSize),
            writeCandidateSums = cms.untracked.bool(writeCandidateSums_flag),
            includeLeptoninIso = cms.untracked.bool(includeLeptonsInCandidatePt_flag)
            )
            producerName = "puppiIso%sR%02i%s" % (testMuonInIsoWeights, int(coneSize*100), testLeptonsInCandidatePt)
            setattr(process,producerName,puppiIsoLepton)
            process.puppiSequence+=getattr(process,producerName)
            iPUPPIMuonIsoInfo_list.append(producerName)
            if (writeCandidateSums_flag):
                for CandidateSum in CandidateSum_list:
                    if ("Charged" in CandidateSum):
                        puppiIsoLepton = cms.EDProducer('PUPPIMuonIsoProducer',
                        muons = cms.InputTag("slimmedMuons"),
                        electrons = cms.InputTag("slimmedElectrons"),
                        pfCands = cms.InputTag("packedPFCandidates"),
                        puppi = cms.InputTag("puppi%s" % testMuonInIsoWeights),
                        dRConeSize = cms.untracked.double(coneSize),
                        writeCandidateSums = cms.untracked.bool(writeCandidateSums_flag),
                        includeLeptoninIso = cms.untracked.bool(includeLeptonsInCandidatePt_flag)
                        )
                        if (includeLeptonsInCandidatePt_flag):
                            producerName = "puppiIso%sR%02isumChargedParticlePt" % (testMuonInIsoWeights, int(coneSize*100))
                        else:
                            producerName = "puppiIso%sR%02isumChargedHadronPt" % (testMuonInIsoWeights, int(coneSize*100))
                        setattr(process,producerName,puppiIsoLepton)
                        process.puppiSequence+=getattr(process,producerName)
                        iPUPPIMuonIsoInfo_list.append(producerName)
                    elif ("Neutral" not in iPUPPIMuonIsoInfo_list or "Photon" not in iPUPPIMuonIsoInfo_list):
                        puppiIsoLepton = cms.EDProducer('PUPPIMuonIsoProducer',
                        muons = cms.InputTag("slimmedMuons"),
                        electrons = cms.InputTag("slimmedElectrons"),
                        pfCands = cms.InputTag("packedPFCandidates"),
                        puppi = cms.InputTag("puppi%s" % testMuonInIsoWeights),
                        dRConeSize = cms.untracked.double(coneSize),
                        writeCandidateSums = cms.untracked.bool(writeCandidateSums_flag),
                        includeLeptoninIso = cms.untracked.bool(includeLeptonsInCandidatePt_flag)
                        )
                        producerName = "puppiIso%sR%02i%s" % (testMuonInIsoWeights, int(coneSize*100), CandidateSum)
                        setattr(process,producerName,puppiIsoLepton)
                        process.puppiSequence+=getattr(process,producerName)
                        iPUPPIMuonIsoInfo_list.append(producerName)


process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('myOutputFile_.root')
)

# Additional output definition for debuggin
process.TFileService = cms.Service("TFileService",
                                  fileName = cms.string('electrons_to_check.root')
                                   )
  
process.p = cms.Path(
    process.PFSequence*
    process.pfDeltaBetaWeightingSequence*
    process.puppiSequence
)

process.e = cms.EndPath(process.out)
