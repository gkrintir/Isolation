import FWCore.ParameterSet.Config as cms

process = cms.Process("OWNPARTICLES")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(2500) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        #'root://xrootd.unl.edu//store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU4bx50_PHYS14_25_V1-v1/00000/003B199E-0F81-E411-8E76-0025905A60B0.root'
        'root://eoscms////eos/cms/store/mc/Phys14DR/TToLeptons_t-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0260CBE1-9F6A-E411-88C8-E0CB4E29C514.root'
        #'root://xrootd.unl.edu//store/mc/Phys14DR/QCD_Pt-20toInf_MuEnrichedPt15_PionKaonDecay_Tune4C_13TeV_pythia8/MINIAODSIM\
#/PU20bx25_PHYS14_25_V1-v3/10000/325BE5B9-AAA6-E411-8371-001E673972E2.root'
    )
)


# PUPPI weights
#test muon inlusion
process.packedPFCandidatesWoMuon  = cms.EDFilter("CandPtrSelector", 
          src = cms.InputTag("packedPFCandidates"), 
          cut = cms.string("fromPV>=2 && abs(pdgId)!=13" ) ) 

##with muon included
from CommonTools.PileupAlgos.Puppi_cff import puppi
process.puppi = puppi.clone()
process.puppi.candName=cms.InputTag("packedPFCandidates") 
process.puppi.vertexName=cms.InputTag("offlineSlimmedPrimaryVertices")

process.puppiSequence = cms.Sequence(process.puppi)

##with muon excluded
process.puppiNoLeptons = puppi.clone()
process.puppiNoLeptons.candName=cms.InputTag("packedPFCandidatesWoMuon")#packedPFCandidates 
process.puppiNoLeptons.vertexName=cms.InputTag("offlineSlimmedPrimaryVertices")

#process.puppiSequence = cms.Sequence(process.packedPFCandidatesWoMuon + process.puppiNoLeptons + process.puppi)

## with muon excluded 
#end of test


process.myProducerLabel = cms.EDProducer('PUPPILeptonIsoProducer',
    electrons = cms.InputTag("slimmedMuons"),
    muons = cms.InputTag("slimmedMuons"),
    pfCands = cms.InputTag("packedPFCandidates"),
    puppi = cms.InputTag("puppi", "PuppiWeights"),
    dRConeSize = cms.untracked.double(0.4),
    writeCandidateSums = cms.untracked.bool(True),
    includeLeptoninIso = cms.untracked.bool(True)
)


process.myProducerLabel1 = cms.EDProducer('PUPPILeptonIsoProducer_Matthias_readNOWeights',
    muons = cms.InputTag("slimmedMuons"),
    electrons = cms.InputTag("slimmedElectrons"),
    pfCands = cms.InputTag("packedPFCandidates"),
    puppi = cms.InputTag("puppi"),
    dRConeSize = cms.untracked.double(0.3),
    writeCandidateSums = cms.untracked.bool(True),
    includeLeptoninIso = cms.untracked.bool(False)
)

process.myProducerLabel2 = cms.EDProducer('PFWeightedLeptonIsoProducer_new',
    leptons = cms.InputTag("slimmedMuons"),
    pfCands = cms.InputTag("packedPFCandidates"),
    pfWeightedHadrons = cms.InputTag("pfWeightedNeutralHadrons"),
    pfWeightedPhotons =cms.InputTag("pfWeightedPhotons"),
    dRConeSize = cms.untracked.double(0.4),
    writeCandidateSums = cms.untracked.bool(True),
    includeLeptoninIso = cms.untracked.bool(True)

)

'''
#PF-Weighted Candidates
## PF ChargedParticles
process.pfAllChargedParticles = cms.EDFilter("CandPtrSelector",
    src = cms.InputTag("packedPFCandidates"),
    cut = cms.string("abs(charge)>0 && fromPV>=3")
)
## PF Pileup ChargedParticles
process.pfPileUpAllChargedParticles = cms.EDFilter("CandPtrSelector",
    src = cms.InputTag("packedPFCandidates"),
    cut = cms.string("abs(charge)>0 && fromPV<=2")
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

process.PFSequence = cms.Sequence(process.pfAllChargedParticles+process.pfPileUpAllChargedParticles+process.pfAllPhotons+process.pfAllNeutralHadrons)

## PF weights
#test muon inlusion
process.pfAllChargedParticlesMuonOut = cms.EDFilter("CandPtrSelector",
    src = cms.InputTag("packedPFCandidates"),
    cut = cms.string(
        "
        abs(charge)>0 && fromPV>=3 && abs(pdgId)!=13
        "
    )
)


from CommonTools.ParticleFlow.deltaBetaWeights_cfi import *
##with muon included
process.pfWeightedPhotonsMuonIn = pfWeightedPhotons.clone()
process.pfWeightedPhotonsMuonIn.src  =  cms.InputTag('pfAllPhotons')
process.pfWeightedPhotonsMuonIn.chargedFromPV  = cms.InputTag('pfAllChargedParticles')
process.pfWeightedPhotonsMuonIn.chargedFromPU  = cms.InputTag("pfPileUpAllChargedParticles")

process.pfWeightedNeutralHadronsMuonIn = pfWeightedNeutralHadrons.clone()
process.pfWeightedNeutralHadronsMuonIn.src  = cms.InputTag("pfAllNeutralHadrons")
process.pfWeightedNeutralHadronsMuonIn.chargedFromPV  = cms.InputTag("pfAllChargedParticles")
process.pfWeightedNeutralHadronsMuonIn.chargedFromPU  = cms.InputTag("pfPileUpAllChargedParticles")

#process.pfDeltaBetaWeightingSequence = cms.Sequence(process.pfWeightedPhotons*process.pfWeightedNeutralHadrons)


##with muon excluded
process.pfWeightedPhotonsMuonOut = pfWeightedPhotons.clone()
process.pfWeightedPhotonsMuonOut.src  =  cms.InputTag('pfAllPhotons')
process.pfWeightedPhotonsMuonOut.chargedFromPV  = cms.InputTag('pfAllChargedParticlesMuonOut')
process.pfWeightedPhotonsMuonOut.chargedFromPU  = cms.InputTag("pfPileUpAllChargedParticles")

process.pfWeightedNeutralHadronsMuonOut = pfWeightedNeutralHadrons.clone()
process.pfWeightedNeutralHadronsMuonOut.src  = cms.InputTag("pfAllNeutralHadrons")
process.pfWeightedNeutralHadronsMuonOut.chargedFromPV  = cms.InputTag("pfAllChargedParticlesMuonOut")
process.pfWeightedNeutralHadronsMuonOut.chargedFromPU  = cms.InputTag("pfPileUpAllChargedParticles")


process.pfDeltaBetaWeightingSequence = cms.Sequence(
    process.pfWeightedPhotonsMuonIn * process.pfWeightedNeutralHadronsMuonIn + 
    process.pfAllChargedParticlesMuonOut * process.pfWeightedPhotonsMuonOut * process.pfWeightedNeutralHadronsMuonOut
    )

pfWeightedIsoLeptonListInfo=[]
CandidateSumList = ["sumChargedHadronPt", "sumChargedParticlePt", "sumNeutralHadronEt", "sumPhotonEt"]
writeCandidateSums_flag = True
includeLeptoninIso_flag = True
Configurations = ["MuonIn", "MuonOut"]

if includeLeptoninIso_flag:
    includeLeptoninIso_flags = ["wLepton", "woLepton"] 
else:
    includeLeptoninIso_flags = ["woLepton"]
for cfg in Configurations:
    for includeLepton_flag in includeLeptoninIso_flags:
        if 'wo' in includeLepton_flag:
            includeLeptoninIso_flag = False
        else:
            includeLeptoninIso_flag = True
        for coneSize in [0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5]:
            pfWeightedLeptonIso = cms.EDProducer('PFWeightedLeptonIsoProducer_new',
            muons = cms.InputTag("slimmedMuons"),
            electrons = cms.InputTag("slimmedElectrons"),
            pfCands = cms.InputTag("packedPFCandidates"),
            pfWeightedHadrons = cms.InputTag("pfWeightedNeutralHadrons%s" % cfg),
            pfWeightedPhotons =cms.InputTag("pfWeightedPhotons%s" % cfg),
            writeCandidateSums = cms.untracked.bool(writeCandidateSums_flag),
            includeLeptoninIso = cms.untracked.bool(includeLeptoninIso_flag),
            dRConeSize = cms.untracked.double(coneSize)                                    
            )
            producerName =  "pfWeightedIso%sR%02i%s" % (cfg, int(coneSize*100), includeLepton_flag)
            pfWeightedIsoLeptonListInfo.append(producerName)
            setattr(process,producerName,pfWeightedLeptonIso)
            process.pfDeltaBetaWeightingSequence+=getattr(process,producerName)
            if (writeCandidateSums_flag):
                for CandidateSum in CandidateSumList:
                    if ("Charged" in CandidateSum):
                        pfWeightedLeptonIso = cms.EDProducer('PFWeightedLeptonIsoProducer_new',
                        muons = cms.InputTag("slimmedMuons"),
                        electrons =cms.InputTag("slimmedElectrons"),
                        pfCands = cms.InputTag("packedPFCandidates"),
                        pfWeightedHadrons = cms.InputTag("pfWeightedNeutralHadrons%s" % cfg),
                        pfWeightedPhotons =cms.InputTag("pfWeightedPhotons%s" % cfg),
                        writeCandidateSums = cms.untracked.bool(writeCandidateSums_flag),
                        includeLeptoninIso = cms.untracked.bool(includeLeptoninIso_flag),
                        dRConeSize = cms.untracked.double(coneSize)                                                  
                        )
                        if (includeLeptoninIso_flag):
                            producerName = "pfWeightedIso%sR%02iChargedParticlePt" % (cfg, int(coneSize*100))
                        else:
                            producerName = "pfWeightedIso%sR%02iChargedHadronPt" % (cfg, int(coneSize*100))
                        pfWeightedIsoLeptonListInfo.append(producerName)
                        setattr(process,producerName,pfWeightedLeptonIso)
                        process.pfDeltaBetaWeightingSequence+=getattr(process,producerName)
                    elif ("Neutral" not in pfWeightedIsoLeptonListInfo or "Photon" not in pfWeightedIsoLeptonListInfo):
                        print '!!!!!', CandidateSum
                        pfWeightedLeptonIso = cms.EDProducer('PFWeightedLeptonIsoProducer_new',
                        muons = cms.InputTag("slimmedMuons"),
                        electrons =cms.InputTag("slimmedElectrons"),
                        pfCands = cms.InputTag("packedPFCandidates"),
                        pfWeightedHadrons = cms.InputTag("pfWeightedNeutralHadrons%s" % cfg),
                        pfWeightedPhotons =cms.InputTag("pfWeightedPhotons%s" % cfg),
                        writeCandidateSums = cms.untracked.bool(writeCandidateSums_flag),
                        includeLeptoninIso = cms.untracked.bool(includeLeptoninIso_flag),
                        dRConeSize = cms.untracked.double(coneSize)
                        )
                        producerName = "pfWeightedIso%sR%02i%s" % (cfg, int(coneSize*100), CandidateSum)
                        pfWeightedIsoLeptonListInfo.append(producerName)
                        setattr(process,producerName,pfWeightedLeptonIso)
                        process.pfDeltaBetaWeightingSequence+=getattr(process,producerName)


'''
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('myOutputFile_.root')
)

# Additional output definition for debuggin
process.TFileService = cms.Service("TFileService",
                                  fileName = cms.string('electrons_to_check.root')
                                   )
  
process.p = cms.Path(
                     #process.PFSequence*
                     #process.pfDeltaBetaWeightingSequence
                     process.puppiSequence*
                     #process.myProducerLabel*
                     process.myProducerLabel1
                     #process.myProducerLabel2
                     )

process.e = cms.EndPath(process.out)
