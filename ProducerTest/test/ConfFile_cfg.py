import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

### standard includes (can't live without)                                                                                            


process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/user/g/gkrintir/github/isolation/CMSSW_7_3_0/src/Isolation/myOutputFile_.root'
    )
)

process.demo = cms.EDAnalyzer('DemoAnalyzer',
                              puppi = cms.InputTag("myProducerLabel1", "sumChargedHadronPtInPuppiIso"),
                              puppiMatthias = cms.InputTag("myProducerLabel1")
                              #puppiMatthias = cms.InputTag("puppi", "PuppiWeights")
                              )


process.TFileService = cms.Service("TFileService",
                                  fileName = cms.string('electrons_to_check_Euler_constant_.root')
                                   )


process.p = cms.Path(process.demo)
