import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

### standard includes (can't live without)                                                                                          \
                                                                                                                                     


process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use                                                                     
    fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/user/g/gkrintir/github/Diff13/CMSSW/CMSSW_7_2_3_patch1/src/Isolation/myOutputFile.root'
    )
)

process.demo = cms.EDAnalyzer('DemoAnalyzer',
                              puppi = cms.InputTag("myProducerLabel", "MuonPuppiIso"),
                              puppiMatthias = cms.InputTag("myProducerLabel1")
                              #puppiMatthias = cms.InputTag("puppi", "PuppiWeights")                                                 
                              )


process.TFileService = cms.Service("TFileService",
                                  fileName = cms.string('valueMap_to_check.root')
                                   )


process.p = cms.Path(process.demo)
