// -*- C++ -*-
//
// Package:    Analyzer/DemoAnalyzer
// Class:      DemoAnalyzer
// 
/**\class DemoAnalyzer DemoAnalyzer.cc Analyzer/DemoAnalyzer/plugins/DemoAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Georgios Krintiras
//         Created:  Fri, 13 Mar 2015 11:02:38 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Common/interface/ValueMap.h"

#include "TH2.h"
#include "TString.h"
#include "math.h"
//
// class declaration
//

class DemoAnalyzer : public edm::EDAnalyzer {
   public:
      explicit DemoAnalyzer(const edm::ParameterSet&);
      ~DemoAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      // ----------member data ---------------------------
    
      edm::EDGetTokenT<edm::ValueMap<double> > _puppiToken;
      edm::EDGetTokenT<edm::ValueMap<double> > _puppiTokenMatthias;

      edm::Service<TFileService> fs_;
      //Funtions

      // Histograms      

      // Detectors
            
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
DemoAnalyzer::DemoAnalyzer(const edm::ParameterSet& iConfig):

  _puppiToken(consumes<edm::ValueMap<double> >(iConfig.getParameter<edm::InputTag>("puppi"))),
  _puppiTokenMatthias(consumes<edm::ValueMap<double> >(iConfig.getParameter<edm::InputTag>("puppiMatthias"))),
  fs_()

{
  //now do what ever initialization is needed
  using namespace edm;

  //h_track_pt_ = fs_->make<TH1F>("Track_Pt", "Track_Pt", 200, 0., 100.); //ok
}


DemoAnalyzer::~DemoAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
DemoAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
#ifdef rrDEBUG
  std::cout << "Famos analysis" << std::endl;
#endif
  // get event and run number
#ifdef rrDEBUG
  int t_Run   = iEvent.id().run();
  int t_Event = iEvent.id().event();
  std::cout
    << " #################################### Run " << t_Run 
    << " Event "                                    << t_Event << " #################################### " 
    << std::endl;
#endif

   using namespace edm;
   
   edm::Handle<edm::ValueMap<double> > puppiIso;
   iEvent.getByToken(_puppiToken, puppiIso);
   //assert(puppiIso.isValid());
   edm::Handle<edm::ValueMap<double> > puppiIsoMatthias;
   iEvent.getByToken(_puppiTokenMatthias, puppiIsoMatthias);
   //assert(puppiIsoMatthias.isValid());

   edm::Handle<edm::ValueMap<float> > puppiweights;
   InputTag token("puppi"); 
   iEvent.getByLabel(token, puppiweights);
   
   edm::Handle<edm::ValueMap<double> > pfweights;
   InputTag token1("myProducerLabel2", "LeptonPFWeightedIso");
   iEvent.getByLabel(token1, pfweights);

   edm::Handle<pat::MuonCollection>  leptons_;
   iEvent.getByLabel("slimmedMuons", leptons_);
   //   assert(leptons.isValid());

   std::vector<const reco::Candidate *> leptons;

   for (const pat::Muon &mu : *leptons_) leptons.push_back(&mu);
  
   if (puppiIso.isValid()!=puppiIsoMatthias.isValid())
     std::cout<<puppiIso.isValid()<<" "<<puppiIsoMatthias.isValid()<<" "<< puppiweights.isValid()<<" "<<pfweights.isValid()<<" "<<leptons.size()<<std::endl;

#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif

}


// ------------ method called once each job just before starting event loop  ------------
void 
DemoAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DemoAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
DemoAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
DemoAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
DemoAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
DemoAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DemoAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


//define this as a plug-in
DEFINE_FWK_MODULE(DemoAnalyzer);
