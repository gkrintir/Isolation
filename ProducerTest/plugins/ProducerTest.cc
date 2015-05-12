// -*- C++ -*-
//
// Package:    test/ProducerTest
// Class:      ProducerTest
// 
/**\class ProducerTest ProducerTest.cc test/ProducerTest/plugins/ProducerTest.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Georgios Krintiras
//         Created:  Mon, 27 Apr 2015 23:04:18 GMT
//
//


// system include files
#include <memory>
#include <vector>
// user include files

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "DataFormats/Common/interface/ValueMap.h"

//
// class declaration
//

class ProducerTest : public edm::EDProducer {
   public:
      explicit ProducerTest(const edm::ParameterSet&);
      ~ProducerTest();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
  
      //template<typename Hand, typename T>
      void storeMuonIso(edm::Event &iEvent,
			const edm::Handle<pat::MuonCollection > & handle, 
			const std::vector<double> & values,
			const std::string    & label) const ;

  /*  
      void storeElectronIso(edm::Event &iEvent,
			    const edm::Handle<pat::ElectronCollection > & handle, 
			    const std::vector<double> & values,
			    const std::string    & label) const ;
  */
      // ----------member data ---------------------------
      edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
      edm::EDGetTokenT<pat::MuonCollection> muonToken_;
      edm::EDGetTokenT<pat::JetCollection> jetToken_;
      edm::EDGetTokenT<pat::PackedCandidateCollection> pfToken_;
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
ProducerTest::ProducerTest(const edm::ParameterSet& iConfig):
    electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
    muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
    jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
    pfToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCands"))) 

{
  //register your products
/* Examples
   produces<ExampleData2>();

   //if do put with a label
   produces<ExampleData2>("label");
 
   //if you want to put into the Run
   produces<ExampleData2,InRun>();
*/
   //now do what ever other initialization is needed

  produces<edm::ValueMap<double> >("MuonIso");
  produces<edm::ValueMap<double> >("ElectronIso");

  
}


ProducerTest::~ProducerTest()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
ProducerTest::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;

    edm::Handle<pat::MuonCollection> muons;
    iEvent.getByToken(muonToken_, muons);
    edm::Handle<pat::ElectronCollection> electrons;
    iEvent.getByToken(electronToken_, electrons);
    edm::Handle<pat::PackedCandidateCollection> pfs;
    iEvent.getByToken(pfToken_, pfs);
    edm::Handle<pat::JetCollection> jets;
    iEvent.getByToken(jetToken_, jets);
    
    edm::Handle<reco::PFCandidateCollection> puppi;
    iEvent.getByLabel("puppi", puppi);
    assert(puppi.isValid());
    //    const reco::PFCandidateCollection *PFCol = puppi.product();
        
    edm::Handle<edm::ValueMap<float> > weights;
    iEvent.getByLabel(edm::InputTag("puppi", "PuppiWeights"), weights);
    assert(weights.isValid());
    const edm::ValueMap<float> puppi_weights = (*weights.product());
    std::cout<<"weight size!! "<<puppi_weights.size()<<std::endl;
	
    /* 1st Way */
    /*
    for (unsigned int i = 0; i <  puppi->size(); ++i) {
      reco::PFCandidateRef myRef(puppi,i);
      std::cout<<"weight!" << (*weights)[myRef]<<std::endl;

    }
    */

    /* 2nd Way */
    edm::Handle< edm::View<reco::PFCandidate> > pfs_handle;
    iEvent.getByLabel("puppi", pfs_handle);
    assert(pfs_handle.isValid());

    for (unsigned int i = 0; i <  pfs_handle->size(); ++i) {
      const reco::PFCandidate &pf =  pfs_handle->at(i); //read candidate e.g pf.pt() etc
      //Make a reference or... 
      edm::RefToBase<reco::PFCandidate>  pf_base_ref;
      pf_base_ref = pfs_handle->refAt(i);
      std::cout<<"weight!" << (*weights)[pf_base_ref]<<std::endl;
      //...a pointer
      //edm::Ptr<reco::PFCandidate> pf_base_Ptr = pfs_handle->ptrAt(i);
      //std::cout<<"weight!" << (*weights)[pf_base_Ptr]<<std::endl;
    }

}

// ------------ method called once each job just before starting event loop  ------------
void 
ProducerTest::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ProducerTest::endJob() {
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ProducerTest::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


//template<typename Hand, typename T>
void
ProducerTest:: storeMuonIso(edm::Event &iEvent,
			    const edm::Handle<pat::MuonCollection > & handle,
			    const std::vector<double> & values,
			    const std::string    & label) const {

  
    using namespace edm; 
    using namespace std;
    std::cout<<handle->size()<<std::endl;
    std::cout<<values.size()<<std::endl;
    auto_ptr<ValueMap<double> > valMap(new ValueMap<double>());
    typename edm::ValueMap<double>::Filler filler(*valMap);

    filler.insert(handle, values.begin(), values.end());
    filler.fill();
    iEvent.put(valMap, label);
}
/*
void
ProducerTest:: storeElectronIso(edm::Event &iEvent,
				const edm::Handle<pat::ElectronCollection > & handle,
				const std::vector<double> & values,
				const std::string    & label) const {
  
    using namespace edm; 
    using namespace std;
    auto_ptr<ValueMap<double> > valMap(new ValueMap<double>());
    typename edm::ValueMap<double>::Filler filler(*valMap);
    filler.insert(handle, values.begin(), values.end());
    filler.fill();
    iEvent.put(valMap, label);
}
*/
//define this as a plug-in
DEFINE_FWK_MODULE(ProducerTest);
