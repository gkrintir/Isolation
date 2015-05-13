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
  
      template<class Hand, typename T>
      void storeMuonIso(edm::Event &iEvent,
			const edm::Handle<Hand > & handle, 
			const std::vector<T> & values,
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
    
  //edm::EDGetTokenT<pat::PackedCandidateCollection> pfToken_;
      typedef edm::View<reco::Candidate> CandidateView;
      edm::EDGetTokenT< CandidateView > tokenPFCandidates_;
      std::vector<edm::EDGetTokenT<edm::ValueMap<float> > > tokenElectronIsoVals_;
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
    muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons")))
    //    pfToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCands"))) 

{
  tokenPFCandidates_= consumes<CandidateView>(iConfig.getParameter<edm::InputTag>("pfCands")); //e.g. 

  produces<edm::ValueMap<double> > ("MuonPuppiIso");
  //produces<edm::ValueMap<double> >("ElectronPuppiIso");

  
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
    //edm::Handle<pat::PackedCandidateCollection> pfs;
    // iEvent.getByToken(pfToken_, pfs);
    
    edm::Handle<CandidateView> pfs_;
    iEvent.getByToken(tokenPFCandidates_,pfs_);
    //    const pat::PackedCandidateCollection *pfs = dynamic_cast<const pat::PackedCandidateCollection*>(pfs_.product());
    const CandidateView *pfs = pfs_.product();
    
    edm::Handle<reco::PFCandidateCollection> puppi;
    iEvent.getByLabel("puppi", puppi);
    assert(puppi.isValid());
    //    const reco::PFCandidateCollection *PFCol = puppi.product();

    
    edm::Handle<edm::ValueMap<float> > weights;
    iEvent.getByLabel(edm::InputTag("puppi", "PuppiWeights"), weights);
    assert(weights.isValid());
    const edm::ValueMap<float> puppi_weights = (*weights.product());
    std::cout<<"weight size!! "<<puppi_weights.size()<<std::endl;
    std::cout<<"weight size!! "<<pfs->size()<<std::endl;	
   

    /* 2nd Way */
    /*
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
    */

    std::vector<const reco::Candidate *> leptons;

    for (const pat::Muon &mu : *muons) {
      std::cout<<"muons!!!!!!!!!!!!! "<<std::endl;
      leptons.push_back(&mu);
    }
    for (const pat::Electron &el : *electrons) leptons.push_back(&el);
    std::vector<double>  muon_isolation_;
    //std::vector<double>  electron_isolation(electrons->size(), -999);

    for (const reco::Candidate *lep : leptons) {
      //  if (lep->pt() < 5) continue;
        // initialize sums
        float charged = 0, neutral = 0, photons  = 0;

        // now get a list of the PF candidates used to build this lepton, so to exclude them
	std::vector<reco::CandidatePtr> footprint;
	for (unsigned int i = 0, n = lep->numberOfSourceCandidatePtrs(); i < n; ++i) {
	  footprint.push_back(lep->sourceCandidatePtr(i));
	}

	// now loop on pf candidates
	/* 1st Way */
    	//
	for (unsigned int i = 0; i <  pfs->size(); ++i) {
	  //for(CandidateView::const_iterator itPF = pfs->begin(); itPF!=pfs->end(); itPF++) {
	  //const pat::PackedCandidate* lPack = dynamic_cast<const pat::PackedCandidate*>(&(*itPF));
	  //if (lPack==0)continue;
	  //const pat::PackedCandidateRef myRef((*PackedCandidateCollection)pfs,i);
	  //std::cout<<lPack->charge()<<std::endl;
	  //std::cout<<"weight!" << (*weights)[myRef]<<std::endl;
	  const pat::PackedCandidate *pf =  dynamic_cast<const pat::PackedCandidate*>(&pfs->at(i));
	  //      const pat::PackedCandidateRef myRef(lPack->ref());
	  edm::RefToBase<reco::Candidate>  pf_base_ref;
	  pf_base_ref = pfs->refAt(i);
      
	  float weight = (*weights)[pf_base_ref];
	  if (deltaR(*pf,*lep) < 0.4) {
	    // pfcandidate-based footprint removal
	    if (std::find(footprint.begin(), footprint.end(), reco::CandidatePtr(pfs,i)) != footprint.end()) {
	      continue;
	    }
	    if (pf->charge() == 0) {
	      if (pf->pdgId() == 22) photons += weight*pf->pt();
	      else
		if (pf->pt() > 0.5) neutral += weight*pf->pt();
	    } else {
	      if (weight==1) charged += weight*pf->pt();
	      
	    }
	  }	
	}
  
	// do deltaBeta
	std::cout<<muons->size()<<std::endl;
	std::cout<<muon_isolation_.size()<<std::endl;
	double rel_iso = (charged + neutral + photons)/lep->pt();
	//std::cout<<"weight!" << iso<<std::endl;
	if (abs(lep->pdgId())==13) {
	  std::cout<<" gemizo!!!!!!!"<<std::endl;
	  muon_isolation_.push_back(rel_iso);//(iso);
	}
	//else if (lep->pdgId()==11) electron_isolation[1]=1;//.push_back(iso);
        printf("%-8s of pt %6.1f, eta %+4.2f: relIso = %5.2f\n",
	       abs(lep->pdgId())==13 ? "muon" : "electron",
	       lep->pt(), lep->eta(), rel_iso);
	
    }
    std::cout<<!muon_isolation_.empty()<<std::endl;
    std::cout<< muon_isolation_.capacity()<<" "<< muon_isolation_.capacity()<<std::endl;
    if(!muon_isolation_.empty())
      storeMuonIso(iEvent, muons, muon_isolation_,  "MuonPuppiIso");

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


template<class Hand, typename T>
void
ProducerTest:: storeMuonIso(edm::Event &iEvent,
			    const edm::Handle<Hand > & handle,
			    const std::vector<T> & values,
			    const std::string    & label) const {

  
    using namespace edm; 
    using namespace std;
    std::cout<<"func"<<handle->size()<<std::endl;
    std::cout<<"func"<<values.size()<<std::endl;
    std::auto_ptr<edm::ValueMap<T> > valMap(new edm::ValueMap<T>());
    typename edm::ValueMap<T>::Filler filler(*valMap);

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
