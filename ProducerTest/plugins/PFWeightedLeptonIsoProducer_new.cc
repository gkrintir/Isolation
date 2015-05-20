// -*- C++ -*-
//
// Package:    test/PFWeightedLeptonIsoProducer
// Class:      PFWeightedLeptonIsoProducer
// 
/**\class PFWeightedLeptonIsoProducer PFWeightedLeptonIsoProducer.cc test/PFWeightedLeptonIsoProducer/plugins/PFWeightedLeptonIsoProducer.cc

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

class PFWeightedLeptonIsoProducer_new : public edm::EDProducer {
   public:
      explicit PFWeightedLeptonIsoProducer_new(const edm::ParameterSet&);
      ~PFWeightedLeptonIsoProducer_new();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
  
      template<class Hand, typename T>
      void storeLeptonIso(edm::Event &iEvent,
			const edm::Handle<Hand > & handle, 
			const std::vector<T> & values,
			const std::string    & label) const ;


      // ----------member data ---------------------------
      edm::EDGetTokenT< edm::View<reco::RecoCandidate> > leptonToken_;
      typedef edm::View<reco::Candidate> candidateView_;
      edm::EDGetTokenT< candidateView_ > pFCandidatesToken_;
      edm::EDGetTokenT< candidateView_ > pfWeightedNeutralHadronsToken_;
      edm::EDGetTokenT< candidateView_ > pfWeightedPhotonsToken_;
   
      float dRConeSize_;
  

};


//
// constructors and destructor
//
PFWeightedLeptonIsoProducer_new::PFWeightedLeptonIsoProducer_new(const edm::ParameterSet& iConfig):
    leptonToken_(consumes<edm::View<reco::RecoCandidate> >(iConfig.getParameter<edm::InputTag>("leptons"))),
    pFCandidatesToken_(consumes<candidateView_>(iConfig.getParameter<edm::InputTag>("pfCands"))),
    pfWeightedNeutralHadronsToken_(consumes<candidateView_>(iConfig.getParameter<edm::InputTag>("pfWeightedHadrons"))),
    pfWeightedPhotonsToken_(consumes<candidateView_>(iConfig.getParameter<edm::InputTag>("pfWeightedPhotons")))

{

  dRConeSize_  = iConfig.getUntrackedParameter<double>("dRConeSize");

  produces<edm::ValueMap<double> > ("LeptonPFWeightedIso");
  
}


PFWeightedLeptonIsoProducer_new::~PFWeightedLeptonIsoProducer_new()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
PFWeightedLeptonIsoProducer_new::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    //Handle Particle Collections for PF-weighted isolation
    edm::Handle<edm::View<reco::RecoCandidate> > leptons;
    iEvent.getByToken(leptonToken_, leptons);
    edm::Handle<candidateView_> pfCharged;
    iEvent.getByToken(pFCandidatesToken_,pfCharged);
    edm::Handle<candidateView_> pfNU;
    iEvent.getByToken(pfWeightedNeutralHadronsToken_,pfNU);
    edm::Handle<candidateView_> pfPH;
    iEvent.getByToken(pfWeightedPhotonsToken_,pfPH);

    std::vector<double>  leptonIsolation(leptons->size());

    for (unsigned int ilepton = 0; ilepton < leptons->size(); ++ilepton)
    {
        const reco::RecoCandidate& lepton = (leptons->at(ilepton));
 
        // initialize sums
        float charged = 0, neutral = 0, photons  = 0;

        // now get a list of the PF candidates used to build this lepton, so to exclude them
	std::vector<reco::CandidatePtr> footprint;
	for (unsigned int i = 0, n = lepton.numberOfSourceCandidatePtrs(); i < n; ++i) {
	  footprint.push_back(lepton.sourceCandidatePtr(i));
	}

	// now loop on pf charged candidates
	
	for (unsigned int ipfCh = 0; ipfCh <  pfCharged->size(); ++ipfCh) {
	  const pat::PackedCandidate *pf =  dynamic_cast<const pat::PackedCandidate*>(&pfCharged->at(ipfCh));
	  if (deltaR(*pf,lepton) < dRConeSize_) { 

	    // pfcandidate-based footprint removal
	    if (std::find(footprint.begin(), footprint.end(), reco::CandidatePtr(pfCharged,ipfCh)) != footprint.end()) {
	      continue;
	    }

	    if (pf->charge() != 0 && pf->fromPV()>2) charged += pf->pt();
	  }
	}
	
	// now loop on PF-weighted neutral candidates
	for (unsigned int ipfNU = 0; ipfNU <  pfNU->size(); ++ipfNU) {
	  const reco::Candidate& pf = pfNU->at(ipfNU);
	  if (deltaR(pf,lepton) < dRConeSize_) { 

	    // pfcandidate-based footprint removal
	    if (std::find(footprint.begin(), footprint.end(), reco::CandidatePtr(pfNU,ipfNU)) != footprint.end()) {
	      continue;
	    }
	    neutral += pf.pt();
	  }
	}

	// now loop on PF-weighted photon candidates
	for (unsigned int ipfPH = 0; ipfPH <  pfPH->size(); ++ipfPH) {
	  const reco::Candidate& pf = pfPH->at(ipfPH);
	  if (deltaR(pf,lepton) < dRConeSize_) { 
	    // pfcandidate-based footprint removal
	    if (std::find(footprint.begin(), footprint.end(), reco::CandidatePtr(pfPH,ipfPH)) != footprint.end()) {
	      continue;
	    }
	    photons += pf.pt();
	  }
	}
	
	
	leptonIsolation[ilepton]=(charged + neutral + photons)/lepton.pt();
    }
    storeLeptonIso(iEvent, leptons, leptonIsolation, "LeptonPFWeightedIso"); 

}

// ------------ method called once each job just before starting event loop  ------------
void 
PFWeightedLeptonIsoProducer_new::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PFWeightedLeptonIsoProducer_new::endJob() {
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PFWeightedLeptonIsoProducer_new::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


template<class Hand, typename T>
void
PFWeightedLeptonIsoProducer_new:: storeLeptonIso(edm::Event &iEvent,
			    const edm::Handle<Hand > & handle,
			    const std::vector<T> & values,
			    const std::string    & label) const {

  
    std::auto_ptr<edm::ValueMap<T> > valMap(new edm::ValueMap<T>());
    typename edm::ValueMap<T>::Filler filler(*valMap);

    filler.insert(handle, values.begin(), values.end());
    filler.fill();
    iEvent.put(valMap, label);

}

//define this as a plug-in
DEFINE_FWK_MODULE(PFWeightedLeptonIsoProducer_new);
