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

// user include files

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
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

    std::vector<const reco::Candidate *> leptons;

    for (const pat::Muon &mu : *muons) leptons.push_back(&mu);
    for (const pat::Electron &el : *electrons) leptons.push_back(&el);
    std::vector<double>  muon_isolation(muons->size(), -999);
    std::vector<double>  electron_isolation(electrons->size(), -999);

    for (const reco::Candidate *lep : leptons) {
        if (lep->pt() < 5) continue;
        // initialize sums
        double charged = 0, neutral = 0, pileup  = 0;
        // now get a list of the PF candidates used to build this lepton, so to exclude them
        std::vector<reco::CandidatePtr> footprint;
        for (unsigned int i = 0, n = lep->numberOfSourceCandidatePtrs(); i < n; ++i) {
            footprint.push_back(lep->sourceCandidatePtr(i));
        }
        // now loop on pf candidates
        for (unsigned int i = 0, n = pfs->size(); i < n; ++i) {
            const pat::PackedCandidate &pf = (*pfs)[i];
            if (deltaR(pf,*lep) < 0.2) {
                // pfcandidate-based footprint removal
                if (std::find(footprint.begin(), footprint.end(), reco::CandidatePtr(pfs,i)) != footprint.end()) {
                    continue;
                }
                if (pf.charge() == 0) {
                    if (pf.pt() > 0.5) neutral += pf.pt();
                } else if (pf.fromPV() >= 2) {
                    charged += pf.pt();
                } else {
                    if (pf.pt() > 0.5) pileup += pf.pt();
                }
            }
        }

        // do deltaBeta
	//std::cout<<muons->size()<<std::endl;
	std::cout<<muon_isolation.size()<<std::endl;
	double iso = charged + std::max(0.0, neutral-0.5*pileup);
	if (lep->pdgId()==13) muon_isolation[1]=iso;//(iso);
	
	//else if (lep->pdgId()==11) electron_isolation[1]=1;//.push_back(iso);
        printf("%-8s of pt %6.1f, eta %+4.2f: relIso = %5.2f\n",
                    abs(lep->pdgId())==13 ? "muon" : "electron",
                    lep->pt(), lep->eta(), iso/lep->pt());

    }
    storeMuonIso(iEvent, muons, muon_isolation,  "MuonIso");
    //    storeElectronIso(iEvent, electrons, electron_isolation,  "ElectronIso");

    // Let's compute the fraction of charged pt from particles with dz < 0.1 cm
    for (const pat::Jet &j :  *jets) {
        if (j.pt() < 40 || fabs(j.eta()) > 2.4) continue;
        double in = 0, out = 0; 
        for (unsigned int id = 0, nd = j.numberOfDaughters(); id < nd; ++id) {
            const pat::PackedCandidate &dau = dynamic_cast<const pat::PackedCandidate &>(*j.daughter(id));
            if (dau.charge() == 0) continue;
            (fabs(dau.dz())<0.1 ? in : out) += dau.pt();
        }
        double sum = in + out;
        printf("Jet with pt %6.1f, eta %+4.2f, beta(0.1) = %+5.3f, pileup mva disc %+.2f\n",
                j.pt(),j.eta(), sum ? in/sum : 0, j.userFloat("pileupJetId:fullDiscriminant"));
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
