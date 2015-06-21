// -*- C++ -*-
//
// Package:    test/PUPPIMuonIsoProducer
// Class:      PUPPIMuonIsoProducer
// 
/**\class PUPPIMuonIsoProducer PUPPIMuonIsoProducer.cc test/PUPPIMuonIsoProducer/plugins/PUPPIMuonIsoProducer.cc
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
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "DataFormats/Common/interface/ValueMap.h"

#include <functional>

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "TH1.h"
//
// class declaration
//

class PUPPIMuonIsoProducer : public edm::EDProducer
{
    public:
        explicit PUPPIMuonIsoProducer(const edm::ParameterSet&);
        ~PUPPIMuonIsoProducer();

        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    private:
        virtual void beginJob() override;
        virtual void produce(edm::Event&, const edm::EventSetup&) override;
        virtual void endJob() override;

        template<class Hand, typename T>
	void storeLeptonIsoInfo(edm::Event &iEvent,
				const edm::Handle<Hand > & handle, 
				const std::vector<T> & values,
				const std::string    & label) const ;

        edm::EDGetTokenT< edm::View<reco::RecoCandidate> > _muonToken;
        edm::EDGetTokenT< edm::View<reco::RecoCandidate> > _electronToken;
        edm::EDGetTokenT<edm::View<reco::Candidate> > _pfCandidatesToken;
        edm::EDGetTokenT<edm::View<reco::Candidate> > _puppiToken;

        float _dRConeSize;
 
        bool _writeCandidateSums;
        bool _includeLeptoninIso;

};


//
// constructors and destructor
//
PUPPIMuonIsoProducer::PUPPIMuonIsoProducer(const edm::ParameterSet& iConfig):
    _muonToken(consumes<edm::View<reco::RecoCandidate> >(iConfig.getParameter<edm::InputTag>("muons"))),
    _electronToken(consumes<edm::View<reco::RecoCandidate> >(iConfig.getParameter<edm::InputTag>("electrons"))),
    _pfCandidatesToken(consumes<edm::View<reco::Candidate>>(iConfig.getParameter<edm::InputTag>("pfCands"))),
    _puppiToken(consumes<edm::View<reco::Candidate> >(iConfig.getParameter<edm::InputTag>("puppi")))

{
    _dRConeSize  = iConfig.getUntrackedParameter<double>("dRConeSize");

    produces<edm::ValueMap<double> > ("LeptonPuppiIso");


    _writeCandidateSums = iConfig.getUntrackedParameter<bool>("writeCandidateSums");
    if (_writeCandidateSums) 
    {
        _includeLeptoninIso = iConfig.getUntrackedParameter<bool>("includeLeptoninIso");
        if (_includeLeptoninIso) 
        {
            produces<edm::ValueMap<double> > ("sumChargedParticlePt");
        }
	else
	{
	    produces<edm::ValueMap<double> > ("sumChargedHadronPt");
	}
        produces<edm::ValueMap<double> > ("sumNeutralHadronEt");
	produces<edm::ValueMap<double> > ("sumPhotonEt");
    }


}


PUPPIMuonIsoProducer::~PUPPIMuonIsoProducer()
{
}


void
PUPPIMuonIsoProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  
    edm::Handle<edm::View<reco::RecoCandidate> > muons;
    iEvent.getByToken(_muonToken, muons);

    edm::Handle<edm::View<reco::RecoCandidate> > electrons;
    iEvent.getByToken(_electronToken, electrons);

    edm::Handle<edm::View<reco::Candidate> > pfCandidates;
    iEvent.getByToken(_pfCandidatesToken,pfCandidates);

    edm::Handle<edm::View<reco::Candidate> > puppipfCandidates;
    iEvent.getByToken(_puppiToken, puppipfCandidates);


    std::vector<double>  leptonIsolation(muons->size());
    std::vector<double>  sumChargedCandidatePtInPuppiIso(muons->size());
    std::vector<double>  sumNeutralHadronPtInPuppiIso(muons->size());
    std::vector<double>  sumPhotonEtInPuppiIso(muons->size());
    
    std::vector<const reco::RecoCandidate *> leptons;
    for (const reco::RecoCandidate &mu : *muons) leptons.push_back(&mu);
    for (const reco::RecoCandidate &el : *electrons) leptons.push_back(&el);

    for (unsigned int ilepton = 0; ilepton < leptons.size(); ++ilepton)
    {
	const reco::RecoCandidate& lepton = (*leptons[ilepton]);

	// initialize sums
        double charged = 0, neutral = 0, photons  = 0;
	
	// now get a list of the PF candidates used to build this lepton, so to exclude them
	std::vector<reco::CandidatePtr> footprint;
	for (unsigned int i = 0, n = lepton.numberOfSourceCandidatePtrs(); i < n; ++i)
	{
	    footprint.push_back(lepton.sourceCandidatePtr(i));
	}
	
	// now loop on PF-charged candidates
	for (unsigned int ipf = 0; ipf <  pfCandidates->size(); ++ipf)
	{
	    const pat::PackedCandidate *pf =  dynamic_cast<const pat::PackedCandidate*>(&pfCandidates->at(ipf));
	    if (deltaR(*pf,lepton) < _dRConeSize)
            { 
	       
	        // pfcandidate-based footprint removal
	        if (std::find(footprint.begin(), footprint.end(), reco::CandidatePtr(pfCandidates,ipf)) != footprint.end())
		{
		    continue;
		}
		  
		if (abs(pf->charge()!=0) && pf->fromPV()>=2)
		{
		    if (!_includeLeptoninIso)
		    {
			if (not (abs(pf->pdgId())==11 or abs(pf->pdgId())==13) ) charged += pf->pt();
		    }
		    else
		    {
			charged += pf->pt();
		    }
		}
	    }
	}
	
        for (unsigned int ipf = 0; ipf <  puppipfCandidates->size(); ++ipf)
        {
	    const reco::Candidate& pf = puppipfCandidates->at(ipf);
	    
            if (deltaR(pf,lepton) < _dRConeSize)
            { 
	      
		if (pf.charge() == 0)
                {
		
 		    if (pf.pdgId() == 22 && pf.pt() > 0.5)
                    {
                        photons += pf.pt();
                    }
                    else if (pf.pdgId() != 22 && pf.pt() > 0.5)
                    {
                        neutral += pf.pt();
                    }
		}
	    }	
        }
	
	if (abs(lepton.pdgId())==13)
	  {
	    const reco::Muon* muon = dynamic_cast<const reco::Muon*>(&lepton); 
	    leptonIsolation[ilepton]=(charged + neutral + photons)/muon->pt();
	    sumChargedCandidatePtInPuppiIso[ilepton]=charged;
	    sumNeutralHadronPtInPuppiIso[ilepton]=neutral;
	    sumPhotonEtInPuppiIso[ilepton]=photons;
	  }

    }
    storeLeptonIsoInfo(iEvent, muons, leptonIsolation, "LeptonPuppiIso"); 
    if (_writeCandidateSums)
    {
        if (_includeLeptoninIso)
        {
	    storeLeptonIsoInfo(iEvent, muons, sumChargedCandidatePtInPuppiIso,  "sumChargedParticlePt");
        }
	else
	{
	    storeLeptonIsoInfo(iEvent, muons, sumChargedCandidatePtInPuppiIso,  "sumChargedHadronPt");
        }
        storeLeptonIsoInfo(iEvent, muons, sumNeutralHadronPtInPuppiIso,  "sumNeutralHadronEt");
        storeLeptonIsoInfo(iEvent, muons, sumPhotonEtInPuppiIso,  "sumPhotonEt");
    }

}

// ------------ method called once each job just before starting event loop  ------------
void 
PUPPIMuonIsoProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PUPPIMuonIsoProducer::endJob() {
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PUPPIMuonIsoProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

template<class Hand, typename T>
void
PUPPIMuonIsoProducer:: storeLeptonIsoInfo(edm::Event &iEvent,
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
DEFINE_FWK_MODULE(PUPPIMuonIsoProducer);
