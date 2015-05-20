// -*- C++ -*-
//
// Package:    test/PUPPILeptonIsoProducer
// Class:      PUPPILeptonIsoProducer
// 
/**\class PUPPILeptonIsoProducer PUPPILeptonIsoProducer.cc test/PUPPILeptonIsoProducer/plugins/PUPPILeptonIsoProducer.cc
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

//
// class declaration
//

class PUPPILeptonIsoProducer_Matthias : public edm::EDProducer
{
    public:
        explicit PUPPILeptonIsoProducer_Matthias(const edm::ParameterSet&);
        ~PUPPILeptonIsoProducer_Matthias();

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

        edm::EDGetTokenT<edm::View<reco::RecoCandidate>> _leptonToken;
        edm::EDGetTokenT<edm::View<reco::Candidate>> _pfCandidatesToken;

        edm::EDGetTokenT<edm::ValueMap<float> > _puppiToken;

        float _dRConeSize;
 
        bool _writeCandidateSums;
        bool _includeLeptoninIso;

};


//
// constructors and destructor
//
PUPPILeptonIsoProducer_Matthias::PUPPILeptonIsoProducer_Matthias(const edm::ParameterSet& iConfig):
    _leptonToken(consumes<edm::View<reco::RecoCandidate>>(iConfig.getParameter<edm::InputTag>("leptons"))),
    _pfCandidatesToken(consumes<edm::View<reco::Candidate>>(iConfig.getParameter<edm::InputTag>("pfCands"))),
    _puppiToken(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("puppi")))

{
    _dRConeSize  = iConfig.getUntrackedParameter<double>("dRConeSize");

    produces<edm::ValueMap<double> > ();
    
    _writeCandidateSums = iConfig.getUntrackedParameter<bool>("writeCandidateSums");
    if (_writeCandidateSums) 
    {
        produces<edm::ValueMap<double> > ("sumChargedHadronPtInPuppiIso");
	produces<edm::ValueMap<double> > ("sumNeutralHadronEtInPuppiIso");
	produces<edm::ValueMap<double> > ("sumPhotonEtInPuppiIso");
    }

    _includeLeptoninIso = iConfig.getUntrackedParameter<bool>("includeLeptoninIso");
}


PUPPILeptonIsoProducer_Matthias::~PUPPILeptonIsoProducer_Matthias()
{
}


void
PUPPILeptonIsoProducer_Matthias::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    edm::Handle<edm::View<reco::RecoCandidate>> leptons;
    iEvent.getByToken(_leptonToken, leptons);
    
    edm::Handle<edm::View<reco::Candidate>> pfCandidates;
    iEvent.getByToken(_pfCandidatesToken,pfCandidates);

    edm::Handle<edm::ValueMap<float> > weights;
    iEvent.getByToken(_puppiToken, weights);


    std::vector<double>  leptonIsolation(leptons->size());
    std::vector<double>  sumChargedHadronPtInPuppiIso(leptons->size());
    std::vector<double>  sumNeutralHadronPtInPuppiIso(leptons->size());
    std::vector<double>  sumPhotonEtInPuppiIso(leptons->size());

    for (unsigned int ilepton = 0; ilepton < leptons->size(); ++ilepton)
    {
        const reco::RecoCandidate& lepton = (leptons->at(ilepton));
        
        std::vector<reco::CandidatePtr> footprint;
	for (unsigned int i = 0, n = lepton.numberOfSourceCandidatePtrs(); i < n; ++i)
	{
	    footprint.push_back(lepton.sourceCandidatePtr(i));
	}
	double charged = 0, neutral = 0, photons  = 0;
        
        for (unsigned int ipf = 0; ipf <  pfCandidates->size(); ++ipf)
        {
            const reco::Candidate& pf = pfCandidates->at(ipf);
            float weight = (*weights)[pfCandidates->refAt(ipf)];

            if (deltaR(pf,lepton) < _dRConeSize)
            { 
	        if (!_includeLeptoninIso)
	        {
                    if (std::find(footprint.begin(), footprint.end(), reco::CandidatePtr(pfCandidates,ipf)) != footprint.end())
		    {
		        continue;
		    }
		}

		if (pf.charge() == 0)
                {
                    if (pf.pdgId() == 22 && pf.pt() > 0.5)
                    {
                        photons += weight*pf.pt();
                    }
                    else if (pf.pdgId() != 22 && pf.pt() > 0.5)
                    {
                        neutral += weight*pf.pt();
                    }
		}
		else
		{
		      if (weight==1)
		      {
		          charged += weight*pf.pt();
		      }
		}
	    }	
        }
		
        leptonIsolation[ilepton]=(charged + neutral + photons)/lepton.pt();
	sumChargedHadronPtInPuppiIso[ilepton]=charged;
	sumNeutralHadronPtInPuppiIso[ilepton]=neutral;
	sumPhotonEtInPuppiIso[ilepton]=photons;
    }
    storeLeptonIsoInfo(iEvent, leptons, leptonIsolation, ""); //name?
    if (_writeCandidateSums)
    {
        storeLeptonIsoInfo(iEvent, leptons, sumChargedHadronPtInPuppiIso,  "sumChargedHadronPtInPuppiIso");
	storeLeptonIsoInfo(iEvent, leptons, sumNeutralHadronPtInPuppiIso,  "sumNeutralHadronEtInPuppiIso");
	storeLeptonIsoInfo(iEvent, leptons, sumPhotonEtInPuppiIso,  "sumPhotonEtInPuppiIso");
    }

}

// ------------ method called once each job just before starting event loop  ------------
void 
PUPPILeptonIsoProducer_Matthias::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PUPPILeptonIsoProducer_Matthias::endJob() {
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PUPPILeptonIsoProducer_Matthias::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

template<class Hand, typename T>
void
PUPPILeptonIsoProducer_Matthias:: storeLeptonIsoInfo(edm::Event &iEvent,
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
DEFINE_FWK_MODULE(PUPPILeptonIsoProducer_Matthias);
