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

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "TH1.h"
//
// class declaration
//

class PUPPILeptonIsoProducer_Matthias_readNOWeights : public edm::EDProducer
{
    public:
        explicit PUPPILeptonIsoProducer_Matthias_readNOWeights(const edm::ParameterSet&);
        ~PUPPILeptonIsoProducer_Matthias_readNOWeights();

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


       //Debugging
       edm::Service<TFileService> fs_;
       TH1F* h_sCHPt;
       TH1F* h_sNHEt_readNOWeights;
       TH1F* h_sPHEt_readNOWeights;
       TH1F* h_pfid ;
       TH1F* h_weights_NH, * h_weights_NH_check;
       TH1F* h_weights_PH, * h_weights_PH_check;
       TH1F* h_sNHEt, * h_sNHEt_weighted;
       TH1F* h_sPHEt, *  h_sPHEt_weighted;
       TH1F* h_sCHPt_MuonPFIsostruct;
       TH1F* h_sCPPt_MuonPFIsostruct;
       TH1F* h_sNHEt_MuonPFIsostruct;
       TH1F* h_sPHEt_MuonPFIsostruct;


};


//
// constructors and destructor
//
PUPPILeptonIsoProducer_Matthias_readNOWeights::PUPPILeptonIsoProducer_Matthias_readNOWeights(const edm::ParameterSet& iConfig):
    _muonToken(consumes<edm::View<reco::RecoCandidate> >(iConfig.getParameter<edm::InputTag>("muons"))),
    _electronToken(consumes<edm::View<reco::RecoCandidate> >(iConfig.getParameter<edm::InputTag>("electrons"))),
    _pfCandidatesToken(consumes<edm::View<reco::Candidate>>(iConfig.getParameter<edm::InputTag>("pfCands"))),
    _puppiToken(consumes<edm::View<reco::Candidate> >(iConfig.getParameter<edm::InputTag>("puppi"))),
    fs_() //Debugging

{
    _dRConeSize  = iConfig.getUntrackedParameter<double>("dRConeSize");

    produces<edm::ValueMap<double> > ();
    
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

     //Debugging                                                                                                                      
    h_sCHPt = fs_->make<TH1F>("sCHPt_wL", "sCHPt",200, 0. , 100.);
    h_sNHEt_readNOWeights = fs_->make<TH1F>("sNHEt_readNOWeights", "sNHEt_readNOWeights",200, 0. , 100.);
    h_sPHEt_readNOWeights = fs_->make<TH1F>("sPHEt_readNOWeights", "sPHEt_readNOWeights",200, 0. , 100.);
    h_pfid  = fs_->make<TH1F>("pfID", "pfId", 2*2212, 0., 2213.);
    h_weights_NH = fs_->make<TH1F>("weights_NH", "weights_NH", 200, 0., 10.);
    h_weights_NH_check = fs_->make<TH1F>("weights_NH_check", "weights_NH_check", 200, 0., 10.);
    h_weights_PH = fs_->make<TH1F>("weights_PH", "weights_PH", 200, 0., 10.);
    h_weights_PH_check = fs_->make<TH1F>("weights_PH_check", "weights_PH_check", 200, 0., 10.);
    h_sNHEt = fs_->make<TH1F>("sNHEt", "sNHEt",200, 0. , 100.);
    h_sNHEt_weighted = fs_->make<TH1F>("sNHEt_weighted", "sNHEt_weighted",200, 0. , 100.);
    h_sPHEt = fs_->make<TH1F>("sPHEt", "sPHEt",200, 0. , 100.);
    h_sPHEt_weighted = fs_->make<TH1F>("sPHEt_weighted", "sPHEt_weighted",200, 0. , 100.);
    h_sCHPt_MuonPFIsostruct = fs_->make<TH1F>("sCHPt_MuonPFIsostruct", "sCHPtMuonPFIsostruct",200, 0. , 100.);
    h_sCPPt_MuonPFIsostruct = fs_->make<TH1F>("sCHPtMuonPFIsostruct", "sCHPtMuonPFIsostruct",200, 0. , 100.);
    h_sNHEt_MuonPFIsostruct = fs_->make<TH1F>("sNHEtMuonPFIsostruct", "sNHEtMuonPFIsostruct",200, 0. , 100.);
    h_sPHEt_MuonPFIsostruct = fs_->make<TH1F>("sPHEtMuonPFIsostruct", "sPHEtMuonPFIsostruct",200, 0. , 100.);

}


PUPPILeptonIsoProducer_Matthias_readNOWeights::~PUPPILeptonIsoProducer_Matthias_readNOWeights()
{
}


void
PUPPILeptonIsoProducer_Matthias_readNOWeights::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    edm::Handle<edm::View<reco::RecoCandidate> > muons;
    iEvent.getByToken(_muonToken, muons);

    edm::Handle<edm::View<reco::RecoCandidate> > electrons;
    iEvent.getByToken(_electronToken, electrons);

    edm::Handle<edm::View<reco::Candidate> > pfCandidates;
    iEvent.getByToken(_pfCandidatesToken,pfCandidates);

    edm::Handle<edm::View<reco::Candidate> > puppipfCandidates;
    iEvent.getByToken(_puppiToken, puppipfCandidates);

    //Debugging
    edm::Handle<edm::ValueMap<float> > weights;
    iEvent.getByLabel("puppi", "PuppiWeights", weights);
    
        
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
        
        std::vector<reco::CandidatePtr> footprint;
	
	for (unsigned int i = 0, n = lepton.numberOfSourceCandidatePtrs(); i < n; ++i)
	{
	    footprint.push_back(lepton.sourceCandidatePtr(i));
	}

	double charged = 0, neutral = 0, photons  = 0;
	
	//Debugging
	double charged_wght = 0, neutral_wght = 0, photons_wght  = 0;
	double neutral_nowght = 0, photons_nowght  = 0;

	for (unsigned int ipf = 0; ipf <  pfCandidates->size(); ++ipf)
	{
	      const pat::PackedCandidate *pf =  dynamic_cast<const pat::PackedCandidate*>(&pfCandidates->at(ipf));
	      if (deltaR(*pf,lepton) < _dRConeSize)
              { 
	       
		  if (std::find(footprint.begin(), footprint.end(), reco::CandidatePtr(pfCandidates,ipf)) != footprint.end())
		  {
		      continue;
		  }
		  
		  //Debugging
		  if (pf->charge() == 0)
		  {
		    
		      float weight = (*weights)[pfCandidates->refAt(ipf)];
		      
		      if (pf->pdgId() == 22 && pf->pt() > 0.5)
		      {
			  h_weights_PH_check->Fill(weight);
			  photons_wght += weight*pf->pt();
			  photons_nowght += pf->pt();
		      }
		      else if (pf->pdgId() != 22 && pf->pt() > 0.5)
		      {
			  h_weights_NH_check->Fill(weight);
			  neutral_wght += weight*pf->pt();
			  neutral_nowght += pf->pt();
		      }
		  }

		  if (abs(pf->charge()!=0) && pf->fromPV()>=2)
		  {
		      if (!_includeLeptoninIso)
		      {
			  if (not (abs(pf->pdgId())==11 or abs(pf->pdgId())==13) ) charged += pf->pt();
		      }
		      else
		      {
			  charged_wght += pf->pt();
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
	      h_sCHPt_MuonPFIsostruct->Fill(muon->pfIsolationR03().sumChargedHadronPt);
	      h_sCPPt_MuonPFIsostruct->Fill(muon->pfIsolationR03().sumChargedParticlePt);
	      h_sNHEt_MuonPFIsostruct->Fill(muon->pfIsolationR03().sumNeutralHadronEt);
              h_sPHEt_MuonPFIsostruct->Fill(muon->pfIsolationR03().sumPhotonEt);
	      h_sCHPt->Fill(charged);
	      h_sNHEt_readNOWeights->Fill(neutral);
	      h_sPHEt_readNOWeights->Fill(photons);
	      h_weights_PH->Fill(photons_wght/photons_nowght);
	      h_weights_NH->Fill(neutral_wght/neutral_nowght);
	      h_sNHEt_weighted->Fill(neutral_wght);
	      h_sPHEt_weighted->Fill(photons_wght);
	      h_sNHEt->Fill(neutral_nowght);
	      h_sPHEt->Fill(photons_nowght);
	  }
    }
    storeLeptonIsoInfo(iEvent, muons, leptonIsolation, ""); //name?
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
PUPPILeptonIsoProducer_Matthias_readNOWeights::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PUPPILeptonIsoProducer_Matthias_readNOWeights::endJob() {
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PUPPILeptonIsoProducer_Matthias_readNOWeights::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

template<class Hand, typename T>
void
PUPPILeptonIsoProducer_Matthias_readNOWeights:: storeLeptonIsoInfo(edm::Event &iEvent,
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
DEFINE_FWK_MODULE(PUPPILeptonIsoProducer_Matthias_readNOWeights);
