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

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "TH1.h"
#include "TH2.h"

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
      void storeLeptonIsoInfo(edm::Event &iEvent,
			      const edm::Handle<Hand > & handle, 
			      const std::vector<T> & values,
			      const std::string    & label) const ;


      // ----------member data ---------------------------
      edm::EDGetTokenT< edm::View<reco::RecoCandidate> > muonToken_;
      edm::EDGetTokenT< edm::View<reco::RecoCandidate> > electronToken_;
      typedef edm::View<reco::Candidate> candidateView_;
      edm::EDGetTokenT< candidateView_ > pFCandidatesToken_;
      edm::EDGetTokenT< candidateView_ > pfWeightedNeutralHadronsToken_;
      edm::EDGetTokenT< candidateView_ > pfWeightedPhotonsToken_;
   
      float _dRConeSize;

      bool _writeCandidateSums;
      bool _includeLeptoninIso;

      //Debugging
      edm::Service<TFileService> fs_;
      TH1F* h_sCHPt;
      TH1F* h_sCHPt_MuonPFIsostruct;
      TH1F* h_sCPPt_MuonPFIsostruct;
      TH1F* h_sNHEt;
      TH1F* h_sPHEt;
      TH2F* h_pfid ;
};


//
// constructors and destructor
//
PFWeightedLeptonIsoProducer_new::PFWeightedLeptonIsoProducer_new(const edm::ParameterSet& iConfig):
    muonToken_(consumes<edm::View<reco::RecoCandidate> >(iConfig.getParameter<edm::InputTag>("muons"))),
    electronToken_(consumes<edm::View<reco::RecoCandidate> >(iConfig.getParameter<edm::InputTag>("electrons"))),
    pFCandidatesToken_(consumes<candidateView_>(iConfig.getParameter<edm::InputTag>("pfCands"))),
    pfWeightedNeutralHadronsToken_(consumes<candidateView_>(iConfig.getParameter<edm::InputTag>("pfWeightedHadrons"))),
    pfWeightedPhotonsToken_(consumes<candidateView_>(iConfig.getParameter<edm::InputTag>("pfWeightedPhotons"))), 
    fs_() //Debugging

{

    _dRConeSize  = iConfig.getUntrackedParameter<double>("dRConeSize");
    
    produces<edm::ValueMap<double> > ("LeptonPFWeightedIso");
    
    _writeCandidateSums = iConfig.getUntrackedParameter<bool>("writeCandidateSums");
    if (_writeCandidateSums)
    {
        produces<edm::ValueMap<double> > ("sumChargedCandidatePt");
        produces<edm::ValueMap<double> > ("sumNeutralHadronEt");
        produces<edm::ValueMap<double> > ("sumPhotonEt");

    }

    _includeLeptoninIso = iConfig.getUntrackedParameter<bool>("includeLeptoninIso");

    //Debugging                                                                                                                              
    h_sCHPt = fs_->make<TH1F>("sCHPt_wL", "sCHPt_wL",200, 0. , 100.);
    h_sCHPt_MuonPFIsostruct = fs_->make<TH1F>("sCHPt_Isostruct", "sCHPt_Isostruct",200, 0. , 100.);
    h_sCPPt_MuonPFIsostruct = fs_->make<TH1F>("sCPPt_Isostruct", "sCPPt_Isostruct",200, 0. , 100.);
    h_sNHEt = fs_->make<TH1F>("sNHEt_wL", "sNHEt_wL",200, 0. , 100.);
    h_sPHEt = fs_->make<TH1F>("sPHEt_wL", "sPHEt_wL",200, 0. , 100.);
    h_pfid  = fs_->make<TH2F>("pfID", "pfId", 2*2212, 0., 2213.,200, 0. , 100. );
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
    edm::Handle<edm::View<reco::RecoCandidate> > muons;
    iEvent.getByToken(muonToken_, muons);
    edm::Handle<edm::View<reco::RecoCandidate> > electrons;
    iEvent.getByToken(electronToken_, electrons);
    edm::Handle<candidateView_> pfCharged;
    iEvent.getByToken(pFCandidatesToken_,pfCharged);
    edm::Handle<candidateView_> pfNU;
    iEvent.getByToken(pfWeightedNeutralHadronsToken_,pfNU);
    edm::Handle<candidateView_> pfPH;
    iEvent.getByToken(pfWeightedPhotonsToken_,pfPH);

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
	//const reco::RecoCandidate& lepton = (leptons->at(ilepton)); //in case 'leptons' is handled as collection 
 
        // initialize sums
        float charged = 0, neutral = 0, photons  = 0;

        // now get a list of the PF candidates used to build this lepton, so to exclude them
	std::vector<reco::CandidatePtr> footprint;
	for (unsigned int i = 0, n = lepton.numberOfSourceCandidatePtrs(); i < n; ++i) 
	{
	    footprint.push_back(lepton.sourceCandidatePtr(i));
	}

	// now loop on PF-charged candidates
	for (unsigned int ipfCh = 0; ipfCh <  pfCharged->size(); ++ipfCh)
        {
            const pat::PackedCandidate *pf =  dynamic_cast<const pat::PackedCandidate*>(&pfCharged->at(ipfCh));
            if (deltaR(*pf,lepton) < _dRConeSize)
            {

	        // pfcandidate-based footprint removal                                                                               
                if (std::find(footprint.begin(), footprint.end(), reco::CandidatePtr(pfCharged,ipfCh)) != footprint.end())
                {
                    continue;
                }
	
                if (abs(pf->charge()) > 0 && pf->fromPV()>=2)
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
	// now loop on PF-weighted neutral candidates
	for (unsigned int ipfNU = 0; ipfNU <  pfNU->size(); ++ipfNU) 
	{
	    const reco::Candidate& pf = pfNU->at(ipfNU);
	    if (deltaR(pf,lepton) < _dRConeSize) 
	    { 
	        neutral += pf.pt();
	    }
	}

	// now loop on PF-weighted photon candidates
	for (unsigned int ipfPH = 0; ipfPH <  pfPH->size(); ++ipfPH) 
	{
	    const reco::Candidate& pf = pfPH->at(ipfPH);
	    if (deltaR(pf,lepton) < _dRConeSize) 
	    { 
	        photons += pf.pt();
	    }
	}
	
	if (abs(lepton.pdgId())==13)
	{
	    const reco::Muon* muon = dynamic_cast<const reco::Muon*>(&lepton); 
	    leptonIsolation[ilepton]=(charged + neutral + photons)/muon->pt();
	    sumChargedCandidatePtInPuppiIso[ilepton]=charged;
	    sumNeutralHadronPtInPuppiIso[ilepton]=neutral;
	    sumPhotonEtInPuppiIso[ilepton]=photons;
	    h_sCHPt->Fill(charged);
	    h_sCHPt_MuonPFIsostruct->Fill(muon->pfIsolationR03().sumChargedHadronPt);
	    h_sCPPt_MuonPFIsostruct->Fill(muon->pfIsolationR03().sumChargedParticlePt);
	    h_sNHEt->Fill(neutral);
	    h_sPHEt->Fill(photons);
	}

    }
    storeLeptonIsoInfo(iEvent, muons, leptonIsolation, "LeptonPFWeightedIso"); 
    if (_writeCandidateSums)
    {
        storeLeptonIsoInfo(iEvent, muons, sumChargedCandidatePtInPuppiIso,  "sumChargedCandidatePt");
        storeLeptonIsoInfo(iEvent, muons, sumNeutralHadronPtInPuppiIso,  "sumNeutralHadronEt");
	storeLeptonIsoInfo(iEvent, muons, sumPhotonEtInPuppiIso,  "sumPhotonEt");

    }


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
PFWeightedLeptonIsoProducer_new:: storeLeptonIsoInfo(edm::Event &iEvent,
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
