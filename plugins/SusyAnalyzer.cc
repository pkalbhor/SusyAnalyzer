// -*- C++ -*-
//
// Package:    SusyAnalyzer/SusyAnalyzer
// Class:      SusyAnalyzer
//
/**\class SusyAnalyzer SusyAnalyzer.cc SusyAnalyzer/SusyAnalyzer/plugins/SusyAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Pritam Kalbhor
//         Created:  Fri, 15 Jun 2018 22:39:50 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
 #include "FWCore/Utilities/interface/InputTag.h"
 #include "DataFormats/TrackReco/interface/Track.h"
 #include "DataFormats/TrackReco/interface/TrackFwd.h"

//Newly added, which are needed for this to work
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include <TTree.h>
#include <vector>

//For Jets and electron
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "PhysicsTools/PatAlgos/plugins/PATJetProducer.h"



//
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

    using namespace edm;
    using namespace std;
    using namespace reco;
    using namespace pat;


using reco::TrackCollection;

class SusyAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit SusyAnalyzer(const edm::ParameterSet&);
      ~SusyAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      //edm::EDGetTokenT<TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
        edm::InputTag jettag;
        edm::EDGetTokenT<edm::View<pat::Jet>>  JetToken_;

        edm::InputTag slimmedElectrons;
        edm::EDGetTokenT<edm::View<pat::Electron>>  slimmedElectronsToken_;

        edm::InputTag slimmedMuons;
        edm::EDGetTokenT<edm::View<pat::Muon>>  slimmedMuonsToken_;

        edm::InputTag prunedGenParticles;
        edm::EDGetTokenT<edm::View<reco::GenParticle>>  prunedGenParticlesToken_;

        TTree *MyTree;
        std::vector<float> PtJet;
        std::vector<float> EtaJet;

        std::vector<float> PtElec;
        std::vector<float> EtaElec;
      
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
SusyAnalyzer::SusyAnalyzer(const edm::ParameterSet& iConfig):
//  tracksToken_(consumes<TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tracks")))

jettag(iConfig.getUntrackedParameter<edm::InputTag>("jettag")),
JetToken_(consumes<edm::View<pat::Jet>>(jettag)),

slimmedElectrons(iConfig.getUntrackedParameter<edm::InputTag>("slimmedElectrons")),
slimmedElectronsToken_(consumes<edm::View<pat::Electron>>(slimmedElectrons)),

slimmedMuons(iConfig.getUntrackedParameter<edm::InputTag>("slimmedMuons")),
slimmedMuonsToken_(consumes<edm::View<pat::Muon>>(slimmedMuons)),

prunedGenParticles(iConfig.getUntrackedParameter<edm::InputTag>("prunedGenParticles")),
prunedGenParticlesToken_(consumes<edm::View<reco::GenParticle>>(prunedGenParticles))

{
   //now do what ever initialization is needed

   usesResource("TFileService");
   edm::Service<TFileService> fs;
   MyTree = fs->make<TTree>("DYJet","DY");

   MyTree->Branch("PtJet", &PtJet);
   MyTree->Branch("EtaJet", &EtaJet);

   MyTree->Branch("PtElec", &PtElec);
   MyTree->Branch("EtaElec", &EtaElec);

}


SusyAnalyzer::~SusyAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
SusyAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
/*
    Handle<TrackCollection> tracks;
    iEvent.getByToken(tracksToken_, tracks);
    for(TrackCollection::const_iterator itTrack = tracks->begin();
        itTrack != tracks->end();
        ++itTrack) {
      // do something with track parameters, e.g, plot the charge.
      // int charge = itTrack->charge();
    }
*/

        PtJet.clear(); EtaJet.clear();

        edm::Handle< View<pat::Jet>> jetCands;
             iEvent.getByToken(JetToken_,jetCands);
             for(View<pat::Jet>::const_iterator iJet = jetCands->begin(); iJet != jetCands->end(); ++iJet){
                std::cout << "jetpt: " << iJet->pt() << std::endl;
                std::cout << "jeteta" << iJet->eta() << std::endl;
                std::cout << "jetphi" << iJet->phi() <<std::endl;

                PtJet.push_back(iJet->pt());
                EtaJet.push_back(iJet->eta());
             }

        PtElec.clear(); EtaElec.clear();
        edm::Handle< View<pat::Electron>> eCand;
             iEvent.getByToken(slimmedElectronsToken_,eCand);
             for(View<pat::Electron>::const_iterator i = eCand->begin(); i != eCand->end(); ++i){
                PtElec.push_back(i->pt());
                EtaElec.push_back(i->eta());
             }

MyTree->Fill();

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
SusyAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
SusyAnalyzer::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SusyAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SusyAnalyzer);
