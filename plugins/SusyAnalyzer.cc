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

#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/METReco/interface/MET.h"
#include <iostream>
#include <DataFormats/PatCandidates/interface/Photon.h>
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/Math/interface/deltaR.h"
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

        edm::InputTag slimmedPhotons;
        edm::EDGetTokenT<edm::View<pat::Photon>>  slimmedPhotonsToken_;

	edm::InputTag PVertices;
	edm::EDGetTokenT<edm::View<reco::Vertex>> PVerticesToken_;

	edm::InputTag MetTag;
	edm::EDGetTokenT<edm::View<pat::MET>> MetTagToken_;

        TTree *MyTree;
        std::vector<float> PtJet;
        std::vector<float> EtaJet;
	std::vector<float> PhiJet;
	std::vector<float> EnJet;

        std::vector<float> PtElec;
        std::vector<float> EtaElec;
        std::vector<float> PhiElec;
 	std::vector<float> sieie;
	std::vector<int> convVeto;
	std::vector<int> mhits;
        std::vector<float> dEtaIn;
        std::vector<float> dPhiIn;
        std::vector<float> hoe;
        std::vector<float> ooemoop;
        std::vector<float> d0vtx;
        std::vector<float> dzvtx;                        
 //   	float hoe2;
	std::vector<float> ElecSieieBarrel;
	std::vector<float> ElecSieieEC;
	std::vector<int> Elec_prompt; 

        std::vector<float> PtMuon;
        std::vector<float> EtaMuon;
        std::vector<float> PhiMuon;
	std::vector<float> ChargeMu;
        float ChgIso;
	float ChgPU;
	float NeuIso;
	double dBIsoMu;

        std::vector<float> PtGenPart;
        std::vector<float> EtaGenPart;
        std::vector<float> PhiGenPart;


        std::vector<float> PtGamma;
        std::vector<float> EtaGamma;
        std::vector<float> PhiGamma;
        std::vector<float> goodPhotons;
        std::vector<float> photon_isEB;
        std::vector<float> photon_genMatched;
        std::vector<float> photon_hadTowOverEM;
        std::vector<float> photon_pfChargedIso;
        std::vector<float> photon_pfGammaIso;
        std::vector<float> photon_pfNeutralIso;
        std::vector<float> photon_hasPixelSeed;
	std::vector<float> PhSieieBarrel; //Barrel SigmaIetaIeta 
	std::vector<float> PhSieieEC;
	std::vector<int> photon_nonPrompt;
	std::vector<int> photon_electronFakes;	

	double metpt;
	double metphi;
	double meteta;
        double metLorentz;
        double metsign;
	double metEn;
	double metSig;
	
	
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
prunedGenParticlesToken_(consumes<edm::View<reco::GenParticle>>(prunedGenParticles)),

slimmedPhotons(iConfig.getUntrackedParameter<edm::InputTag>("slimmedPhotons")),
slimmedPhotonsToken_(consumes<edm::View<pat::Photon>>(slimmedPhotons)),

PVertices(iConfig.getUntrackedParameter<edm::InputTag>("PVertices")),
PVerticesToken_(consumes<edm::View<reco::Vertex>>(PVertices)),

MetTag(iConfig.getUntrackedParameter<edm::InputTag>("MetTag")),
MetTagToken_(consumes<edm::View<pat::MET>>(MetTag))
{
   //now do what ever initialization is needed

   usesResource("TFileService");
   edm::Service<TFileService> fs;
   MyTree = fs->make<TTree>("DYJet","DY");

   MyTree->Branch("PtJet", &PtJet);
   MyTree->Branch("EtaJet", &EtaJet);
   MyTree->Branch("PhiJet", &PhiJet);
   MyTree->Branch("EnJet", &EnJet);

   MyTree->Branch("PtElec", &PtElec);
   MyTree->Branch("EtaElec", &EtaElec);
   MyTree->Branch("PhiElec", &PhiElec);
   MyTree->Branch("sieie", &sieie);
   MyTree->Branch("convVeto",&convVeto);
   MyTree->Branch("mhits", &mhits);
   MyTree->Branch("dEtaIn", &dEtaIn);
   MyTree->Branch("dPhiIn", &dPhiIn);
   MyTree->Branch("hoe", &hoe);
   MyTree->Branch("ooemoop", &ooemoop);
   MyTree->Branch("d0vtx", &d0vtx);
   MyTree->Branch("dzvtx", &dzvtx);
//   MyTree->Branch("hoe2", &hoe2);
   MyTree->Branch("ElecSieieBarrel", &ElecSieieBarrel);
   MyTree->Branch("ElecSieieEC", &ElecSieieEC);
   MyTree->Branch("Prompt_Electron", &Elec_prompt);

   MyTree->Branch("PtMuon", &PtMuon);
   MyTree->Branch("EtaMuon", &EtaMuon);
   MyTree->Branch("PhiMuon", &PhiMuon);
   MyTree->Branch("Charge_of_muon", &ChargeMu);
   MyTree->Branch("ChargeIso", &ChgIso);
   MyTree->Branch("ChgPU_SumPUPt", &ChgPU);
   MyTree->Branch("NeuIso", &NeuIso);
   MyTree->Branch("dBIsoMu", &dBIsoMu);

   MyTree->Branch("PtGenPart", &PtGenPart);
   MyTree->Branch("EtaGenPart", &EtaGenPart);
   MyTree->Branch("PhiGenPart", &PhiGenPart);

   MyTree->Branch("PtGamma", &PtGamma);
   MyTree->Branch("EtaGamma", &EtaGamma);
   MyTree->Branch("PhiGamma", &PhiGamma);
   MyTree->Branch("goodPhotons", &goodPhotons);
   MyTree->Branch("photon_isEB", &photon_isEB);
   MyTree->Branch("photon_genMatched", &photon_genMatched);
   MyTree->Branch("photon_hadTowOverEM", &photon_hadTowOverEM);
   MyTree->Branch("photon_pfChargedIso", &photon_pfChargedIso);
   MyTree->Branch("photon_pfGammaIso", &photon_pfGammaIso);
   MyTree->Branch("photon_pfNeutralIso", &photon_pfNeutralIso);
   MyTree->Branch("photon_hasPixelSeed", &photon_hasPixelSeed);
   MyTree->Branch("Photon_Sigma_ieie_Barrel", &PhSieieBarrel);
   MyTree->Branch("Photon_Sigma_ieie_EP", &PhSieieEC);
   MyTree->Branch("photon_nonPrompt", &photon_nonPrompt);
   MyTree->Branch("photon_electronFakes", &photon_electronFakes);

   MyTree->Branch("metpt", &metpt);
   MyTree->Branch("meteta", &meteta);
   MyTree->Branch("metphi", &metphi);
   MyTree->Branch("metEn", &metEn);
   MyTree->Branch("metSignificance", &metSig);
 

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

	edm::Handle< View<reco::Vertex>> vertx;
	iEvent.getByToken(PVerticesToken_, vertx); 
      
         edm::Handle< View<reco::GenParticle>> iGen;
         iEvent.getByToken(prunedGenParticlesToken_,iGen);
 
        PtJet.clear(); EtaJet.clear(); PhiJet.clear(); EnJet.clear();

	//Jets
        edm::Handle< View<pat::Jet>> jetCands;
             iEvent.getByToken(JetToken_,jetCands);
             for(View<pat::Jet>::const_iterator iJet = jetCands->begin(); iJet != jetCands->end(); ++iJet){
                PtJet.push_back(iJet->pt());
                EtaJet.push_back(iJet->eta());
                PhiJet.push_back(iJet->phi());
		EnJet.push_back(iJet->energy());
             }

	//Electrons
        PtElec.clear(); EtaElec.clear(); PhiElec.clear(); sieie.clear(); convVeto.clear(); mhits.clear(); dEtaIn.clear(); dPhiIn.clear(); hoe.clear(); ooemoop.clear(); d0vtx.clear(); dzvtx.clear(); ElecSieieBarrel.clear(); ElecSieieEC.clear(); Elec_prompt.clear();
        edm::Handle< View<pat::Electron>> eCand;
             iEvent.getByToken(slimmedElectronsToken_,eCand);
             for(View<pat::Electron>::const_iterator i = eCand->begin(); i != eCand->end(); ++i){


                double ElecEta=i->eta();
                bool isBarrelElec=false; //to distinguish barrel and endcap electrons
                bool isEndcapElec=false;

                if(fabs(ElecEta) < 1.4442)
                        isBarrelElec=true;
                else if(fabs(ElecEta)>1.566 && fabs(ElecEta)<2.5)
                        isEndcapElec=true;
                else{
                      isBarrelElec=false;
                      isEndcapElec=false;
                }

                if(isBarrelElec){
                        ElecSieieBarrel.push_back(i->full5x5_sigmaIetaIeta());
                }else if(isEndcapElec){
                        ElecSieieEC.push_back(i->full5x5_sigmaIetaIeta());
                }
		
	      if (iGen.isValid()){//genLevel Stuff
        	// loop over gen particles and find nonprompt ELECTRONS
	                 int matchedGenPrompt = 0;
                         int matchedGenNonPrompt = 0;
			for(View<reco::GenParticle>::const_iterator j = iGen->begin(); j!=iGen->end(); ++j){//genparticle loop
				if( j->pdgId() == 11 && ( ( j->status() / 10 ) == 2 || j->status() == 1 || j->status() == 2) ){
					if(deltaR(j->p4(),i->p4()) < 0.1 && abs(j->mother()->pdgId())==23)matchedGenPrompt++;
					else matchedGenNonPrompt++;
				}//conditions for prompt not-prompt electron
			}//End of genparticle loop

		
		if(matchedGenPrompt > 0) Elec_prompt.push_back(true);
		else Elec_prompt.push_back(false);

        	}//end of genLevel Stuff

                PtElec.push_back(i->pt());
                EtaElec.push_back(i->eta());
                PhiElec.push_back(i->phi());
		sieie.push_back(i->full5x5_sigmaIetaIeta());
		convVeto.push_back(i->passConversionVeto());
		mhits.push_back(i->gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS));
		dEtaIn.push_back(i->deltaEtaSuperClusterTrackAtVtx());
		dPhiIn.push_back(i->deltaPhiSuperClusterTrackAtVtx());
		hoe.push_back(i->hadronicOverEm());
		ooemoop.push_back(fabs(1.0/i->ecalEnergy() - i->eSuperClusterOverP()/i->ecalEnergy()));
		d0vtx.push_back(i->gsfTrack()->dxy(vertx->at(0).position()));
		dzvtx.push_back(i->gsfTrack()->dz(vertx->at(0).position()));			
	      }
	//End for Electrons
	//
	//Muons
       PtMuon.clear(); EtaMuon.clear(); PhiMuon.clear(); ChargeMu.clear(); 
        edm::Handle< View<pat::Muon>> mCand;
             iEvent.getByToken(slimmedMuonsToken_,mCand);
             for(View<pat::Muon>::const_iterator i = mCand->begin(); i!=mCand->end(); ++i){
                PtMuon.push_back(i->pt());
                EtaMuon.push_back(i->eta());
                PhiMuon.push_back(i->phi());
		ChargeMu.push_back(i->charge());
		ChgIso=i->pfIsolationR04().sumChargedHadronPt;
		ChgPU=i->pfIsolationR04().sumPUPt;
		NeuIso=i->pfIsolationR04().sumNeutralHadronEt+i->pfIsolationR04().sumPhotonEt;
		dBIsoMu=(ChgIso+std::max(0., NeuIso-0.5*ChgPU))/i->pt();
             }

	//GenParticles
       PtGenPart.clear(); EtaGenPart.clear(); PhiGenPart.clear();
            for(View<reco::GenParticle>::const_iterator j = iGen->begin(); j!=iGen->end(); ++j){
                PtGenPart.push_back(j->pt());
                EtaGenPart.push_back(j->eta());
                PhiGenPart.push_back(j->phi());
             } //End for Genparticles
	
	//Photons
        PtGamma.clear(); EtaGamma.clear(); PhiGamma.clear(); goodPhotons.clear(); photon_isEB.clear(); photon_genMatched.clear(); photon_hadTowOverEM.clear(); photon_pfChargedIso.clear(); photon_pfGammaIso.clear(); photon_pfNeutralIso.clear(); photon_hasPixelSeed.clear(); PhSieieBarrel.clear(); PhSieieEC.clear(); photon_nonPrompt.clear(); photon_electronFakes.clear();

        edm::Handle< View<pat::Photon>> pCand;
             iEvent.getByToken(slimmedPhotonsToken_,pCand);
             for(View<pat::Photon>::const_iterator p = pCand->begin(); p != pCand->end(); ++p){

		double PhEta=p->eta();	
		bool isBarrelPhoton=false; //to distinguish barrel and endcap photons
		bool isEndcapPhoton=false;

		if(fabs(PhEta) < 1.4442)
			isBarrelPhoton=true;
		else if(fabs(PhEta)>1.566 && fabs(PhEta)<2.5)
			isEndcapPhoton=true;
		else{
		      isBarrelPhoton=false;
		      isEndcapPhoton=false;
		}

		if(isBarrelPhoton){
			PhSieieBarrel.push_back(p->full5x5_sigmaIetaIeta());		
		}else if(isEndcapPhoton){
			PhSieieEC.push_back(p->full5x5_sigmaIetaIeta());
		}


      if (iGen.isValid()){//genLevel Stuff
        // loop over gen particles and find nonprompt and hadronization photons
	        int matchedGenPrompt = 0;
                int matchedGenNonPrompt = 0 ;
                bool photonMatchGenE = false;	

		//pdgId()==set PDG identifier
		//deltaR==functions to compute deltaR

		for(View<reco::GenParticle>::const_iterator j = iGen->begin(); j!=iGen->end(); ++j){
			if( j->pdgId() == 22 && ( ( j->status() / 10 ) == 2 || j->status() == 1 || j->status() == 2 ) ){
			if( deltaR(j->p4(), p->p4()) < 0.2 ){
			if( abs(j->mother()->pdgId()) > 100 && abs(j->mother()->pdgId()) < 1000000 && abs(j->mother()->pdgId()) != 2212 ) matchedGenNonPrompt++ ;
			if( abs(j->mother()->pdgId()) <= 100 || abs(j->mother()->pdgId()) == 2212 ){
			if( j->pt()/p->pt() > 0.5 && j->pt()/p->pt() < 1.5 )
			  matchedGenPrompt++ ;
			}//for prompt photons
			}//gen matching
		}

		  if( abs(j->pdgId()) == 11 && j->status() == 1 && abs(j->mother()->pdgId()) <=25 ){
		    if( deltaR(j->p4(),p->p4()) < 0.2 && (j->pt()/p->pt() > 0.9 && j->pt()/j->pt() < 1.1) ){
		      		photonMatchGenE = true;
			    }
			}

	}//end of loop for GenParticle

	if( matchedGenPrompt > 0 || matchedGenNonPrompt == 0) photon_nonPrompt.push_back(false);
	else if( matchedGenNonPrompt > 0 ) photon_nonPrompt.push_back(true);
	else photon_nonPrompt.push_back(false);
	
	//check if photon is fake or not.
	if( photonMatchGenE )//make sure that photon is matched to gen electron and has similar pT as that of gen e.
	       photon_electronFakes.push_back(true);
	else
	       photon_electronFakes.push_back(false);

	}//GenLevelStuff

                PtGamma.push_back(p->pt());
                EtaGamma.push_back(p->eta());
		PhiGamma.push_back(p->phi());
	      goodPhotons.push_back( p->pt() );
	      photon_isEB.push_back( p->isEB() );
	     // photon_genMatchedi.push_back( p->genPhoton() );
	      photon_hadTowOverEM.push_back( p->hadTowOverEm() ) ;
	      photon_pfChargedIso.push_back(      p->chargedHadronIso() );
	      photon_pfGammaIso.push_back(        p->photonIso() );
	      photon_pfNeutralIso.push_back(      p->neutralHadronIso() );
	      photon_hasPixelSeed.push_back( p->hasPixelSeed() );
            }// End of loop for Photon

	//MET
        metpt=0; metphi=0; meteta=0; metEn=0; metSig=0; 
        edm::Handle< View<pat::MET>> MetCand;
	     iEvent.getByToken(MetTagToken_,MetCand);
             for(View<pat::MET>::const_iterator iMet = MetCand->begin(); iMet != MetCand->end(); ++iMet){
                metpt=	iMet->pt();
                meteta=iMet->eta();
                metphi=iMet->phi();
		metEn=iMet->energy();
		metSig=iMet->metSignificance();
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
