// system include files
#include <memory>
#include <string>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/HcalRecHit/interface/HBHERecHit.h"
#include "DataFormats/HcalRecHit/interface/HFRecHit.h"
#include "DataFormats/HcalRecHit/interface/HORecHit.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
/////#include <map>
#include<iostream>

//
// class declaration
//

class AlCaDiJetsProducer : public edm::EDProducer {
 public:
  explicit AlCaDiJetsProducer(const edm::ParameterSet&);
  ~AlCaDiJetsProducer();
  virtual void beginJob() ;
  virtual void produce(edm::Event &, const edm::EventSetup&);
  virtual void endJob();
 private:
  bool select (edm::Handle<reco::PFJetCollection> jt);
  // ----------member data ---------------------------
  
  edm::InputTag   labelPhoton_, labelPFJet_, labelHBHE_, labelHF_, labelHO_, labelTrigger_, labelPFCandidate_, labelVertex_, labelPFMET_, 
    labelPFMETtype1_, labelGsfEle_, labelRho_, labelConv_, labelBeamSpot_, labelLoosePhot_, labelTightPhot_;
  double          minPtJet_, minPtPhoton_;
  int             nAll_, nSelect_;
  
  //edm::EDGetTokenT<reco::PhotonCollection>                                                tok_Photon_; 
  edm::EDGetTokenT<reco::PFJetCollection>                                                 tok_PFJet_;
  edm::EDGetTokenT<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit>>> tok_HBHE_;
  edm::EDGetTokenT<edm::SortedCollection<HFRecHit,edm::StrictWeakOrdering<HFRecHit>>>     tok_HF_;
  edm::EDGetTokenT<edm::SortedCollection<HORecHit,edm::StrictWeakOrdering<HORecHit>>>     tok_HO_;
  //edm::EDGetTokenT<edm::TriggerResults>                                                   tok_TrigRes_;
  edm::EDGetTokenT<reco::PFCandidateCollection>                                           tok_PFCand_;
  edm::EDGetTokenT<reco::VertexCollection>                                                tok_Vertex_;
  //edm::EDGetTokenT<reco::PFMETCollection>                                                 tok_PFMET_;
  //edm::EDGetTokenT<reco::PFMETCollection>                                                 tok_PFType1MET_;
  //edm::EDGetTokenT<reco::GsfElectronCollection>                                           tok_GsfElec_;
  //edm::EDGetTokenT<double>                                                                tok_Rho_;
  //edm::EDGetTokenT<reco::ConversionCollection>                                            tok_Conv_;
  //edm::EDGetTokenT<reco::BeamSpot>                                                        tok_BS_;
  //edm::EDGetTokenT<edm::ValueMap<Bool_t> >                                                tok_loosePhoton_;
  //edm::EDGetTokenT<edm::ValueMap<Bool_t> >                                                tok_tightPhoton_;
};

AlCaDiJetsProducer::AlCaDiJetsProducer(const edm::ParameterSet& iConfig) : nAll_(0), nSelect_(0) {
  //std::cout << "Before labelPhoton_ = iConfig.getParameter<edm::InputTag>" << std::endl;
  // Take input 
  //labelPhoton_     = iConfig.getParameter<edm::InputTag>("PhoInput");
  labelPFJet_      = iConfig.getParameter<edm::InputTag>("PFjetInput");
  labelHBHE_       = iConfig.getParameter<edm::InputTag>("HBHEInput");
  labelHF_         = iConfig.getParameter<edm::InputTag>("HFInput");
  labelHO_         = iConfig.getParameter<edm::InputTag>("HOInput");
  //std::cout << "Before labelTrigger_    = edm::InputTag" << std::endl;
  //labelTrigger_    = edm::InputTag("TriggerResults::HLT");
  labelPFCandidate_= iConfig.getParameter<edm::InputTag>("particleFlowInput");
  labelVertex_     = iConfig.getParameter<edm::InputTag>("VertexInput");
  //labelPFMET_      = iConfig.getParameter<edm::InputTag>("METInput");
  //labelPFMETtype1_ = iConfig.getParameter<edm::InputTag>("Type1METInput");
  //labelGsfEle_     = iConfig.getParameter<edm::InputTag>("gsfeleInput");
  //labelRho_        = iConfig.getParameter<edm::InputTag>("rhoInput");
  //labelConv_       = iConfig.getParameter<edm::InputTag>("ConversionsInput");
  //labelBeamSpot_   = iConfig.getParameter<edm::InputTag>("BeamSpotInput");
  //std::cout <<"Before labelLoosePhot_  = edm::InputTag" << std::endl;
  //labelLoosePhot_  = edm::InputTag("PhotonIDProd", "PhotonCutBasedIDLoose");
  //labelTightPhot_  = edm::InputTag("PhotonIDProd", "PhotonCutBasedIDTight");
  minPtJet_        = iConfig.getParameter<double>("MinPtJet");
  //minPtPhoton_     = iConfig.getParameter<double>("MinPtPhoton");

  //tok_Photon_ = consumes<reco::PhotonCollection>(labelPhoton_);
  tok_PFJet_  = consumes<reco::PFJetCollection>(labelPFJet_);
  tok_HBHE_   = consumes<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit>>>(labelHBHE_);
  tok_HF_     = consumes<edm::SortedCollection<HFRecHit,edm::StrictWeakOrdering<HFRecHit>>>(labelHF_);
  tok_HO_     = consumes<edm::SortedCollection<HORecHit,edm::StrictWeakOrdering<HORecHit>>>(labelHO_);
  //tok_TrigRes_= consumes<edm::TriggerResults>(labelTrigger_);
  tok_PFCand_ = consumes<reco::PFCandidateCollection>(labelPFCandidate_);
  tok_Vertex_ = consumes<reco::VertexCollection>(labelVertex_);
  //tok_PFMET_  = consumes<reco::PFMETCollection>(labelPFMET_);
  //tok_PFType1MET_ = consumes<reco::PFMETCollection>(labelPFMETtype1_);
  //tok_loosePhoton_ = consumes<edm::ValueMap<Bool_t> >(labelLoosePhot_);
  //tok_tightPhoton_ = consumes<edm::ValueMap<Bool_t> >(labelTightPhot_);
  //tok_GsfElec_ = consumes<reco::GsfElectronCollection>(labelGsfEle_);
  //tok_Rho_ = consumes<double>(labelRho_);
  //tok_Conv_        = consumes<reco::ConversionCollection>(labelConv_);
  //tok_BS_          = consumes<reco::BeamSpot>(labelBeamSpot_);

  // register your products
  //std::cout << "Before produces<reco::PhotonCollection>(labelPhoton_.label())" << std::endl;
  //produces<reco::PhotonCollection>(labelPhoton_.label());
  produces<reco::PFJetCollection>(labelPFJet_.label());
  produces<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit>>>(labelHBHE_.label());
  produces<edm::SortedCollection<HFRecHit,edm::StrictWeakOrdering<HFRecHit>>>(labelHF_.label());
  produces<edm::SortedCollection<HORecHit,edm::StrictWeakOrdering<HORecHit>>>(labelHO_.label());
  //std::cout << "before produces<edm::TriggerResults>(labelTrigger_.label())" << std::endl; 
  //produces<edm::TriggerResults>(labelTrigger_.label());
  //produces<edm::ValueMap<Bool_t>>(labelLoosePhot_.label());
  //produces<edm::ValueMap<Bool_t>>(labelTightPhot_.label());
  //produces<double>(labelRho_.label());
  produces<reco::PFCandidateCollection>(labelPFCandidate_.label());
  produces<reco::VertexCollection>(labelVertex_.label());
  //produces<reco::PFMETCollection>(labelPFMET_.label());
  //produces<reco::PFMETCollection>(labelPFMETtype1_.label());
  //produces<reco::GsfElectronCollection>(labelGsfEle_.label());
  //produces<reco::ConversionCollection>(labelConv_.label());
  //std::cout << "before produces<reco::BeamSpot>(labelBeamSpot_.label())" << std::endl;
  //produces<reco::BeamSpot>(labelBeamSpot_.label());
  //std::cout << "All produces command completed" << std::endl;
}

AlCaDiJetsProducer::~AlCaDiJetsProducer()
{
  //std::cout << "Inside Destructor" << std::endl;
}

void AlCaDiJetsProducer::beginJob() {
  //std::cout << "Inside beginJob()" << std::endl;
}

void AlCaDiJetsProducer::endJob() {
  //std::cout << "Inside endJob()" << std::endl;
  edm::LogInfo("AlcaDiJets") << "Accepts " << nSelect_ << " events from a total of " << nAll_ << " events";
}

bool AlCaDiJetsProducer::select (edm::Handle<reco::PFJetCollection> jt) {
  //std::cout << "Inside select()" << std::endl;
  if ((((*jt)[0]).pt()) < minPtJet_ || (((*jt)[1]).pt()) < minPtJet_)    return false;
  return true;
}

// ------------ method called to produce the data  ------------
void AlCaDiJetsProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //std::cout << "Inside produce(---)" << std::endl;
  ++nAll_;
  // Access the collections from iEvent
  //std::cout << "Access PhotonCollection" << std::endl;

  //edm::Handle<reco::PhotonCollection> pho;
  //iEvent.getByToken(tok_Photon_,pho);
  //if (!pho.isValid()){
  //edm::LogWarning("AlCaDiJets") << "AlCaDiJetsProducer: Error! can't get photon product!";
  //return ;
  //}
  //const reco::PhotonCollection photon = *(pho.product());

  edm::Handle<reco::PFJetCollection> pfjet;
  iEvent.getByToken(tok_PFJet_,pfjet);
  if (!pfjet.isValid()){
    edm::LogWarning("AlCaDiJets") << "AlCaDiJetsProducer: Error! can't get pfjet product!";
    return ;
  }
  const reco::PFJetCollection pfjets = *(pfjet.product());

  edm::Handle<reco::PFCandidateCollection> pfc;
  iEvent.getByToken(tok_PFCand_,pfc);
  if (!pfc.isValid()){
    edm::LogWarning("AlCaDiJets") << "AlCaDiJetsProducer: Error! can't get pfcandidate product!";
    return ;
  }
  const reco::PFCandidateCollection pfcand = *(pfc.product());

  edm::Handle<reco::VertexCollection> vt;
  iEvent.getByToken(tok_Vertex_,vt);
  if (!vt.isValid()){
    edm::LogWarning("AlCaDiJets") << "AlCaDiJetsProducer: Error! can't get vertex product!";
    return ;
  }
  const reco::VertexCollection vtx = *(vt.product());

  /*edm::Handle<reco::PFMETCollection> pfmt;
  iEvent.getByToken(tok_PFMET_,pfmt);
  if (!pfmt.isValid()){
    edm::LogWarning("AlCaDiJets") << "AlCaDiJetsProducer: Error! can't get pfmet product!";
    return ;
  }
  const reco::PFMETCollection pfmet = *(pfmt.product());

  edm::Handle<reco::PFMETCollection> pfmt1;
  iEvent.getByToken(tok_PFType1MET_,pfmt1);
  if (!pfmt1.isValid()){
    edm::LogWarning("AlCaDiJets") << "AlCaDiJetsProducer: Error! can't get pftype1met product!";
    return ;
  }
  const reco::PFMETCollection pfmet1 = *(pfmt1.product());*/

  edm::Handle<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> > > hbhe;
  iEvent.getByToken(tok_HBHE_,hbhe);
  if (!hbhe.isValid()){
    edm::LogWarning("AlCaDiJets") << "AlCaDiJetsProducer: Error! can't get hbhe product!";
    return ;
  }
  const edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> > Hithbhe = *(hbhe.product());

  edm::Handle<edm::SortedCollection<HORecHit,edm::StrictWeakOrdering<HORecHit> > > ho;
  iEvent.getByToken(tok_HO_,ho);
  if(!ho.isValid()) {
    edm::LogWarning("AlCaDiJets") << "AlCaDiJetsProducer: Error! can't get ho product!";
    return ;
  }
  const edm::SortedCollection<HORecHit,edm::StrictWeakOrdering<HORecHit> > Hitho = *(ho.product());
    
  edm::Handle<edm::SortedCollection<HFRecHit,edm::StrictWeakOrdering<HFRecHit> > > hf;
  iEvent.getByToken(tok_HF_,hf);
  if(!hf.isValid()) {
    edm::LogWarning("AlCaDiJets") << "AlCaDiJetsProducer: Error! can't get hf product!";
    return ;
  }
  const edm::SortedCollection<HFRecHit,edm::StrictWeakOrdering<HFRecHit> > Hithf = *(hf.product());

  /*edm::Handle<edm::TriggerResults> trig;
  iEvent.getByToken(tok_TrigRes_,trig);
  if (!trig.isValid()){
    edm::LogWarning("AlCaDiJets") << "AlCaDiJetsProducer: Error! can't get trigger product!";
    return ;
  }
  const edm::TriggerResults trigres = *(trig.product());

  edm::Handle<double> rh;
  iEvent.getByToken(tok_Rho_,rh);
  if (!rh.isValid()){
    edm::LogWarning("AlCaDiJets") << "AlCaDiJetsProducer: Error! can't get rho product!";
    return ;
  }
  const double rho_val = *(rh.product());

  edm::Handle<reco::GsfElectronCollection> gsf;
  iEvent.getByToken(tok_GsfElec_,gsf);
  if (!gsf.isValid()){
    edm::LogWarning("AlCaDiJets") << "AlCaDiJetsProducer: Error! can't get gsfelectron product!";
    return ;
  }
  const reco::GsfElectronCollection gsfele = *(gsf.product());

  edm::Handle<reco::ConversionCollection> con;
  iEvent.getByToken(tok_Conv_,con);
  if (!con.isValid()){
    edm::LogWarning("AlCaDiJets") << "AlCaDiJetsProducer: Error! can't get conversion product!";
    return ;
  }
  const reco::ConversionCollection conv = *(con.product());

  edm::Handle<reco::BeamSpot> bs;
  iEvent.getByToken(tok_BS_,bs);
  if (!bs.isValid()){
    edm::LogWarning("AlCaDiJets") << "AlCaDiJetsProducer: Error! can't get beamspot product!";
    return ;
  }
  const reco::BeamSpot beam = *(bs.product());*/

  // See if this event is useful
  bool accept = select(pfjet);

  if (accept) {
    ++nSelect_;

    //Copy from standard place
    std::auto_ptr<reco::PFJetCollection>  miniPFjetCollection(new reco::PFJetCollection);
    for(reco::PFJetCollection::const_iterator pfjetItr=pfjets.begin(); pfjetItr!=pfjets.end(); ++pfjetItr) {
      miniPFjetCollection->push_back(*pfjetItr);
    }

    /*std::auto_ptr<reco::PhotonCollection> miniPhotonCollection(new reco::PhotonCollection);
    for(reco::PhotonCollection::const_iterator phoItr=photon.begin();
        phoItr!=photon.end(); phoItr++) {
	miniPhotonCollection->push_back(*phoItr);
	}*/

    std::auto_ptr<reco::PFCandidateCollection> miniPFCandCollection(new reco::PFCandidateCollection);
    for(reco::PFCandidateCollection::const_iterator pfcItr=pfcand.begin(); pfcItr!=pfcand.end(); ++pfcItr) {
      miniPFCandCollection->push_back(*pfcItr);
    }

    std::auto_ptr<reco::VertexCollection> miniVtxCollection(new reco::VertexCollection);
    for(reco::VertexCollection::const_iterator vtxItr=vtx.begin(); vtxItr!=vtx.end(); ++vtxItr) {
      miniVtxCollection->push_back(*vtxItr);
    }

    /*std::auto_ptr<reco::PFMETCollection> miniPFMETCollection(new reco::PFMETCollection);
    for(reco::PFMETCollection::const_iterator pfmetItr=pfmet.begin();
        pfmetItr!=pfmet.end(); pfmetItr++) {
      miniPFMETCollection->push_back(*pfmetItr);
    }

    std::auto_ptr<reco::PFMETCollection> miniPFMET1Collection(new reco::PFMETCollection);
    for(reco::PFMETCollection::const_iterator pfmet1Itr=pfmet1.begin();
        pfmet1Itr!=pfmet1.end(); pfmet1Itr++) {
      miniPFMET1Collection->push_back(*pfmet1Itr);
      }*/

    std::auto_ptr<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit>>>  miniHBHECollection(new edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit>>);
    for(edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> >::const_iterator hbheItr=Hithbhe.begin(); hbheItr!=Hithbhe.end(); ++hbheItr) {
      miniHBHECollection->push_back(*hbheItr);
    }

    std::auto_ptr<edm::SortedCollection<HORecHit,edm::StrictWeakOrdering<HORecHit>>>  miniHOCollection(new edm::SortedCollection<HORecHit,edm::StrictWeakOrdering<HORecHit>>);
    for(edm::SortedCollection<HORecHit,edm::StrictWeakOrdering<HORecHit> >::const_iterator hoItr=Hitho.begin(); hoItr!=Hitho.end(); ++hoItr) {
      miniHOCollection->push_back(*hoItr);
    }

    std::auto_ptr<edm::SortedCollection<HFRecHit,edm::StrictWeakOrdering<HFRecHit>>>  miniHFCollection(new edm::SortedCollection<HFRecHit,edm::StrictWeakOrdering<HFRecHit>>);
    for(edm::SortedCollection<HFRecHit,edm::StrictWeakOrdering<HFRecHit> >::const_iterator hfItr=Hithf.begin(); hfItr!=Hithf.end(); ++hfItr) {
      miniHFCollection->push_back(*hfItr);
    }

    /*std::auto_ptr<reco::GsfElectronCollection> miniGSFeleCollection(new reco::GsfElectronCollection);
    for(reco::GsfElectronCollection::const_iterator gsfItr=gsfele.begin();
        gsfItr!=gsfele.end(); gsfItr++) {
      miniGSFeleCollection->push_back(*gsfItr);
    }

    std::auto_ptr<reco::ConversionCollection> miniConversionCollection(new reco::ConversionCollection);
    for(reco::ConversionCollection::const_iterator convItr=conv.begin();
        convItr!=conv.end(); convItr++) {
      miniConversionCollection->push_back(*convItr);
      }
    
    std::auto_ptr<reco::BeamSpot> miniBeamSpotCollection(new reco::BeamSpot(beam.position(),beam.sigmaZ(),
									    beam.dxdz(),beam.dydz(),beam.BeamWidthX(),
									    beam.covariance(),beam.type()));
    
    std::auto_ptr<edm::TriggerResults> miniTriggerCollection(new edm::TriggerResults);
    *miniTriggerCollection = trigres;

    std::auto_ptr<double> miniRhoCollection(new double);
    *miniRhoCollection = rho_val;

    edm::Handle<edm::ValueMap<Bool_t> > loosePhotonQual;
    iEvent.getByToken(tok_loosePhoton_, loosePhotonQual);
    std::auto_ptr<edm::ValueMap<Bool_t> > miniLoosePhoton;
    if (loosePhotonQual.isValid()) *miniLoosePhoton = *(loosePhotonQual.product());

    edm::Handle<edm::ValueMap<Bool_t> > tightPhotonQual;
    iEvent.getByToken(tok_tightPhoton_, tightPhotonQual);
    std::auto_ptr<edm::ValueMap<Bool_t> > miniTightPhoton;
    if (tightPhotonQual.isValid()) *miniTightPhoton = *(tightPhotonQual.product());*/

    //Put them in the event
    //iEvent.put( miniPhotonCollection,      labelPhoton_.label());
    iEvent.put( miniPFjetCollection,       labelPFJet_.label());
    iEvent.put( miniHBHECollection,        labelHBHE_.label());
    iEvent.put( miniHFCollection,          labelHF_.label());
    iEvent.put( miniHOCollection,          labelHO_.label());
    //iEvent.put( miniTriggerCollection,     labelTrigger_.label());
    iEvent.put( miniPFCandCollection,      labelPFCandidate_.label());
    iEvent.put( miniVtxCollection,         labelVertex_.label());
    //iEvent.put( miniPFMETCollection,       labelPFMET_.label());
    //iEvent.put( miniPFMET1Collection,      labelPFMETtype1_.label());
    //iEvent.put( miniGSFeleCollection,      labelGsfEle_.label());
    //iEvent.put( miniRhoCollection,         labelRho_.label());
    //iEvent.put( miniConversionCollection,  labelConv_.label());
    //iEvent.put( miniBeamSpotCollection,    labelBeamSpot_.label());
    //iEvent.put( miniLoosePhoton,           labelLoosePhot_.label());
    //iEvent.put( miniTightPhoton,           labelTightPhot_.label());
  }
  return;

}

DEFINE_FWK_MODULE(AlCaDiJetsProducer); 
