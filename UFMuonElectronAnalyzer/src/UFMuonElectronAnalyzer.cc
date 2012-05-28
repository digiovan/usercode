// -*- C++ -*-
//
// Package:    UFMuonElectronAnalyzer
// Class:      UFMuonElectronAnalyzer
// 
/**\class UFMuonElectronAnalyzer UFMuonElectronAnalyzer/src/UFMuonElectronAnalyzer.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Gian Piero Di Giovanni,32 4-B08,+41227674961,
//         Created:  Thur Oct 21 10:44:13 CEST 2010
// $Id: UFMuonElectronAnalyzer.cc,v 1.1.2.1 2012/05/27 19:57:01 digiovan Exp $
//
//


// system include files
#include <memory>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>

#include "TLorentzVector.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TBranch.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonIsolation.h"
#include "DataFormats/MuonReco/interface/MuonChamberMatch.h"
#include "DataFormats/MuonReco/interface/MuonSegmentMatch.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include <Geometry/Records/interface/MuonGeometryRecord.h>
#include <Geometry/CSCGeometry/interface/CSCGeometry.h>
#include <Geometry/CSCGeometry/interface/CSCLayer.h>
#include <Geometry/CSCGeometry/interface/CSCLayerGeometry.h>

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include "DataFormats/MuonDetId/interface/DTWireId.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"

#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "DataFormats/CSCRecHit/interface/CSCRecHit2DCollection.h"


#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"

#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
//
// math classes
//
#include "DataFormats/Math/interface/deltaPhi.h"
//
// trigger
// 
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "DataFormats/Common/interface/TriggerResults.h" 
#include "FWCore/Common/interface/TriggerNames.h" 
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

//
// vertexing
//
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
//
// gen particles
//
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

// HEEP
#include "SHarper/HEEPAnalyzer/interface/HEEPEventHelper.h"
#include "SHarper/HEEPAnalyzer/interface/HEEPEvent.h"

// 2010.11.21 Adding the Muon Cocktail
#include "DataFormats/MuonReco/interface/MuonCocktails.h"
#include "DataFormats/Math/interface/deltaR.h"

// PU Info
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

// Add the data formats
#include "ZPtToMuMu/UFMuonElectronAnalyzer/interface/DataFormats.h"

//
// class declaration
//


class UFMuonElectronAnalyzer : public edm::EDAnalyzer {

public:

  explicit UFMuonElectronAnalyzer(const edm::ParameterSet&);
  ~UFMuonElectronAnalyzer();
  
  int _numEvents;
  edm::ESHandle<TransientTrackBuilder> transientTrackBuilder;

  edm::ESHandle<GlobalTrackingGeometry> globalTrackingGeometry;
  //  MuonServiceProxy* theService;

  // electrons
  typedef std::vector<heep::Ele> electronsColl;

  // ele-muon pairs
  typedef std::pair<reco::Muon,heep::Ele> MuElePair;
  typedef std::vector<MuElePair>          MuElePairs;

  // ele-track, e.g. ele-cocktail
  typedef std::pair<reco::Track,heep::Ele> TrackElePair;
  typedef std::vector<TrackElePair> TrackElePairs;

  // general event information	
  _EventInfo eventInfo;
    
  // vertex information
  _VertexInfo vertexInfo;

  _MuonInfo  _muon;
  _EleInfo   _ele;
  _TrackInfo _cocktail;

  // combined muon-ele info
  float _recoCandMass;
  float _recoCandPt;

  // combined muon-ele info
  float _recoCandMassCktl;
  float _recoCandPtCktl;


  // where to save all the info  
  TFile* _outFile;	
  TTree* _outTree;


private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  void initMuon(_MuonInfo& muon);
  void initEle ( _EleInfo& ele );

  virtual void beginRun(edm::Run const &, edm::EventSetup const&);

  // method to select the muons
  bool isHltPassed (const edm::Event&, const edm::EventSetup&, const std::string& triggerName);
  bool isHltMatched(const edm::Event&, const edm::EventSetup&, 
                    const std::string& triggerName, const std::string& filterName,
                    const reco::Muon&);
  void addVersion(const HLTConfigProvider hltConfig_,
                  std::string& triggerBaseName,
                  std::string& triggerName);    
  std::string findFilterName(const HLTConfigProvider hltConfig_, const std::string& triggerName);
  double DR(double eta1, double eta2, double phi1, double phi2);

  bool passKinCuts(const reco::Muon& muon,
                   edm::Handle<reco::BeamSpot> beamSpotHandle);

  void displaySelection();

  // helpers for the electrons selections
  heep::EventHelper evtHelper_; //this is our magic class where all the nastyness is contained
  heep::Event heepEvt_;
  edm::ParameterSet heepPset;

  // muons
  edm::InputTag _muonColl;
  edm::InputTag _beamSpot;		
  std::string _getFilename;	

  // variable to cuts over
  double _isGlobal;
  double _isTracker;  

  double _ptMin;
  double _etaMax;
  double _normChiSquareMax;
  double _d0Max;

  int _numPixelLayersMin;   
  int _numTrackerLayersMin; 
  int _numStripLayersMin;   

  double _validFracTrackerMin; 

  int _numValidMuonHitsMin;
  int _numValidPixelHitsMin;
  int _numValidStripHitsMin;
  int _numValidTrackerHitsMin;
  int _numSegmentMatchesMin;
  int _numOfMatchedStationsMin;


  // module config parameters
  std::string   processName_;
  std::vector < std::string > triggerNames_;
  std::vector < std::string > triggerBaseNames_;
  std::vector < std::string > filterNames_;
  edm::InputTag triggerResultsTag_;
  edm::InputTag triggerEventTag_;


  // additional class data memebers
  bool _checkTrigger; // activate or not the trigger checking   
  edm::Handle<edm::TriggerResults>   triggerResultsHandle_;
  edm::Handle<trigger::TriggerEvent> triggerEventHandle_;
  HLTConfigProvider hltConfig_;


  MuElePairs const GetMuElePairs(reco::MuonCollection const* muons,
                                 electronsColl const* electrons) const;
  
  // placeholder: function not yet implemented
  TrackElePairs const GetTrackElePairs(reco::TrackCollection const* tracks,
                                       electronsColl const* electrons) const;
  

  TLorentzVector const GetLorentzVector(UFMuonElectronAnalyzer::MuElePair const* pair);// const; 
  TLorentzVector const GetLorentzVector(UFMuonElectronAnalyzer::TrackElePair const* pair);// const; 

 
  // useful switches
  bool _isVerbose;    // debug mode
  bool _isMonteCarlo;
};


//
// constructors and destructor
//
UFMuonElectronAnalyzer::UFMuonElectronAnalyzer(const edm::ParameterSet& iConfig):
  _numEvents(0),evtHelper_(),heepEvt_()
{

  //now do what ever initialization is needed
  heepPset = iConfig.getParameter<edm::ParameterSet>("heepPSet");
  evtHelper_.setup(heepPset);

  _getFilename  = iConfig.getUntrackedParameter<std::string>("getFilename", "badger.root");
  _muonColl	= iConfig.getParameter<edm::InputTag>("muonColl");

  _beamSpot	= iConfig.getUntrackedParameter<edm::InputTag>("beamSpotTag",
                                                               edm::InputTag("offlineBeamSpot") );

  _isVerbose	= iConfig.getUntrackedParameter<bool>("isVerbose",   false);
  _isMonteCarlo	= iConfig.getParameter<bool>("isMonteCarlo");

  //edm::ParameterSet serviceParameters = iConfig.getParameter<edm::ParameterSet>("ServiceParameters");
  //theService = new MuonServiceProxy(serviceParameters);

  // selection cuts
  _ptMin                  = iConfig.getParameter<double>("ptMin");
  _etaMax                 = iConfig.getParameter<double>("etaMax");

  _normChiSquareMax	  = iConfig.getParameter<double>("normChiSquareMax");

  _numPixelLayersMin      = iConfig.getParameter<int>("minPixelLayers");  
  _numTrackerLayersMin    = iConfig.getParameter<int>("minTrackerLayers");
  _numStripLayersMin      = iConfig.getParameter<int>("minStripLayers");  
                                                 
  _validFracTrackerMin    = iConfig.getParameter<int>("validFracTrackerMin");

  _numValidMuonHitsMin	  = iConfig.getParameter<int>("minMuonHits");
  _numValidPixelHitsMin	  = iConfig.getParameter<int>("minPixelHits");
  _numValidStripHitsMin	  = iConfig.getParameter<int>("minStripHits");
  _numValidTrackerHitsMin = iConfig.getParameter<int>("minTrackerHits");
  _numSegmentMatchesMin	  = iConfig.getParameter<int>("minSegmentMatches");
  _numOfMatchedStationsMin= iConfig.getParameter<int>("minNumOfMatchedStations");
  
  _d0Max                  = iConfig.getParameter<double>("d0Max");

  //HLT trigger initialization
  _checkTrigger	     = iConfig.getParameter<bool>("checkTrigger");

  processName_       = iConfig.getParameter<std::string>("processName");
  triggerBaseNames_  = iConfig.getParameter<std::vector <std::string> >("triggerNames");

  if (triggerBaseNames_[0]=="") {
    triggerBaseNames_.clear();
    triggerNames_.clear();
  }

  if (triggerNames_.size() > 2){
    std::cout << "ERROR: triggerNames has its maximum size to 2! -> change the code in case\n";
    assert(0);
  }
  
  for (unsigned int i=0; i<triggerBaseNames_.size();i++) {
    triggerNames_.push_back("HLT_placeHolder");
    filterNames_.push_back("filter_placeHolder");
  }

  if (triggerNames_.size() == 0 && _checkTrigger){
    std::cout << "ERROR: You need to pass at least one valid HLT trigger\n";
    assert(0);
  }

  triggerResultsTag_ = iConfig.getParameter<edm::InputTag>("triggerResults");
  triggerEventTag_   = iConfig.getParameter<edm::InputTag>("triggerEvent");

}



UFMuonElectronAnalyzer::~UFMuonElectronAnalyzer() {}


// ------------ method called to for each event  ------------
void UFMuonElectronAnalyzer::analyze(const edm::Event& iEvent, 
                                     const edm::EventSetup& iSetup)
{

  // H L T   H A N D L E S
  iEvent.getByLabel(triggerResultsTag_,triggerResultsHandle_);
  if (!triggerResultsHandle_.isValid()) {
    std::cout << "UFMuonElectronAnalyzer::analyze: Error in getting TriggerResults product from Event!" << std::endl;
    return;
  }
  iEvent.getByLabel(triggerEventTag_,triggerEventHandle_);
  if (!triggerEventHandle_.isValid()) {
    std::cout << "UFMuonElectronAnalyzer::analyze: Error in getting TriggerEvent product from Event!" << std::endl;
    return;
  }
  // sanity check
  assert(triggerResultsHandle_->size()==hltConfig_.size());
  
  // if all the HLT paths are not fired, will discard the event immediately
  if (_checkTrigger) {
    bool isPassed=false;
     
    for (unsigned int iTrigger=0; iTrigger<triggerNames_.size(); iTrigger++) 
      if ( isHltPassed(iEvent,iSetup,triggerNames_[iTrigger]) ) isPassed=true;
     
    if (!isPassed) {
      if (_isVerbose) std::cout << "None of the HLT paths fired -> discard the event\n";
      return;
    }
     
  }

  int theRun   = iEvent.id().run();
  int theLumi  = iEvent.luminosityBlock();
  int theEvent = iEvent.id().event();
  int theBx    = iEvent.bunchCrossing();
  int theOrbit = iEvent.orbitNumber();
 

  eventInfo.run   = theRun;
  eventInfo.lumi  = theLumi;
  eventInfo.event = theEvent;
  eventInfo.bx    = theBx;
  eventInfo.orbit = theOrbit;
  

  //vertices
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByLabel("offlinePrimaryVertices", vertices);
 
  for (int i=0;i<20;i++) {
    vertexInfo.isValid[i]  = 0;
    vertexInfo.x[i]        = -999;     
    vertexInfo.y[i]        = -999;     
    vertexInfo.z[i]        = -999;     
    vertexInfo.xErr[i]     = -999;
    vertexInfo.yErr[i]     = -999;
    vertexInfo.zErr[i]     = -999;
    vertexInfo.chi2[i]     = -999;
    vertexInfo.ndf[i]      = -999;
    vertexInfo.normChi2[i] = -999;
  }      
  vertexInfo.nVertices   = 0;
  
  // init (vertices)
  if (vertices.isValid()) {
    
    //std::cout << "vertex->size():"<< vertices->size() << std::endl;
    int iVertex=0;
    for (reco::VertexCollection::const_iterator vtx = vertices->begin(); vtx!=vertices->end(); ++vtx){

      if (!vtx->isValid()) {
        vertexInfo.isValid[iVertex] = 0;
        continue;
      }
    
      vertexInfo.isValid[iVertex] = 1;
      
      vertexInfo.x[iVertex]        = vtx->position().X();	
      vertexInfo.y[iVertex]        = vtx->position().Y();	
      vertexInfo.z[iVertex]        = vtx->position().Z();	
      vertexInfo.xErr[iVertex]     = vtx->xError();	
      vertexInfo.yErr[iVertex]     = vtx->yError();	
      vertexInfo.zErr[iVertex]     = vtx->zError();	
      vertexInfo.chi2[iVertex]     = vtx->chi2();	
      vertexInfo.ndf[iVertex]      = vtx->ndof();	
      vertexInfo.normChi2[iVertex] = vtx->normalizedChi2();
      
      iVertex++;
      vertexInfo.nVertices++;
    }
  }
  else std::cout << "VertexCollection is NOT valid -> vertex Info NOT filled!\n";



  // ***************************************************************************
  // E L E C T R O N S
  //
  evtHelper_.makeHeepEvent(iEvent,iSetup,heepEvt_);
  const std::vector<heep::Ele>& eles = heepEvt_.heepEles();
 
  // Select all the electrons which pass the HEEP Id
  electronsColl elesSelected;
  
   // Select all the electrons which pass the HEEP Id
  for(size_t eleNr=0;eleNr<eles.size();eleNr++){
          
    if( eles[eleNr].cutCode()==0x0 ) elesSelected.push_back(eles[eleNr]);
  }
  // ***************************************************************************


  // ===========================================================================
  // M U O N S
  //
  edm::Handle<reco::MuonCollection> muons;
  iEvent.getByLabel(_muonColl, muons);
  
  // B E A M S P O T
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByLabel(_beamSpot, beamSpotHandle);

  // loop over muons
  reco::MuonCollection muonsSelected; muonsSelected.clear(); 

  for (reco::MuonCollection::const_iterator muon = muons->begin(), 
         muonsEnd = muons->end(); muon !=muonsEnd; ++muon){

    // relative tracker isolation cut    
    double isovar = muon->isolationR03().sumPt;
    isovar /= muon->pt(); // relative combine isolation

    if (_isVerbose) {
      std::cout << "Checking isolation...\n";
      std::cout << "Rel. Tracker Iso. = " << isovar << std::endl;
    }
    
    // kinematic cuts
    if ( !passKinCuts(*muon, beamSpotHandle) ) {
      if (_isVerbose) std::cout << "Muon NOT passing kinematic selections\n"; 
      continue;
    }
  
    // trigger cut
    if (_checkTrigger) {
      
      bool passSelection=!true;

      for (unsigned int iTrigger=0; iTrigger<triggerNames_.size(); iTrigger++) {

        if ( isHltMatched(iEvent,iSetup,triggerNames_[iTrigger], filterNames_[iTrigger], *muon) ) 
          passSelection=true;
      }
      
      if ( !passSelection ) continue;
    }
    

    if (_isVerbose) std::cout << "Muon passing the desired selections\n";     

    // put this muons in the collection
    muonsSelected.push_back(*muon);
  
  }
  // ===========================================================================
 
  // at least one muon and one electron
  if (_isVerbose) {
    std::cout<< "muons: "     << muonsSelected.size() << "  "
             << "electrons: " << elesSelected.size()  << std::endl;
  }
  
  if (elesSelected.size()  < 1 ||
      muonsSelected.size() < 1  ) {
    return;
  }

  // grab the candidates
  MuElePairs mueleCands = GetMuElePairs(&muonsSelected, &elesSelected);
  
  // loop over the candidates
  for (MuElePairs::const_iterator pair= mueleCands.begin(); 
       pair != mueleCands.end(); ++pair){
  
    reco::Muon mu  = pair->first;
    heep::Ele  ele = pair->second;

    initMuon(_muon);
    initEle (_ele );
    _recoCandMass = -999;
    _recoCandPt   = -999;
    
 
    // store all the info
    // muon
    _muon.isGlobal     = mu.isGlobalMuon(); 
    _muon.isTracker    = mu.isTrackerMuon(); 
    _muon.isStandAlone = mu.isStandAloneMuon(); 

    reco::TrackRef global;
    if      (mu.isGlobalMuon())     global = mu.globalTrack();
    else if (mu.isTrackerMuon())    global = mu.innerTrack();
    else {
      std::cout << "ERROR: The muon is NOT global NOR tracker ?!?\n";
      return ;
    }

    _muon.charge = mu.charge();
    _muon.pt     = mu.pt();  
    _muon.ptErr  = global->ptError(); 
    _muon.eta    = mu.eta(); 
    _muon.phi    = mu.phi();
   
    if (mu.isTrackerMuon()) {
      _muon.trkPt   = mu.innerTrack()->pt();                    
      _muon.trkPtErr= mu.innerTrack()->ptError();
      _muon.trketa  = mu.innerTrack()->eta();                   
      _muon.trkPhi  = mu.innerTrack()->phi();                   
    }

    _muon.normChiSquare=global->normalizedChi2();
    _muon.d0= global->dxy( beamSpotHandle->position() );
    _muon.dz= global->dz ( beamSpotHandle->position() );

    _muon.numPixelLayers   = global->hitPattern().pixelLayersWithMeasurement();   
    _muon.numTrackerLayers = global->hitPattern().trackerLayersWithMeasurement(); 
    _muon.numStripLayers   = global->hitPattern().stripLayersWithMeasurement();   
    
    _muon.validFracTracker = global->validFraction(); 

    _muon.numValidMuonHits    = global->hitPattern().numberOfValidMuonHits();
    _muon.numValidPixelHits   = global->hitPattern().numberOfValidPixelHits();
    _muon.numValidTrackerHits = global->hitPattern().numberOfValidTrackerHits();
    _muon.numValidStripHits   = global->hitPattern().numberOfValidStripHits();
    _muon.numSegmentMatches   = mu.numberOfMatches();
    _muon.numOfMatchedStations= mu.numberOfMatchedStations();

    // save trigger informations
    for (unsigned int iTrigger=0;iTrigger<triggerNames_.size();iTrigger++)
      _muon.isHltMatched[iTrigger] = (isHltMatched(iEvent, iSetup, triggerNames_[iTrigger], filterNames_[iTrigger], mu) );

    
    // - relative combined isolation -
    double isovar = mu.isolationR03().sumPt;
    isovar       += mu.isolationR03().hadEt; //HCAL
    isovar       /= mu.pt(); // relative combine isolation

    _muon.relCombIso=isovar;
    // - end relative combined isolation -

    
    //isolation info
    _muon.trackIsoSumPt = mu.isolationR03().sumPt ;
    _muon.ecalIso       = mu.isolationR03().emEt ;
    _muon.hcalIso       = mu.isolationR03().hadEt ;
    
    // PF Isolation
    _muon.isPFMuon = mu.isPFMuon();
    
    if ( mu.isPFMuon() ) {
      
      reco::Candidate::LorentzVector pfmuon = mu.pfP4();
      
      _muon.pfPt  = pfmuon.Pt();
      _muon.pfEta = pfmuon.Eta();
      _muon.pfPhi = pfmuon.Phi();

      _muon.sumChargedHadronPtR03   = mu.pfIsolationR03().sumChargedHadronPt  ;
      _muon.sumChargedParticlePtR03 = mu.pfIsolationR03().sumChargedParticlePt;
      _muon.sumNeutralHadronEtR03   = mu.pfIsolationR03().sumNeutralHadronEt  ;
      _muon.sumPhotonEtR03          = mu.pfIsolationR03().sumPhotonEt         ;
      _muon.sumPUPtR03              = mu.pfIsolationR03().sumPUPt             ;
                                        
      _muon.sumChargedHadronPtR04   = mu.pfIsolationR04().sumChargedHadronPt  ;
      _muon.sumChargedParticlePtR04 = mu.pfIsolationR04().sumChargedParticlePt;
      _muon.sumNeutralHadronEtR04   = mu.pfIsolationR04().sumNeutralHadronEt  ;
      _muon.sumPhotonEtR04          = mu.pfIsolationR04().sumPhotonEt         ;
      _muon.sumPUPtR04              = mu.pfIsolationR04().sumPUPt             ;
    }


    // ele
    _ele.charge = ele.charge(); 
    _ele.et     = ele.et(); 
    _ele.eta    = ele.eta(); 
    _ele.phi    = ele.phi();

    _ele.isEB    = ele.isEB();
    _ele.isEE    = ele.isEE();



    // combine the info
    TLorentzVector mother=GetLorentzVector(&*pair); 
    _recoCandMass = mother.M();
    _recoCandPt   = mother.Pt();
    
    // store them in a ntuple
    _outTree->Fill();
    
    if (_isVerbose) {
      std::cout<<"\t"<<theRun    <<":"<<theLumi<<":"<<theEvent<<"\t"<<mother.M()
               <<"\t"<<_muon.eta<<":"<<_ele.eta
               <<"\t"<<_muon.pt <<":"<<_ele.et
               <<std::endl;
    }

  }
  
  return;
  
}
	
// ------------ method called once each job just before starting event loop  ------------
void UFMuonElectronAnalyzer::beginJob()
{
  
  displaySelection();

  // create the file
  _outFile	= new TFile(_getFilename.c_str(), "recreate");	
  _outFile->cd();

  // create the tree
  _outTree	= new TTree("tree", "myTree");


  // eventInfo;
  _outTree->Branch("eventInfo",  &eventInfo, "run/I:lumi/I:event/I:bx/I:orbit/I");

  _outTree->Branch("vertexInfo", &vertexInfo, "nVertices/I:isValid[20]/I:"
		   "x[20]/F:y[20]/F:z[20]/F:xErr[20]/F:yErr[20]/F:zErr[20]/F:"
		   "chi2[20]/F:ndf[20]/I:normChi2[20]/F");

   _outTree->Branch("muon"     , &_muon    	, 
                   "isTracker/I:isStandAlone/I:isGlobal/I:"
                   "charge/I:pt/F:ptErr/F:eta/F:phi/F:"
                   "trkPt/F:trkPtErr/F:trkEta/F:trkPhi/F:"
                   "normChiSquare/F:d0/F:dz/F:"
                   "numPixelLayers/I:"
                   "numTrackerLayers/I:"
                   "numStripLayers/I:"
                   "validFracTracker/F:"
                   "numValidMuonHits/I:"
                   "numValidPixelHits/I:"    
                   "numValidTrackerHits/I:"  
                   "numValidStripHits/I:"    
                   "numSegmentMatches/I:"    
                   "numOfMatchedStations/I:"
                   "trackIsoSumPt/F:"      
                   "hcalIso/F:"
                   "ecalIso/F:"
                   "relCombIso/F:"
                   "isPFMuon/I:"
                   "pfPt/F:"
                   "pfEta/F:"
                   "pfPhi/F:"
                   "sumChargedHadronPtR03/F:"
                   "sumChargedParticlePtR03/F:"
                   "sumNeutralHadronEtR03/F:"
                   "sumPhotonEtR03/F:"
                   "sumPUPtR03/F:"
                   "sumChargedHadronPtR04/F:"
                   "sumChargedParticlePtR04/F:"
                   "sumNeutralHadronEtR04/F:"
                   "sumPhotonEtR04/F:"
                   "sumPUPtR04/F:"
                   "isHltMatched[2]/I");

  

   _outTree->Branch("ele", &_ele,         
                    "charge/I:et/F:eta/F:phi/F:isEB/I:isEE");

  // mass and pt for the fir
  _outTree->Branch("recoCandMass", &_recoCandMass, "recoCandMass/F");
  _outTree->Branch("recoCandPt" ,  &_recoCandPt  , "recoCandPt/F");

}


// ------------ method called once each job just after ending the event loop  ------------
void UFMuonElectronAnalyzer::endJob() {

  std::cout<<"number of candidate muon-ele events: "
           <<_outTree->GetEntries()<<std::endl;
  _outFile->cd();
  _outTree->Write();

}

TLorentzVector const UFMuonElectronAnalyzer::GetLorentzVector(UFMuonElectronAnalyzer::MuElePair const* pair) /*const*/ {

  TLorentzVector muon, ele, sum;
  double const MASS_MUON = 0.105658367;    //GeV/c2
  double const MASS_ELE  = 0.000510998910; //GeV/c2

  reco::TrackRef const muonTrack = pair->first.innerTrack();
  heep::Ele const electron = pair->second;

  muon.SetPtEtaPhiM(muonTrack->pt(), muonTrack->eta(), muonTrack->phi(), MASS_MUON);
  ele.SetPtEtaPhiM (electron.et()  , electron.eta()  , electron.phi()  , MASS_ELE);

  sum = muon+ele;
  return sum;
  
}

TLorentzVector const UFMuonElectronAnalyzer::GetLorentzVector(UFMuonElectronAnalyzer::TrackElePair const* pair) /*const*/ {
  
  TLorentzVector muon, ele, sum;
  double const MASS_MUON = 0.105658367;    //GeV/c2
  double const MASS_ELE  = 0.000510998910; //GeV/c2
  
  reco::Track const muonTrack = pair->first;
  heep::Ele const electron = pair->second;

  muon.SetPtEtaPhiM(muonTrack.pt(), muonTrack.eta(), muonTrack.phi(), MASS_MUON);
  ele.SetPtEtaPhiM (electron.et() , electron.eta() , electron.phi() , MASS_ELE);

  sum = muon+ele;
  return sum;
  
}


// check the HLT configuration for each run. It may change you know ;-)
void UFMuonElectronAnalyzer::beginRun(edm::Run const& iRun, 
                                      edm::EventSetup const& iSetup)
{
  using namespace std;
  using namespace edm;

  bool changed(true);
  if (hltConfig_.init(iRun,iSetup,processName_,changed)) {
    
    // if you set the trigger to look for the lowest unprescaled single mu trigger 
    // only there is no need to check its existence as it automatically will be
    // found later in the code
    if (changed) {
      
      // check if trigger name in (new) config
      const unsigned int n(hltConfig_.size());
      
      // loop over all the trigger names provided
      unsigned int triggerSize = triggerNames_.size();
      for (unsigned int iTrigger=0; iTrigger<triggerSize; iTrigger++) {

        addVersion(hltConfig_, triggerBaseNames_[iTrigger], triggerNames_[iTrigger]);
        //std::cout << "The trigger after  is " << triggerNames_[iTrigger]  << std::endl;

	const unsigned int triggerIndex(hltConfig_.triggerIndex(triggerNames_[iTrigger]));
	if (triggerIndex>=n) {
          std::cout << "\n\nHLTEventAnalyzerAOD::analyze:"
                    << " TriggerName \"" << triggerNames_[iTrigger] 
                    << "\" is NOT available in (new) config!" << std::endl << std::endl;
          std::cout << " The available TriggerNames are: " << std::endl;
	  hltConfig_.dump("Triggers");
          
          throw cms::Exception("UFDiMuonsAnalyzer")<< "Throwing an exception because "
                                                   << "the trigger path name you want to check DOES NOT EXIST";
	}
        else filterNames_[iTrigger] = findFilterName ( hltConfig_, triggerNames_[iTrigger] ) ;
        

      }
      // dear god you do not want to uncomment them... but one day you could be
      // interested so I leave them there as a potential reference.
      //hltConfig_.dump("Streams");
      //hltConfig_.dump("Datasets");
      //hltConfig_.dump("PrescaleTable");
      //hltConfig_.dump("ProcessPSet");
    }
  } else {
    cout << "HLTEventAnalyzerAOD::analyze:"
	 << " config extraction failure with process name "
	 << processName_ << endl;
    
    throw cms::Exception("UFDiMuonsAnalyzer")<<"Wrong processName_(\""<<processName_
                                             <<"\"): please double check what you passed "
                                             <<"in the python file... ";
    
  }
}

// this method will simply check is the selected HLT path (via triggerName)
// is run and accepted and no error are found
//
// bool true  if (run && accept && !error)
//      false if any other combination
bool UFMuonElectronAnalyzer::isHltPassed(const edm::Event& iEvent, 
                                         const edm::EventSetup& iSetup, 
                                         const std::string& triggerName) {
  
  using namespace std;
  using namespace edm;
  using namespace reco;
  using namespace trigger;

  if (_isVerbose)
    std::cout << "\nANALYZING THE TRIGGER "<< triggerName << "...\n";
      
  const unsigned int n(hltConfig_.size());
  const unsigned int triggerIndex(hltConfig_.triggerIndex(triggerName));
  assert(triggerIndex==iEvent.triggerNames(*triggerResultsHandle_).triggerIndex(triggerName));

  // abort on invalid trigger name
  if (triggerIndex>=n) {
    cout << "UFMuonElectronAnalyzer::isHltPassed: path "
	 << triggerName << " - not found!" << endl;
    return false;
  }

  // Results from TriggerResults product
  if (_isVerbose) {
    std::cout << " Trigger path status:"
              << " WasRun=" << triggerResultsHandle_->wasrun(triggerIndex)
              << " Accept=" << triggerResultsHandle_->accept(triggerIndex)
              << " Error =" << triggerResultsHandle_->error(triggerIndex)
              << std::endl;
  }

  bool wasRun = triggerResultsHandle_->wasrun(triggerIndex);
  bool accept = triggerResultsHandle_->accept(triggerIndex);
  bool error  = triggerResultsHandle_->error (triggerIndex);

  bool isPassed=true;  

  if (!wasRun) isPassed=false;
  if (!accept) isPassed=false;
  if ( error ) isPassed=false;

  return isPassed;
}

// same check for isHltPassed +
// check if the muon is the one firing the HLT path
bool UFMuonElectronAnalyzer::isHltMatched(const edm::Event& iEvent, 
                                          const edm::EventSetup& iSetup, 
                                          const std::string& triggerName, 
                                          const std::string& filterName, 
                                          const reco::Muon& muon){

  if (_isVerbose)  std::cout << "analyzing hlt matching" << std::endl;

  bool isMatched=false;

  using namespace std;
  using namespace edm;
  using namespace reco;
  using namespace trigger;

  //cout << endl;
  //std::cout << "ANALYZING THE TRIGGER "<< triggerName << "..."<<std::endl;

  const unsigned int n(hltConfig_.size());
  const unsigned int triggerIndex(hltConfig_.triggerIndex(triggerName));
  assert(triggerIndex==iEvent.triggerNames(*triggerResultsHandle_).triggerIndex(triggerName));

  // abort on invalid trigger name
  if (triggerIndex>=n) {
    cout << "UFDiMuonsAnalyzer::isHltMatched: path "
	 << triggerName << " - not found!" << endl;
    return isMatched;
  }

  // Results from TriggerResults product
  //std::cout << " Trigger path status:"
  //          << " WasRun=" << triggerResultsHandle_->wasrun(triggerIndex)
  //          << " Accept=" << triggerResultsHandle_->accept(triggerIndex)
  //          << " Error =" << triggerResultsHandle_->error(triggerIndex)
  //          << std::endl;
  
  bool wasRun = triggerResultsHandle_->wasrun(triggerIndex);
  bool accept = triggerResultsHandle_->accept(triggerIndex);
  bool error  = triggerResultsHandle_->error (triggerIndex);

  if (!wasRun) return isMatched;
  if (!accept) return isMatched;
  if ( error ) return isMatched;
  
  // modules on this trigger path
  const unsigned int m(hltConfig_.size(triggerIndex));
  const vector<string>& moduleLabels(hltConfig_.moduleLabels(triggerIndex));

  // Results from TriggerResults product
  //std::cout << " Trigger path status:"
  //          << " WasRun=" << triggerResultsHandle_->wasrun(triggerIndex)
  //          << " Accept=" << triggerResultsHandle_->accept(triggerIndex)
  //          << " Error =" << triggerResultsHandle_->error(triggerIndex)
  //          << std::endl;
  
  const unsigned int moduleIndex(triggerResultsHandle_->index(triggerIndex));
  //std::cout << " Last active module - label/type: "
  //          << moduleLabels[moduleIndex] << "/" << hltConfig_.moduleType(moduleLabels[moduleIndex])
  //          << " [" << moduleIndex << " out of 0-" << (m-1) << " on this path]"
  //          << std::endl;
  assert (moduleIndex<m);

  // Results from TriggerEvent product - Attention: must look only for
  // modules actually run in this path for this event!
  for (unsigned int j=0; j<=moduleIndex; ++j) {
    
    const string& moduleLabel(moduleLabels[j]);
    const string  moduleType(hltConfig_.moduleType(moduleLabel));
  
    // check only the last module
    if (moduleLabel != filterName) continue;

    // check whether the module is packed up in TriggerEvent product
    const unsigned int filterIndex(triggerEventHandle_->filterIndex(InputTag(moduleLabel,"",processName_)));
  
    if (filterIndex<triggerEventHandle_->sizeFilters()) {
      //std::cout << " 'L3' filter in slot " << j 
      //          << " - label/type "        << moduleLabel 
      //          << "/" << moduleType << std::endl;
    
      const Vids& VIDS (triggerEventHandle_->filterIds (filterIndex));
      const Keys& KEYS (triggerEventHandle_->filterKeys(filterIndex));
      const size_type nI(VIDS.size());
      const size_type nK(KEYS.size());
      assert(nI==nK);

      const size_type n(max(nI,nK));
      if (_isVerbose)
        std::cout << "   " << n  << " accepted 'L3' objects found: " << std::endl;

      const TriggerObjectCollection& TOC(triggerEventHandle_->getObjects());
    
      for (size_type i=0; i!=n; ++i) {
        const TriggerObject& TO(TOC[KEYS[i]]);
        
   
        double DPT=fabs(TO.pt()-muon.pt());

        if (_isVerbose) {
          std::cout << "   " << i << " " << VIDS[i] << "/" << KEYS[i] << ": "
                    << TO.id()  << " " << TO.pt() << " " << TO.eta() << " " 
                    << TO.phi() << " " << TO.mass()
                    << std::endl;

          std::cout << " muon: eta=" << muon.eta()
                    << ", phi=" << muon.phi()
                    << std::endl;

          std::cout << " DPT=" << DPT
                    << " [max is 0.2]"  
                    << std::endl;


        }
        
        if ( DR( TO.eta(),muon.eta(), TO.phi(),muon.phi() ) < 0.2) isMatched=true;      
	//&& fabs(TO.eta())<2.1  ) isMatched=true;      
        // removed the DPT/PT used in the 2010 analyses
        //&& (DPT/muon.pt() < 1) ) isMatched=true;      
      }
    }
    
  }
 
  if (_isVerbose) std::cout << "is HTL Matched? " << isMatched << std::endl;
  return isMatched;
}



double UFMuonElectronAnalyzer::DR(double eta1, double eta2,
                                  double phi1, double phi2){
  
  double diffEta = eta1 - eta2;
  double diffPhi = phi1 - phi2;
  double dr = sqrt(diffEta*diffEta + diffPhi*diffPhi);

  //std::cout << "diffEta: " << diffEta << std::endl;
  //std::cout << "diffPhi: " << diffPhi << std::endl;
  //std::cout << "dr: " << dr << std::endl;

  return dr;
  
}



UFMuonElectronAnalyzer::MuElePairs const UFMuonElectronAnalyzer::GetMuElePairs(reco::MuonCollection const* muons,
                                                                           electronsColl const* electrons) const {
                                                                           
  
  MuElePairs muelepairs; muelepairs.clear();
  
  for (reco::MuonCollection::const_iterator muon = muons->begin(), 
         muonsEnd = muons->end(); muon != muonsEnd; ++muon){

    for (electronsColl::const_iterator ele = electrons->begin(),
           eleEnd = electrons->end(); ele != eleEnd; ++ele){

      muelepairs.push_back( UFMuonElectronAnalyzer::MuElePair(*muon,*ele) );
    }
  }
  
  return muelepairs;
}

void UFMuonElectronAnalyzer::displaySelection() {

  std::cout << "\n\n*** UFMuonElectronAnalyzer Configuration ***\n";
  std::cout << " - Events saved in file: " << _getFilename << std::endl;	

  // variable to cuts over
  std::cout << " - _ptMin:            " << _ptMin << std::endl;
  std::cout << " - _etaMax:           " << _etaMax<< std::endl;
  std::cout << " - _normChiSquareMax: " << _normChiSquareMax << std::endl;
  std::cout << " - _d0Max:            " << _d0Max << std::endl;            

  std::cout << " - _numPixelLayersMin:   " << _numPixelLayersMin  << std::endl;   
  std::cout << " - _numTrackerLayersMin: " << _numTrackerLayersMin<< std::endl; 
  std::cout << " - _numStripLayersMin:   " << _numStripLayersMin  << std::endl;   

  std::cout << " - _validFracTrackerMin: " << _validFracTrackerMin<< std::endl; 

  std::cout << " - _numValidMuonHitsMin:    " << _numValidMuonHitsMin     << std::endl;
  std::cout << " - _numValidPixelHitsMin:   " << _numValidPixelHitsMin    << std::endl;
  std::cout << " - _numValidStripHitsMin:   " << _numValidStripHitsMin    << std::endl;
  std::cout << " - _numValidTrackerHitsMin: " << _numValidTrackerHitsMin  << std::endl;
  std::cout << " - _numSegmentMatchesMin:   " << _numSegmentMatchesMin    << std::endl;  
  std::cout << " - _numOfMatchedStationsMin:" << _numOfMatchedStationsMin << std::endl;  

  // module config parameters
  std::cout << " - _checkTrigger: " << _checkTrigger << std::endl;
  std::cout << " - Triggers To Probe:\n";
  unsigned int triggerSize = triggerNames_.size();
  for (unsigned int i=0; i < triggerSize; i++) 
    std::cout << "    * triggerNames["<<i<<"]: " << triggerBaseNames_[i] << std::endl;
  
  std::cout << std::endl << std::endl;

}

void UFMuonElectronAnalyzer::initMuon(_MuonInfo& muon) {

  muon.isTracker    = -999;
  muon.isStandAlone = -999;
  muon.isGlobal     = -999;

  muon.charge = -999;
  muon.pt     = -999;
  muon.eta    = -999; 
  muon.phi    = -999;
  
  muon.normChiSquare=-999;
  muon.d0= -999;
  muon.dz= -999;
  
  muon.numPixelLayers = -999; 
  muon.numTrackerLayers = -999;
  muon.numStripLayers = -999;  
  
  muon.validFracTracker = -999;

  muon.numValidMuonHits    = -999;
  muon.numValidPixelHits   = -999;
  muon.numValidTrackerHits = -999;
  muon.numValidStripHits   = -999;
  muon.numSegmentMatches   = -999;
  muon.numOfMatchedStations= -999;
  
  muon.trackIsoSumPt = -999;
  muon.hcalIso       = -999;
  muon.ecalIso       = -999;
  muon.relCombIso    = -999;

  muon.isPFMuon = -999;

  muon.pfPt  = -999;
  muon.pfEta = -999;
  muon.pfPhi = -999;

  muon.sumChargedHadronPtR03   = -999;
  muon.sumChargedParticlePtR03 = -999;
  muon.sumNeutralHadronEtR03   = -999;
  muon.sumPhotonEtR03          = -999;
  muon.sumPUPtR03              = -999;
  
  muon.sumChargedHadronPtR04   = -999;
  muon.sumChargedParticlePtR04 = -999;
  muon.sumNeutralHadronEtR04   = -999;
  muon.sumPhotonEtR04          = -999;
  muon.sumPUPtR04              = -999;

  muon.isHltMatched[0] = -999;
  muon.isHltMatched[1] = -999;

}

void UFMuonElectronAnalyzer::initEle(_EleInfo& ele) {

  ele.charge = -999;
  ele.et = -999;
  ele.eta = -999;
  ele.phi = -999;

  ele.isEB = -999;
  ele.isEE = -999;

}

void UFMuonElectronAnalyzer::addVersion(const HLTConfigProvider hltConfig_,
                                        std::string& triggerBaseName,
                                        std::string& triggerName) {

  //std::cout << "The trigger is " << triggerBaseName  << std::endl;
  
  // This function looks gets a trigger name, e.g. HLT_Mu15 and add to it
  // its version number, e.g. HLT_Mu15_v2, if the trigger is present in 
  // the hltConfig  
  const std::vector<std::string>& triggerNames = hltConfig_.triggerNames();
    
  for (size_t ts = 0; ts< triggerNames.size() ; ts++){
    std::string trig = triggerNames[ts];
    size_t f = trig.find(triggerBaseName);

    if (f != std::string::npos)  {

      //adding version extension
      for (unsigned int iVersion = 1; iVersion<30; iVersion++ ){
        std::string trigWithVersion= triggerBaseName;
        trigWithVersion.append("_v");
        std::stringstream ss;
        ss << iVersion;
        std::string version = ss.str(); 
        trigWithVersion.append(version);
        
        if (trig==trigWithVersion) {
          triggerName.replace(triggerName.begin(),    triggerName.end(),
                              trigWithVersion.begin(),trigWithVersion.end());
          return;
        }
        
      }
    }
    
  }
  
  std::cout << "TriggerBaseName \"" << triggerBaseName 
            << "\" was already \"versionized\" OR "
            << " the root is not found in the list of trigger which"
            << " means cmssw is going to crash providing the list"
            << " of available triggers\n";

  triggerName.replace(triggerName.begin(),    triggerName.end(),
                      triggerBaseName.begin(),triggerBaseName.end());

  return;

}


std::string UFMuonElectronAnalyzer::findFilterName(const HLTConfigProvider hltConfig_, 
                                                   const std::string& triggerName){
  
  std::string L3FilterName_;
  
  const std::vector<std::string>& moduleLabs = hltConfig_.moduleLabels(triggerName); 
    
  // the l3 filter name is just the last module.... 
  size_t moduleLabsSizeMinus2 = moduleLabs.size() - 2 ;
  
 
  L3FilterName_ = moduleLabs[moduleLabsSizeMinus2];

  if (_isVerbose) std::cout<<"triggerName="    << triggerName
                           <<" -> filterName=" << L3FilterName_ << std::endl;        
  
  // *********************************** //
  return L3FilterName_ ;
  //return moduleLabs[moduleLabsSizeMinus2];
  // *********************************** //
}

bool UFMuonElectronAnalyzer::passKinCuts(const reco::Muon& muon,
                                    edm::Handle<reco::BeamSpot> beamSpotHandle){
  
  bool passKinCuts=false;
  
  // =========================================================== //
  // What was corresponding to the old Loose VBTF
  // =========================================================== //
  if (_isVerbose) {
    std::cout<< "is Global?"     << muon.isGlobalMuon()     << std::endl;
    std::cout<< "is Tracker?"    << muon.isTrackerMuon()    << std::endl;
    std::cout<< "is StandAlone?" << muon.isStandAloneMuon() << std::endl;
  }

  // reconstruction cuts
  if (!muon.isGlobalMuon()    && _isGlobal )    return passKinCuts; // gbl muon
  if (!muon.isTrackerMuon()   && _isTracker)    return passKinCuts; // trk muon

  // do not accept muons which are standalone only
  if(!muon.isGlobalMuon() && !muon.isTrackerMuon()) return passKinCuts;

  reco::TrackRef globalTrack;
  if (muon.isGlobalMuon()) globalTrack = muon.globalTrack();
  else if (muon.isTrackerMuon()) globalTrack = muon.innerTrack();
  else {
    // redundant: just in case...
    std::cout << "ERROR: The muon is NOT global NOR tracker ?!?\n";
    return false;
  }

  if (_isVerbose) {
    std::cout<< "muon.pt(): "  << muon.pt() << " [ptMin="  << _ptMin  <<"]" << std::endl;
    std::cout<< "fabs(muon.eta()): " << fabs(muon.eta())<< " [etaMax=" << _etaMax <<"]" << std::endl;
    std::cout<< "numberOfValidTrackerHits(): " << globalTrack->hitPattern().numberOfValidTrackerHits() 
             << " [min=" << _numValidTrackerHitsMin << "]" << std::endl;
    
    std::cout<<"d0: " << fabs(globalTrack->dxy(beamSpotHandle->position()))
             << " [max=" << _d0Max << "]" << std::endl;
    
  }

    
  // kinematic cuts
  if (muon.pt()        < _ptMin ) return passKinCuts; // pt cut
  if (fabs(muon.eta()) > _etaMax) return passKinCuts; // eta cut
  if (globalTrack->hitPattern().numberOfValidTrackerHits() < _numValidTrackerHitsMin) return passKinCuts; // # hits in tracker

  // beam spot cut
  if (fabs(globalTrack->dxy(beamSpotHandle->position())) > _d0Max) return passKinCuts;


  // =========================================================== //
  // What was corresponding to the old Tight VBTF
  // + some additions
  // =========================================================== //
  if (_isVerbose) {
    std::cout << "numberOfValidMuonHits: " 
              << globalTrack->hitPattern().numberOfValidMuonHits() 
              << " [min=" << _numValidMuonHitsMin << "]" 
              << std::endl;
  
    std::cout << "numberOfValidPixelHits(): " 
              << globalTrack->hitPattern().numberOfValidPixelHits() 
              << " [min=" << _numValidPixelHitsMin << "]" 
              << std::endl;
     
    std::cout << "numberOfValidStripHits(): " 
              << globalTrack->hitPattern().numberOfValidStripHits() 
              << " [min=" << _numValidStripHitsMin << "]" 
              << std::endl;
    
    std::cout << "pixelLayersWithMeasurement(): " 
              << globalTrack->hitPattern().pixelLayersWithMeasurement() 
              << " [min=" << _numPixelLayersMin << "]" 
              << std::endl;

    std::cout << "trackerLayersWithMeasurement(): " 
              << globalTrack->hitPattern().trackerLayersWithMeasurement() 
              << " [min=" << _numTrackerLayersMin << "]" 
              << std::endl;

    std::cout << "stripLayersWithMeasurement(): " 
              << globalTrack->hitPattern().stripLayersWithMeasurement() 
              << " [min=" << _numStripLayersMin << "]" 
              << std::endl;

    std::cout << "validFraction(): " 
              << globalTrack-> validFraction()
              << " [min=" << _validFracTrackerMin << "]" 
              << std::endl;

    std::cout << "muon.numberOfMatches(): " 
              <<  muon.numberOfMatches() 
              << " [min=" << _numSegmentMatchesMin << "]" 
              << std::endl;
  
    std::cout << "muon.numberOfMatchedStations(): " 
              <<  muon.numberOfMatchedStations() 
              << " [min=" << _numOfMatchedStationsMin << "]" 
              << std::endl;
  
    std::cout<<"globalTrack->normalizedChi2(): " 
             << globalTrack->normalizedChi2()
             << " [max=" << _normChiSquareMax <<"]" 
             << std::endl;
   
  }

  if (globalTrack->hitPattern().pixelLayersWithMeasurement()   < _numPixelLayersMin)   return passKinCuts;   
  if (globalTrack->hitPattern().trackerLayersWithMeasurement() < _numTrackerLayersMin) return passKinCuts; 
  if (globalTrack->hitPattern().stripLayersWithMeasurement()   < _numTrackerLayersMin) return passKinCuts;   
  		 
  if (globalTrack->validFraction() < _validFracTrackerMin) return passKinCuts; 

  if ( globalTrack->hitPattern().numberOfValidMuonHits() < _numValidMuonHitsMin  )  return passKinCuts;
  if ( globalTrack->hitPattern().numberOfValidPixelHits() < _numValidPixelHitsMin ) return passKinCuts;
  if ( globalTrack->hitPattern().numberOfValidStripHits() < _numValidStripHitsMin ) return passKinCuts;
  if ( muon.numberOfMatches() < _numSegmentMatchesMin   )         return passKinCuts;
  if ( muon.numberOfMatchedStations() < _numOfMatchedStationsMin) return passKinCuts;
  if ( globalTrack->normalizedChi2() > _normChiSquareMax)         return passKinCuts;

  if (_isVerbose) std::cout << "passing kinematic cuts\n"; 

  passKinCuts=true;
  return passKinCuts;
}

DEFINE_FWK_MODULE(UFMuonElectronAnalyzer);


