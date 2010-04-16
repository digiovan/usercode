#include "TrigEff.h"
#include <map>
#include <vector>
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Framework/interface/TriggerNames.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticleFwd.h"
#include "DataFormats/L1Trigger/interface/L1ParticleMap.h"
#include "DataFormats/L1Trigger/interface/L1ParticleMapFwd.h"
//#include "DataFormats/HepMCCandidate/interface/GenParticleCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"


// from HLTEventAnalyzerRaw
// need access to class objects being referenced to get their content!
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/HcalIsolatedTrack/interface/IsolatedPixelTrackCandidate.h"
#include "DataFormats/L1Trigger/interface/L1HFRings.h"
#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticle.h"

#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include <cassert>

using namespace trigger;
using namespace reco;
using namespace std;
using namespace edm;

TrigEff::TrigEff(const edm::ParameterSet& pset):edm::EDAnalyzer(){
  // not all these tags are useful: need some clean-up
  genTag    = pset.getParameter<InputTag>("genTag");
  L1extraTag= pset.getParameter<InputTag>("L1extraTag");
  Level1Tag = pset.getParameter<InputTag>("Level1Tag");
  Level2Tag = pset.getParameter<InputTag>("Level2Tag");
  Level3Tag = pset.getParameter<InputTag>("Level3Tag");
  tracksTag = pset.getParameter<InputTag>("tracksTag");
  muonsTag  = pset.getParameter<InputTag>("muonsTag");
  saMuonsTag= pset.getParameter<InputTag>("saMuonsTag");
  eta_cut   = pset.getParameter<double>("eta_cut");
  pt_cut    = pset.getParameter<double>("pt_cut");
  outputFile= pset.getParameter<string>("outputFile");
  matching_cone  = pset.getParameter<double>("matching_cone");
  TriggerEventTag= pset.getParameter<InputTag>("TriggerEventTag");
  HLTriggerTag = pset.getParameter<InputTag>("HLTriggerTag");

  minvLow  = pset.getParameter<double>("minvLow");
  minvHigh = pset.getParameter<double>("minvHigh");
   
  // snippet from HLTEventAnalyzerRaw
  // unfortunately HLTEventAnalyzerRaw is not used in the production datasets
  // but only in RelVal or you need to produce them yourself
  // why CMS does not want to store l2 info?
  processName_ = pset.getParameter<std::string>("processName");
  triggerName_ = pset.getParameter<std::string>("triggerName");
  triggerResultsTag_ = pset.getParameter<edm::InputTag>("triggerResults");
  triggerEventWithRefsTag_ = pset.getParameter<edm::InputTag>("triggerEventWithRefs");
  //
 
  // L1GTUtils
  m_nameAlgTechTrig = pset.getParameter<std::string> ("AlgorithmName");
  //
  
  // get the module name you are interested in looking into
  level2module = pset.getParameter<string>("level2module");
  level3module = pset.getParameter<string>("level3module");
  
  // printLevel
  // -1 - No printout
  //  0 - only skipped collection
  //  1 - verbose
  printLevel = pset.getUntrackedParameter<int>("printLevel",0);
  
  file  = new TFile(outputFile.c_str(),"RECREATE");

  // These collections are added after Ivan's email about we need
  // to do trigger efficiency. They contain all the info about 
  // reco'd muons: global, standalone and tracker ones
  // reco Muons Collection
  recoMuons = new TTree("recoMuons", "recoMuons");
  
  //---------------------------------------------------------------------
  // general information Booking
  //--------------------------------------------------------------------- 
  recoMuons->Branch("Run",   &Run,   "Run/I"  );
  recoMuons->Branch("Event", &Event, "Event/I");
  recoMuons->Branch("Lumi",  &Lumi,  "Lumi/I" );
  recoMuons->Branch("Bx",    &Bx,    "Bx/I"   );
  recoMuons->Branch("Orbit", &Orbit, "Orbit/I");
  
  recoMuons->Branch("muonSize", &muonSize, "muonSize/I");
   
  recoMuons->Branch("isGlobalMuon"        , &isGlobalMuon       );
  recoMuons->Branch("isTrackerMuon"       , &isTrackerMuon      );
  recoMuons->Branch("isStandAloneMuon"    , &isStandAloneMuon   );
  recoMuons->Branch("isMuonAllArbitrated" , &isMuonAllArbitrated);

  recoMuons->Branch("isEnergyValid", &isEnergyValid);
  recoMuons->Branch("caloCompatibility", &caloCompatibility);
  recoMuons->Branch("em",    &em    ); 
  recoMuons->Branch("emS9",  &emS9  ); 
  recoMuons->Branch("emS25", &emS25 );
  recoMuons->Branch("emMax", &emMax );
  recoMuons->Branch("had",   &had   );
  recoMuons->Branch("hadS9", &hadS9 );
  recoMuons->Branch("hadMax",&hadMax);
  recoMuons->Branch("ho",    &ho    );
  recoMuons->Branch("hoS9",  &hoS9  );
  //---------------------------------------------------------------------

  //---------------------------------------------------------------------
  // Global Muon Block Booking
  //---------------------------------------------------------------------
  recoMuons->Branch("gmrPt"               , &gmrPt               );
  recoMuons->Branch("gmrEta"              , &gmrEta              );
  recoMuons->Branch("gmrPhi"              , &gmrPhi              );
  recoMuons->Branch("gmrP"                , &gmrP                );
  recoMuons->Branch("gmrPx"               , &gmrPx               );
  recoMuons->Branch("gmrPy"               , &gmrPy               );
  recoMuons->Branch("gmrPz"               , &gmrPz               );
  recoMuons->Branch("gmrTheta"            , &gmrTheta            );
  recoMuons->Branch("gmrVx"               , &gmrVx               );
  recoMuons->Branch("gmrVy"               , &gmrVy               );
  recoMuons->Branch("gmrVz"               , &gmrVz               );
  recoMuons->Branch("gmrCharge"           , &gmrCharge           );
  recoMuons->Branch("gmrNDoF"             , &gmrNDoF             );
  recoMuons->Branch("gmrChi2"             , &gmrChi2             );
  recoMuons->Branch("gmrChi2Norm"         , &gmrChi2Norm         );
  recoMuons->Branch("gmrDXY"              , &gmrDXY              );
  recoMuons->Branch("gmrDTheta"           , &gmrDTheta           );
  recoMuons->Branch("gmrDPt"              , &gmrDPt              );
  recoMuons->Branch("gmrDEta"             , &gmrDEta             );
  recoMuons->Branch("gmrDPhi"             , &gmrDPhi             );
  recoMuons->Branch("gmrDDXY"             , &gmrDDXY             );
  recoMuons->Branch("gmrIso03nTracks"     , &gmrIso03nTracks     );
  recoMuons->Branch("gmrIso03sumPt"       , &gmrIso03sumPt       );
  recoMuons->Branch("gmrDz"               , &gmrDz               );
  recoMuons->Branch("gmrD0"               , &gmrD0               );
  recoMuons->Branch("gmrDsz"              , &gmrDsz              );
  recoMuons->Branch("gmrDDz"              , &gmrDDz              );
  recoMuons->Branch("gmrDD0"              , &gmrDD0              );
  recoMuons->Branch("gmrDDsz"             , &gmrDDsz             );
  recoMuons->Branch("gmrInnerX"           , &gmrInnerX           );
  recoMuons->Branch("gmrInnerY"           , &gmrInnerY           );
  recoMuons->Branch("gmrInnerZ"           , &gmrInnerZ           );
  recoMuons->Branch("gmrOuterX"           , &gmrOuterX           );
  recoMuons->Branch("gmrOuterY"           , &gmrOuterY           );
  recoMuons->Branch("gmrOuterZ"           , &gmrOuterZ           );
  recoMuons->Branch("gmrValHits"          , &gmrValHits          );

  //---------------------------------------------------------------------
  // Standalone Muon Block Booking
  //---------------------------------------------------------------------
  recoMuons->Branch("stdEnergy"   , &stdEnergy   );
  recoMuons->Branch("stdPt"       , &stdPt       );
  recoMuons->Branch("stdEta"      , &stdEta      );
  recoMuons->Branch("stdPhi"      , &stdPhi      );
  recoMuons->Branch("stdPx"       , &stdPx       );
  recoMuons->Branch("stdPy"       , &stdPy       );
  recoMuons->Branch("stdPz"       , &stdPz       );
  recoMuons->Branch("stdVx"       , &stdVx       );
  recoMuons->Branch("stdVy"       , &stdVy       );
  recoMuons->Branch("stdVz"       , &stdVz       );
  recoMuons->Branch("stdCharge"   , &stdCharge   );
  recoMuons->Branch("stdDPt"      , &stdDPt      );
  recoMuons->Branch("stdDEta"     , &stdDEta     );
  recoMuons->Branch("stdDPhi"     , &stdDPhi     );
  recoMuons->Branch("stdDz"       , &stdDz       );
  recoMuons->Branch("stdD0"       , &stdD0       );
  recoMuons->Branch("stdNDoF"     , &stdNDoF     ); 
  recoMuons->Branch("stdChi2"     , &stdChi2     );
  recoMuons->Branch("stdChi2Norm" , &stdChi2Norm );
  recoMuons->Branch("stdDXY"      , &stdDXY      );
  recoMuons->Branch("stdTheta"    , &stdTheta    );
  recoMuons->Branch("stdDTheta"   , &stdDTheta   );
  recoMuons->Branch("stdDDz"      , &stdDDz      );
  recoMuons->Branch("stdDD0"      , &stdDD0      );     
  recoMuons->Branch("stdValHits"  , &stdValHits  );
   
  //---------------------------------------------------------------------
  // Tracker Muon Block Booking
  //---------------------------------------------------------------------
  recoMuons->Branch("trkEnergy"   , &trkEnergy   );
  recoMuons->Branch("trkPt"       , &trkPt       );
  recoMuons->Branch("trkEta"      , &trkEta      );
  recoMuons->Branch("trkPhi"      , &trkPhi      );
  recoMuons->Branch("trkPx"       , &trkPx       );
  recoMuons->Branch("trkPy"       , &trkPy       );
  recoMuons->Branch("trkPz"       , &trkPz       );
  recoMuons->Branch("trkVx"       , &trkVx       );
  recoMuons->Branch("trkVy"       , &trkVy       );
  recoMuons->Branch("trkVz"       , &trkVz       );
  recoMuons->Branch("trkCharge"   , &trkCharge   );
  recoMuons->Branch("trkDPt"      , &trkDPt      );
  recoMuons->Branch("trkDEta"     , &trkDEta     );
  recoMuons->Branch("trkDPhi"     , &trkDPhi     );
  recoMuons->Branch("trkDz"       , &trkDz       );
  recoMuons->Branch("trkD0"       , &trkD0       );
  recoMuons->Branch("trkNDoF"     , &trkNDoF     ); 
  recoMuons->Branch("trkChi2"     , &trkChi2     );
  recoMuons->Branch("trkChi2Norm" , &trkChi2Norm );
  recoMuons->Branch("trkDXY"      , &trkDXY      );
  recoMuons->Branch("trkTheta"    , &trkTheta    );
  recoMuons->Branch("trkDTheta"   , &trkDTheta   );
  recoMuons->Branch("trkDDz"      , &trkDDz      );
  recoMuons->Branch("trkDD0"      , &trkDD0      );     
  recoMuons->Branch("trkValHits"  , &trkValHits  );

  // CSC segment for the tracker muon
  //chamber
  recoMuons->Branch("trkNchambers"    , &trkNchambers    );
  recoMuons->Branch("trkNofMatches"   , &trkNofMatches   );
  recoMuons->Branch("trkIsMatchValid" , &trkIsMatchValid );

  recoMuons->Branch("trkNSegs"        , &trkNSegs           );
  recoMuons->Branch("trkChamberId"    , &trkSegChamberId    );
  recoMuons->Branch("trkRing"         , &trkSegRing         );
  recoMuons->Branch("trkStation"      , &trkSegStation      );
  recoMuons->Branch("trkEndcap"       , &trkSegEndcap       );
  recoMuons->Branch("trkTriggerSector", &trkSegTriggerSector);
  recoMuons->Branch("trkTriggerCscId" , &trkSegTriggerCscId );
  recoMuons->Branch("trkXfromMatch"   , &trkSegXfromMatch   );
  recoMuons->Branch("trkYfromMatch"   , &trkSegYfromMatch   );
  recoMuons->Branch("trkPhifromMatch" , &trkSegPhifromMatch );

  //segment
  recoMuons->Branch("trkSegIsArb", &trkSegIsArb);
  recoMuons->Branch("trkSegX"    , &trkSegX    );
  recoMuons->Branch("trkSegY"    , &trkSegY    );
  recoMuons->Branch("trkSegR"    , &trkSegR    );
  recoMuons->Branch("trkSegPhi"  , &trkSegPhi  );

  //---------------------------------------------------------------------
  // RECHIT information: only for standalone/global muons!
  //---------------------------------------------------------------------
  recoMuons->Branch("rchCSCtype" , &rchCSCtype ); 
  recoMuons->Branch("rchEta"     , &rchEta     ); 
  recoMuons->Branch("rchPhi"     , &rchPhi     ); 
  recoMuons->Branch("rchPhi_02PI", &rchPhi_02PI);

  recoMuons->Branch("rchStation", &rchStation);
  recoMuons->Branch("rchChamber", &rchChamber);
  recoMuons->Branch("rchRing"   , &rchRing   );
  recoMuons->Branch("rchLayer"  , &rchLayer  );
   
  recoMuons->Branch("rchMuonSize",      &rchMuonSize     );
  recoMuons->Branch("rchEtaMatrix"    , &rchEtaMatrix    );
  recoMuons->Branch("rchPhiMatrix"    , &rchPhiMatrix    );
  recoMuons->Branch("rchPhi02PIMatrix", &rchPhi02PIMatrix);
  recoMuons->Branch("rchStationMatrix", &rchStationMatrix);
  recoMuons->Branch("rchChamberMatrix", &rchChamberMatrix);
  recoMuons->Branch("rchRingMatrix"   , &rchRingMatrix   );
  recoMuons->Branch("rchLayerMatrix"  , &rchLayerMatrix  );
  recoMuons->Branch("rchTypeMatrix"   , &rchTypeMatrix   );
   
  //--------------------------------------------------------------------- 
  // old format: keep it until you are sure the TMatrixF works
  recoMuons->Branch("nMu_nCscRchHits", &nMu_nCscRchHits, "nMu_nCscRchHits/I"); 
  recoMuons->Branch("rchEtaList",       rchEtaList,      "rchEtaList[nMu_nCscRchHits]/D"); 
  recoMuons->Branch("rchPhiList",       rchPhiList,      "rchPhiList[nMu_nCscRchHits]/D");
  recoMuons->Branch("rchPhiList_02PI",  rchPhiList_02PI, "rchPhiList_02PI[nMu_nCscRchHits]/D");
   
  resizeRchHits(1);  
  //---------------------------------------------------------------------
   

  //---------------------------------------------------------------------
  // Propagation block booking
  //---------------------------------------------------------------------
  // propagation to ME1/1
  recoMuons->Branch(  "muons_x_mep11",  &muons_x_mep11);
  recoMuons->Branch(  "muons_y_mep11",  &muons_y_mep11);
  recoMuons->Branch(  "muons_z_mep11",  &muons_z_mep11);
  recoMuons->Branch("muons_phi_mep11",&muons_phi_mep11);
  recoMuons->Branch("muons_eta_mep11",&muons_eta_mep11);

  recoMuons->Branch(  "muons_x_mem11",  &muons_x_mem11);
  recoMuons->Branch(  "muons_y_mem11",  &muons_y_mem11);
  recoMuons->Branch(  "muons_z_mem11",  &muons_z_mem11);
  recoMuons->Branch("muons_phi_mem11",&muons_phi_mem11);
  recoMuons->Branch("muons_eta_mem11",&muons_eta_mem11);

  // propagation to ME1
  recoMuons->Branch(  "muons_x_mep1",  &muons_x_mep1);
  recoMuons->Branch(  "muons_y_mep1",  &muons_y_mep1);
  recoMuons->Branch(  "muons_z_mep1",  &muons_z_mep1);
  recoMuons->Branch("muons_phi_mep1",&muons_phi_mep1);
  recoMuons->Branch("muons_eta_mep1",&muons_eta_mep1);

  recoMuons->Branch(  "muons_x_mem1",  &muons_x_mem1);
  recoMuons->Branch(  "muons_y_mem1",  &muons_y_mem1);
  recoMuons->Branch(  "muons_z_mem1",  &muons_z_mem1);
  recoMuons->Branch("muons_phi_mem1",&muons_phi_mem1);
  recoMuons->Branch("muons_eta_mem1",&muons_eta_mem1);

  // propagation to ME2
  recoMuons->Branch(  "muons_x_mep2",  &muons_x_mep2);
  recoMuons->Branch(  "muons_y_mep2",  &muons_y_mep2);
  recoMuons->Branch(  "muons_z_mep2",  &muons_z_mep2);
  recoMuons->Branch("muons_phi_mep2",&muons_phi_mep2);
  recoMuons->Branch("muons_eta_mep2",&muons_eta_mep2);

  recoMuons->Branch(  "muons_x_mem2",  &muons_x_mem2);
  recoMuons->Branch(  "muons_y_mem2",  &muons_y_mem2);
  recoMuons->Branch(  "muons_z_mem2",  &muons_z_mem2);
  recoMuons->Branch("muons_phi_mem2",&muons_phi_mem2);
  recoMuons->Branch("muons_eta_mem2",&muons_eta_mem2);
   
  // propagation to ME3
  recoMuons->Branch(  "muons_x_mep3",  &muons_x_mep3);
  recoMuons->Branch(  "muons_y_mep3",  &muons_y_mep3);
  recoMuons->Branch(  "muons_z_mep3",  &muons_z_mep3);
  recoMuons->Branch("muons_phi_mep3",&muons_phi_mep3);
  recoMuons->Branch("muons_eta_mep3",&muons_eta_mep3);
   
  recoMuons->Branch(  "muons_x_mem3",  &muons_x_mem3);
  recoMuons->Branch(  "muons_y_mem3",  &muons_y_mem3);
  recoMuons->Branch(  "muons_z_mem3",  &muons_z_mem3);
  recoMuons->Branch("muons_phi_mem3",&muons_phi_mem3);
  recoMuons->Branch("muons_eta_mem3",&muons_eta_mem3);
  //---------------------------------------------------------------------

   
  //---------------------------------------------------------------------
  // l1 extra muon collection booking
  //---------------------------------------------------------------------
  l1extraMuons = new TTree("l1extraMuons", "l1extraMuons");
  l1extraMuons->Branch("l1Size", &l1Size, "l1Size/I");
    
  l1extraMuons->Branch("l1Eta", &l1Eta);
  l1extraMuons->Branch("l1Pt" , &l1Pt );
  l1extraMuons->Branch("l1Phi", &l1Phi);
   
  l1extraMuons->Branch("isIsolated"  , &isIsolated  );
  l1extraMuons->Branch("isMip"       , &isMip       );
  l1extraMuons->Branch("isForward"   , &isForward   );
  l1extraMuons->Branch("isRPC"       , &isRPC       );
  l1extraMuons->Branch("detectorType", &detectorType);
  l1extraMuons->Branch("rank"        , &rank        );
  //---------------------------------------------------------------------

  //// csctf
  //csctfTTree->Branch("EndcapTrk"    , &EndcapTrk     );
  //csctfTTree->Branch("SectorTrk"    , &SectorTrk     );
  //csctfTTree->Branch("BxTrk"  	    , &BxTrk         );
  //csctfTTree->Branch("me1ID"  	    , &me1ID         );
  //csctfTTree->Branch("me2ID"        , &me2ID 	     );
  //csctfTTree->Branch("me3ID"  	    , &me3ID         );
  //csctfTTree->Branch("me4ID"  	    , &me4ID         );
  //csctfTTree->Branch("mb1ID"  	    , &mb1ID         );
  //csctfTTree->Branch("ModeTrk"      , &ModeTrk       );
  //csctfTTree->Branch("EtaTrk"  	    , &EtaTrk        );
  //csctfTTree->Branch("PhiTrk"  	    , &PhiTrk        );
  //csctfTTree->Branch("PhiTrk_02PI"  , &PhiTrk_02PI   );
  //csctfTTree->Branch("PtTrk"  	    , &PtTrk         );
  //csctfTTree->Branch("ChargeTrk"    , &ChargeTrk     );
  //csctfTTree->Branch("ChargeValidTrk, &ChargeValidTrk);
  //csctfTTree->Branch("QualityTrk"   , &QualityTrk    );
  //csctfTTree->Branch("ForRTrk"      , &ForRTrk       );
  //csctfTTree->Branch("Phi23Trk"     , &Phi23Trk      );
  //csctfTTree->Branch("Phi12Trk"     , &Phi12Trk      );
  //csctfTTree->Branch("PhiSignTrk"   , &PhiSignTrk    );
  //csctfTTree->Branch("EtaBit"  	    , &EtaBit  	     );
  //csctfTTree->Branch("PhiBit"  	    , &PhiBit  	     );
  //csctfTTree->Branch("LocalPhi"     , &LocalPhi      );
  //csctfTTree->Branch("PtBit"  	    , &PtBit   	     );
  //csctfTTree->Branch("PtTrkData"    , &PtTrkData     );
  //csctfTTree->Branch("QualityTrkData, &QualityTrkData);
  //csctfTTree->Branch("OutputLinkTrk", &OutputLinkTrk );
   
}

// destructor
TrigEff::~TrigEff(void){ file->Write(); file->Close(); }

// analyze
void TrigEff::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){


//   Handle<GenParticleCollection> genParticles;
//   if( genTag.label() != "null" ) iEvent.getByLabel(genTag,genParticles);

//   // Find generated muons
//   vector<GenParticleRef> genMuons;
//   if( genParticles.isValid() ) {
//     for(size_t index=0; index<genParticles->size(); index++){
//       const Candidate & particle = (*genParticles)[index];
//       // Look for stable muons:
//       if( abs(particle.pdgId())==13 && particle.status()==1 ) {
//         genMuons.push_back( GenParticleRef(genParticles, index) );
//       }
//     }
//   }else {
//     if (printLevel > -1)
//       cout<<"Empty generator collection ... skipping"<<endl; 
//   }

  
  // Second: get global muons, stand alone muons, and track
  Handle<MuonCollection>  muons;
  if( muonsTag.label() != "null" ) iEvent.getByLabel(muonsTag, muons);

  //Get the Magnetic field from the setup
  iSetup.get<IdealMagneticFieldRecord>().get(theBField);
  // Get the GlobalTrackingGeometry from the setup
  iSetup.get<GlobalTrackingGeometryRecord>().get(theTrackingGeometry);

  // Get the propagators
  iSetup.get<TrackingComponentsRecord>().get("SmartPropagatorAnyRK", propagatorAlong   );
  iSetup.get<TrackingComponentsRecord>().get("SmartPropagatorAnyOpposite", propagatorOpposite);


  if (printLevel > 0){
    cout<<"**************************************************************"<<endl; 
    cout<<"**************************************************************"<<endl; 
    cout<<"New event"<<endl;
  }

  Run   = iEvent.id().run();
  Event = iEvent.id().event();
  Bx    = iEvent.bunchCrossing();   //overwritten by EVM info until fixed by fw
  Lumi  = iEvent.luminosityBlock();
  Orbit = iEvent.orbitNumber();     //overwritten by EVM info until fixed by fw
  
  // ============================================================================
  // Now start with creating array of di-muon pairs on the 
  // reconstruction level (called true)

  if (printLevel > 0) cout<<"============ DI MUONS PAIRS ================"<<endl; 
  map< MuonRef, vector< pair< MuonRef, double > > > trueDiMuons;
  
  if( muons.isValid() ){

    muonsInit();

    //GP start
    // Get the CSC Geometry 
    iSetup.get<MuonGeometryRecord>().get(cscGeom);
    //GP end  

    //fill the muons
    fillMuons(muons);
      
    recoMuons->Fill();
    muonsDel();
  }

 if (printLevel > 0) {
    cout<<"======================================================"<<endl; 
    cout<<"Level-1: "<<endl;
  }

  Handle<l1extra::L1MuonParticleCollection> l1muons;
  if( L1extraTag.label() != "null" ) iEvent.getByLabel(L1extraTag, l1muons);

  // Assign similar container, mapping muon to pairs of track and 
  // matching Level-1 candidate, if the last exists
  map< MuonRef, vector< pair<TrackRef,l1extra::L1MuonParticleCollection::const_iterator> > > diMuonCandidates_l1;
  
  // Fetch all Level-1 muon candidates
  if( l1muons.isValid() ){


    // ============================================================================
    // fill all the l1 information per event in a ttree
    if (printLevel > 0) {
      printf("============================================================\n");
      printf("%-20s: %d\n\n", "l1muons->size()"     , l1muons->size()        );
      printf("  %-30s: %d\n", "isl1MuonsValid()"    , l1muons.isValid()      );  
    }

    l1Size = l1muons->size();
    l1extraInit();
    
    //fill the l1 extra muons
    for(l1extra::L1MuonParticleCollection::const_iterator iter=l1muons->begin(); iter!=l1muons->end(); iter++) {
      // just print out in case
      if (printLevel > 0)
        cout<<"   l1 candidate: eta=" <<iter->eta()<<" phi=" <<iter->phi()<<endl;
      
      // where the filling is happening
      fillExtra(iter);
    }
    
    l1extraMuons->Fill();
    l1extraDel();
  } 
  
}


//-------------------------------------------------------------------------
// methods for the reco muons collection
void TrigEff::fillMuons(const edm::Handle<reco::MuonCollection> muons)
{
  
  if (printLevel > 0) {
    printf("============================================================\n");
    printf("%-20s: %d\n\n", "muons->size()"     , muons->size()            );
    printf("  %-30s: %d\n", "isMuonsValid()"    , muons.isValid()          );  
  }
  
  muonSize = muons->size();
  
  // Once you know how many muons we have in the event,
  // we resize the rchEtaList and rchPhiList
  
  
  nMu_nCscRchHits = muonSize * MAX_CSC_RECHIT;
  if (printLevel > 0) 
    cout << "nMu_nCscRchHits: " << nMu_nCscRchHits << endl;
  resizeRchHits(nMu_nCscRchHits);
  
  
  int whichMuon =-1;
  for (MuonCollection::const_iterator muon=muons->begin();
       muon!=muons->end(); muon++) {
    
    whichMuon++;

    if (printLevel > 0) {
      printf("  %-30s\n"    , "--------------------------------");
      printf("  %-30s: %d\n", "isGlobalMuon    ()", muon->isGlobalMuon    ());
      printf("  %-30s: %d\n", "isTrackerMuon   ()", muon->isTrackerMuon   ());
      printf("  %-30s: %d\n", "isStandAloneMuon()", muon->isStandAloneMuon());
      printf("  %-30s: %d\n", "combinedMuon    ().isNonnull()", muon->combinedMuon  ().isNonnull());
      printf("  %-30s: %d\n", "track           ().isNonnull()", muon->track         ().isNonnull());
      printf("  %-30s: %d\n", "standAloneMuon  ().isNonnull()", muon->standAloneMuon().isNonnull());
      printf("  %-30s\n\n"  , "--------------------------------");
    }
    
    // old left over from a previous ntuplizer
    if (muon->isIsolationValid()) {
      MuonIsolation isoR03 = muon->isolationR03();
      gmrIso03nTracks->push_back(isoR03.nTracks);
      gmrIso03sumPt  ->push_back(isoR03.sumPt  );
    }
    else {
      gmrIso03nTracks->push_back(-999);
      gmrIso03sumPt  ->push_back(-999);
    }
    //---------------------------------------------

    isGlobalMuon        -> push_back(muon->isGlobalMuon    ()); 
    isTrackerMuon       -> push_back(muon->isTrackerMuon   ()); 
    isStandAloneMuon    -> push_back(muon->isStandAloneMuon());
    isMuonAllArbitrated -> push_back(muon::isGoodMuon(*muon, muon::AllArbitrated));

    isEnergyValid       -> push_back(muon->isEnergyValid());
    caloCompatibility   -> push_back(muon->caloCompatibility());
    em                  -> push_back(muon->calEnergy().em    ); 
    emS9                -> push_back(muon->calEnergy().emS9  );
    emS25               -> push_back(muon->calEnergy().emS25 );
    emMax               -> push_back(muon->calEnergy().emMax );
    had                 -> push_back(muon->calEnergy().had   );
    hadS9               -> push_back(muon->calEnergy().hadS9 );
    hadMax              -> push_back(muon->calEnergy().hadMax);
    ho                  -> push_back(muon->calEnergy().ho    );
    hoS9                -> push_back(muon->calEnergy().hoS9  );

    // global muon
    if (muon->combinedMuon().isNonnull()) {

      TrackRef trackRef = muon->combinedMuon();

      if (printLevel>0) 
        cout << "(GBL) muon->pt(): " << trackRef->pt() 
             << ", muon->eta(): "    << trackRef->eta() 
             << ", muon->phi(): "    << trackRef->phi() << endl;
              
      gmrPt      -> push_back(trackRef->pt    ());
      gmrPhi     -> push_back(trackRef->phi   ());
      gmrEta     -> push_back(trackRef->eta   ());
      gmrP       -> push_back(trackRef->p     ());
      gmrPx      -> push_back(trackRef->px    ());
      gmrPy      -> push_back(trackRef->py    ());
      gmrPz      -> push_back(trackRef->pz    ());
      gmrTheta   -> push_back(trackRef->theta ());
      gmrCharge  -> push_back(trackRef->charge());
      gmrVx      -> push_back(trackRef->vx    ());
      gmrVy      -> push_back(trackRef->vy    ());
      gmrVz      -> push_back(trackRef->vz    ());
      gmrDXY     -> push_back(trackRef->dxy              ());
      gmrDDXY    -> push_back(trackRef->dxyError         ());
      gmrDPhi    -> push_back(trackRef->phiError         ());
      gmrDEta    -> push_back(trackRef->etaError         ());
      gmrDPt     -> push_back(trackRef->ptError          ());
      gmrChi2    -> push_back(trackRef->chi2             ());
      gmrNDoF    -> push_back(trackRef->ndof             ());
      gmrChi2Norm-> push_back(trackRef->normalizedChi2   ());
      gmrDTheta  -> push_back(trackRef->thetaError       ());
      gmrDz      -> push_back(trackRef->dz               ());
      gmrD0      -> push_back(trackRef->d0               ());
      gmrDsz     -> push_back(trackRef->dsz              ());
      gmrDDz     -> push_back(trackRef->dzError          ());
      gmrDD0     -> push_back(trackRef->d0Error          ());
      gmrDDsz    -> push_back(trackRef->dszError         ());
      gmrValHits -> push_back(trackRef->numberOfValidHits());

      gmrInnerX  -> push_back(trackRef->innerPosition().X());
      gmrInnerY  -> push_back(trackRef->innerPosition().Y());
      gmrInnerZ  -> push_back(trackRef->innerPosition().Z());   
  
      gmrOuterX  -> push_back(trackRef->outerPosition().X());
      gmrOuterY  -> push_back(trackRef->outerPosition().Y());
      gmrOuterZ  -> push_back(trackRef->outerPosition().Z());   

    }
    else {
      gmrPt      -> push_back(-999);
      gmrPhi     -> push_back(-999);
      gmrEta     -> push_back(-999);
      gmrP       -> push_back(-999);
      gmrPx      -> push_back(-999);
      gmrPy      -> push_back(-999);
      gmrPz      -> push_back(-999);
      gmrTheta   -> push_back(-999);
      gmrCharge  -> push_back(-999);
      gmrVx      -> push_back(-999);
      gmrVy      -> push_back(-999);
      gmrVz      -> push_back(-999);
      gmrDXY     -> push_back(-999);
      gmrDDXY    -> push_back(-999);
      gmrDPhi    -> push_back(-999);
      gmrDEta    -> push_back(-999);
      gmrDPt     -> push_back(-999);
      gmrChi2    -> push_back(-999);
      gmrNDoF    -> push_back(-999);
      gmrChi2Norm-> push_back(-999);
      gmrDTheta  -> push_back(-999);
      gmrDz      -> push_back(-999);
      gmrD0      -> push_back(-999);
      gmrDsz     -> push_back(-999);
      gmrDDz     -> push_back(-999);
      gmrDD0     -> push_back(-999);
      gmrDDsz    -> push_back(-999);
      gmrValHits -> push_back(-999);

      gmrInnerX  -> push_back(-999);
      gmrInnerY  -> push_back(-999);
      gmrInnerZ  -> push_back(-999);   
                                  
      gmrOuterX  -> push_back(-999);
      gmrOuterY  -> push_back(-999);
      gmrOuterZ  -> push_back(-999);   
    }

    
    // standalone muon
    if (muon->standAloneMuon().isNonnull()) {
      TrackRef trackRefStd = muon->standAloneMuon();

      if (printLevel>0) 
        cout << "(STA) muon->pt(): " << trackRefStd->pt() 
             << ", muon->eta(): "    << trackRefStd->eta() 
             << ", muon->phi(): "    << trackRefStd->phi() << endl;

      stdPt      -> push_back(trackRefStd->pt       ());
      stdEta     -> push_back(trackRefStd->eta      ());
      stdPhi     -> push_back(trackRefStd->phi      ());
      stdPx      -> push_back(trackRefStd->px       ());
      stdPy      -> push_back(trackRefStd->py       ());
      stdPz      -> push_back(trackRefStd->pz       ());
      stdVx      -> push_back(trackRefStd->vx       ());
      stdVy      -> push_back(trackRefStd->vy       ());
      stdVz      -> push_back(trackRefStd->vz       ());
      stdCharge  -> push_back(trackRefStd->charge   ());
      stdDPt     -> push_back(trackRefStd->ptError  ());
      stdDEta    -> push_back(trackRefStd->etaError ());
      stdDPhi    -> push_back(trackRefStd->phiError ());
      stdDz      -> push_back(trackRefStd->dz       ());
      stdD0      -> push_back(trackRefStd->d0       ());

      stdChi2    -> push_back(trackRefStd->chi2             ());
      stdNDoF    -> push_back(trackRefStd->ndof             ());
      stdChi2Norm-> push_back(trackRefStd->normalizedChi2   ());
      stdTheta   -> push_back(trackRefStd->theta            ());
      stdDTheta  -> push_back(trackRefStd->thetaError       ());
      stdDDz     -> push_back(trackRefStd->dzError          ());
      stdDD0     -> push_back(trackRefStd->d0Error          ());
      stdValHits -> push_back(trackRefStd->numberOfValidHits());

      stdInnerX  -> push_back(trackRefStd->innerPosition().X());
      stdInnerY  -> push_back(trackRefStd->innerPosition().Y());
      stdInnerZ  -> push_back(trackRefStd->innerPosition().Z());   

      stdOuterX  -> push_back(trackRefStd->outerPosition().X());
      stdOuterY  -> push_back(trackRefStd->outerPosition().Y());
      stdOuterZ  -> push_back(trackRefStd->outerPosition().Z());   
    }
    else {
      stdPt      -> push_back(-999);
      stdEta     -> push_back(-999);
      stdPhi     -> push_back(-999);
      stdPx      -> push_back(-999);
      stdPy      -> push_back(-999);
      stdPz      -> push_back(-999);
      stdVx      -> push_back(-999);
      stdVy      -> push_back(-999);
      stdVz      -> push_back(-999);
      stdCharge  -> push_back(-999);
      stdDPt     -> push_back(-999);
      stdDEta    -> push_back(-999);
      stdDPhi    -> push_back(-999);
      stdDz      -> push_back(-999);
      stdD0      -> push_back(-999);

      stdChi2    -> push_back(-999);
      stdNDoF    -> push_back(-999);
      stdChi2Norm-> push_back(-999);
      stdTheta   -> push_back(-999);
      stdDTheta  -> push_back(-999);
      stdDDz     -> push_back(-999);
      stdDD0     -> push_back(-999);
      stdValHits -> push_back(-999);

      stdInnerX  -> push_back(-999);
      stdInnerY  -> push_back(-999);
      stdInnerZ  -> push_back(-999);   
                              
      stdOuterX  -> push_back(-999);
      stdOuterY  -> push_back(-999);
      stdOuterZ  -> push_back(-999);   
      
    }
    
    //tracker muon    
    if (muon->track().isNonnull()) {

      TrackRef trackRefTrk = muon->track();

      if (printLevel>0) 
        cout << "(TRK) muon->pt(): " << trackRefTrk->pt() 
             << ", muon->eta(): "    << trackRefTrk->eta() 
             << ", muon->phi(): "    << trackRefTrk->phi() << endl;

      trkPt     ->push_back(trackRefTrk->pt       ());
      trkEta    ->push_back(trackRefTrk->eta      ());
      trkPhi    ->push_back(trackRefTrk->phi      ());
      trkPx     ->push_back(trackRefTrk->px       ());
      trkPy     ->push_back(trackRefTrk->py       ());
      trkPz     ->push_back(trackRefTrk->pz       ());
      trkVx     ->push_back(trackRefTrk->vx       ());
      trkVy     ->push_back(trackRefTrk->vy       ());
      trkVz     ->push_back(trackRefTrk->vz       ());
      trkCharge ->push_back(trackRefTrk->charge   ());
      trkDPt    ->push_back(trackRefTrk->ptError  ());
      trkDEta   ->push_back(trackRefTrk->etaError ());
      trkDPhi   ->push_back(trackRefTrk->phiError ());
      trkDz     ->push_back(trackRefTrk->dz       ());
      trkD0     ->push_back(trackRefTrk->d0       ());
      
      trkChi2    ->push_back(trackRefTrk->chi2             ());
      trkNDoF    ->push_back(trackRefTrk->ndof             ());
      trkChi2Norm->push_back(trackRefTrk->normalizedChi2   ());
      trkTheta   ->push_back(trackRefTrk->theta            ());
      trkDTheta  ->push_back(trackRefTrk->thetaError       ());
      trkDDz     ->push_back(trackRefTrk->dzError          ());
      trkDD0     ->push_back(trackRefTrk->d0Error          ());
      trkValHits ->push_back(trackRefTrk->numberOfValidHits());
    }
    else {
      trkPt     ->push_back(-999);
      trkEta    ->push_back(-999);
      trkPhi    ->push_back(-999);
      trkPx     ->push_back(-999);
      trkPy     ->push_back(-999);
      trkPz     ->push_back(-999);
      trkVx     ->push_back(-999);
      trkVy     ->push_back(-999);
      trkVz     ->push_back(-999);
      trkCharge ->push_back(-999);
      trkDPt    ->push_back(-999);
      trkDEta   ->push_back(-999);
      trkDPhi   ->push_back(-999);
      trkDz     ->push_back(-999);
      trkD0     ->push_back(-999);

      trkChi2    ->push_back(-999);
      trkNDoF    ->push_back(-999);
      trkChi2Norm->push_back(-999);
      trkTheta   ->push_back(-999);
      trkDTheta  ->push_back(-999);
      trkDDz     ->push_back(-999);
      trkDD0     ->push_back(-999);
      trkValHits ->push_back(-999);
    }


    // --------------------------------------------------------------------------
    // propagation...
    // --------------------------------------------------------------------------
    TrajectoryStateOnSurface tsos;

    // which TrackRef to use?
    TrackRef trackRef;
    // get the tracker muon
    if (muon->combinedMuon  ().isNonnull() ||
        muon->track         ().isNonnull()  )    trackRef = muon->track();
    else if (muon->standAloneMuon().isNonnull()) trackRef = muon->standAloneMuon();


    // ... to ME1/1
    // track at ME+1/1 surface, 5.9 m - extrapolation
    tsos = surfExtrapTrkSam(trackRef, 590);  
    if (tsos.isValid()) {
      double xx = tsos.globalPosition().x();
      double yy = tsos.globalPosition().y();
      double zz = tsos.globalPosition().z();
      
      muons_x_mep11->push_back(xx);	
      muons_y_mep11->push_back(yy);	
      muons_z_mep11->push_back(zz);	

      double rr = sqrt(xx*xx + yy*yy);
      double cosphi = xx/rr;
      
      if (yy>=0) muons_phi_mep11->push_back(acos(cosphi));
      else       muons_phi_mep11->push_back(2*PI-acos(cosphi));

    
      double abspseta = -log( tan( atan(fabs(rr/zz))/2.0 ) );
      if (zz>=0) muons_eta_mep11->push_back(abspseta);
      else       muons_eta_mep11->push_back(-abspseta);

      if (printLevel>0) {
        cout<<"I am projection the track to the ME+1/1 surface" << endl;
        cout << "muons_x_mep11:" << xx << endl;
        cout << "muons_y_mep11:" << yy << endl;
        cout << "muons_z_mep11:" << zz << endl;      
        if (yy>=0) cout << "muons_phi_mep11:" << acos(cosphi) << endl;
        else       cout << "muons_phi_mep11:" << 2*PI-acos(cosphi) << endl;        
        if (zz>=0) cout << "muons_eta_mep11:" <<  abspseta << endl;
        else       cout << "muons_eta_mep11:" << -abspseta << endl;
      }
      
    }

    // track at ME-1/1 surface, -5.9 m - extrapolation
    tsos = surfExtrapTrkSam(trackRef, -590);  
    if (tsos.isValid()) {
      double xx = tsos.globalPosition().x();
      double yy = tsos.globalPosition().y();
      double zz = tsos.globalPosition().z();
      
      muons_x_mem11->push_back(xx);	
      muons_y_mem11->push_back(yy);	
      muons_z_mem11->push_back(zz);	

      double rr = sqrt(xx*xx + yy*yy);
      double cosphi = xx/rr;
      
      if (yy>=0) muons_phi_mem11->push_back(acos(cosphi));
      else       muons_phi_mem11->push_back(2*PI-acos(cosphi));

    
      double abspseta = -log( tan( atan(fabs(rr/zz))/2.0 ) );
      if (zz>=0) muons_eta_mem11->push_back(abspseta);
      else       muons_eta_mem11->push_back(-abspseta);

      if (printLevel>0) {
        cout<<"I am projection the track to the ME-1/1 surface" << endl;
        cout << "muons_x_mem11:" << xx << endl;
        cout << "muons_y_mem11:" << yy << endl;
        cout << "muons_z_mem11:" << zz << endl;      
        if (yy>=0) cout << "muons_phi_mem11:" << acos(cosphi) << endl;
        else       cout << "muons_phi_mem11:" << 2*PI-acos(cosphi) << endl;        
        if (zz>=0) cout << "muons_eta_mem11:" <<  abspseta << endl;
        else       cout << "muons_eta_mem11:" << -abspseta << endl;
      }
      
    }
    
    // ... to ME1
    // track at ME+1 surface, 7.1 m - extrapolation
    tsos = surfExtrapTrkSam(trackRef, 710);  
    if (tsos.isValid()) {
      double xx = tsos.globalPosition().x();
      double yy = tsos.globalPosition().y();
      double zz = tsos.globalPosition().z();
      
      muons_x_mep1->push_back(xx);	
      muons_y_mep1->push_back(yy);	
      muons_z_mep1->push_back(zz);	

      double rr = sqrt(xx*xx + yy*yy);
      double cosphi = xx/rr;
      
      if (yy>=0) muons_phi_mep1->push_back(acos(cosphi));
      else       muons_phi_mep1->push_back(2*PI-acos(cosphi));

    
      double abspseta = -log( tan( atan(fabs(rr/zz))/2.0 ) );
      if (zz>=0) muons_eta_mep1->push_back(abspseta);
      else       muons_eta_mep1->push_back(-abspseta);

      if (printLevel>0) {
        cout<<"I am projection the track to the ME+1/(2-3) surface" << endl;
        cout << "muons_x_mep1:" << xx << endl;
        cout << "muons_y_mep1:" << yy << endl;
        cout << "muons_z_mep1:" << zz << endl;      
        if (yy>=0) cout << "muons_phi_mep1:" << acos(cosphi) << endl;
        else       cout << "muons_phi_mep1:" << 2*PI-acos(cosphi) << endl;        
        if (zz>=0) cout << "muons_eta_mep1:" <<  abspseta << endl;
        else       cout << "muons_eta_mep1:" << -abspseta << endl;
      }
      
    }

    // track at ME-1 surface, -7.1 m - extrapolation
    tsos = surfExtrapTrkSam(trackRef, -710);  
    if (tsos.isValid()) {
      double xx = tsos.globalPosition().x();
      double yy = tsos.globalPosition().y();
      double zz = tsos.globalPosition().z();
      
      muons_x_mem1->push_back(xx);	
      muons_y_mem1->push_back(yy);	
      muons_z_mem1->push_back(zz);	

      double rr = sqrt(xx*xx + yy*yy);
      double cosphi = xx/rr;
      
      if (yy>=0) muons_phi_mem1->push_back(acos(cosphi));
      else       muons_phi_mem1->push_back(2*PI-acos(cosphi));

    
      double abspseta = -log( tan( atan(fabs(rr/zz))/2.0 ) );
      if (zz>=0) muons_eta_mem1->push_back(abspseta);
      else       muons_eta_mem1->push_back(-abspseta);

      if (printLevel>0) {
        cout<<"I am projection the track to the ME-1/(2-3) surface" << endl;
        cout << "muons_x_mem1:" << xx << endl;
        cout << "muons_y_mem1:" << yy << endl;
        cout << "muons_z_mem1:" << zz << endl;      
        if (yy>=0) cout << "muons_phi_mem1:" << acos(cosphi) << endl;
        else       cout << "muons_phi_mem1:" << 2*PI-acos(cosphi) << endl;        
        if (zz>=0) cout << "muons_eta_mem1:" <<  abspseta << endl;
        else       cout << "muons_eta_mem1:" << -abspseta << endl;
      }
      
    }


    // ... to ME2
    // track at ME+2 surface, 8.5 m - extrapolation
    tsos = surfExtrapTrkSam(trackRef, 850);  
    if (tsos.isValid()) {
      double xx = tsos.globalPosition().x();
      double yy = tsos.globalPosition().y();
      double zz = tsos.globalPosition().z();
      
      muons_x_mep2->push_back(xx);	
      muons_y_mep2->push_back(yy);	
      muons_z_mep2->push_back(zz);	

      double rr = sqrt(xx*xx + yy*yy);
      double cosphi = xx/rr;
      
      if (yy>=0) muons_phi_mep2->push_back(acos(cosphi));
      else       muons_phi_mep2->push_back(2*PI-acos(cosphi));

    
      double abspseta = -log( tan( atan(fabs(rr/zz))/2.0 ) );
      if (zz>=0) muons_eta_mep2->push_back(abspseta);
      else       muons_eta_mep2->push_back(-abspseta);

      if (printLevel>0) {
        cout<<"I am projection the track to the ME+2 surface" << endl;
        cout << "muons_x_mep2:" << xx << endl;
        cout << "muons_y_mep2:" << yy << endl;
        cout << "muons_z_mep2:" << zz << endl;      
        if (yy>=0) cout << "muons_phi_mep2:" << acos(cosphi) << endl;
        else       cout << "muons_phi_mep2:" << 2*PI-acos(cosphi) << endl;        
        if (zz>=0) cout << "muons_eta_mep2:" << abspseta << endl;
        else       cout << "muons_eta_mep2:" << -abspseta << endl;
      }
      
    }

    // track at ME-2 surface, -8.5 m - extrapolation
    tsos = surfExtrapTrkSam(trackRef, -850);  
    if (tsos.isValid()) {
      double xx = tsos.globalPosition().x();
      double yy = tsos.globalPosition().y();
      double zz = tsos.globalPosition().z();
      
      muons_x_mem2->push_back(xx);	
      muons_y_mem2->push_back(yy);	
      muons_z_mem2->push_back(zz);	

      double rr = sqrt(xx*xx + yy*yy);
      double cosphi = xx/rr;
      
      if (yy>=0) muons_phi_mem2->push_back(acos(cosphi));
      else       muons_phi_mem2->push_back(2*PI-acos(cosphi));

    
      double abspseta = -log( tan( atan(fabs(rr/zz))/2.0 ) );
      if (zz>=0) muons_eta_mem2->push_back(abspseta);
      else       muons_eta_mem2->push_back(-abspseta);

      if (printLevel>0) {
        cout<<"I am projection the track to the ME+2 surface" << endl;
        cout << "muons_x_mem2:" << xx << endl;
        cout << "muons_y_mem2:" << yy << endl;
        cout << "muons_z_mem2:" << zz << endl;      
        if (yy>=0) cout << "muons_phi_mem2:" << acos(cosphi) << endl;
        else       cout << "muons_phi_mem2:" << 2*PI-acos(cosphi) << endl;        
        if (zz>=0) cout << "muons_eta_mem2:" << abspseta << endl;
        else       cout << "muons_eta_mem2:" << -abspseta << endl;
      }
      
    }


    // ... to ME3
    // track at ME+3 surface, 9.7 m - extrapolation
    tsos = surfExtrapTrkSam(trackRef, 970);  
    if (tsos.isValid()) {
      double xx = tsos.globalPosition().x();
      double yy = tsos.globalPosition().y();
      double zz = tsos.globalPosition().z();
      
      muons_x_mep3->push_back(xx);	
      muons_y_mep3->push_back(yy);	
      muons_z_mep3->push_back(zz);	

      double rr = sqrt(xx*xx + yy*yy);
      double cosphi = xx/rr;
      
      if (yy>=0) muons_phi_mep3->push_back(acos(cosphi));
      else       muons_phi_mep3->push_back(2*PI-acos(cosphi));

    
      double abspseta = -log( tan( atan(fabs(rr/zz))/2.0 ) );
      if (zz>=0) muons_eta_mep3->push_back(abspseta);
      else       muons_eta_mep3->push_back(-abspseta);

      if (printLevel>0) {
        cout<<"I am projection the track to the ME+2 surface" << endl;
        cout << "muons_x_mep3:" << xx << endl;
        cout << "muons_y_mep3:" << yy << endl;
        cout << "muons_z_mep3:" << zz << endl;      
        if (yy>=0) cout << "muons_phi_mep3:" << acos(cosphi) << endl;
        else       cout << "muons_phi_mep3:" << 2*PI-acos(cosphi) << endl;        
        if (zz>=0) cout << "muons_eta_mep3:" << abspseta << endl;
        else       cout << "muons_eta_mep3:" << -abspseta << endl;
      }
      
    }

    // track at ME-3 surface, -9.7 m - extrapolation
    tsos = surfExtrapTrkSam(trackRef, -970);  
    if (tsos.isValid()) {
      double xx = tsos.globalPosition().x();
      double yy = tsos.globalPosition().y();
      double zz = tsos.globalPosition().z();
      
      muons_x_mem3->push_back(xx);	
      muons_y_mem3->push_back(yy);	
      muons_z_mem3->push_back(zz);	

      double rr = sqrt(xx*xx + yy*yy);
      double cosphi = xx/rr;
      
      if (yy>=0) muons_phi_mem3->push_back(acos(cosphi));
      else       muons_phi_mem3->push_back(2*PI-acos(cosphi));

    
      double abspseta = -log( tan( atan(fabs(rr/zz))/2.0 ) );
      if (zz>=0) muons_eta_mem3->push_back(abspseta);
      else       muons_eta_mem3->push_back(-abspseta);

      if (printLevel>0) {
        cout<<"I am projection the track to the ME+2 surface" << endl;
        cout << "muons_x_mem3:" << xx << endl;
        cout << "muons_y_mem3:" << yy << endl;
        cout << "muons_z_mem3:" << zz << endl;      
        if (yy>=0) cout << "muons_phi_mem3:" << acos(cosphi) << endl;
        else       cout << "muons_phi_mem3:" << 2*PI-acos(cosphi) << endl;        
        if (zz>=0) cout << "muons_eta_mem3:" << abspseta << endl;
        else       cout << "muons_eta_mem3:" << -abspseta << endl;
      }
      
    }
    
    //---------------------------------------------------------------------
    // RECHIT information in CSC: only for standalone/global muons!
    //---------------------------------------------------------------------
    // An artificial rank to sort RecHits 
    // the closer RecHit to key station/layer -> the smaller the rank
    // KK's original idea
    const int lutCSC[4][6] = { {26,24,22,21,23,25} , 
			       { 6, 4, 2, 1, 3, 5} , 
			       {16,14,12,11,13,15} , 
			       {36,34,32,31,33,35} };
    
    if( muon->isGlobalMuon()     ||
        muon->isStandAloneMuon()  ){

      // var needed to register the better ranked muon
       int   globalTypeRCH = -999;
       double globalEtaRCH = -999; 
       double globalPhiRCH = -999;
      int globalStationRCH = -999;
      int globalChamberRCH = -999;
         int globalRingRCH = -999;
        int globalLayerRCH = -999;
      
      int iRechit = 0;
      for(trackingRecHit_iterator hit  = muon->outerTrack()->recHitsBegin(); 
 	  hit != muon->outerTrack()->recHitsEnd(); 
 	  hit++){
        
        // some basic protection
	if ( !((*hit)->isValid()) ) continue;
        
	// Hardware ID of the RecHit (in terms of wire/strips/chambers)
	DetId detid = (*hit)->geographicalId();
        
	// Interested in muon systems only
	if( detid.det() != DetId::Muon ) continue;
        //Look only at CSC Hits (CSC id is 2)
        if (detid.subdetId() != MuonSubdetId::CSC) continue;
          
        CSCDetId id(detid.rawId());
          
        // another sanity check
        if  (id.station() < 1) continue;
        if  (id.layer()   < 1) continue;
          
        // Look up some stuff specific to CSCRecHit2D
        const CSCRecHit2D* CSChit =dynamic_cast<const CSCRecHit2D*>(&**hit);
        
        LocalPoint rhitlocal = CSChit->localPosition();
        GlobalPoint gp = GlobalPoint(0.0, 0.0, 0.0);
        
        const CSCChamber* cscchamber = cscGeom->chamber(id);
        
        if (!cscchamber) continue;
          
        gp = cscchamber->toGlobal(rhitlocal);
        
        // identify the rechit position
        int pos = ( ((id.station()-1)*6) + (id.layer()-1) ) + (MAX_CSC_RECHIT*whichMuon);
        
        // --------------------------------------------------
        // this part has to be deprecated once we are sure of 
        // the TMatrixF usage ;)
        // fill the rechits array
        rchEtaList[pos]      = gp.eta();
        rchPhiList[pos]      = gp.phi();
        float phi02PI = gp.phi();
        if (gp.phi() < 0) phi02PI += (2*PI); 
        
        rchPhiList_02PI[pos] = phi02PI;
        
        if (printLevel > 0) {
          cout << "iRechit: " << iRechit+1 
               << " -> pos: " << pos;
          cout << "   - RCH Type: "         << lutCSC[id.station()-1][id.layer()-1]
               <<      " RCH Eta: "         << gp.eta()	    
               <<      " RCH Phi: "         << gp.phi()
               <<      " RCH Phi [0,2PI]: " << phi02PI 
               << endl;
        }
        
        // -------------------------------------------------- 
        // new Matrix block
        if (whichMuon < MAX_MUONS && iRechit < MAX_CSC_RECHIT) {
          rchEtaMatrix[whichMuon][iRechit]     = gp.eta(); 
          rchPhiMatrix[whichMuon][iRechit]     = gp.phi();
          rchPhi02PIMatrix[whichMuon][iRechit] = phi02PI;
          rchStationMatrix[whichMuon][iRechit] = id.station();
          rchChamberMatrix[whichMuon][iRechit] = id.chamber();
          rchRingMatrix[whichMuon][iRechit]    = id.ring();
          rchLayerMatrix[whichMuon][iRechit]   = id.layer();
          rchTypeMatrix[whichMuon][iRechit]    = lutCSC[id.station()-1][id.layer()-1];
        }
        else
          cout << "ERROR: too many segment or muons "
               << "MAX_MUONS is currently set to "      << MAX_MUONS      << endl
               << "whichMuon loop is currently at "     << whichMuon      << endl
               << "MAX_CSC_RECHIT is currently set to " << MAX_CSC_RECHIT << endl
               << "iRechit loop is currently at "       << iRechit        << endl;
        
        iRechit++;
        // --------------------------------------------------
        
        
        // --------------------------------------------------
        // See if this hit is closer to the "key" position
        if( lutCSC[id.station()-1][id.layer()-1]<globalTypeRCH || globalTypeRCH<0 ){
          globalTypeRCH    = lutCSC[id.station()-1][id.layer()-1];
          globalEtaRCH     = gp.eta();
          globalPhiRCH     = gp.phi();
          globalStationRCH = id.station();
          globalChamberRCH = id.chamber();
          globalRingRCH    = id.ring();
          globalLayerRCH   = id.layer();
        }
        // -------------------------------------------------- 
        
      }//end loop rechit
    
      rchMuonSize->push_back(iRechit);

      // at the end of the loop, write only the best rechit
      rchCSCtype  -> push_back(globalTypeRCH);       
      rchEta      -> push_back(globalEtaRCH);     
      rchPhi      -> push_back(globalPhiRCH);     
      rchStation  -> push_back(globalStationRCH);
      rchChamber  -> push_back(globalChamberRCH);
      rchRing     -> push_back(globalRingRCH);
      rchLayer    -> push_back(globalLayerRCH);
    
      if (printLevel > 0) {
        cout << "\n######### CLOSER #########";
        cout << "\n - RCH Type: " << globalTypeRCH
             <<     " RCH Eta: "  << globalEtaRCH	    
             <<     " RCH Phi: "  << globalPhiRCH;
      }
      
      if ( globalPhiRCH < 0) {
        rchPhi_02PI -> push_back(globalPhiRCH + (2*PI));
 	if (printLevel > 0) 
 	  cout << " RCH Phi [0,2PI]: " << globalPhiRCH + (2*PI) << endl;
      }
      else {
        rchPhi_02PI -> push_back(globalPhiRCH         );
 	if (printLevel > 0) 
 	  cout << " RCH Phi [0,2PI]: " << globalPhiRCH          << endl;
      }
      
      if (printLevel > 0) 
	cout << "##########################\n\n";
      
    }//isGlobal || isStandAlone
    else {
      // if only tracker muons 
      rchCSCtype  -> push_back(-999);       
      rchEta      -> push_back(-999);     
      rchPhi      -> push_back(-999);     
      rchPhi_02PI -> push_back(-999);
      rchStation  -> push_back(-999);
      rchChamber  -> push_back(-999);
      rchRing     -> push_back(-999);
      rchLayer    -> push_back(-999);
    }
    
    
    // only global and tracker muons 
    if( (muon->isGlobalMuon()     == 0) ||
        (muon->isTrackerMuon()    == 1)  ){
      
      if (printLevel > 0) {
        cout << "Number of Chambers:" << muon->numberOfChambers() << endl;
        cout << "Number of Matches:"  << muon->numberOfMatches()  << endl;
        cout << "StationMask: "       << muon->stationMask()      << endl;
        cout << "Is Matches Valid?:"  << muon->isMatchesValid()   << endl;
      }
      
      trkNchambers->push_back(muon->numberOfChambers());
      trkNofMatches->push_back(muon->numberOfMatches());
      trkIsMatchValid->push_back(muon->isMatchesValid());

      if (printLevel > 0)
        cout << "IsGoodMuon? " 
             << muon::isGoodMuon(*muon, muon::AllArbitrated) << endl;
      
      // fill the vectors only if the muon is arbitrated
      if (muon::isGoodMuon(*muon, muon::AllArbitrated)) {

        int iChamber=0;
        int iSegment=0;
        for( std::vector<MuonChamberMatch>::const_iterator chamber = muon->matches().begin();
             chamber != muon->matches().end(); ++chamber) {
          
          if (printLevel > 0) cout << "Chamber:" << iChamber+1 << endl;
          
          // we are interested only in CSC chamber
          if( chamber->detector() != MuonSubdetId::CSC ) continue; 

          //get the handle to the info, CSCDetId            
          CSCDetId cscId(chamber->id.rawId());
            
          // extract the information
          int chamberId     = cscId.chamber();
          int ring          = cscId.ring();
          int station       = cscId.station();
          int endcap        = cscId.endcap();
          int triggerSector = cscId.triggerSector();
          int triggerCscId  = cscId.triggerCscId();
          
          float xx = chamber->x;
          float yy = chamber->y;
          
          float rr = sqrt(xx*xx + yy*yy);
          float cosphi = xx/rr;
          
          //look at the segments
          for ( std::vector<reco::MuonSegmentMatch>::const_iterator segment = chamber->segmentMatches.begin();
                segment != chamber->segmentMatches.end(); ++segment ){
            
            if (printLevel > 0) cout << "Segment:" << iSegment+1 << endl;

            int   segIsArb=(segment->Arbitrated)>>8;

            if (!segIsArb) continue;

            float segX=segment->x;
            float segY=segment->y;
            
            float segR = sqrt(segX*segX + segY*segY);
            float segCosPhi = segX/segR;
           
            if (printLevel>0) {
              cout << "###### IS CSC ########"             << endl;
              cout << "trkSegChamberId:"     << chamberId     << endl;
              cout << "trkSegRing:"          << ring          << endl;
              cout << "trkSegStation:"       << station       << endl;
              cout << "trkSegEndcap:"        << endcap        << endl;
              cout << "trkSegTriggerSector:" << triggerSector << endl;
              cout << "trkSegTriggerCscId :" << triggerCscId  << endl;
            
              cout<< "trkSegXfromMatch:"     << xx   << endl;
              cout<< "trkSegYfromMatch:"     << yy    << endl;
              
              if (yy>=0) cout << "trkSegPhifromMatch:" <<      acos(cosphi) << endl;
              else       cout << "trkSegPhifromMatch:" << 2*PI-acos(cosphi) << endl;
              
              cout << "SEGMENT:"                  << endl;
              cout << "trkSegIsArb: " << segIsArb << endl; 
              cout << "trkSegX: "     << segX     << endl; 
              cout << "trkSegY: "     << segY     << endl; 
              
              cout << "segCosPhi:"       << segCosPhi       << endl;
              cout << "acos(segCosPhi):" << acos(segCosPhi) << endl;
              
              if (segY>=0) cout << "trkSegPhi::" <<      acos(segCosPhi) << endl;
              else         cout << "trkSegPhi::" << 2*PI-acos(segCosPhi) << endl;        
            }
            
            
            // fill the matrices
            if (whichMuon < MAX_MUONS && iSegment<MAX_TRK_SEGS) {
              trkSegChamberId    [whichMuon][iSegment] = chamberId;
              trkSegRing         [whichMuon][iSegment] = ring;
              trkSegStation      [whichMuon][iSegment] = station;
              trkSegEndcap       [whichMuon][iSegment] = endcap;
              trkSegTriggerSector[whichMuon][iSegment] = triggerSector;
              trkSegTriggerCscId [whichMuon][iSegment] = triggerCscId;
            
              trkSegXfromMatch[whichMuon][iSegment] = xx;
              trkSegYfromMatch[whichMuon][iSegment] = yy;
              
              if (yy>=0) trkSegPhifromMatch[whichMuon][iSegment] =      acos(cosphi);
              else       trkSegPhifromMatch[whichMuon][iSegment] = 2*PI-acos(cosphi);
              
              trkSegIsArb[whichMuon][iSegment]=segIsArb;
              trkSegX[whichMuon][iSegment]=segX;
              trkSegY[whichMuon][iSegment]=segY;
            
              if (segY>0)trkSegPhi[whichMuon][iSegment]=     acos(segCosPhi);
              else       trkSegPhi[whichMuon][iSegment]=2*PI-acos(segCosPhi);
            }         
            else
              cout << "ERROR: too many segment or muons "
                   << "MAX_MUONS is currently set to "    << MAX_MUONS    << endl
                   << "whichMuon loop is currently at "   << whichMuon    << endl
                   << "MAX_TRK_SEGS is currently set to " << MAX_TRK_SEGS << endl
                   << "iSegment loop is currently at "    << iSegment     << endl;

            iSegment++;        
          }//end of segment loop 
          iChamber++;        
        }//end loop on the chambers
        trkNSegs->push_back(iSegment);
      }// is muon good?
      else trkNSegs->push_back(-999);
    }//isTrackerMuon
    else trkNSegs->push_back(-999);
  }
  
  
}


void TrigEff::muonsInit()
{
  
  isGlobalMuon         = new vector<int>;	  
  isTrackerMuon	       = new vector<int>;
  isStandAloneMuon     = new vector<int>; 
  isMuonAllArbitrated  = new vector<int>; 

  isEnergyValid  = new vector<int>; 
  caloCompatibility  = new vector<float>; 
  em     = new vector<float>;
  emS9   = new vector<float>;
  emS25  = new vector<float>;
  emMax  = new vector<float>;
  had    = new vector<float>;
  hadS9  = new vector<float>;
  hadMax = new vector<float>;
  ho     = new vector<float>;
  hoS9   = new vector<float>;

  gmrPt	               = new vector<float>;
  gmrEta	       = new vector<float>;
  gmrPhi	       = new vector<float>;
  gmrP	               = new vector<float>;
  gmrPx	               = new vector<float>;
  gmrPy	               = new vector<float>;
  gmrPz	               = new vector<float>;
  gmrTheta	       = new vector<float>;
  gmrVx	               = new vector<float>;
  gmrVy	               = new vector<float>;
  gmrVz	               = new vector<float>;
  gmrCharge            = new vector<float>;
  gmrNDoF	       = new vector<float>;
  gmrChi2	       = new vector<float>;
  gmrChi2Norm          = new vector<float>;
  gmrDXY	       = new vector<float>;
  gmrDTheta            = new vector<float>;
  gmrDPt	       = new vector<float>;
  gmrDEta	       = new vector<float>;
  gmrDPhi	       = new vector<float>;
  gmrDDXY	       = new vector<float>;
  gmrIso03nTracks      = new vector<float>;
  gmrIso03sumPt        = new vector<float>;
  gmrDz	               = new vector<float>;
  gmrD0	               = new vector<float>;
  gmrDsz	       = new vector<float>;
  gmrDDz	       = new vector<float>;
  gmrDD0	       = new vector<float>;
  gmrDDsz	       = new vector<float>;
  gmrInnerX            = new vector<float>;
  gmrInnerY            = new vector<float>;
  gmrInnerZ            = new vector<float>;
  gmrOuterX            = new vector<float>;
  gmrOuterY            = new vector<float>;
  gmrOuterZ            = new vector<float>;
  gmrValHits           = new vector<int>;

  stdEnergy    = new vector<float>;
  stdPt        = new vector<float>;
  stdEta       = new vector<float>;
  stdPhi       = new vector<float>;
  stdPx        = new vector<float>;
  stdPy        = new vector<float>;
  stdPz        = new vector<float>;
  stdVx        = new vector<float>;
  stdVy        = new vector<float>;
  stdVz        = new vector<float>;
  stdCharge    = new vector<float>;
  stdDPt       = new vector<float>;
  stdDEta      = new vector<float>;
  stdDPhi      = new vector<float>;
  stdDz        = new vector<float>;
  stdD0        = new vector<float>;
  stdNDoF      = new vector<float>;  
  stdChi2      = new vector<float>;
  stdChi2Norm  = new vector<float>;
  stdDXY       = new vector<float>;
  stdTheta     = new vector<float>;
  stdDTheta    = new vector<float>;
  stdDDz       = new vector<float>;
  stdDD0       = new vector<float>;     
  stdValHits   = new vector<int>;
  stdInnerX    = new vector<float>;
  stdInnerY    = new vector<float>;
  stdInnerZ    = new vector<float>;
  stdOuterX    = new vector<float>;
  stdOuterY    = new vector<float>;
  stdOuterZ    = new vector<float>;
   
  trkEnergy            = new vector<float>;
  trkPt                = new vector<float>;
  trkEta               = new vector<float>;
  trkPhi               = new vector<float>;
  trkPx                = new vector<float>;
  trkPy                = new vector<float>;
  trkPz                = new vector<float>;
  trkVx                = new vector<float>;
  trkVy                = new vector<float>;
  trkVz                = new vector<float>;
  trkCharge            = new vector<float>;
  trkDPt               = new vector<float>;
  trkDEta              = new vector<float>;
  trkDPhi              = new vector<float>;
  trkDz                = new vector<float>;
  trkD0                = new vector<float>;
  trkNDoF              = new vector<float>; 
  trkChi2              = new vector<float>;
  trkChi2Norm          = new vector<float>;
  trkDXY               = new vector<float>;
  trkTheta             = new vector<float>;
  trkDTheta            = new vector<float>;
  trkDDz               = new vector<float>;
  trkDD0               = new vector<float>;     
  trkValHits           = new vector<int>;

  // ------------------------------------------------------
  //segment
  trkNchambers        = new vector<int>;
  trkNofMatches       = new vector<int>;
  trkIsMatchValid     = new vector<int>;

  trkNSegs            = new vector<int>;
  trkSegChamberId.ResizeTo(MAX_MUONS,MAX_TRK_SEGS);
  trkSegRing.ResizeTo(MAX_MUONS,MAX_TRK_SEGS);
  trkSegStation.ResizeTo(MAX_MUONS,MAX_TRK_SEGS);
  trkSegEndcap.ResizeTo(MAX_MUONS,MAX_TRK_SEGS);
  trkSegTriggerSector.ResizeTo(MAX_MUONS,MAX_TRK_SEGS);
  trkSegTriggerCscId.ResizeTo(MAX_MUONS,MAX_TRK_SEGS);

  trkSegXfromMatch.ResizeTo(MAX_MUONS,MAX_TRK_SEGS);
  trkSegYfromMatch.ResizeTo(MAX_MUONS,MAX_TRK_SEGS);
  trkSegPhifromMatch.ResizeTo(MAX_MUONS,MAX_TRK_SEGS);
 
  trkSegIsArb.ResizeTo(MAX_MUONS,MAX_TRK_SEGS);
  trkSegX.ResizeTo(MAX_MUONS,MAX_TRK_SEGS);
  trkSegY.ResizeTo(MAX_MUONS,MAX_TRK_SEGS);
  trkSegR.ResizeTo(MAX_MUONS,MAX_TRK_SEGS);
  trkSegPhi.ResizeTo(MAX_MUONS,MAX_TRK_SEGS);


  for (int row=0; row < MAX_MUONS; row++) 
    for (int col=0; col < MAX_TRK_SEGS; col++) {
      trkSegChamberId[row][col] = -999;   
      trkSegRing[row][col] = -999;        
      trkSegStation[row][col] = -999;     
      trkSegEndcap[row][col] = -999;      
      trkSegTriggerSector[row][col] = -999;
      trkSegTriggerCscId[row][col] = -999;
         
      trkSegXfromMatch[row][col] = -999;  
      trkSegYfromMatch[row][col] = -999;  
      trkSegPhifromMatch[row][col] = -999;

      trkSegIsArb[row][col] = -999;
      trkSegX[row][col]     = -999;
      trkSegY[row][col]     = -999;
      trkSegR[row][col]     = -999;
      trkSegPhi[row][col]   = -999;
    }
  // ------------------------------------------------------  


  rchCSCtype  = new vector<int>  ; 
  rchEta      = new vector<float>;     
  rchPhi      = new vector<float>;     
  rchPhi_02PI = new vector<float>;
  rchStation  = new vector<int>;
  rchChamber  = new vector<int>;
  rchRing     = new vector<int>;
  rchLayer    = new vector<int>;

  rchMuonSize = new vector<int>;
  rchEtaMatrix.ResizeTo(MAX_MUONS,MAX_CSC_RECHIT);
  rchPhiMatrix.ResizeTo(MAX_MUONS,MAX_CSC_RECHIT);
  rchPhi02PIMatrix.ResizeTo(MAX_MUONS,MAX_CSC_RECHIT);
  rchStationMatrix.ResizeTo(MAX_MUONS,MAX_CSC_RECHIT);
  rchChamberMatrix.ResizeTo(MAX_MUONS,MAX_CSC_RECHIT);
  rchRingMatrix.ResizeTo(MAX_MUONS,MAX_CSC_RECHIT);
  rchLayerMatrix.ResizeTo(MAX_MUONS,MAX_CSC_RECHIT);
  rchTypeMatrix.ResizeTo(MAX_MUONS,MAX_CSC_RECHIT);

  for (int row=0; row < MAX_MUONS; row++) 
    for (int col=0; col < MAX_CSC_RECHIT; col++) {
      rchEtaMatrix[row][col]     = -999;
      rchPhiMatrix[row][col]     = -999;
      rchPhi02PIMatrix[row][col] = -999;
      rchStationMatrix[row][col] = -999;
      rchChamberMatrix[row][col] = -999;
      rchRingMatrix[row][col]    = -999;
      rchLayerMatrix[row][col]   = -999;
      rchTypeMatrix[row][col]    = -999;
    }
  

  // propagation to ME1/1
    muons_x_mep11 = new vector<float>;
    muons_y_mep11 = new vector<float>;
    muons_z_mep11 = new vector<float>;
  muons_phi_mep11 = new vector<float>;
  muons_eta_mep11 = new vector<float>;
                
    muons_x_mem11 = new vector<float>;
    muons_y_mem11 = new vector<float>;
    muons_z_mem11 = new vector<float>;
  muons_phi_mem11 = new vector<float>;
  muons_eta_mem11 = new vector<float>;

  // propagation to ME1
    muons_x_mep1 = new vector<float>;
    muons_y_mep1 = new vector<float>;
    muons_z_mep1 = new vector<float>;
  muons_phi_mep1 = new vector<float>;
  muons_eta_mep1 = new vector<float>;
                 
    muons_x_mem1 = new vector<float>;
    muons_y_mem1 = new vector<float>;
    muons_z_mem1 = new vector<float>;
  muons_phi_mem1 = new vector<float>;
  muons_eta_mem1 = new vector<float>;

  // propagation to ME2
    muons_x_mep2 = new vector<float>;
    muons_y_mep2 = new vector<float>;
    muons_z_mep2 = new vector<float>;
  muons_phi_mep2 = new vector<float>;
  muons_eta_mep2 = new vector<float>;
                  
    muons_x_mem2 = new vector<float>;
    muons_y_mem2 = new vector<float>;
    muons_z_mem2 = new vector<float>;
  muons_phi_mem2 = new vector<float>;
  muons_eta_mem2 = new vector<float>;

  // propagation to ME3                 
    muons_x_mep3 = new vector<float>;
    muons_y_mep3 = new vector<float>;
    muons_z_mep3 = new vector<float>;
  muons_phi_mep3 = new vector<float>;
  muons_eta_mep3 = new vector<float>;
                 
    muons_x_mem3 = new vector<float>;
    muons_y_mem3 = new vector<float>;
    muons_z_mem3 = new vector<float>;
  muons_phi_mem3 = new vector<float>;
  muons_eta_mem3 = new vector<float>;
  
}


void TrigEff::muonsDel() {
  
  vector<int>().swap(*isGlobalMuon);	  
  vector<int>().swap(*isTrackerMuon);
  vector<int>().swap(*isStandAloneMuon); 
  vector<int>().swap(*isMuonAllArbitrated); 

  vector<int>().swap(*isEnergyValid); 
  vector<float>().swap(*caloCompatibility); 
  vector<float>().swap(*em); 
  vector<float>().swap(*emS9);
  vector<float>().swap(*emS25);
  vector<float>().swap(*emMax);
  vector<float>().swap(*had);
  vector<float>().swap(*hadS9);
  vector<float>().swap(*hadMax);
  vector<float>().swap(*ho);
  vector<float>().swap(*hoS9);

  vector<float>().swap(*gmrPt);
  vector<float>().swap(*gmrEta);
  vector<float>().swap(*gmrPhi);
  vector<float>().swap(*gmrP);
  vector<float>().swap(*gmrPx);
  vector<float>().swap(*gmrPy);
  vector<float>().swap(*gmrPz);
  vector<float>().swap(*gmrTheta);
  vector<float>().swap(*gmrVx);
  vector<float>().swap(*gmrVy);
  vector<float>().swap(*gmrVz);
  vector<float>().swap(*gmrCharge);
  vector<float>().swap(*gmrNDoF);
  vector<float>().swap(*gmrChi2);
  vector<float>().swap(*gmrChi2Norm);
  vector<float>().swap(*gmrDXY);
  vector<float>().swap(*gmrDTheta);
  vector<float>().swap(*gmrDPt);
  vector<float>().swap(*gmrDEta);
  vector<float>().swap(*gmrDPhi);
  vector<float>().swap(*gmrDDXY);
  vector<float>().swap(*gmrIso03nTracks);
  vector<float>().swap(*gmrIso03sumPt);
  vector<float>().swap(*gmrChi2Norm);
  vector<float>().swap(*gmrDz);
  vector<float>().swap(*gmrD0);
  vector<float>().swap(*gmrDsz);
  vector<float>().swap(*gmrDDz);
  vector<float>().swap(*gmrDD0);
  vector<float>().swap(*gmrDDsz);
  vector<float>().swap(*gmrInnerX);
  vector<float>().swap(*gmrInnerY);
  vector<float>().swap(*gmrInnerZ);
  vector<float>().swap(*gmrOuterX);
  vector<float>().swap(*gmrOuterY);
  vector<float>().swap(*gmrOuterZ);
  vector<int>().swap(*gmrValHits);


  vector<float>().swap(*stdEnergy);
  vector<float>().swap(*stdPt);
  vector<float>().swap(*stdEta);
  vector<float>().swap(*stdPhi);
  vector<float>().swap(*stdPx);
  vector<float>().swap(*stdPy);
  vector<float>().swap(*stdPz);
  vector<float>().swap(*stdVx);
  vector<float>().swap(*stdVy);
  vector<float>().swap(*stdVz);
  vector<float>().swap(*stdCharge);
  vector<float>().swap(*stdDPt);
  vector<float>().swap(*stdDEta);
  vector<float>().swap(*stdDPhi);
  vector<float>().swap(*stdDz);
  vector<float>().swap(*stdD0);
  vector<float>().swap(*stdNDoF); 
  vector<float>().swap(*stdChi2);
  vector<float>().swap(*stdChi2Norm);
  vector<float>().swap(*stdDXY);
  vector<float>().swap(*stdTheta);
  vector<float>().swap(*stdDTheta);
  vector<float>().swap(*stdDDz);
  vector<float>().swap(*stdDD0);     
  vector<int>().swap(*stdValHits);
  vector<float>().swap(*stdInnerX);
  vector<float>().swap(*stdInnerY);
  vector<float>().swap(*stdInnerZ);
  vector<float>().swap(*stdOuterX);
  vector<float>().swap(*stdOuterY);
  vector<float>().swap(*stdOuterZ);
   

  vector<float>().swap(*trkEnergy);
  vector<float>().swap(*trkPt);
  vector<float>().swap(*trkEta);
  vector<float>().swap(*trkPhi);
  vector<float>().swap(*trkPx);
  vector<float>().swap(*trkPy);
  vector<float>().swap(*trkPz);
  vector<float>().swap(*trkVx);
  vector<float>().swap(*trkVy);
  vector<float>().swap(*trkVz);
  vector<float>().swap(*trkCharge);
  vector<float>().swap(*trkDPt);
  vector<float>().swap(*trkDEta);
  vector<float>().swap(*trkDPhi);
  vector<float>().swap(*trkDz);
  vector<float>().swap(*trkD0);
  vector<float>().swap(*trkNDoF); 
  vector<float>().swap(*trkChi2);
  vector<float>().swap(*trkChi2Norm );
  vector<float>().swap(*trkDXY);
  vector<float>().swap(*trkTheta);
  vector<float>().swap(*trkDTheta);
  vector<float>().swap(*trkDDz);
  vector<float>().swap(*trkDD0);     
  vector<int>().swap(*trkValHits);

  // CSC segment for the tracker muon
  //chamber
  vector<int>().swap(*trkNchambers); 
  vector<int>().swap(*trkNofMatches);
  vector<int>().swap(*trkIsMatchValid);

  vector<int>().swap(*trkNSegs);

  vector<int>  ().swap(*rchCSCtype ); 
  vector<float>().swap(*rchEta     );     
  vector<float>().swap(*rchPhi     );     
  vector<float>().swap(*rchPhi_02PI);

  vector<int>().swap(*rchStation);
  vector<int>().swap(*rchChamber);
  vector<int>().swap(*rchRing);
  vector<int>().swap(*rchLayer);
  
  vector<int>().swap(*rchMuonSize);

  delete [] rchEtaList;
  delete [] rchPhiList;
  delete [] rchPhiList_02PI;


    vector<float>().swap(*muons_x_mep11);
    vector<float>().swap(*muons_y_mep11);
    vector<float>().swap(*muons_z_mep11);
  vector<float>().swap(*muons_phi_mep11);
  vector<float>().swap(*muons_eta_mep11);
                                     
    vector<float>().swap(*muons_x_mem11);
    vector<float>().swap(*muons_y_mem11);
    vector<float>().swap(*muons_z_mem11);
  vector<float>().swap(*muons_phi_mem11);
  vector<float>().swap(*muons_eta_mem11);

    vector<float>().swap(*muons_x_mep1);
    vector<float>().swap(*muons_y_mep1);
    vector<float>().swap(*muons_z_mep1);
  vector<float>().swap(*muons_phi_mep1);
  vector<float>().swap(*muons_eta_mep1);
                                     
    vector<float>().swap(*muons_x_mem1);
    vector<float>().swap(*muons_y_mem1);
    vector<float>().swap(*muons_z_mem1);
  vector<float>().swap(*muons_phi_mem1);
  vector<float>().swap(*muons_eta_mem1);

    vector<float>().swap(*muons_x_mep2);
    vector<float>().swap(*muons_y_mep2);
    vector<float>().swap(*muons_z_mep2);
  vector<float>().swap(*muons_phi_mep2);
  vector<float>().swap(*muons_eta_mep2);
                       
    vector<float>().swap(*muons_x_mem2);
    vector<float>().swap(*muons_y_mem2);
    vector<float>().swap(*muons_z_mem2);
  vector<float>().swap(*muons_phi_mem2);
  vector<float>().swap(*muons_eta_mem2);
                       
    vector<float>().swap(*muons_x_mep3);
    vector<float>().swap(*muons_y_mep3);
    vector<float>().swap(*muons_z_mep3);
  vector<float>().swap(*muons_phi_mep3);
  vector<float>().swap(*muons_eta_mep3);
                       
    vector<float>().swap(*muons_x_mem3);
    vector<float>().swap(*muons_y_mem3);
    vector<float>().swap(*muons_z_mem3);
  vector<float>().swap(*muons_phi_mem3);
  vector<float>().swap(*muons_eta_mem3);

}


//-------------------------------------------------------------------------
// fill the l1 extra collection
//-------------------------------------------------------------------------
void TrigEff::fillExtra(l1extra::L1MuonParticleCollection::const_iterator l1muon) {
  if (printLevel > 0) {
    // placeholder for some printout 
  }
  
  l1Eta->push_back(l1muon->eta());
  l1Pt ->push_back(l1muon->pt());
  l1Phi->push_back(l1muon->phi());
          
  isIsolated->push_back(l1muon->isIsolated());
  isMip     ->push_back(l1muon->isMip()     );
  isForward ->push_back(l1muon->isForward() );
  isRPC     ->push_back(l1muon->isRPC()     );
  detectorType ->push_back(l1muon->gmtMuonCand().detector());
  rank ->push_back(l1muon->gmtMuonCand().rank());
}
  

void TrigEff::l1extraInit() {
  
  l1Eta = new vector<float>; 
  l1Pt  = new vector<float>; 
  l1Phi = new vector<float>; 
                               
  isIsolated = new vector<int>;   
  isMip      = new vector<int>;   
  isForward  = new vector<int>;   
  isRPC      = new vector<int>;   

  detectorType = new vector<int>;  
  rank         = new vector<int>;
}

void TrigEff::l1extraDel() {
  
  vector<float>().swap(*l1Eta); 
  vector<float>().swap(*l1Pt ); 
  vector<float>().swap(*l1Phi); 
                               
  vector<int>().swap(*isIsolated);   
  vector<int>().swap(*isMip     );   
  vector<int>().swap(*isForward );   
  vector<int>().swap(*isRPC     );   

  vector<int>().swap(*detectorType);  
  vector<int>().swap(*rank);
}


void TrigEff::resizeRchHits(int nMu_nMaxCscRchHits){
  
  if (!rchEtaList)      delete [] rchEtaList;
  if (!rchPhiList)      delete [] rchPhiList;
  if (!rchPhiList_02PI) delete [] rchPhiList_02PI;

  rchEtaList      = new Double_t[nMu_nMaxCscRchHits];
  rchPhiList      = new Double_t[nMu_nMaxCscRchHits];
  rchPhiList_02PI = new Double_t[nMu_nMaxCscRchHits];
  
  recoMuons->SetBranchAddress("rchEtaList"     , rchEtaList); 
  recoMuons->SetBranchAddress("rchPhiList"     , rchPhiList);
  recoMuons->SetBranchAddress("rchPhiList_02PI", rchPhiList_02PI);
  
  for (int i = 0; i < nMu_nMaxCscRchHits; i++) {
    rchEtaList[i]      = -999;
    rchPhiList[i]      = -999;
    rchPhiList_02PI[i] = -999;
  }
}

// void TrigEff::csctfInit() {	
    
//   EndcapTrk     = new vector<int>;
//   SectorTrk     = new vector<int>;  
//   BxTrk         = new vector<int>;  
       	    
//   me1ID         = new vector<int>;  
//   me2ID         = new vector<int>;
//   me3ID         = new vector<int>;  
//   me4ID         = new vector<int>;  
//   mb1ID         = new vector<int>;  

//   ModeTrk       = new vector<int>  ;
//   EtaTrk        = new vector<float>;  
//   PhiTrk        = new vector<float>;  
//   PhiTrk_02PI   = new vector<float>;  
//   PtTrk         = new vector<float>;  
    
//   ChargeValidTrk = new vector<int>;
//   ChargeTrk      = new vector<int>;
//   QualityTrk     = new vector<int>;
//   ForRTrk        = new vector<int>;
//   Phi23Trk       = new vector<int>;
//   Phi12Trk       = new vector<int>;
//   PhiSignTrk     = new vector<int>;

//   EtaBit        = new vector<int>;
//   PhiBit        = new vector<int>;
//   LocalPhi      = new vector<int>;
//   PtBit         = new vector<int>;

//   PtTrkData     = new vector<float>;
//   QualityTrkData= new vector<int>;   
//   OutputLinkTrk = new vector<int>;      

// }

// void TrigEff::csctfDel() {	

//   vector<int>().swap(*EndcapTrk); 
//   vector<int>().swap(*SectorTrk); 
//   vector<int>().swap(*BxTrk    ); 
  
//   vector<int>().swap(*me1ID); 
//   vector<int>().swap(*me2ID); 
//   vector<int>().swap(*me3ID); 
//   vector<int>().swap(*me4ID); 
//   vector<int>().swap(*mb1ID); 

//   vector<int>  ().swap(*ModeTrk    ); 
//   vector<float>().swap(*EtaTrk     ); 
//   vector<float>().swap(*PhiTrk     ); 
//   vector<float>().swap(*PhiTrk_02PI); 
//   vector<float>().swap(*PtTrk      ); 

//   vector<int>().swap(*ChargeTrk      ); 
//   vector<int>().swap(*ChargeValidTrk ); 
//   vector<int>().swap(*QualityTrk     ); 
//   vector<int>().swap(*ForRTrk        ); 
//   vector<int>().swap(*Phi23Trk       ); 
//   vector<int>().swap(*Phi12Trk       ); 
//   vector<int>().swap(*PhiSignTrk     ); 

//   vector<int>().swap(*EtaBit);
//   vector<int>().swap(*PhiBit);
//   vector<int>().swap(*LocalPhi);
//   vector<int>().swap(*PtBit );
  
//   vector<float>().swap(*PtTrkData);  
//   vector<int>().swap(*QualityTrkData);
//   vector<int>().swap(*OutputLinkTrk);
// }


// this code snippet was taken from 
// http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/L1TriggerDPG/src/L1MuonRecoTreeProducer.cc?revision=1.6&view=markup
// to get track position at a particular (xy) plane given its z
TrajectoryStateOnSurface TrigEff::surfExtrapTrkSam(reco::TrackRef track, double z)
{
  Plane::PositionType pos(0, 0, z);
  Plane::RotationType rot;
  Plane::PlanePointer myPlane = Plane::build(pos, rot);

  FreeTrajectoryState recoStart = freeTrajStateMuon(track);
  TrajectoryStateOnSurface recoProp;
  recoProp = propagatorAlong->propagate(recoStart, *myPlane);
  if (!recoProp.isValid()) {
    recoProp = propagatorOpposite->propagate(recoStart, *myPlane);
  }
  return recoProp;
}

FreeTrajectoryState TrigEff::freeTrajStateMuon(reco::TrackRef track)
{
  GlobalPoint  innerPoint(track->innerPosition().x(),  track->innerPosition().y(),  track->innerPosition().z());
  GlobalVector innerVec  (track->innerMomentum().x(),  track->innerMomentum().y(),  track->innerMomentum().z());  
  
  FreeTrajectoryState recoStart(innerPoint, innerVec, track->charge(), &*theBField);
  
  return recoStart;
}


