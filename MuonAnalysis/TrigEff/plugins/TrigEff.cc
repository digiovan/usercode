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
  L1extraTag= pset.getParameter<InputTag>("L1extraTag");
  muonsTag  = pset.getParameter<InputTag>("muonsTag");
  outputFile= pset.getParameter<string>("outputFile");

  csctfTag      = pset.getParameter<InputTag>("csctfTag");
  csctfLctsTag  = pset.getParameter<InputTag>("csctfLctsTag");

  cscSegTag = pset.getParameter<InputTag>("cscSegTag");

  // L1GTUtils
  m_nameAlgTechTrig = pset.getParameter<std::string> ("AlgorithmName");
    
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
   
  recoMuons->Branch("isGlobalMuon"           , &isGlobalMuon       );
  recoMuons->Branch("isTrackerMuon"          , &isTrackerMuon      );
  recoMuons->Branch("isStandAloneMuon"       , &isStandAloneMuon   );
  recoMuons->Branch("isMuonAllArbitrated"    , &isMuonAllArbitrated);
  recoMuons->Branch("isTMLastStationAngTight", &isTMLastStationAngTight);
  recoMuons->Branch("isGlobalMuonPromptTight", &isGlobalMuonPromptTight);

  recoMuons->Branch("isEnergyValid"    , &isEnergyValid);
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
  //recoMuons->Branch("gmrEnergy"           , &gmrEnergy           );
  //recoMuons->Branch("gmrDEnergy"          , &gmrDEnergy          );
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
  //recoMuons->Branch("stdEnergy"   , &stdEnergy   );
  //recoMuons->Branch("stdDEnergy"  , &stdDEnergy  );
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
  //recoMuons->Branch("trkEnergy"   , &trkEnergy   );
  //recoMuons->Branch("trkDEnergy"  , &trkDEnergy  );
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
  
  recoMuons->Branch("trkNSegs"           , &trkNSegs);

  recoMuons->Branch("trkSegChamberId"    , trkSegChamberId    ,"trkSegChamberId[muonSize][100]/I");
  recoMuons->Branch("trkSegRing"         , trkSegRing         ,"trkSegRing[muonSize][100]/I");
  recoMuons->Branch("trkSegStation"      , trkSegStation      ,"trkSegStation[muonSize][100]/I");
  recoMuons->Branch("trkSegEndcap"       , trkSegEndcap       ,"trkSegEndcap[muonSize][100]/I");
  recoMuons->Branch("trkSegTriggerSector", trkSegTriggerSector,"trkSegTriggerSector[muonSize][100]/I");
  recoMuons->Branch("trkSegTriggerCscId" , trkSegTriggerCscId ,"trkSegTriggerCscId[muonSize][100]/I");
  recoMuons->Branch("trkSegXfromMatch"   , trkSegXfromMatch   ,"trkSegXfromMatch[muonSize][100]/F");
  recoMuons->Branch("trkSegYfromMatch"   , trkSegYfromMatch   ,"trkSegYfromMatch[muonSize][100]/F");
  recoMuons->Branch("trkSegZfromMatch"   , trkSegZfromMatch   ,"trkSegZfromMatch[muonSize][100]/F");
  recoMuons->Branch("trkSegRfromMatch"   , trkSegRfromMatch   ,"trkSegRfromMatch[muonSize][100]/F");
  recoMuons->Branch("trkSegPhifromMatch" , trkSegPhifromMatch ,"trkSegPhifromMatch[muonSize][100]/F");
  recoMuons->Branch("trkSegEtafromMatch" , trkSegEtafromMatch ,"trkSegEtafromMatch[muonSize][100]/F");

  //segment
  recoMuons->Branch("trkSegIsArb", &trkSegIsArb, "trkSegIsArb[muonSize][100]/I");
  recoMuons->Branch("trkSegX"    , &trkSegX    , "trkSegX[muonSize][100]/F");
  recoMuons->Branch("trkSegY"    , &trkSegY    , "trkSegY[muonSize][100]/F");
  recoMuons->Branch("trkSegZ"    , &trkSegZ    , "trkSegZ[muonSize][100]/F");
  recoMuons->Branch("trkSegR"    , &trkSegR    , "trkSegR[muonSize][100]/F");
  recoMuons->Branch("trkSegPhi"  , &trkSegPhi  , "trkSegPhi[muonSize][100]/F");
  recoMuons->Branch("trkSegEta"  , &trkSegEta  , "trkSegEta[muonSize][100]/F");
  recoMuons->Branch("trkSegDxDz"    , &trkSegDxDz    , "trkSegDxDz[muonSize][100]/F");
  recoMuons->Branch("trkSegDyDz"    , &trkSegDyDz    , "trkSegDyDz[muonSize][100]/F");
  recoMuons->Branch("trkSegDxDzErr" , &trkSegDxDzErr , "trkSegDxDzErr[muonSize][100]/F");
  recoMuons->Branch("trkSegDyDzErr" , &trkSegDyDzErr , "trkSegDyDzErr[muonSize][100]/F");

  //---------------------------------------------------------------------
  // RECHIT information: only for standalone/global muons!
  //---------------------------------------------------------------------
  recoMuons->Branch("rchCSCtype" , &rchCSCtype ); 
  recoMuons->Branch("rchEtaLocal", &rchEtaLocal); 
  recoMuons->Branch("rchPhiLocal", &rchPhiLocal); 
  recoMuons->Branch("rchEta"     , &rchEta     ); 
  recoMuons->Branch("rchPhi"     , &rchPhi     ); 
  recoMuons->Branch("rchPhi_02PI", &rchPhi_02PI);

  recoMuons->Branch("rchStation", &rchStation);
  recoMuons->Branch("rchChamber", &rchChamber);
  recoMuons->Branch("rchRing"   , &rchRing   );
  recoMuons->Branch("rchLayer"  , &rchLayer  );

  recoMuons->Branch("rchMuonSize",      &rchMuonSize     );
 
  recoMuons->Branch("rchEtaMatrixLocal", rchEtaMatrixLocal, "rchEtaMatrixLocal[muonSize][35]/F");
  recoMuons->Branch("rchPhiMatrixLocal", rchPhiMatrixLocal, "rchPhiMatrixLocal[muonSize][35]/F");
  recoMuons->Branch("rchEtaMatrix"    , rchEtaMatrix    , "rchEtaMatrix[muonSize][35]/F");
  recoMuons->Branch("rchPhiMatrix"    , rchPhiMatrix    , "rchPhiMatrix[muonSize][35]/F");
  recoMuons->Branch("rchPhi02PIMatrix", rchPhi02PIMatrix, "rchPhi02PIMatrix[muonSize][35]/F");
  recoMuons->Branch("rchStationMatrix", rchStationMatrix, "rchStationMatrix[muonSize][35]/I");
  recoMuons->Branch("rchChamberMatrix", rchChamberMatrix, "rchChamberMatrix[muonSize][35]/I");
  recoMuons->Branch("rchRingMatrix"   , rchRingMatrix   , "rchRingMatrix[muonSize][35]/I");
  recoMuons->Branch("rchLayerMatrix"  , rchLayerMatrix  , "rchLayerMatrix[muonSize][35]/I");
  recoMuons->Branch("rchTypeMatrix"   , rchTypeMatrix   , "rchTypeMatrix[muonSize][35]/I");
   
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
  recoMuons->Branch(  "muons_x_me11",  &muons_x_me11);
  recoMuons->Branch(  "muons_y_me11",  &muons_y_me11);
  recoMuons->Branch(  "muons_z_me11",  &muons_z_me11);
  recoMuons->Branch("muons_phi_me11",&muons_phi_me11);
  recoMuons->Branch("muons_eta_me11",&muons_eta_me11);

  // propagation to ME1
  recoMuons->Branch(  "muons_x_me1",  &muons_x_me1);
  recoMuons->Branch(  "muons_y_me1",  &muons_y_me1);
  recoMuons->Branch(  "muons_z_me1",  &muons_z_me1);
  recoMuons->Branch("muons_phi_me1",&muons_phi_me1);
  recoMuons->Branch("muons_eta_me1",&muons_eta_me1);

  // propagation to ME2
  recoMuons->Branch(  "muons_x_me2",  &muons_x_me2);
  recoMuons->Branch(  "muons_y_me2",  &muons_y_me2);
  recoMuons->Branch(  "muons_z_me2",  &muons_z_me2);
  recoMuons->Branch("muons_phi_me2",&muons_phi_me2);
  recoMuons->Branch("muons_eta_me2",&muons_eta_me2);

  // propagation to ME3
  recoMuons->Branch(  "muons_x_me3",  &muons_x_me3);
  recoMuons->Branch(  "muons_y_me3",  &muons_y_me3);
  recoMuons->Branch(  "muons_z_me3",  &muons_z_me3);
  recoMuons->Branch("muons_phi_me3",&muons_phi_me3);
  recoMuons->Branch("muons_eta_me3",&muons_eta_me3);

  //---------------------------------------------------------------------
  // all segment information
  //---------------------------------------------------------------------
  recoMuons->Branch("segsSize"     ,&segsSize     ,"segsSize/I");
  recoMuons->Branch("cscsegs_loc_x",&cscsegs_loc_x);
  recoMuons->Branch("cscsegs_loc_y",&cscsegs_loc_y);
  recoMuons->Branch("cscsegs_loc_z",&cscsegs_loc_z);

  recoMuons->Branch("cscsegs_loc_theta",&cscsegs_loc_theta);
  recoMuons->Branch("cscsegs_loc_eta",  &cscsegs_loc_eta);
  recoMuons->Branch("cscsegs_loc_phi",  &cscsegs_loc_phi);

  recoMuons->Branch("cscsegs_loc_dir_theta",&cscsegs_loc_dir_theta);
  recoMuons->Branch("cscsegs_loc_dir_eta",  &cscsegs_loc_dir_eta);
  recoMuons->Branch("cscsegs_loc_dir_phi",  &cscsegs_loc_dir_phi);

  recoMuons->Branch("cscsegs_gbl_x",&cscsegs_gbl_x);
  recoMuons->Branch("cscsegs_gbl_y",&cscsegs_gbl_y);
  recoMuons->Branch("cscsegs_gbl_z",&cscsegs_gbl_z);

  recoMuons->Branch("cscsegs_gbl_theta",&cscsegs_gbl_theta);
  recoMuons->Branch("cscsegs_gbl_eta",  &cscsegs_gbl_eta);
  recoMuons->Branch("cscsegs_gbl_phi",  &cscsegs_gbl_phi);

  recoMuons->Branch("cscsegs_gbl_dir_theta",&cscsegs_gbl_dir_theta);
  recoMuons->Branch("cscsegs_gbl_dir_eta",  &cscsegs_gbl_dir_eta);
  recoMuons->Branch("cscsegs_gbl_dir_phi",  &cscsegs_gbl_dir_phi);

  recoMuons->Branch("cscsegs_endcap" ,&cscsegs_endcap );
  recoMuons->Branch("cscsegs_station",&cscsegs_station);
  recoMuons->Branch("cscsegs_ring"   ,&cscsegs_ring   );
  recoMuons->Branch("cscsegs_chamber",&cscsegs_chamber);
  
  //---------------------------------------------------------------------
  // segments belonging to the muon
  //---------------------------------------------------------------------
  recoMuons->Branch("muonNsegs",&muonNsegs);

  recoMuons->Branch("muon_cscsegs_loc_x"      , muon_cscsegs_loc_x      ,"muon_cscsegs_loc_x[muonSize][16]/F");
  recoMuons->Branch("muon_cscsegs_loc_y"      , muon_cscsegs_loc_y      ,"muon_cscsegs_loc_y[muonSize][16]/F");
  recoMuons->Branch("muon_cscsegs_loc_eta"    , muon_cscsegs_loc_eta    ,"muon_cscsegs_loc_eta[muonSize][16]/F");
  recoMuons->Branch("muon_cscsegs_loc_phi"    , muon_cscsegs_loc_phi    ,"muon_cscsegs_loc_phi[muonSize][16]/F");
  recoMuons->Branch("muon_cscsegs_loc_dir_eta", muon_cscsegs_loc_dir_eta,"muon_cscsegs_loc_dir_eta[muonSize][16]/F");
  recoMuons->Branch("muon_cscsegs_loc_dir_phi", muon_cscsegs_loc_dir_phi,"muon_cscsegs_loc_dir_phi[muonSize][16]/F");

  recoMuons->Branch("muon_cscsegs_gbl_x"      , muon_cscsegs_gbl_x      ,"muon_cscsegs_gbl_x[muonSize][16]/F");
  recoMuons->Branch("muon_cscsegs_gbl_y"      , muon_cscsegs_gbl_y      ,"muon_cscsegs_gbl_y[muonSize][16]/F");
  recoMuons->Branch("muon_cscsegs_gbl_eta"    , muon_cscsegs_gbl_eta    ,"muon_cscsegs_gbl_eta[muonSize][16]/F");
  recoMuons->Branch("muon_cscsegs_gbl_phi"    , muon_cscsegs_gbl_phi    ,"muon_cscsegs_gbl_phi[muonSize][16]/F");
  recoMuons->Branch("muon_cscsegs_gbl_dir_eta", muon_cscsegs_gbl_dir_eta,"muon_cscsegs_gbl_dir_eta[muonSize][16]/F");
  recoMuons->Branch("muon_cscsegs_gbl_dir_phi", muon_cscsegs_gbl_dir_phi,"muon_cscsegs_gbl_dir_phi[muonSize][16]/F");

  recoMuons->Branch("muon_cscsegs_dxdz"   , muon_cscsegs_dxdz   ,"muon_cscsegs_dxdz[muonSize][16]/F");
  recoMuons->Branch("muon_cscsegs_dydz"   , muon_cscsegs_dydz   ,"muon_cscsegs_dydz[muonSize][16]/F");
  recoMuons->Branch("muon_cscsegs_dxdzErr", muon_cscsegs_dxdzErr,"muon_cscsegs_dxdzErr[muonSize][16]/F");
  recoMuons->Branch("muon_cscsegs_dydzErr", muon_cscsegs_dydzErr,"muon_cscsegs_dydzErr[muonSize][16]/F");

  recoMuons->Branch("muon_cscsegs_endcap" , muon_cscsegs_endcap ,"muon_cscsegs_endcap[muonSize][16]/I");
  recoMuons->Branch("muon_cscsegs_station", muon_cscsegs_station,"muon_cscsegs_station[muonSize][16]/I");
  recoMuons->Branch("muon_cscsegs_ring"   , muon_cscsegs_ring   ,"muon_cscsegs_ring[muonSize][16]/I");
  recoMuons->Branch("muon_cscsegs_chamber", muon_cscsegs_chamber,"muon_cscsegs_chamber[muonSize][16]/I");
  recoMuons->Branch("muon_cscsegs_nhits"  , muon_cscsegs_nhits  ,"muon_cscsegs_nhits[muonSize][16]/I");
  
  recoMuons->Branch("muon_cscsegs_islctable", muon_cscsegs_islctable,"muon_cscsegs_islctable[muonSize][16]/I");
  recoMuons->Branch("muon_cscsegs_ismatched", muon_cscsegs_ismatched,"muon_cscsegs_ismatched[muonSize][16]/I");
  recoMuons->Branch("muon_cscsegs_lctId"    , muon_cscsegs_lctId    ,"muon_cscsegs_lctId[muonSize][16]/I");
  recoMuons->Branch("muon_cscsegs_nmatched" , muon_cscsegs_nmatched ,"muon_cscsegs_nmatched[muonSize][16]/I");
  

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

  csctfTTree = new TTree("csctfTTree","csctfTTree");
  // csctf
  csctfTTree->Branch("SizeTrk"       , &SizeTrk,      "SizeTrk/I");
  csctfTTree->Branch("EndcapTrk"     , &EndcapTrk     );
  csctfTTree->Branch("SectorTrk"     , &SectorTrk     );
  csctfTTree->Branch("BxTrk"  	     , &BxTrk         );
  csctfTTree->Branch("me1ID"  	     , &me1ID         );
  csctfTTree->Branch("me2ID"         , &me2ID 	      );
  csctfTTree->Branch("me3ID"  	     , &me3ID         );
  csctfTTree->Branch("me4ID"  	     , &me4ID         );
  csctfTTree->Branch("mb1ID"  	     , &mb1ID         );
  csctfTTree->Branch("OutputLinkTrk" , &OutputLinkTrk );
  csctfTTree->Branch("ModeTrk"       , &ModeTrk       );
  csctfTTree->Branch("EtaTrk"  	     , &EtaTrk        );
  csctfTTree->Branch("PhiTrk"  	     , &PhiTrk        );
  csctfTTree->Branch("PhiTrk_02PI"   , &PhiTrk_02PI   );
  csctfTTree->Branch("PtTrk"  	     , &PtTrk         );
  csctfTTree->Branch("ChargeTrk"     , &ChargeTrk     );
  csctfTTree->Branch("ChargeValidTrk", &ChargeValidTrk);
  csctfTTree->Branch("QualityTrk"    , &QualityTrk    );
  csctfTTree->Branch("ForRTrk"       , &ForRTrk       );
  csctfTTree->Branch("Phi23Trk"      , &Phi23Trk      );
  csctfTTree->Branch("Phi12Trk"      , &Phi12Trk      );
  csctfTTree->Branch("PhiSignTrk"    , &PhiSignTrk    );
  csctfTTree->Branch("EtaBitTrk"     , &EtaBitTrk     );
  csctfTTree->Branch("PhiBitTrk"     , &PhiBitTrk     );
  csctfTTree->Branch("PtBitTrk"      , &PtBitTrk      );

  csctfTTree->Branch("NumLCTsTrk"       , &NumLCTsTrk       );

  csctfTTree->Branch("trLctEndcap"      , trLctEndcap      , "trLctEndcap[SizeTrk][4]");
  csctfTTree->Branch("trLctSector"      , trLctSector      , "trLctSector[SizeTrk][4]");
  csctfTTree->Branch("trLctSubSector"   , trLctSubSector   , "trLctSubSector[SizeTrk][4]");
  csctfTTree->Branch("trLctBx"          , trLctBx          , "trLctBx[SizeTrk][4]");
  csctfTTree->Branch("trLctBx0"         , trLctBx0         , "trLctBx0[SizeTrk][4]");
                                                                                     
  csctfTTree->Branch("trLctStation"     , trLctStation     , "trLctStation[SizeTrk][4]");
  csctfTTree->Branch("trLctRing"        , trLctRing        , "trLctRing[SizeTrk][4]");
  csctfTTree->Branch("trLctChamber"     , trLctChamber     , "trLctChamber[SizeTrk][4]");
  csctfTTree->Branch("trLctTriggerCSCID", trLctTriggerCSCID, "trLctTriggerCSCID[SizeTrk][4]");
  csctfTTree->Branch("trLctFpga"        , trLctFpga        , "trLctFpga[SizeTrk][4]");
                                                                                     
  csctfTTree->Branch("trLctlocalPhi"    , trLctlocalPhi    , "trLctlocalPhi[SizeTrk][4]");
  csctfTTree->Branch("trLctglobalPhi"   , trLctglobalPhi   , "trLctglobalPhi[SizeTrk][4]");
  csctfTTree->Branch("trLctglobalEta"   , trLctglobalEta   , "trLctglobalEta[SizeTrk][4]");
                                                                                  
  csctfTTree->Branch("trLctstripNum"    , trLctstripNum    , "trLctstripNum[SizeTrk][4]");
  csctfTTree->Branch("trLctwireGroup"   , trLctwireGroup   , "trLctwireGroup[SizeTrk][4]");

  // all lcts
  csctfTTree->Branch("SizeLCTs"       , &SizeLCTs       ,"SizeLCTs/I");
  csctfTTree->Branch("lctEndcap"      , &lctEndcap      );
  csctfTTree->Branch("lctSector"      , &lctSector      );
  csctfTTree->Branch("lctSubSector"   , &lctSubSector   );
  csctfTTree->Branch("lctBx"          , &lctBx          );
  csctfTTree->Branch("lctBx0"         , &lctBx0         );
  csctfTTree->Branch("lctStation"     , &lctStation     );
  csctfTTree->Branch("lctRing"        , &lctRing        );
  csctfTTree->Branch("lctChamber"     , &lctChamber     );
  csctfTTree->Branch("lctTriggerCSCID", &lctTriggerCSCID);
  csctfTTree->Branch("lctFpga"        , &lctFpga        );
  csctfTTree->Branch("lctlocalPhi"    , &lctlocalPhi    );
  csctfTTree->Branch("lctglobalPhi"   , &lctglobalPhi   );
  csctfTTree->Branch("lctglobalEta"   , &lctglobalEta   );
  csctfTTree->Branch("lctstripNum"    , &lctstripNum    );
  csctfTTree->Branch("lctwireGroup"   , &lctwireGroup   );

  // -----------------------------------------------------------------------------
  // csctf does not preserve information about the LCT (stubs) which forms
  // the track so we need to retrieve this information. In order to do so
  // we need to initialize the Sector Receiver LUTs in the software
  
  bzero(srLUTs_ , sizeof(srLUTs_));
  int sector=1;    // assume SR LUTs are all same for every sector
  bool TMB07=true; // specific TMB firmware
  // Create a pset for SR/PT LUTs: if you do not change the value in the 
  // configuration file, it will load the default minitLUTs
  edm::ParameterSet srLUTset;
  srLUTset.addUntrackedParameter<bool>("ReadLUTs", false);
  srLUTset.addUntrackedParameter<bool>("Binary",   false);
  srLUTset.addUntrackedParameter<std::string>("LUTPath", "./");
 
   // positive endcap
  int endcap = 1; 
  for(int station=1,fpga=0; station<=4 && fpga<5; station++)
    {
      if(station==1)
	for(int subSector=0; subSector<2 && fpga<5; subSector++)
	  srLUTs_[fpga++][1] = new CSCSectorReceiverLUT(endcap,sector,subSector+1,
							station, srLUTset, TMB07);
      else
	srLUTs_[fpga++][1]   = new CSCSectorReceiverLUT(endcap,  sector,   0, 
							station, srLUTset, TMB07);
    }

  // negative endcap
  endcap = 2; 
  for(int station=1,fpga=0; station<=4 && fpga<5; station++)
    {
      if(station==1)
	for(int subSector=0; subSector<2 && fpga<5; subSector++)
	  srLUTs_[fpga++][0] = new CSCSectorReceiverLUT(endcap,sector,subSector+1,
							station, srLUTset, TMB07);
      else
	srLUTs_[fpga++][0]   = new CSCSectorReceiverLUT(endcap,  sector,   0, 
							station, srLUTset, TMB07);
    }
  // -----------------------------------------------------------------------------
  
  edm::ParameterSet serviceParameters = pset.getParameter<edm::ParameterSet>("ServiceParameters");
  theService = new MuonServiceProxy(serviceParameters);

}

void TrigEff::endJob() {
  //free the CSCTF array of pointers
  for(int j=0; j<2; j++) 
    for(int i=0; i<5; i++) 
      delete srLUTs_[i][j]; 
  
  // delete ts;
  // delete tpts;
}


// destructor
TrigEff::~TrigEff(void){ 
  file->Write(); 
  file->Close(); 
}

// analyze
void TrigEff::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  //Get the Magnetic field from the setup
  iSetup.get<IdealMagneticFieldRecord>().get(theBField);
  // Get the GlobalTrackingGeometry from the setup
  iSetup.get<GlobalTrackingGeometryRecord>().get(theTrackingGeometry);
  theService->update(iSetup);

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
 
  //GP start
  // Get the CSC Geometry 
  iSetup.get<MuonGeometryRecord>().get(cscGeom);
  //GP end  
 
  // ============================================================================
  // fill the reco muon block
  // Second: get global muons, stand alone muons, and track
  Handle<MuonCollection>  muons;
  if( muonsTag.label() != "null" ) iEvent.getByLabel(muonsTag, muons);

  if (printLevel > 0) cout<<"============ FILLING MUONS  ================"<<endl; 
  if( muons.isValid() ){

    muonsInit();

    //fill the muons
    fillMuons(muons);
  
    // ==========================================================================
    //
    // look at SEGMENTs (from the CSC Validation, A. Kubik)
    //
    // ==========================================================================
    // get CSC segment collection
    edm::Handle<CSCSegmentCollection> cscSegments;
    //cscSegTag = cms.InputTag("cscSegments"),
    if( cscSegTag.label() != "null" ) iEvent.getByLabel(cscSegTag, cscSegments);
    if( cscSegments.isValid()) fillSegments(cscSegments, cscGeom);
   
    Handle<CSCCorrelatedLCTDigiCollection> CSCTFlcts;
    if( csctfLctsTag.label() != "null" ) 
      iEvent.getByLabel(csctfLctsTag, CSCTFlcts);

    if( cscSegments.isValid() && CSCTFlcts.isValid() ) 
      fillSegmentsMuons(muons, cscSegments, cscGeom, CSCTFlcts);
      
    
    recoMuons->Fill();
    muonsDel();
  }
  

  
  if (printLevel > 0) {
    cout<<"======================================================"<<endl; 
    cout<<"Level-1 extra: "<<endl;
  }

  // ============================================================================
  // fill all the l1 information per event in a ttree
  Handle<l1extra::L1MuonParticleCollection> l1muons;
  if( L1extraTag.label() != "null" ) iEvent.getByLabel(L1extraTag, l1muons);

  // Fetch all Level-1 muon candidates
  if( l1muons.isValid() ){
    
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

  // ============================================================================
  // fill the csctf raw block
  Handle<L1CSCTrackCollection> CSCTFtracks;
  if( csctfTag.label() != "null" ) {
    iEvent.getByLabel(csctfTag, CSCTFtracks);
    if (printLevel > 0) cout<<"============ FILLING CSCTF RAW  ============\n";
  }

  if( CSCTFtracks.isValid() ){
    
    //initialize
    csctfInit();

    if( iSetup.get< L1MuTriggerScalesRcd > ().cacheIdentifier() != m_scalesCacheID ||
        iSetup.get< L1MuTriggerPtScaleRcd >().cacheIdentifier() != m_ptScaleCacheID ){
      
      ESHandle< L1MuTriggerScales > scales;
      iSetup.get< L1MuTriggerScalesRcd >().get(scales);
      ts = scales.product();
      ESHandle< L1MuTriggerPtScale > ptscales;
      iSetup.get< L1MuTriggerPtScaleRcd >().get(ptscales);
      tpts = ptscales.product();
      m_scalesCacheID  = iSetup.get< L1MuTriggerScalesRcd  >().cacheIdentifier();
      m_ptScaleCacheID = iSetup.get< L1MuTriggerPtScaleRcd >().cacheIdentifier();
      
      std::cout  << "Changing triggerscales and triggerptscales...\n";
    }    
    
    //fill the csctf information
    fillCSCTF(CSCTFtracks,
              ts,
              tpts,
              srLUTs_);

    Handle<CSCCorrelatedLCTDigiCollection> CSCTFlcts;
    if( csctfLctsTag.label() != "null" ) {
      iEvent.getByLabel(csctfLctsTag, CSCTFlcts);
      if (printLevel > 0) cout<<"========== FILLING ALL LCTS RAW  ==========\n";
    }

    if( CSCTFlcts.isValid() ) {
      //fill all the LCTs csctf information
      fillAllLCTs(CSCTFlcts, srLUTs_);
    }
    else
      cout << "Invalid CSCCorrelatedLCTDigiCollection... skipping it\n";

    
    // fill the ttree
    csctfTTree->Fill();

    //clean-up the pointers
    csctfDel();
  } 
  else 
    cout << "Invalid L1CSCTrackCollection... skipping it\n";
  
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
    int isAllArb = muon::isGoodMuon(*muon, muon::AllArbitrated);
    isMuonAllArbitrated -> push_back(isAllArb);
    int isTMLSAT = muon::isGoodMuon(*muon, muon::TMLastStationAngTight);
    isTMLastStationAngTight -> push_back(isTMLSAT);
    int isGBLPT = muon::isGoodMuon(*muon, muon::GlobalMuonPromptTight);
    isGlobalMuonPromptTight -> push_back(isGBLPT);

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

    
    // z planes
    int endcapPlane=0;
    if (trackRef->eta() > 0) endcapPlane=1; 
    if (trackRef->eta() < 0) endcapPlane=-1; 

    float zzPlaneME11 = endcapPlane*585;  
    float zzPlaneME1  = endcapPlane*615;  
    float zzPlaneME2  = endcapPlane*830;  
    float zzPlaneME3  = endcapPlane*935;  

    // ... to ME1/1
    // track at ME1/1 surface, +/-5.85 m - extrapolation
    tsos = surfExtrapTrkSam(trackRef, zzPlaneME11);  
    if (tsos.isValid()) {
      double xx = tsos.globalPosition().x();
      double yy = tsos.globalPosition().y();
      double zz = tsos.globalPosition().z();
      
      muons_x_me11->push_back(xx);	
      muons_y_me11->push_back(yy);	
      muons_z_me11->push_back(zz);	

      double rr = sqrt(xx*xx + yy*yy);
      double cosphi = xx/rr;
      
      if (yy>=0) muons_phi_me11->push_back(acos(cosphi));
      else       muons_phi_me11->push_back(2*PI-acos(cosphi));

    
      double abspseta = -log( tan( atan(fabs(rr/zz))/2.0 ) );
      if (zz>=0) muons_eta_me11->push_back(abspseta);
      else       muons_eta_me11->push_back(-abspseta);

      if (printLevel>0) {
        cout<<"I am projection the track to the ME";
        if (endcapPlane>0) cout << "+";
        if (endcapPlane<0) cout << "-";
        cout<< "1/1 surface" << endl;
        cout<< "zzPlaneM11=" << zzPlaneME11 << endl;

        cout << "muons_x_me11:" << xx << endl;
        cout << "muons_y_me11:" << yy << endl;
        cout << "muons_z_me11:" << zz << endl;      

        if (yy>=0) cout << "muons_phi_me11:" << acos(cosphi) << endl;
        else       cout << "muons_phi_me11:" << 2*PI-acos(cosphi) << endl;        
        if (zz>=0) cout << "muons_eta_me11:" <<  abspseta << endl;
        else       cout << "muons_eta_me11:" << -abspseta << endl;
      }
    }
    else {
      if (printLevel>0) cout << "extrapolation to ME1/1 NOT valid\n";
      muons_x_me11->push_back(-999);	
      muons_y_me11->push_back(-999);	
      muons_z_me11->push_back(-999);	
      muons_phi_me11->push_back(-999);
      muons_eta_me11->push_back(-999);
    }


    // ... to ME1
    // track at ME1 surface, +/-6.15 m - extrapolation
    tsos = surfExtrapTrkSam(trackRef, zzPlaneME1);  
    if (tsos.isValid()) {
      double xx = tsos.globalPosition().x();
      double yy = tsos.globalPosition().y();
      double zz = tsos.globalPosition().z();
      
      muons_x_me1->push_back(xx);	
      muons_y_me1->push_back(yy);	
      muons_z_me1->push_back(zz);	

      double rr = sqrt(xx*xx + yy*yy);
      double cosphi = xx/rr;
      
      if (yy>=0) muons_phi_me1->push_back(acos(cosphi));
      else       muons_phi_me1->push_back(2*PI-acos(cosphi));

    
      double abspseta = -log( tan( atan(fabs(rr/zz))/2.0 ) );
      if (zz>=0) muons_eta_me1->push_back(abspseta);
      else       muons_eta_me1->push_back(-abspseta);

      if (printLevel>0) {
          cout<<"I am projection the track to the ME";
        if (endcapPlane>0) cout << "+";
        if (endcapPlane<0) cout << "-";
        cout<< "1 surface" << endl;
        cout<< "zzPlaneM1=" << zzPlaneME1 << endl;

        cout << "muons_x_me1:" << xx << endl;
        cout << "muons_y_me1:" << yy << endl;
        cout << "muons_z_me1:" << zz << endl;      
        if (yy>=0) cout << "muons_phi_me1:" << acos(cosphi) << endl;
        else       cout << "muons_phi_me1:" << 2*PI-acos(cosphi) << endl;        
        if (zz>=0) cout << "muons_eta_me1:" <<  abspseta << endl;
        else       cout << "muons_eta_me1:" << -abspseta << endl;
      }
    }
    else {
      if (printLevel>0) cout << "extrapolation to ME1 NOT valid\n";
      muons_x_me1->push_back(-999);	
      muons_y_me1->push_back(-999);	
      muons_z_me1->push_back(-999);	
      muons_phi_me1->push_back(-999);
      muons_eta_me1->push_back(-999);
    }


    // ... to ME2
    // track at ME2 surface, +/-8.30 m - extrapolation
    tsos = surfExtrapTrkSam(trackRef, zzPlaneME2);  
    if (tsos.isValid()) {
      double xx = tsos.globalPosition().x();
      double yy = tsos.globalPosition().y();
      double zz = tsos.globalPosition().z();
      
      muons_x_me2->push_back(xx);	
      muons_y_me2->push_back(yy);	
      muons_z_me2->push_back(zz);	

      double rr = sqrt(xx*xx + yy*yy);
      double cosphi = xx/rr;
      
      if (yy>=0) muons_phi_me2->push_back(acos(cosphi));
      else       muons_phi_me2->push_back(2*PI-acos(cosphi));

    
      double abspseta = -log( tan( atan(fabs(rr/zz))/2.0 ) );
      if (zz>=0) muons_eta_me2->push_back(abspseta);
      else       muons_eta_me2->push_back(-abspseta);

      if (printLevel>0) {
        cout<<"I am projection the track to the ME";
        if (endcapPlane>0) cout << "+";
        if (endcapPlane<0) cout << "-";
        cout<< "2 surface" << endl;
        cout<< "zzPlaneM2=" << zzPlaneME2 << endl;
        
        cout << "muons_x_me2:" << xx << endl;
        cout << "muons_y_me2:" << yy << endl;
        cout << "muons_z_me2:" << zz << endl;      
        if (yy>=0) cout << "muons_phi_me2:" << acos(cosphi) << endl;
        else       cout << "muons_phi_me2:" << 2*PI-acos(cosphi) << endl;        
        if (zz>=0) cout << "muons_eta_me2:" << abspseta << endl;
        else       cout << "muons_eta_me2:" << -abspseta << endl;
      }
    }
    else {
      if (printLevel>0) cout << "extrapolation to ME2 NOT valid\n";
      muons_x_me2->push_back(-999);	
      muons_y_me2->push_back(-999);	
      muons_z_me2->push_back(-999);	
      muons_phi_me2->push_back(-999);
      muons_eta_me2->push_back(-999);
    }
    
    // ... to ME3
    // track at ME3 surface, +/-9.35 m - extrapolation
    tsos = surfExtrapTrkSam(trackRef, zzPlaneME3);  
    if (tsos.isValid()) {
      double xx = tsos.globalPosition().x();
      double yy = tsos.globalPosition().y();
      double zz = tsos.globalPosition().z();
      
      muons_x_me3->push_back(xx);	
      muons_y_me3->push_back(yy);	
      muons_z_me3->push_back(zz);	
      
      double rr = sqrt(xx*xx + yy*yy);
      double cosphi = xx/rr;
      
      if (yy>=0) muons_phi_me3->push_back(acos(cosphi));
      else       muons_phi_me3->push_back(2*PI-acos(cosphi));
      
      
      double abspseta = -log( tan( atan(fabs(rr/zz))/2.0 ) );
      if (zz>=0) muons_eta_me3->push_back(abspseta);
      else       muons_eta_me3->push_back(-abspseta);
      
      if (printLevel>0) {
        cout<<"I am projection the track to the ME";
        if (endcapPlane>0) cout << "+";
        if (endcapPlane<0) cout << "-";
        cout<< "3 surface" << endl;
        cout<< "zzPlaneM3=" << zzPlaneME3 << endl;
                
        cout << "muons_x_me3:" << xx << endl;
        cout << "muons_y_me3:" << yy << endl;
        cout << "muons_z_me3:" << zz << endl;      
        if (yy>=0) cout << "muons_phi_me3:" << acos(cosphi) << endl;
        else       cout << "muons_phi_me3:" << 2*PI-acos(cosphi) << endl;        
        if (zz>=0) cout << "muons_eta_me3:" << abspseta << endl;
        else       cout << "muons_eta_me3:" << -abspseta << endl;
      }
    }
    else {
      if (printLevel>0) cout << "extrapolation to ME3 NOT valid\n";
      muons_x_me3->push_back(-999);	
      muons_y_me3->push_back(-999);	
      muons_z_me3->push_back(-999);	
      muons_phi_me3->push_back(-999);
      muons_eta_me3->push_back(-999);
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
        float  localEtaRCH = -999; 
        float  localPhiRCH = -999;
        float globalEtaRCH = -999; 
        float globalPhiRCH = -999;
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

        // for debugging purpouses
        if (printLevel > 0) {
          cout << "rhitlocal.x="<<rhitlocal.x();
          cout << "rhitlocal.y="<<rhitlocal.y();
          cout << "rhitlocal.z="<<rhitlocal.z();
        }

        GlobalPoint gp = GlobalPoint(0.0, 0.0, 0.0);
        
        const CSCChamber* cscchamber = cscGeom->chamber(id);
        
        if (!cscchamber) continue;
          
        gp = cscchamber->toGlobal(rhitlocal);
        
        // identify the rechit position
        int pos = ( ((id.station()-1)*6) + (id.layer()-1) ) + (MAX_CSC_RECHIT*whichMuon);
        
        // --------------------------------------------------
        // this part has to be deprecated once we are sure of 
        // the Matrix usage ;)
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
               <<      " RCH X: " << gp.x()  
               <<      " RCH Y: " << gp.y()  
               <<      " RCH Z: " << gp.z()  
               << endl;
        }
        
        // -------------------------------------------------- 
        // new Matrix block
        if (whichMuon < MAX_MUONS && iRechit < MAX_CSC_RECHIT) {
          rchEtaMatrixLocal[whichMuon][iRechit] = rhitlocal.eta(); 
          rchPhiMatrixLocal[whichMuon][iRechit] = rhitlocal.phi();
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
          localEtaRCH      = rhitlocal.eta();
          localPhiRCH      = rhitlocal.phi();
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
      rchEtaLocal -> push_back(localEtaRCH);     
      rchPhiLocal -> push_back(localPhiRCH);     
      rchEta      -> push_back(globalEtaRCH);     
      rchPhi      -> push_back(globalPhiRCH);     
      rchStation  -> push_back(globalStationRCH);
      rchChamber  -> push_back(globalChamberRCH);
      rchRing     -> push_back(globalRingRCH);
      rchLayer    -> push_back(globalLayerRCH);
    
      if (printLevel > 0) {
        cout << "\n######### CLOSER #########";
        cout << "\n - RCH Type: " << globalTypeRCH
             <<     " RCH Eta Local: " << localEtaRCH	    
             <<     " RCH Phi Local: " << localPhiRCH
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
        cout << "Is Muon All Arbitrated? " 
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
          
          /* original coordinates are local to the chamber
            float xx = chamber->x;
            float yy = chamber->y;
          */

          GlobalPoint gpChb = theService->trackingGeometry()->idToDet(cscId)->surface().toGlobal(LocalPoint(chamber->x,chamber->y,0));
          
          float xx = gpChb.x();
          float yy = gpChb.y();
          float zz = gpChb.z();
        
          float rr = sqrt(xx*xx + yy*yy);
          float phi = gpChb.phi();
          if (phi < 0   ) phi += 2*PI;
          if (phi > 2*PI) phi -= 2*PI;
          float eta = gpChb.eta();


          //look at the segments
          for ( std::vector<reco::MuonSegmentMatch>::const_iterator segment = chamber->segmentMatches.begin();
                segment != chamber->segmentMatches.end(); ++segment ){
            
            if (printLevel > 0) cout << "Segment:" << iSegment+1 << endl;

            int   segIsArb=(segment->Arbitrated)>>8;

            if (!segIsArb) continue;

            /* original coordinates are local to the chamber
               float segX=segment->x;
               float segY=segment->y;
            */
            
            // global coordinates
            GlobalPoint gpSeg = theService->trackingGeometry()->idToDet(cscId)->surface().toGlobal(LocalPoint(segment->x,segment->y,0));

            float segX=gpSeg.x();
            float segY=gpSeg.y();
            float segZ=gpSeg.z();
            
            float segR = sqrt(segX*segX + segY*segY);
            float segPhi = gpSeg.phi();
            if (segPhi < 0   ) segPhi += 2*PI;
            if (segPhi > 2*PI) segPhi -= 2*PI;
            float segEta = gpSeg.eta();

            float segDxDz = segment->dXdZ;
            float segDyDz = segment->dYdZ;
            float segDxDzErr = segment->dXdZErr;
            float segDyDzErr = segment->dYdZErr;


            if (printLevel>0) {
              cout << "###### IS CSC ########"                << endl;
              cout << "trkSegChamberId:"     << chamberId     << endl;
              cout << "trkSegRing:"          << ring          << endl;
              cout << "trkSegStation:"       << station       << endl;
              cout << "trkSegEndcap:"        << endcap        << endl;
              cout << "trkSegTriggerSector:" << triggerSector << endl;
              cout << "trkSegTriggerCscId :" << triggerCscId  << endl;
            
              cout<< "trkSegXfromMatch:"     << xx   << endl;
              cout<< "trkSegYfromMatch:"     << yy   << endl;
              cout<< "trkSegZfromMatch:"     << zz   << endl;
              cout<< "trkSegRfromMatch:"     << rr   << endl;
              
              cout << "trkSegPhifromMatch:" << phi << endl;
              cout << "trkSegEtafromMatch:" << eta << endl;

              cout << "SEGMENT:"                  << endl;
              cout << "trkSegIsArb: " << segIsArb << endl; 
              cout << "trkSegX: "     << segX     << endl; 
              cout << "trkSegY: "     << segY     << endl; 
              cout << "trkSegZ: "     << segZ     << endl; 
              cout << "trkSegR: "     << segR     << endl; 
              
              cout << "segPhi:"       << segPhi       << endl;
              cout << "segEta:"       << segEta       << endl;

              cout << "segDxDz:"      << segDxDz      << endl;
              cout << "segDyDz:"      << segDyDz      << endl;
              cout << "segDxDzErr:"   << segDxDzErr   << endl;
              cout << "segDyDzErr:"   << segDyDzErr   << endl;
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
              trkSegZfromMatch[whichMuon][iSegment] = zz;
              trkSegRfromMatch[whichMuon][iSegment] = rr;
              
              trkSegPhifromMatch[whichMuon][iSegment] = phi;
              trkSegEtafromMatch[whichMuon][iSegment] = eta;

              
              trkSegIsArb[whichMuon][iSegment]=segIsArb;
              trkSegX[whichMuon][iSegment]=segX;
              trkSegY[whichMuon][iSegment]=segY;
              trkSegZ[whichMuon][iSegment]=segZ;
              trkSegR[whichMuon][iSegment]=segR;
            
              trkSegPhi[whichMuon][iSegment]= segPhi;
              trkSegEta[whichMuon][iSegment]= segEta;

              trkSegDxDz[whichMuon][iSegment]   = segDxDz;
              trkSegDyDz[whichMuon][iSegment]   = segDyDz;
              trkSegDxDzErr[whichMuon][iSegment]= segDxDzErr;
              trkSegDyDzErr[whichMuon][iSegment]= segDyDzErr;


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
        if (printLevel>0) cout << "trkNSegs="<<iSegment<<endl;
        trkNSegs->push_back(iSegment);
      }// is muon good?
      else trkNSegs->push_back(-999);
    }//isTrackerMuon
    else trkNSegs->push_back(-999);
  }

}


void TrigEff::muonsInit()
{
  
  isGlobalMuon        = new vector<int>;	  
  isTrackerMuon	      = new vector<int>;
  isStandAloneMuon    = new vector<int>; 
  isMuonAllArbitrated = new vector<int>; 
  isTMLastStationAngTight = new vector<int>; 
  isGlobalMuonPromptTight = new vector<int>; 

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

  //gmrEnergy            = new vector<float>;
  //gmrDEnergy           = new vector<float>;
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

  //stdEnergy    = new vector<float>;
  //stdDEnergy   = new vector<float>;
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
   
  //trkEnergy            = new vector<float>;
  //trkDEnergy           = new vector<float>;
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
      trkSegZfromMatch[row][col] = -999;  
      trkSegRfromMatch[row][col] = -999;  
      trkSegPhifromMatch[row][col] = -999;
      trkSegEtafromMatch[row][col] = -999;

      trkSegIsArb[row][col] = -999;
      trkSegX[row][col]     = -999;
      trkSegY[row][col]     = -999;
      trkSegZ[row][col]     = -999;
      trkSegR[row][col]     = -999;
      trkSegPhi[row][col]   = -999;
      trkSegEta[row][col]   = -999;

      trkSegDxDz[row][col]     = -999;
      trkSegDyDz[row][col]     = -999;
      trkSegDxDzErr[row][col]  = -999;
      trkSegDyDzErr[row][col]  = -999;
    }
  // ------------------------------------------------------  


  rchCSCtype  = new vector<int>  ; 
  rchEtaLocal = new vector<float>;     
  rchPhiLocal = new vector<float>;     
  rchEta      = new vector<float>;     
  rchPhi      = new vector<float>;     
  rchPhi_02PI = new vector<float>;
  rchStation  = new vector<int>;
  rchChamber  = new vector<int>;
  rchRing     = new vector<int>;
  rchLayer    = new vector<int>;

  rchMuonSize = new vector<int>;
  for (int row=0; row < MAX_MUONS; row++) 
    for (int col=0; col < MAX_CSC_RECHIT; col++) {
      rchEtaMatrixLocal[row][col] = -999;
      rchPhiMatrixLocal[row][col] = -999;
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
    muons_x_me11 = new vector<float>;
    muons_y_me11 = new vector<float>;
    muons_z_me11 = new vector<float>;
  muons_phi_me11 = new vector<float>;
  muons_eta_me11 = new vector<float>;

  // propagation to ME1
    muons_x_me1 = new vector<float>;
    muons_y_me1 = new vector<float>;
    muons_z_me1 = new vector<float>;
  muons_phi_me1 = new vector<float>;
  muons_eta_me1 = new vector<float>;
                 
  // propagation to ME2
    muons_x_me2 = new vector<float>;
    muons_y_me2 = new vector<float>;
    muons_z_me2 = new vector<float>;
  muons_phi_me2 = new vector<float>;
  muons_eta_me2 = new vector<float>;

  // propagation to ME3                 
    muons_x_me3 = new vector<float>;
    muons_y_me3 = new vector<float>;
    muons_z_me3 = new vector<float>;
  muons_phi_me3 = new vector<float>;
  muons_eta_me3 = new vector<float>;
                 
  // csc segments
  cscsegs_loc_x = new vector<float>;
  cscsegs_loc_y = new vector<float>;
  cscsegs_loc_z = new vector<float>;

  cscsegs_loc_theta = new vector<float>;
  cscsegs_loc_eta   = new vector<float>;
  cscsegs_loc_phi   = new vector<float>;

  cscsegs_loc_dir_theta = new vector<float>;
  cscsegs_loc_dir_eta   = new vector<float>;
  cscsegs_loc_dir_phi   = new vector<float>;

  cscsegs_gbl_x = new vector<float>;
  cscsegs_gbl_y = new vector<float>;
  cscsegs_gbl_z = new vector<float>;

  cscsegs_gbl_theta = new vector<float>;
  cscsegs_gbl_eta   = new vector<float>;
  cscsegs_gbl_phi   = new vector<float>;

  cscsegs_gbl_dir_theta = new vector<float>;
  cscsegs_gbl_dir_eta   = new vector<float>;
  cscsegs_gbl_dir_phi   = new vector<float>;

  cscsegs_endcap  = new vector<int>;
  cscsegs_station = new vector<int>;
  cscsegs_ring    = new vector<int>;
  cscsegs_chamber = new vector<int>;

  muonNsegs = new std::vector<int>; 

  for (int row=0; row < MAX_MUONS; row++) 
    for (int col=0; col < MAX_SEGS_STD; col++) {
 
      muon_cscsegs_loc_x[row][col] = -999;
      muon_cscsegs_loc_y[row][col] = -999;
      muon_cscsegs_loc_eta[row][col] = -999;
      muon_cscsegs_loc_phi[row][col] = -999;
      muon_cscsegs_loc_dir_eta[row][col] = -999;
      muon_cscsegs_loc_dir_phi[row][col] = -999;

      muon_cscsegs_gbl_x[row][col] = -999;
      muon_cscsegs_gbl_y[row][col] = -999;
      muon_cscsegs_gbl_eta[row][col] = -999;
      muon_cscsegs_gbl_phi[row][col] = -999;
      muon_cscsegs_gbl_dir_eta[row][col] = -999;
      muon_cscsegs_gbl_dir_phi[row][col] = -999;

      muon_cscsegs_dxdz[row][col] = -999;
      muon_cscsegs_dydz[row][col] = -999;
      muon_cscsegs_dxdzErr[row][col] = -999;
      muon_cscsegs_dydzErr[row][col] = -999;

      muon_cscsegs_endcap[row][col] = -999;
      muon_cscsegs_station[row][col] = -999;
      muon_cscsegs_ring[row][col] = -999;
      muon_cscsegs_chamber[row][col] = -999;
      muon_cscsegs_nhits[row][col] = -999;

      muon_cscsegs_islctable[row][col] = -999;
      muon_cscsegs_ismatched[row][col] = -999;
      muon_cscsegs_lctId[row][col] = -999;
      muon_cscsegs_nmatched[row][col] = -999;
    }
}


void TrigEff::muonsDel() {
  
  vector<int>().swap(*isGlobalMuon);	  
  vector<int>().swap(*isTrackerMuon);
  vector<int>().swap(*isStandAloneMuon); 
  vector<int>().swap(*isMuonAllArbitrated); 
  vector<int>().swap(*isTMLastStationAngTight); 
  vector<int>().swap(*isGlobalMuonPromptTight); 

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

  //vector<float>().swap(*gmrEnergy);
  //vector<float>().swap(*gmrDEnergy);
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


  //vector<float>().swap(*stdEnergy);
  //vector<float>().swap(*stdDEnergy);
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
   

  //vector<float>().swap(*trkEnergy);
  //vector<float>().swap(*trkDEnergy);
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
  vector<float>().swap(*rchEtaLocal);     
  vector<float>().swap(*rchPhiLocal);     
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

    vector<float>().swap(*muons_x_me11);
    vector<float>().swap(*muons_y_me11);
    vector<float>().swap(*muons_z_me11);
  vector<float>().swap(*muons_phi_me11);
  vector<float>().swap(*muons_eta_me11);
                                    
    vector<float>().swap(*muons_x_me1);
    vector<float>().swap(*muons_y_me1);
    vector<float>().swap(*muons_z_me1);
  vector<float>().swap(*muons_phi_me1);
  vector<float>().swap(*muons_eta_me1);
                                    
    vector<float>().swap(*muons_x_me2);
    vector<float>().swap(*muons_y_me2);
    vector<float>().swap(*muons_z_me2);
  vector<float>().swap(*muons_phi_me2);
  vector<float>().swap(*muons_eta_me2);
                       
    vector<float>().swap(*muons_x_me3);
    vector<float>().swap(*muons_y_me3);
    vector<float>().swap(*muons_z_me3);
  vector<float>().swap(*muons_phi_me3);
  vector<float>().swap(*muons_eta_me3);

  vector<float>().swap(*cscsegs_loc_x);
  vector<float>().swap(*cscsegs_loc_y);
  vector<float>().swap(*cscsegs_loc_z);

  vector<float>().swap(*cscsegs_loc_theta);
  vector<float>().swap(*cscsegs_loc_eta);
  vector<float>().swap(*cscsegs_loc_phi);

  vector<float>().swap(*cscsegs_loc_dir_theta);
  vector<float>().swap(*cscsegs_loc_dir_eta);
  vector<float>().swap(*cscsegs_loc_dir_phi);

  vector<float>().swap(*cscsegs_gbl_x);
  vector<float>().swap(*cscsegs_gbl_y);
  vector<float>().swap(*cscsegs_gbl_z);

  vector<float>().swap(*cscsegs_gbl_theta);
  vector<float>().swap(*cscsegs_gbl_eta);
  vector<float>().swap(*cscsegs_gbl_phi);

  vector<float>().swap(*cscsegs_gbl_dir_theta);
  vector<float>().swap(*cscsegs_gbl_dir_eta);
  vector<float>().swap(*cscsegs_gbl_dir_phi);

  vector<int>().swap(*cscsegs_endcap);
  vector<int>().swap(*cscsegs_station);
  vector<int>().swap(*cscsegs_ring);
  vector<int>().swap(*cscsegs_chamber);

  std::vector<int>().swap(*muonNsegs); 
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

void TrigEff::csctfInit() {	
  
  EndcapTrk     = new vector<int>;
  SectorTrk     = new vector<int>;  
  BxTrk         = new vector<int>;  
     	    
  me1ID         = new vector<int>;  
  me2ID         = new vector<int>;
  me3ID         = new vector<int>;  
  me4ID         = new vector<int>;  
  mb1ID         = new vector<int>;  

  OutputLinkTrk = new vector<int>;      

  ModeTrk       = new vector<int>  ;
  EtaTrk        = new vector<float>;  
  PhiTrk        = new vector<float>;  
  PhiTrk_02PI   = new vector<float>;  
  PtTrk         = new vector<float>;  
  
  ChargeValidTrk = new vector<int>;
  ChargeTrk      = new vector<int>;
  QualityTrk     = new vector<int>;
  ForRTrk        = new vector<int>;
  Phi23Trk       = new vector<int>;
  Phi12Trk       = new vector<int>;
  PhiSignTrk     = new vector<int>;

  EtaBitTrk      = new vector<int>;
  PhiBitTrk      = new vector<int>;
  PtBitTrk       = new vector<int>;

  NumLCTsTrk    = new vector<int>; 
    
  /*
    trLctEndcap.ResizeTo(MAX_CSCTF_TRK,MAX_LCTS_PER_TRK); 
    trLctSector.ResizeTo(MAX_CSCTF_TRK,MAX_LCTS_PER_TRK); 
    trLctSubSector.ResizeTo(MAX_CSCTF_TRK,MAX_LCTS_PER_TRK); 
    trLctBx.ResizeTo(MAX_CSCTF_TRK,MAX_LCTS_PER_TRK); 
    trLctBx0.ResizeTo(MAX_CSCTF_TRK,MAX_LCTS_PER_TRK); 
  
    trLctStation.ResizeTo(MAX_CSCTF_TRK,MAX_LCTS_PER_TRK); 
    trLctRing.ResizeTo(MAX_CSCTF_TRK,MAX_LCTS_PER_TRK); 
    trLctChamber.ResizeTo(MAX_CSCTF_TRK,MAX_LCTS_PER_TRK); 
    trLctTriggerCSCID.ResizeTo(MAX_CSCTF_TRK,MAX_LCTS_PER_TRK); 
    trLctFpga.ResizeTo(MAX_CSCTF_TRK,MAX_LCTS_PER_TRK);	  
    
    trLctlocalPhi.ResizeTo(MAX_CSCTF_TRK,MAX_LCTS_PER_TRK); 
    trLctglobalPhi.ResizeTo(MAX_CSCTF_TRK,MAX_LCTS_PER_TRK);   
    trLctglobalEta.ResizeTo(MAX_CSCTF_TRK,MAX_LCTS_PER_TRK); 
  
    trLctstripNum.ResizeTo(MAX_CSCTF_TRK,MAX_LCTS_PER_TRK);   
    trLctwireGroup.ResizeTo(MAX_CSCTF_TRK,MAX_LCTS_PER_TRK);     
  */

  // all LCTs
  lctEndcap       = new vector<int>; 
  lctSector       = new vector<int>; 
  lctSubSector    = new vector<int>; 
  lctBx           = new vector<int>; 
  lctBx0          = new vector<int>; 
  lctStation      = new vector<int>; 
  lctRing         = new vector<int>; 
  lctChamber      = new vector<int>; 
  lctTriggerCSCID = new vector<int>; 
  lctFpga         = new vector<int>; 
  lctlocalPhi     = new vector<int>; 
  lctglobalPhi    = new vector<int>; 
  lctglobalEta    = new vector<int>; 
  lctstripNum     = new vector<int>; 
  lctwireGroup    = new vector<int>;   


}

void TrigEff::csctfDel() {	

  vector<int>().swap(*EndcapTrk); 
  vector<int>().swap(*SectorTrk); 
  vector<int>().swap(*BxTrk    ); 

  vector<int>().swap(*me1ID); 
  vector<int>().swap(*me2ID); 
  vector<int>().swap(*me3ID); 
  vector<int>().swap(*me4ID); 
  vector<int>().swap(*mb1ID); 

  vector<int>().swap(*OutputLinkTrk);

  vector<int>  ().swap(*ModeTrk    ); 
  vector<float>().swap(*EtaTrk     ); 
  vector<float>().swap(*PhiTrk     ); 
  vector<float>().swap(*PhiTrk_02PI); 
  vector<float>().swap(*PtTrk      ); 

  vector<int>().swap(*ChargeTrk      ); 
  vector<int>().swap(*ChargeValidTrk ); 
  vector<int>().swap(*QualityTrk     ); 
  vector<int>().swap(*ForRTrk        ); 
  vector<int>().swap(*Phi23Trk       ); 
  vector<int>().swap(*Phi12Trk       ); 
  vector<int>().swap(*PhiSignTrk     ); 

  vector<int>().swap(*EtaBitTrk);
  vector<int>().swap(*PhiBitTrk);
  vector<int>().swap(*PtBitTrk );

  vector<int>().swap(*NumLCTsTrk );

  vector<int>().swap(*lctEndcap); 
  vector<int>().swap(*lctSector); 
  vector<int>().swap(*lctSubSector); 
  vector<int>().swap(*lctBx); 
  vector<int>().swap(*lctBx0); 
  vector<int>().swap(*lctStation); 
  vector<int>().swap(*lctRing); 
  vector<int>().swap(*lctChamber); 
  vector<int>().swap(*lctTriggerCSCID); 
  vector<int>().swap(*lctFpga);     
  vector<int>().swap(*lctlocalPhi); 
  vector<int>().swap(*lctglobalPhi);   
  vector<int>().swap(*lctglobalEta); 
  vector<int>().swap(*lctstripNum);   
  vector<int>().swap(*lctwireGroup);   

}


void TrigEff::fillCSCTF(const edm::Handle<L1CSCTrackCollection> tracks,
                        const L1MuTriggerScales  *ts, 
                        const L1MuTriggerPtScale *tpts, 
                        CSCSectorReceiverLUT* srLUTs_[5][2]) {
  
  int nTrk=0; 
  // loop over CSCTF tracks
  for(L1CSCTrackCollection::const_iterator trk=tracks->begin(); 
      trk<tracks->end(); trk++){
    
    nTrk++;

    // Standard Pt LUTs	  
    edm::ParameterSet ptLUTset;
    ptLUTset.addUntrackedParameter<bool>("ReadLUTs", false);
    ptLUTset.addUntrackedParameter<bool>("Binary",   false);
    ptLUTset.addUntrackedParameter<std::string>("LUTPath", "./");

    /*      
    // Reading Pt LUTs from file
    edm::ParameterSet ptLUTset;	    
    ptLUTset.addUntrackedParameter<bool>("ReadPtLUT", true);
    ptLUTset.addUntrackedParameter<bool>("isBinary",  true);
    ptLUTset.addUntrackedParameter<bool>("isBeamStartConf", true);
    
    edm::FileInPath pt_lut_file;  
    pt_lut_file = edm::FileInPath("../results/cfg/L1CSCPtLUT.bin");
    ptLUTset.addParameter<edm::FileInPath>("PtLUTFile",pt_lut_file);
    */
    
    //CSCTFPtLUT* ptLUT = new CSCTFPtLUT(ptLUTset, ts, tpts);
    CSCTFPtLUT ptLUT(ptLUTset, ts, tpts);
    
    //define the variables
    // trk->first.endcap() = 2 for - endcap
    //                     = 1 for + endcap
    int    trEndcap = (trk->first.endcap()==2 ? trk->first.endcap()-3 : trk->first.endcap());
    int    trSector = 6*(trk->first.endcap()-1)+trk->first.sector();
    int    trBX     = trk->first.BX();
    int    trMe1ID = trk->first.me1ID();
    int    trMe2ID = trk->first.me2ID();
    int    trMe3ID = trk->first.me3ID();
    int    trMe4ID = trk->first.me4ID();
    int    trMb1ID = trk->first.mb1ID();

    int    trEtaBit = trk->first.eta_packed();
    int    trPhiBit = trk->first.localPhi();
    int    trCharge = trk->first.chargeValue();

    int    trRank   = trk->first.rank();    

    // PtAddress gives an handle on other parameters
    ptadd thePtAddress(trk->first.ptLUTAddress());
	
    int trPhiSign = thePtAddress.delta_phi_sign;
    int trPhi12   = thePtAddress.delta_phi_12;
    int trPhi23   = thePtAddress.delta_phi_23;
    int trMode    = thePtAddress.track_mode;
    int trForR    = thePtAddress.track_fr;
	
    int trOutputLink = trk->first.outputLink();

    //Pt needs some more workaround since it is not in the unpacked data
    //ptdat thePtData  = ptLUT->Pt(thePtAddress);
    ptdat thePtData  = ptLUT.Pt(thePtAddress);
	
    int trPtBit       = 0;
    int trQuality     = 0;
    int trChargeValid = 0;

    // front or rear bit? 
    if (thePtAddress.track_fr) {
      trPtBit = (thePtData.front_rank&0x1f);
      trQuality = ((thePtData.front_rank>>5)&0x3);
      trChargeValid = thePtData.charge_valid_front;
    } else {
      trPtBit = (thePtData.rear_rank&0x1f);
      trQuality = ((thePtData.rear_rank>>5)&0x3);
      trChargeValid = thePtData.charge_valid_rear;
    }
	  
    if (printLevel > 0) {
      cout << " ===== CSCTF BIT VALUES ====\n"
           << " Track #"        << nTrk
           << "\n Endcap: "     << trEndcap
           << "\n Sector: "     << trSector 
           << "\n BX: "         << trBX 
           << "\n me1ID: "      << trMe1ID 			
	   << "\n me2ID: "      << trMe2ID
	   << "\n me3ID: "      << trMe3ID
	   << "\n me4ID: "      << trMe4ID 
           << "\n mb1ID: "      << trMb1ID 
           << "\n OutputLink: " << trOutputLink
        
           << "\n Charge: "     << trCharge
           << "\n DPhiSign: "   << trPhiSign
	   << "\n DPhi12: "     << trPhi12
	   << "\n DPhi23: "     << trPhi23
	   << "\n ForR: "       << trForR
                                
           << "\n Mode: "       << trMode 
           << "\n Quality: "    << trQuality 
           << "\n Rank : "      << trRank 
           << "\n ChargeValid: "<< trChargeValid
           << "\n Eta: "        << trEtaBit 
           << "\n Phi: "        << trPhiBit
           << "\n Pt: "         << trPtBit
           << endl;         
    }

    //... in radians
    // Type 2 is CSC
    double trEta = ts->getRegionalEtaScale(2)->getCenter(trk->first.eta_packed());
    double trPhi = ts->getPhiScale()->getLowEdge(trk->first.localPhi());
    //Phi in one sector varies from [0,62] degrees -> Rescale manually to global coords.
    double trPhi02PI = fmod(trPhi + 
                            ((trSector-1)*TMath::Pi()/3) + //sector 1 starts at 15 degrees 
                            (TMath::Pi()/12) , 2*TMath::Pi());
    
    // convert the Pt in human readable values (GeV/c)
    double trPt = tpts->getPtScale()->getLowEdge(trPtBit); 


    if (printLevel > 0) {
      cout << " ===== CSCTF TRK SCALES ====\n"
           << "\n Eta(scales): " << trEta 
           << "\n Phi(scales): " << trPhi02PI
           << "\n Pt(scales): "  << trPt
           << endl;
    }

    // fill the track vectors
    EndcapTrk      -> push_back(trEndcap);
    SectorTrk      -> push_back(trSector); 
    BxTrk          -> push_back(trBX    ); 
    me1ID          -> push_back(trMe1ID);  
    me2ID          -> push_back(trMe2ID);
    me3ID          -> push_back(trMe3ID); 
    me4ID          -> push_back(trMe4ID); 
    mb1ID          -> push_back(trMb1ID);  
    OutputLinkTrk  -> push_back(trOutputLink);  

    ModeTrk        -> push_back(trMode   );
    EtaTrk         -> push_back(trEta    ); 
    PhiTrk         -> push_back(trPhi    ); 
    PhiTrk_02PI    -> push_back(trPhi02PI); 
    PtTrk          -> push_back(trPt     );  

    ChargeValidTrk -> push_back(trChargeValid);
    ChargeTrk      -> push_back(trCharge     );
    QualityTrk     -> push_back(trQuality    );
    ForRTrk        -> push_back(trForR       );
    Phi23Trk       -> push_back(trPhi23      );
    Phi12Trk       -> push_back(trPhi12      );
    PhiSignTrk     -> push_back(trPhiSign    );
      	       
    EtaBitTrk      -> push_back(trEtaBit); 
    PhiBitTrk      -> push_back(trPhiBit); 
    PtBitTrk       -> push_back(trPtBit );  


    // For each trk, get the list of its LCTs
    CSCCorrelatedLCTDigiCollection lctsOfTracks = trk -> second;
  
    int LctTrkId_ = 0;

    for(CSCCorrelatedLCTDigiCollection::DigiRangeIterator lctOfTrks = lctsOfTracks.begin(); 
        lctOfTrks  != lctsOfTracks.end()  ; lctOfTrks++){
	  
      int lctTrkId = 0;	
			
      CSCCorrelatedLCTDigiCollection::Range lctRange = 
        lctsOfTracks.get((*lctOfTrks).first);
	  
      for(CSCCorrelatedLCTDigiCollection::const_iterator 
            lctTrk = lctRange.first ; 
          lctTrk  != lctRange.second; lctTrk++, lctTrkId++){

        trLctEndcap[nTrk-1][LctTrkId_] = (*lctOfTrks).first.zendcap();
        // sector (pos: 1->6, neg: 7 -> 12)
        if ((*lctOfTrks).first.zendcap() > 0)
          trLctSector[nTrk-1][LctTrkId_] = (*lctOfTrks).first.triggerSector();
        else
          trLctSector[nTrk-1][LctTrkId_] = 6+(*lctOfTrks).first.triggerSector();
   	  
        trLctSubSector[nTrk-1][LctTrkId_] = CSCTriggerNumbering::triggerSubSectorFromLabels((*lctOfTrks).first);;
        trLctBx[nTrk-1][LctTrkId_] = lctTrk -> getBX();
        trLctBx0[nTrk-1][LctTrkId_] = lctTrk -> getBX0();
	      
        trLctStation[nTrk-1][LctTrkId_] = (*lctOfTrks).first.station();
        trLctRing[nTrk-1][LctTrkId_] = (*lctOfTrks).first.ring();
        trLctChamber[nTrk-1][LctTrkId_] = (*lctOfTrks).first.chamber();
        trLctTriggerCSCID[nTrk-1][LctTrkId_] = (*lctOfTrks).first.triggerCscId();
        trLctFpga[nTrk-1][LctTrkId_] = 
          ( trLctSubSector[nTrk-1][LctTrkId_] ? trLctSubSector[nTrk-1][LctTrkId_] : (*lctOfTrks).first.station()+1);

        // Check if DetId is within range
        if( trLctSector[nTrk-1][LctTrkId_] < 1       || trLctSector[nTrk-1][LctTrkId_] > 12 || 
            trLctStation[nTrk-1][LctTrkId_] < 1      || trLctStation[nTrk-1][LctTrkId_] >  4 || 
            trLctTriggerCSCID[nTrk-1][LctTrkId_] < 1 || trLctTriggerCSCID[nTrk-1][LctTrkId_] >  9 || 
            lctTrkId < 0 || lctTrkId >  1 ){
          cout <<"  TRACK ERROR: CSC digi are out of range: ";
          continue;
        }

        // handles not to overload the method: mostly for readability	      
        int endcap = (*lctOfTrks).first.zendcap();
        if (endcap < 0) endcap = 0; 

        int StationLctTrk  = (*lctOfTrks).first.station();
        int CscIdLctTrk    = (*lctOfTrks).first.triggerCscId();
        int SubSectorLctTrk = 
          CSCTriggerNumbering::triggerSubSectorFromLabels((*lctOfTrks).first);
        
        int FPGALctTrk    = 
          ( SubSectorLctTrk ? SubSectorLctTrk-1 : StationLctTrk );
        
	      
        // local Phi
        lclphidat lclPhi;
	
        try {
          
          trLctstripNum[nTrk-1][LctTrkId_] = lctTrk->getStrip();
          lclPhi = srLUTs_[FPGALctTrk][endcap] -> localPhi(lctTrk->getStrip(), 
                                                           lctTrk->getPattern(), 
                                                           lctTrk->getQuality(), 
                                                           lctTrk->getBend() );
          
          trLctlocalPhi[nTrk-1][LctTrkId_] = lclPhi.phi_local;
        } 
        catch(...) { 
          bzero(&lclPhi,sizeof(lclPhi)); 
          trLctlocalPhi[nTrk-1][LctTrkId_] = -999;
        }
		
        
        // Global Phi
        gblphidat gblPhi;
	
        try {
          
          trLctwireGroup[nTrk-1][LctTrkId_] = lctTrk->getKeyWG();
          gblPhi = srLUTs_[FPGALctTrk][endcap] -> globalPhiME(lclPhi.phi_local  , 
                                                              lctTrk->getKeyWG(), 
                                                              CscIdLctTrk);
          
          trLctglobalPhi[nTrk-1][LctTrkId_] = gblPhi.global_phi;
          
        } catch(...) { 
          bzero(&gblPhi,sizeof(gblPhi)); 
          trLctglobalPhi[nTrk-1][LctTrkId_] = -999;
        }
		
        // Global Eta
        gbletadat gblEta;
	
        try {
          gblEta = srLUTs_[FPGALctTrk][endcap] -> globalEtaME(lclPhi.phi_bend_local, 
                                                              lclPhi.phi_local     , 
                                                              lctTrk->getKeyWG()   , 
                                                              CscIdLctTrk);
          trLctglobalEta[nTrk-1][LctTrkId_] = gblEta.global_eta;
        } 	  
        catch(...) { 
          bzero(&gblEta,sizeof(gblEta)); 
          trLctglobalEta[nTrk-1][LctTrkId_] = -999;
        } 
	
        ++LctTrkId_;
	
      } // for(CSCCorrelatedLCTDigiCollection::const_iterator lctTrk 
    } // for(CSCCorrelatedLCTDigiCollection::DigiRangeIterator lctOfTrks
	  
    NumLCTsTrk->push_back(LctTrkId_);
  	     
  } //for(L1CSCTrackCollection::const_iterator trk=csctfTrks->begin(); trk<csctfTrks->end(); trk++,nTrk++){

  SizeTrk = nTrk;

}

void TrigEff::fillAllLCTs(const edm::Handle<CSCCorrelatedLCTDigiCollection> corrlcts, 
                          CSCSectorReceiverLUT* srLUTs_[5][2]) {
  
  // ALL LCTs
  int nLCT=0;
  for(CSCCorrelatedLCTDigiCollection::DigiRangeIterator 
        corrLct = corrlcts.product()->begin(); 
      corrLct != corrlcts.product()->end()  ; corrLct++){
      
    nLCT++;
 
    int lctId = 0;	
      
    CSCCorrelatedLCTDigiCollection::Range lctRange = 
      corrlcts.product()->get((*corrLct).first);
			
    for(CSCCorrelatedLCTDigiCollection::const_iterator 
          lct = lctRange.first ; 
        lct != lctRange.second; lct++, lctId++){
		
      lctEndcap->push_back((*corrLct).first.zendcap());
      if ((*corrLct).first.zendcap() > 0)
        lctSector->push_back((*corrLct).first.triggerSector());
      else
        lctSector->push_back(6+(*corrLct).first.triggerSector());
	
      lctSubSector->push_back(CSCTriggerNumbering::triggerSubSectorFromLabels((*corrLct).first));
      lctBx->push_back(lct->getBX());
      lctBx0->push_back(lct->getBX0());
	
      lctStation->push_back((*corrLct).first.station());
      lctRing->push_back((*corrLct).first.ring());
      lctChamber->push_back((*corrLct).first.chamber());
      lctTriggerCSCID->push_back((*corrLct).first.triggerCscId());
      lctFpga->push_back((lctSubSector->back() ? lctSubSector->back() : (*corrLct).first.station()+1));
	

      // Check if DetId is within range
      if( lctSector->back() < 1 || lctSector->back() > 12 || 
          lctStation->back() < 1 || lctStation->back() >  4 || 
          lctTriggerCSCID->back() < 1 || lctTriggerCSCID->back() >  9 || 
          lctId < 0 || lctId >  1 ){
	  
        cout<<"  LCT ERROR: CSC digi are out of range: ";
	  
        continue;
      }

      // handles not to overload the method: mostly for readability	      
      int endcap = (*corrLct).first.zendcap();
      if (endcap < 0) endcap = 0; 
	
      int StationLct  = (*corrLct).first.station();
      int CscIdLct    = (*corrLct).first.triggerCscId();
      int SubSectorLct = 
        CSCTriggerNumbering::triggerSubSectorFromLabels((*corrLct).first);
	      
      int FPGALct    = ( SubSectorLct ? SubSectorLct-1 : StationLct );
	
	
      // local Phi
      lclphidat lclPhi;
	
      try {
	  
        lctstripNum->push_back(lct->getStrip());
        lclPhi = srLUTs_[FPGALct][endcap] -> localPhi(lct->getStrip(), 
                                                      lct->getPattern(), 
                                                      lct->getQuality(), 
                                                      lct->getBend() );
	  
        lctlocalPhi->push_back(lclPhi.phi_local);
      } 
      catch(...) { 
        bzero(&lclPhi,sizeof(lclPhi)); 
        lctlocalPhi->push_back(-999);
      }
		
	
      // Global Phi
      gblphidat gblPhi;
	
      try {
	  
        lctwireGroup->push_back(lct->getKeyWG());

        gblPhi = srLUTs_[FPGALct][endcap] -> globalPhiME(lclPhi.phi_local  , 
                                                         lct->getKeyWG(), 
                                                         CscIdLct);
	  
        lctglobalPhi->push_back(gblPhi.global_phi);
	  
      } catch(...) { 
        bzero(&gblPhi,sizeof(gblPhi)); 
        lctglobalPhi->push_back(-999);
      }
	
      // Global Eta
      gbletadat gblEta;
	
      try {
        gblEta = srLUTs_[FPGALct][endcap] -> globalEtaME(lclPhi.phi_bend_local, 
                                                         lclPhi.phi_local     , 
                                                         lct->getKeyWG()   , 
                                                         CscIdLct);
        lctglobalEta->push_back(gblEta.global_eta);
      } 	  
      catch(...) { 
        bzero(&gblEta,sizeof(gblEta)); 
        lctglobalEta->push_back(-999);
      } 
      
    } // for(CSCCorrelatedLCTDigiCollection::const_iterator lct 
  } // for(CSCCorrelatedLCTDigiCollection::DigiRangeIterator lct

  SizeLCTs = nLCT;

}


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


void TrigEff::fillSegments(edm::Handle<CSCSegmentCollection> cscSegments, 
                           edm::ESHandle<CSCGeometry> cscGeom){

  // get CSC segment collection
  int segsSize = cscSegments->size();
  if (printLevel>0) cout << "segsSize:"<< segsSize << endl;

  // -----------------------
  // loop over segments
  // -----------------------
  int iSegment = 0;
  for(CSCSegmentCollection::const_iterator dSiter=cscSegments->begin(); 
      dSiter != cscSegments->end(); 
      dSiter++) {
    
    iSegment++;
    //
    CSCDetId id  = (CSCDetId)(*dSiter).cscDetId();
    int kEndcap  = id.endcap();
    int kRing    = id.ring();
    int kStation = id.station();
    int kChamber = id.chamber();
    
    cscsegs_endcap->push_back(kEndcap);
    cscsegs_station->push_back(kRing);
    cscsegs_ring->push_back(kStation);
    cscsegs_chamber->push_back(kChamber);

    //
    //float chisq    = (*dSiter).chi2();
    //int nhits      = (*dSiter).nRecHits();
    //int nDOF       = 2*nhits-4;
    //double chisqProb = ChiSquaredProbability( (double)chisq, nDOF );
    LocalPoint localPos = (*dSiter).localPosition();
    float segX = localPos.x();
    float segY = localPos.y();
    float segZ = localPos.z();
    float segPhi = localPos.phi();
    float segEta = localPos.eta();
    float segTheta = localPos.theta();
      
    cscsegs_loc_x->push_back(segX);
    cscsegs_loc_y->push_back(segY);
    cscsegs_loc_z->push_back(segZ);

    cscsegs_loc_theta->push_back(segTheta);
    cscsegs_loc_eta  ->push_back(segEta);
    cscsegs_loc_phi  ->push_back(segPhi);
      
    LocalVector segDir = (*dSiter).localDirection();
    double dirTheta = segDir.theta();
    double dirEta   = segDir.eta();
    double dirPhi   = segDir.phi();

    cscsegs_loc_dir_theta->push_back(dirTheta);
    cscsegs_loc_dir_eta  ->push_back(dirEta);
    cscsegs_loc_dir_phi  ->push_back(dirPhi);


    if (printLevel>0) {
      cout << "iSegment: " << iSegment << endl;
      cout << "kEndcap : "<< kEndcap  << endl;
      cout << "kRing   : "<< kRing    << endl;
      cout << "kStation: "<< kStation << endl;
      cout << "kChamber: "<< kChamber << endl;
      cout << "segX: "<< segX << endl;
      cout << "segY: "<< segY << endl;
      cout << "segZ: "<< segZ << endl;
      cout << "segPhi: "<< segPhi << endl;
      cout << "segEta: "<< segEta << endl;
      cout << "segTheta: "<< segTheta << endl;
      cout << "dirTheta: "<< dirTheta << endl;
      cout << "dirEta: "  << dirEta   << endl;
      cout << "dirPhi: "  << dirPhi   << endl;
    }

    // global transformation
    float globX = 0.;
    float globY = 0.;
    float globZ = 0.;
    float globSegPhi  = 0.;
    float globSegTheta = 0.;
    float globSegEta   = 0.;
    float globDirPhi  = 0.;
    float globDirTheta = 0.;
    float globDirEta   = 0.;

    
    const CSCChamber* cscchamber = cscGeom->chamber(id);
    if (cscchamber) {
      GlobalPoint globalPosition = cscchamber->toGlobal(localPos);
      globX   = globalPosition.x();
      globY   = globalPosition.y();
      globZ   = globalPosition.z();
      globSegPhi   = globalPosition.phi();
      globSegEta   = globalPosition.eta();
      globSegTheta = globalPosition.theta();

      if (printLevel>0) {
        cout << "globX  : " << globX   << endl;
        cout << "globY  : " << globY   << endl;
        cout << "globZ  : " << globZ   << endl;
        cout << "globSegPhi: " << globSegPhi << endl;
        cout << "globSegEta: " << globSegEta << endl;
      }

      GlobalVector globalDirection = cscchamber->toGlobal(segDir);
      globDirTheta = globalDirection.theta();
      globDirEta   = globalDirection.eta();
      globDirPhi   = globalDirection.phi();

      cscsegs_gbl_x->push_back(globX);
      cscsegs_gbl_y->push_back(globY);
      cscsegs_gbl_z->push_back(globZ);

      cscsegs_gbl_theta->push_back(globSegTheta);
      cscsegs_gbl_eta  ->push_back(globSegEta  );
      cscsegs_gbl_phi  ->push_back(globSegPhi  );

      cscsegs_gbl_dir_theta->push_back(globDirTheta);
      cscsegs_gbl_dir_eta  ->push_back(globDirEta  );
      cscsegs_gbl_dir_phi  ->push_back(globDirPhi  );
    } 
    else {
      cscsegs_gbl_x->push_back(-999);
      cscsegs_gbl_y->push_back(-999);
      cscsegs_gbl_z->push_back(-999);

      cscsegs_gbl_theta->push_back(-999);
      cscsegs_gbl_eta  ->push_back(-999);
      cscsegs_gbl_phi  ->push_back(-999);

      cscsegs_gbl_dir_theta->push_back(-999);
      cscsegs_gbl_dir_eta  ->push_back(-999);
      cscsegs_gbl_dir_phi  ->push_back(-999);
    }
  }
}
  

// From Ivan
bool TrigEff::isLCTAble ( const CSCSegment &segment ){

  if (segment . cscDetId().ring() == 4) cout << "IsLCTAble FOUND ME11a\n";

  if (segment . cscDetId().ring() == 4 ) return false;
  if (segment . nRecHits() < 4 )         return false;

  int thisStation = segment.cscDetId().station();
  
  int keyStrip = 999, keyWireg = 999;

  const std::vector<CSCRecHit2D>& theHits = segment . specificRecHits();
	    
  std::vector<CSCRecHit2D>::const_iterator hitIter;

  for (hitIter = theHits.begin(); hitIter!= theHits.end(); hitIter++){

    if (hitIter -> cscDetId(). layer() != 3)
      continue;
    
    keyStrip = halfStrip(hitIter -> channels()[1], hitIter -> positionWithinStrip() );
    
    keyWireg =  hitIter -> wgroups()[0];
    
  }	    

  int alctEnvelopes[6] = { 2, 1, 0, 1, 2, 2 };
  int clctEnvelopes[6] = { 5, 2, 0, 2, 4, 5 };


  int hitsFidAlct = 0;
  int hitsFidClct = 0;

  for (hitIter = theHits.begin(); hitIter!= theHits.end(); hitIter++){

    int thisLayer = hitIter -> cscDetId() . layer();

    int delWgroup = hitIter -> wgroups()[0] - keyWireg;
    int delStrip  = halfStrip(hitIter -> channels()[1],
			      hitIter -> positionWithinStrip()) - keyStrip;
    
    if (thisLayer <=3)
      delWgroup = -delWgroup;

    if (thisStation == 3)
      delWgroup = -delWgroup;

    if (thisStation == 4)
      delWgroup = -delWgroup;

    if ( delWgroup >=0 )
      if ( delWgroup <= alctEnvelopes[ thisLayer - 1 ] ) hitsFidAlct++;

    if ( abs(delStrip)  <= clctEnvelopes[ thisLayer - 1 ] ) hitsFidClct++;

  }

  if ( hitsFidAlct < 4) return false;
  if ( hitsFidClct < 4) return false;
  
  return true;

}




/*
  if( segmentTag.label() != "null" ) {
  iEvent.getByLabel(segmentTag, cscSegments);
  if (printLevel >= 0) std::cout<<"========== Reading segments  ==========" << std::endl;
  }    
  
  if( !cscSegments.isValid() ) { 
  std::cout << "Invalid CSCSegmentCollection... skipping it" << std::endl;
  return;
  }

  if( cscLctsTag.label() != "null" ) {
  iEvent.getByLabel(cscLctsTag, CSCTFlcts);
  if (printLevel > 0) std::cout<<"========== FILLING ALL LCTS RAW  ==========" << std::endl;
  }
  
  if( !CSCTFlcts.isValid() ){
  std::cout << "bad LCT collection... skipping" << std::endl;
  return;
  }

*/

void  TrigEff::fillSegmentsMuons ( const edm::Handle<reco::MuonCollection> muons,
                                   edm::Handle<CSCSegmentCollection> cscSegments, 
                                   edm::ESHandle<CSCGeometry> cscGeom,
                                   const edm::Handle<CSCCorrelatedLCTDigiCollection> CSCTFlcts) {


  if (printLevel >= 0){
    std::cout << "\n============= fillSegmentsMuons =============\n";
   }
  
  // loop over the muons
  int whichMuon =-1;
  for (MuonCollection::const_iterator muon=muons->begin();
       muon!=muons->end(); muon++) {
    
    whichMuon++;
    if (printLevel) cout << "whichMuon="<<whichMuon<<endl;
    // given a muon candidate, loop over the segments and find the segments belonging
    // to the muon and see if they are "LCTAble"

    if (!muon->combinedMuon())   continue;
    if (!muon->standAloneMuon()) continue;

    // get the segments which match the muon candidate
    SegmentVector *segVect = SegmentsInMuon( &(*muon), &(*cscSegments) );
    
    // debugging
    if (printLevel>0)
      std::cout << "segVect -> size()=" << segVect -> size() << endl;

    muonNsegs->push_back( segVect -> size() );
  
    // sanity check
    if (segVect -> size() == 0){
      delete segVect; 
      continue;
    }
    
    int iSegment=-1;
    // loop over the segments
    for(SegmentVector::const_iterator segmentCSC = segVect -> begin();
        segmentCSC != segVect -> end(); segmentCSC++){

      iSegment++;

      if ( iSegment > (MAX_SEGS_STD-1) ) {
        std::cout << "the muon has " << iSegment+1 << ", but the MAX allowed is "
                  << MAX_SEGS_STD << " -> Skipping the segment... " << std::endl;
      }

      LocalPoint localPos = (*segmentCSC)->cscsegcand.localPosition();
      LocalVector  segDir = (*segmentCSC)->cscsegcand.localDirection();    
      
      muon_cscsegs_loc_x      [whichMuon][iSegment] = localPos.x();
      muon_cscsegs_loc_y      [whichMuon][iSegment] = localPos.y();
      muon_cscsegs_loc_eta    [whichMuon][iSegment] = localPos.eta();
      muon_cscsegs_loc_phi    [whichMuon][iSegment] = localPos.phi();
      muon_cscsegs_loc_dir_eta[whichMuon][iSegment] = segDir.eta();
      muon_cscsegs_loc_dir_phi[whichMuon][iSegment] = segDir.phi();


      CSCDetId id  = (CSCDetId) (*segmentCSC)->cscsegcand.cscDetId();
        
      const CSCChamber* cscchamber = cscGeom->chamber(id);
   
      if (!cscchamber) { 
        std::cout << "cscchamber not valid" << std::endl;
        continue;
      }
      
      GlobalPoint globalPosition   = cscchamber->toGlobal(localPos);
      GlobalVector globalDirection = cscchamber->toGlobal(segDir);      


      muon_cscsegs_gbl_x      [whichMuon][iSegment] = globalPosition.x();
      muon_cscsegs_gbl_y      [whichMuon][iSegment] = globalPosition.y();
      muon_cscsegs_gbl_eta    [whichMuon][iSegment] = globalPosition.eta();
      muon_cscsegs_gbl_phi    [whichMuon][iSegment] = globalPosition.phi();
      muon_cscsegs_gbl_dir_eta[whichMuon][iSegment] = globalDirection.eta();
      muon_cscsegs_gbl_dir_phi[whichMuon][iSegment] = globalDirection.phi();

      cout << "GP debugging...\n";
      cout << "kEndcap =" << id.endcap()  << endl;
      cout << "kRing   =" << id.ring()    << endl;
      cout << "kStation=" << id.station() << endl;
      cout << "kChamber=" << id.chamber() << endl;
    
      muon_cscsegs_endcap [whichMuon][iSegment] = id.endcap(); 
      muon_cscsegs_station[whichMuon][iSegment] = id.station();
      muon_cscsegs_ring   [whichMuon][iSegment] = id.ring();   
      muon_cscsegs_chamber[whichMuon][iSegment] = id.chamber();
      muon_cscsegs_nhits  [whichMuon][iSegment] = (*segmentCSC)->cscsegcand.nRecHits();

      muon_cscsegs_dxdz   [whichMuon][iSegment] = -999; //(*segmentCSC)->cscsegcand.dXdZ;
      muon_cscsegs_dydz   [whichMuon][iSegment] = -999; //(*segmentCSC)->cscsegcand.dYdZ;
      muon_cscsegs_dxdzErr[whichMuon][iSegment] = -999; //(*segmentCSC)->cscsegcand.dXdZErr;
      muon_cscsegs_dydzErr[whichMuon][iSegment] = -999; //(*segmentCSC)->cscsegcand.dYdZErr;

      bool isTriggerAble = isLCTAble( (*segmentCSC)->cscsegcand );
      bool isLCTMatched  = isMatched( (*segmentCSC)->cscsegcand, CSCTFlcts );
          
      if (printLevel>0) {
        std::cout << std::endl;
        std::cout <<"segmentCSC->nMatchedHits=" 
                  << (*segmentCSC)->nMatchedHits << std::endl; 
        std::cout <<"isLCTAble?=" << isTriggerAble << std::endl; 
        std::cout <<"isMatched?=" << isLCTMatched  << std::endl; 
      }

      muon_cscsegs_islctable[whichMuon][iSegment] = isTriggerAble;
      muon_cscsegs_ismatched[whichMuon][iSegment] = isLCTMatched;

      if (printLevel>0) {
        std::cout << "muon_cscsegs_islctable[" << whichMuon << "][" << iSegment<< "]:" << muon_cscsegs_islctable[whichMuon][iSegment] << std::endl;
        std::cout << "muon_cscsegs_ismatched[" << whichMuon << "][" << iSegment<< "]:" << muon_cscsegs_ismatched[whichMuon][iSegment] << std::endl;
      }

      int whichLCT=-999;// find the corresponding LCT in the list
      if (isLCTMatched) {
        
        int iLCT=-1;
        for(CSCCorrelatedLCTDigiCollection::DigiRangeIterator 
              corrLct = CSCTFlcts.product()->begin(); 
            corrLct != CSCTFlcts.product()->end()  ; corrLct++){
          
          iLCT++;
      
          CSCCorrelatedLCTDigiCollection::Range lctRange1 = 
            CSCTFlcts.product()->get( (*segmentCSC)->cscsegcand.cscDetId() );
          
          CSCCorrelatedLCTDigiCollection::Range lctRange2 = 
            CSCTFlcts.product()->get((*corrLct).first);
        	
          // find the matching one
          if (lctRange1 == lctRange2) {
            
            whichLCT=iLCT;
            
            if (printLevel > 0 )
              std::cout<< "Corresponds to LCT number:" << whichLCT << std::endl;
          }


        }
      }

      muon_cscsegs_lctId    [whichMuon][iSegment] = whichLCT;
      muon_cscsegs_nmatched [whichMuon][iSegment] = (*segmentCSC)->nMatchedHits;

    }

    delete segVect;
  }
}



TrigEff::SegmentVector* TrigEff::SegmentsInMuon(const reco::Muon* muon, 
                                                const CSCSegmentCollection* segments ){

  
  TrigEff::SegmentVector *retVal = new TrigEff::SegmentVector();

  
  bool isMuonStd=false; // has the muon a standalone component
  if (muon->combinedMuon().isNonnull() || muon->standAloneMuon().isNonnull()) 
    isMuonStd=true;
 
  // return empty vector if the muon is not interesting
  if (!isMuonStd) return retVal;

  
  int iSegment=0;
  //loop over the segments
  for (CSCSegmentCollection::const_iterator segIter = segments -> begin();
       segIter != segments -> end(); segIter++){

    int nHits=segIter -> nRecHits();

    if (printLevel>0) {
      std::cout << " ======================================== " << std::endl;
      std::cout << "Segment in CSC:" << iSegment++ << std::endl;
      std::cout << "# segs hits:"<< nHits << std::endl;
    }

    // no need to verify the segment is in CSC. It is the default assumption in
    // the collection ;-)
    const std::vector<CSCRecHit2D>& theHits = segIter -> specificRecHits();
        
    std::vector<CSCRecHit2D>::const_iterator hitIter;
    
    int nMuonMatchedHits=0;
    int iHit=0;
    // loop over the segments hits
    for (hitIter = theHits.begin(); hitIter!= theHits.end(); hitIter++){

      if (printLevel>0) std::cout << "iHit:" << iHit++ << " , ";
      // check it the hit will match the standalone muon component
      bool isHitMatched=false;

      
      LocalPoint seghitlocal = hitIter -> localPosition();

      double segHitX = seghitlocal.x();      
      double segHitY = seghitlocal.y();      

      if (printLevel > 0)
        std::cout << "segHitX="<<segHitX <<  ", "
                  << "segHitY="<<segHitY;
        

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

        if (segHitX==rhitlocal.x() &&
            segHitY==rhitlocal.y()  )
          isHitMatched=true;
      
      }// end loop trackingRecHit_iterator


      if (printLevel>0){
        if (isHitMatched) std::cout<< " -> Matched" <<std::endl;
        else              std::cout<< " -> NOT Matched" <<std::endl;
      }

      if (isHitMatched) nMuonMatchedHits++;

    }
  

    if (printLevel>0) std::cout<< "segment has "  << nMuonMatchedHits
                               << " hits out of " << nHits 
                               << " matched"      << std::endl;
  

    // fill the the vector with the matching segments
    if (nMuonMatchedHits!=0) {
      Segment* segment = new Segment(*segIter,nMuonMatchedHits);
      retVal->push_back(segment);
    }

  }// end loop on the segments
  
  return retVal;
}

// From Ivan, match the LCTs to segment by looking at strip and WG
// no phi or eta... 
bool TrigEff::isMatched ( const CSCSegment &segment, 
                          edm::Handle<CSCCorrelatedLCTDigiCollection> CSCTFlcts ){
  
  bool retVal=false;

  int LCT_key_strip = -999;
  int LCT_key_wireg = -999;


  CSCCorrelatedLCTDigiCollection::Range lctRange = CSCTFlcts.product()->get( segment.cscDetId() );
  
  if ( printLevel > 0) {
    std::cout << " ============================================ " << std::endl;
    std::cout << " segment CSCDetId " << segment.cscDetId() << std::endl;
  }

  if (segment.cscDetId().ring() == 4) cout << "IsMatched: FOUND ME11a\n";
    
  for(CSCCorrelatedLCTDigiCollection::const_iterator 
        lct = lctRange.first ; lct != lctRange.second; lct++ ){
    
    int closestStrip = 999;
    int closestWireg = 999;
    
    if ( printLevel > 0)
      std::cout << (*lct) << std::endl;
  
    LCT_key_wireg = lct -> getKeyWG();
    LCT_key_strip = lct -> getStrip();
    
    bool foundKeyLayer = false;
    bool matchedWireg  = false;
    bool matchedStrip  = false;
    
    closestStrip = 999;
    closestWireg = 999;
    
    const std::vector<CSCRecHit2D>& theHits = segment.specificRecHits();
    
    std::vector<CSCRecHit2D>::const_iterator hitIter;
    
    for (hitIter = theHits.begin(); hitIter!= theHits.end(); hitIter++){
      
      bool thisHitStripMatch = false;
      bool thisHitWiregMatch = false;
      
      if ( (hitIter -> cscDetId() . layer() != 3) && 
           (hitIter -> cscDetId() . layer() != 4) ) continue;
      
      foundKeyLayer = true;
      
      CSCRecHit2D::ChannelContainer::const_iterator channIter;
      
      for (channIter = hitIter -> channels().begin(); channIter != hitIter -> channels().end();
           channIter++){
        
        if ( printLevel >0)
          std::cout << (*channIter) << ", ";
        
        int deltaStrip = LCT_key_strip - 2*(*channIter) + 1;
        
        //std::cout << "deltaStrip:" << deltaStrip << std::endl;
        if ( fabs(deltaStrip) < fabs(closestStrip) )
          closestStrip = deltaStrip;
        
        // if ( abs(deltaStrip) <= 2 ){
        // suggestion from Vadim...
        
        if ( abs(deltaStrip) <= 10 ){
          matchedStrip = true;
          thisHitStripMatch = true;
        }
       
      }
    
      if ( printLevel > 0)
        std::cout << std::endl;
    
      for (channIter = hitIter -> wgroups().begin(); channIter != hitIter -> wgroups().end();
           channIter++){
      
        if ( printLevel > 0)
          std::cout << (*channIter) << ", ";
       
        int deltaWire = LCT_key_wireg - (*channIter) + 1;
      
        //std::cout << "deltaWire" << deltaWire << std::endl;
      
        if ( abs(deltaWire) < abs(closestWireg) )
          closestWireg = deltaWire;
       
        if ( abs(deltaWire) <= 2 ){
          matchedWireg = true;
          thisHitWiregMatch = true;
        }
      
      }
     
      if ( printLevel > 0)
        std::cout << std::endl;
   
      if (thisHitWiregMatch && thisHitStripMatch){
        retVal = true;
      } 
    }
  }

  return retVal;
}


// GP's test
//int TrigEff::whichLCTMatches ( const CSCSegment &segment, edm::Handle<CSCCorrelatedLCTDigiCollection> CSCTFlcts ){
  
//   int  retValDefault = -999;
//   bool isMatch = false;

//   int LCT_key_strip = -999;
//   int LCT_key_wireg = -999;


//   //CSCCorrelatedLCTDigiCollection::Range lctRange = CSCTFlcts.product()->get( segment.cscDetId() );
  
//   if ( printLevel > 0) {
//     std::cout << " ============================================ " << std::endl;
//     std::cout << " segment CSCDetId " << segment.cscDetId() << std::endl;
//   }

//   int iLCT=-1;
//   if (printLevel>0) std::cout << "iLCT:" << iLCT << std::endl; 

//   for(CSCCorrelatedLCTDigiCollection::DigiRangeIterator 
//         corrLct = CSCTFlcts.product()->begin(); 
//       corrLct != CSCTFlcts.product()->end()  ; corrLct++){
      
//     iLCT++;
    
//     CSCCorrelatedLCTDigiCollection::Range lctRange = 
//       CSCTFlcts.product()->get((*corrLct).first);

//     for(CSCCorrelatedLCTDigiCollection::const_iterator 
//           lct = lctRange.first ; lct != lctRange.second; lct++ ){
    
//       if (printLevel>0) std::cout << "iLCT:" << iLCT << std::endl; 
    
      
//       int closestStrip = 999;
//       int closestWireg = 999;
      
//       if ( printLevel > 0)
//         std::cout << (*lct) << std::endl;
    
//       LCT_key_wireg = lct -> getKeyWG();
//       LCT_key_strip = lct -> getStrip();
      
//       bool foundKeyLayer = false;
//       bool matchedWireg  = false;
//       bool matchedStrip  = false;
      
//       closestStrip = 999;
//       closestWireg = 999;
      
//       const std::vector<CSCRecHit2D>& theHits = segment.specificRecHits();
      
//       std::vector<CSCRecHit2D>::const_iterator hitIter;
      
//       for (hitIter = theHits.begin(); hitIter!= theHits.end(); hitIter++){
        
//         bool thisHitStripMatch = false;
//         bool thisHitWiregMatch = false;
        
//         if ( (hitIter -> cscDetId() . layer() != 3) && 
//              (hitIter -> cscDetId() . layer() != 4) ) continue;
        
//         foundKeyLayer = true;
        
//         CSCRecHit2D::ChannelContainer::const_iterator channIter;
        
//         for (channIter = hitIter -> channels().begin(); channIter != hitIter -> channels().end();
//              channIter++){
          
//           if ( printLevel >0)
//             std::cout << (*channIter) << ", ";
          
//           int deltaStrip = LCT_key_strip - 2*(*channIter) + 1;
          
//           std::cout << "deltaStrip:" << deltaStrip << std::endl;
//           if ( fabs(deltaStrip) < fabs(closestStrip) )
//             closestStrip = deltaStrip;
          
//           // if ( abs(deltaStrip) <= 2 ){
//           // suggestion from Vadim...
          
//           if ( abs(deltaStrip) <= 10 ){
//             matchedStrip = true;
//             thisHitStripMatch = true;
//           }
	 
//         }
      
//         if ( printLevel > 0)
//           std::cout << std::endl;
      
//         for (channIter = hitIter -> wgroups().begin(); channIter != hitIter -> wgroups().end();
//              channIter++){
        
//           if ( printLevel > 0)
//             std::cout << (*channIter) << ", ";
	 
//           int deltaWire = LCT_key_wireg - (*channIter) + 1;
        
//           std::cout << "deltaWire" << deltaWire << std::endl;
	
//           if ( abs(deltaWire) < abs(closestWireg) )
//             closestWireg = deltaWire;
	 
//           if ( abs(deltaWire) <= 2 ){
//             matchedWireg = true;
//             thisHitWiregMatch = true;
//           }
        
//         }
       
//         if ( printLevel > 0)
//           std::cout << std::endl;
     
//         if (thisHitWiregMatch && thisHitStripMatch){
//           isMatch = true;
//         } 
//       }

//       if (isMatch){
//         if (printLevel>0) std::cout << " -> MATCHED" << std::endl; 
//         return iLCT;
//       }
//     }
//   }
  
//   return retValDefault;
  
//}
