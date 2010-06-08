//gSystem->CompileMacro("scripts/Collisions2010/trigEff/noTriggersId.C","Ok")

// ------------------------------------------------------------------------------
// TTree
// ------------------------------------------------------------------------------
#include<TTree.h>
#include<TH1.h>
#include<TH2.h>
#include<TGraphErrors.h>
#include<TCanvas.h>
#include<TLegend.h>
#include<vector>
#include<iostream>
#include<TMath.h>
#include<TROOT.h>
#include<TInterpreter.h>
#include<TStyle.h>
#include<TChain.h>
#include<TString.h>
#include<TPaveStats.h>
#include<TPad.h>
#include<TLatex.h>
// ---- CP error ----
#include <TGraphAsymmErrors.h>
#include <cmath>

#include "PhysicsTools/RooStatsCms/interface/ClopperPearsonBinomialInterval.h"
//--------

#include <stdio.h>
#include <string.h>

using namespace std;

#include "/afs/cern.ch/user/d/digiovan/scripts/Init/modifiedStyle.C"

// helper to identify the missing triggers
struct csctftrack {
  // general info
  int endcap;
  int sector;
  int lctEta[4];
  int lctPhi[4];

  // initialize
  csctftrack(){ 
    endcap=-999; 
    sector=-999; 
    for (int i=0; i<4; i++ ){
      lctEta[i]=-999;
      lctPhi[i]=-999;
    }
  }

  // compute the mode of the track
  int mode() {
    
    int retVal = -999;
    if (lctPhi[0]!=-999 && lctPhi[1]!=-999 && lctPhi[2]!=-999)                    return  2; //ME1-2-3-(4)
    if (lctPhi[0]!=-999 && lctPhi[1]!=-999 && lctPhi[2]==-999 && lctPhi[3]!=-999) return  3;//ME1-2-4
    if (lctPhi[0]!=-999 && lctPhi[1]==-999 && lctPhi[2]!=-999 && lctPhi[3]!=-999) return  4;//ME1-3-4
    if (lctPhi[0]==-999 && lctPhi[1]!=-999 && lctPhi[2]!=-999 && lctPhi[3]!=-999) return  5;//ME2-3-4
    if (lctPhi[0]!=-999 && lctPhi[1]!=-999 && lctPhi[2]==-999 && lctPhi[3]==-999) return  6;//ME1-2
    if (lctPhi[0]!=-999 && lctPhi[1]==-999 && lctPhi[2]!=-999 && lctPhi[3]==-999) return  7;//ME1-3
    if (lctPhi[0]==-999 && lctPhi[1]!=-999 && lctPhi[2]!=-999 && lctPhi[3]==-999) return  8;//ME2-3
    if (lctPhi[0]==-999 && lctPhi[1]!=-999 && lctPhi[2]==-999 && lctPhi[3]!=-999) return  9;//ME2-4
    if (lctPhi[0]==-999 && lctPhi[1]==-999 && lctPhi[2]!=-999 && lctPhi[3]!=-999) return 10;//ME3-4
    if (lctPhi[0]!=-999 && lctPhi[1]==-999 && lctPhi[2]==-999 && lctPhi[3]==-999) return 11;//ME1 singles
    if (lctPhi[0]!=-999 && lctPhi[1]==-999 && lctPhi[2]==-999 && lctPhi[3]!=-999) return 13;//ME1 singles

    return retVal;
  }
   
  // compute phi
  int phi() {
    if (mode() == 11 || 
        mode() == 13  ) 
      return lctPhi[0]; // valid if only singles in ME1
    else if (mode() ==  4 || 
             mode() ==  7 || 
             mode() == 10 ) return lctPhi[2];
    else return lctPhi[1];
  } 


  // compute eta
  int eta() {
    if (mode() == 11 || 
        mode() == 13  ) 
      return lctEta[0]; // valid if only singles in ME1
    else if (mode() ==  4 || 
             mode() ==  7 || 
             mode() == 10 ) return lctEta[2];
    else return lctEta[1];
  } 
};

// CSCTF nuances
bool IsCloseToEdge(csctftrack track);

// DR cut
double DR(double diffeta, double diffphi);


TChain* recoMuons;
TChain* csctfTTree;

// ------------------------------------------------------------------------------
// variables
// ------------------------------------------------------------------------------
double Pi  = TMath::Pi();
double PhiStep = (62*Pi/180)/4096;
double EtaStep = 1.6/128; //(2.5-0.9)/(2^7)

TLegend *tl;
// ------------------------------------------------------------------------------
// switches
// ------------------------------------------------------------------------------

const int MAX_MUONS = 100; 
const int MAX_CSC_RECHIT = 48;
const int MAX_TRK_SEGS = 100;
const int MAX_CSCTF_TRK = 36;
const int MAX_LCTS_PER_TRK = 4;
const int MAX_SEGS_STD=16; 


// main method
void noTriggersId(int printLevel = 0,
                     bool isMC   = !true) {
    

  gROOT->Clear();
  gStyle->SetOptStat(111111);  

  
  // ----------------------------------------------------------------------------
  // Construct the chain of input files
  // ----------------------------------------------------------------------------
  cout << "Loading chain reco... \n";
    
  // get the Nutple  
  #include "scripts/Collisions2010/trigEff/chains/collisionsChain_MinBiasMC_Reco.C"
  #include "scripts/Collisions2010/trigEff/chains/collisionsChain_MinimumBias_Reco.C"
  
  if (isMC) recoMuons  = collisionsChainRecoMC;
  else      recoMuons  = collisionsChainReco;
  cout << "recoMuons->GetEntries(): " << recoMuons->GetEntries() << endl;

  // RAW
  #include "scripts/Collisions2010/trigEff/chains/collisionsChain_MinimumBias_Raw.C"
  #include "scripts/Collisions2010/trigEff/chains/collisionsChain_MinBiasMC_Raw.C"
  cout << "Loading chain raw... \n";
  
  // get the Nutple  
  if (isMC) csctfTTree  = collisionsChainRawMC;
  else      csctfTTree  = collisionsChainRaw;
  cout << "csctfTTree->GetEntries(): " << csctfTTree->GetEntries() << endl;


  // ===========================================================================
  TString eps = "eps/Collisions2010/trigEff/";
  TString png = "png/Collisions2010/trigEff/";
  TString rootPlot = "rootPlot/Collisions2010/trigEff/";
  
  if (isMC) {
    eps+="MC-"; png+="MC-"; rootPlot+="MC-";
  }
  // ===========================================================================

  //----------------------------------------------------------------------------
  // Access the needed variables
  //----------------------------------------------------------------------------
  int    Run, Event, Bx, Lumi, muonSize;

  vector<int>*            isGlobalMuon = new vector<int>();
  vector<int>*        isStandAloneMuon = new vector<int>();
  vector<int>*           isTrackerMuon = new vector<int>();
  vector<int>* isTMLastStationAngTight = new vector<int>();
  vector<int>* isGlobalMuonPromptTight = new vector<int>();

  vector<float>*     ptReco  = new vector<float>();
  vector<float>*     etaReco = new vector<float>();
  vector<float>*     phiReco = new vector<float>();
  vector<float>* gmrChi2Norm = new vector<float>();
  vector<float>*       gmrDz = new vector<float>();
  vector<float>*       gmrD0 = new vector<float>();
  
  vector<float>*      stdEta = new vector<float>();
  vector<float>*      stdPt  = new vector<float>();
  vector<float>*      trkEta = new vector<float>();
  vector<float>*      trkPhi = new vector<float>();
  vector<float>*       trkPt = new vector<float>();
  vector<int>*    trkValHits = new vector<int>();
  vector<float>* trkChi2Norm = new vector<float>();
  vector<float>*       trkDz = new vector<float>();
  vector<float>*       trkD0 = new vector<float>();
  
  vector<float>*      rchEta = new vector<float>();
  vector<float>*      rchPhi = new vector<float>();
  vector<int>*    rchCSCtype = new vector<int>();
  
  // tracker muon variables
  vector<int>*    trkNSegs = new vector<int>();
  
  int      trkSegChamberId[MAX_MUONS][MAX_TRK_SEGS]; 
  int           trkSegRing[MAX_MUONS][MAX_TRK_SEGS];    
  int        trkSegStation[MAX_MUONS][MAX_TRK_SEGS]; 
  int         trkSegEndcap[MAX_MUONS][MAX_TRK_SEGS];  
  int  trkSegTriggerSector[MAX_MUONS][MAX_TRK_SEGS];
  int   trkSegTriggerCscId[MAX_MUONS][MAX_TRK_SEGS]; 
  float   trkSegXfromMatch[MAX_MUONS][MAX_TRK_SEGS];
  float   trkSegYfromMatch[MAX_MUONS][MAX_TRK_SEGS];
  float trkSegPhifromMatch[MAX_MUONS][MAX_TRK_SEGS];

  int trkSegIsArb[MAX_MUONS][MAX_TRK_SEGS];
  float   trkSegX[MAX_MUONS][MAX_TRK_SEGS];
  float   trkSegY[MAX_MUONS][MAX_TRK_SEGS];
  float   trkSegR[MAX_MUONS][MAX_TRK_SEGS];
  float trkSegPhi[MAX_MUONS][MAX_TRK_SEGS];
  float trkSegEta[MAX_MUONS][MAX_TRK_SEGS];

  std::vector<float>*    muons_x_me11 = new vector<float>();
  std::vector<float>*    muons_y_me11 = new vector<float>();
  std::vector<float>*    muons_z_me11 = new vector<float>();
  std::vector<float>*  muons_phi_me11 = new vector<float>();
  std::vector<float>*  muons_eta_me11 = new vector<float>();

  std::vector<float>*    muons_x_me1 = new vector<float>();
  std::vector<float>*    muons_y_me1 = new vector<float>();
  std::vector<float>*    muons_z_me1 = new vector<float>();
  std::vector<float>*  muons_phi_me1 = new vector<float>();
  std::vector<float>*  muons_eta_me1 = new vector<float>();
    
  //--------------------------------------------------------------------------
  // Record information about segments belonging to the STD muon component
  //--------------------------------------------------------------------------
  // how many segments are associated to the muon candidate
  std::vector<int>* muonNsegs = new vector<int>(); 

  // segment position information, local
  float       muon_cscsegs_loc_x[MAX_MUONS][MAX_SEGS_STD];
  float       muon_cscsegs_loc_y[MAX_MUONS][MAX_SEGS_STD];
  float     muon_cscsegs_loc_eta[MAX_MUONS][MAX_SEGS_STD];
  float     muon_cscsegs_loc_phi[MAX_MUONS][MAX_SEGS_STD];
  float muon_cscsegs_loc_dir_eta[MAX_MUONS][MAX_SEGS_STD];
  float muon_cscsegs_loc_dir_phi[MAX_MUONS][MAX_SEGS_STD];

  // segment position information, global
  float       muon_cscsegs_gbl_x[MAX_MUONS][MAX_SEGS_STD];
  float       muon_cscsegs_gbl_y[MAX_MUONS][MAX_SEGS_STD];
  float     muon_cscsegs_gbl_eta[MAX_MUONS][MAX_SEGS_STD];
  float     muon_cscsegs_gbl_phi[MAX_MUONS][MAX_SEGS_STD];
  float muon_cscsegs_gbl_dir_eta[MAX_MUONS][MAX_SEGS_STD];
  float muon_cscsegs_gbl_dir_phi[MAX_MUONS][MAX_SEGS_STD];

  // more on segment direction
  float    muon_cscsegs_dxdz[MAX_MUONS][MAX_SEGS_STD];
  float    muon_cscsegs_dydz[MAX_MUONS][MAX_SEGS_STD];
  float muon_cscsegs_dxdzErr[MAX_MUONS][MAX_SEGS_STD];
  float muon_cscsegs_dydzErr[MAX_MUONS][MAX_SEGS_STD];

  // general segment information
  int  muon_cscsegs_endcap[MAX_MUONS][MAX_SEGS_STD];
  int muon_cscsegs_station[MAX_MUONS][MAX_SEGS_STD];
  int    muon_cscsegs_ring[MAX_MUONS][MAX_SEGS_STD];
  int muon_cscsegs_chamber[MAX_MUONS][MAX_SEGS_STD];
  int   muon_cscsegs_nhits[MAX_MUONS][MAX_SEGS_STD];

  // isLCTAble
  int muon_cscsegs_islctable[MAX_MUONS][MAX_SEGS_STD];
  int muon_cscsegs_ismatched[MAX_MUONS][MAX_SEGS_STD];
  // lctId is the position of the lct in the all LCT collection
  // look also in fillAllLCTs
  int muon_cscsegs_lctId[MAX_MUONS][MAX_SEGS_STD];

  // number of hits belonging to the STD fit
  int muon_cscsegs_nmatched[MAX_MUONS][MAX_SEGS_STD];
  //--------------------------------------------------------------------------

  recoMuons->SetBranchAddress("Run"  , &Run  );
  recoMuons->SetBranchAddress("Event", &Event);
  recoMuons->SetBranchAddress("Bx"   , &Bx   );
  recoMuons->SetBranchAddress("Lumi" , &Lumi );

  recoMuons->SetBranchAddress("muonSize"               , &muonSize        );
  recoMuons->SetBranchAddress("isGlobalMuon"           , &isGlobalMuon    );
  recoMuons->SetBranchAddress("isStandAloneMuon"       , &isStandAloneMuon);
  recoMuons->SetBranchAddress("isTrackerMuon"          , &isTrackerMuon   );
  recoMuons->SetBranchAddress("isTMLastStationAngTight", &isTMLastStationAngTight);
  recoMuons->SetBranchAddress("isGlobalMuonPromptTight", &isGlobalMuonPromptTight);
  
  recoMuons->SetBranchAddress("gmrPt"      , &ptReco );
  recoMuons->SetBranchAddress("gmrEta"     , &etaReco);
  recoMuons->SetBranchAddress("gmrPhi"     , &phiReco);
  recoMuons->SetBranchAddress("gmrChi2Norm", &gmrChi2Norm);
  recoMuons->SetBranchAddress("gmrD0"      ,  &gmrD0);
  recoMuons->SetBranchAddress("gmrDz"      ,  &gmrDz);

  recoMuons->SetBranchAddress("stdEta", &stdEta);
  recoMuons->SetBranchAddress("stdPt" , &stdPt);
  recoMuons->SetBranchAddress("trkEta", &trkEta);
  recoMuons->SetBranchAddress("trkPhi", &trkPhi);
  recoMuons->SetBranchAddress("trkPt",  &trkPt);
  recoMuons->SetBranchAddress("trkValHits",  &trkValHits);
  recoMuons->SetBranchAddress("trkChi2Norm", &trkChi2Norm);
  recoMuons->SetBranchAddress("trkD0", &trkD0);
  recoMuons->SetBranchAddress("trkDz", &trkDz);
   
  recoMuons->SetBranchAddress("rchEta", &rchEta);
  recoMuons->SetBranchAddress("rchPhi", &rchPhi);
  recoMuons->SetBranchAddress("rchCSCtype", &rchCSCtype);
   
  recoMuons->SetBranchAddress("trkNSegs"           , &trkNSegs          );
  recoMuons->SetBranchAddress("trkSegChamberId"    , trkSegChamberId    );
  recoMuons->SetBranchAddress("trkSegRing"         , trkSegRing         );
  recoMuons->SetBranchAddress("trkSegStation"      , trkSegStation      );
  recoMuons->SetBranchAddress("trkSegEndcap"       , trkSegEndcap       );
  recoMuons->SetBranchAddress("trkSegTriggerSector", trkSegTriggerSector);
  recoMuons->SetBranchAddress("trkSegTriggerCscId" , trkSegTriggerCscId );
  recoMuons->SetBranchAddress("trkSegXfromMatch"   , trkSegXfromMatch   );
  recoMuons->SetBranchAddress("trkSegYfromMatch"   , trkSegYfromMatch   );
  recoMuons->SetBranchAddress("trkSegPhifromMatch" , trkSegPhifromMatch );
   
  //segment
  recoMuons->SetBranchAddress("trkSegIsArb", trkSegIsArb);
  recoMuons->SetBranchAddress("trkSegX"    , trkSegX    );
  recoMuons->SetBranchAddress("trkSegY"    , trkSegY    );
  recoMuons->SetBranchAddress("trkSegR"    , trkSegR    );
  recoMuons->SetBranchAddress("trkSegPhi"  , trkSegPhi  );
  recoMuons->SetBranchAddress("trkSegEta"  , trkSegEta  );
   
  // propagation to ME1/1
  recoMuons->SetBranchAddress(  "muons_x_me11",  &muons_x_me11);
  recoMuons->SetBranchAddress(  "muons_y_me11",  &muons_y_me11);
  recoMuons->SetBranchAddress(  "muons_z_me11",  &muons_z_me11);
  recoMuons->SetBranchAddress("muons_phi_me11",&muons_phi_me11);
  recoMuons->SetBranchAddress("muons_eta_me11",&muons_eta_me11);
   
  recoMuons->SetBranchAddress(  "muons_x_me1",  &muons_x_me1);
  recoMuons->SetBranchAddress(  "muons_y_me1",  &muons_y_me1);
  recoMuons->SetBranchAddress(  "muons_z_me1",  &muons_z_me1);
  recoMuons->SetBranchAddress("muons_phi_me1",&muons_phi_me1);
  recoMuons->SetBranchAddress("muons_eta_me1",&muons_eta_me1);

  //---------------------------------------------------------------------
  // segments belonging to the muon
  //---------------------------------------------------------------------
  recoMuons->SetBranchAddress("muonNsegs",&muonNsegs);

  recoMuons->SetBranchAddress("muon_cscsegs_loc_x"      , muon_cscsegs_loc_x      );
  recoMuons->SetBranchAddress("muon_cscsegs_loc_y"      , muon_cscsegs_loc_y      );
  recoMuons->SetBranchAddress("muon_cscsegs_loc_eta"    , muon_cscsegs_loc_eta    );
  recoMuons->SetBranchAddress("muon_cscsegs_loc_phi"    , muon_cscsegs_loc_phi    );
  recoMuons->SetBranchAddress("muon_cscsegs_loc_dir_eta", muon_cscsegs_loc_dir_eta);
  recoMuons->SetBranchAddress("muon_cscsegs_loc_dir_phi", muon_cscsegs_loc_dir_phi);

  recoMuons->SetBranchAddress("muon_cscsegs_gbl_x"      , muon_cscsegs_gbl_x      );
  recoMuons->SetBranchAddress("muon_cscsegs_gbl_y"      , muon_cscsegs_gbl_y      );
  recoMuons->SetBranchAddress("muon_cscsegs_gbl_eta"    , muon_cscsegs_gbl_eta    );
  recoMuons->SetBranchAddress("muon_cscsegs_gbl_phi"    , muon_cscsegs_gbl_phi    );
  recoMuons->SetBranchAddress("muon_cscsegs_gbl_dir_eta", muon_cscsegs_gbl_dir_eta);
  recoMuons->SetBranchAddress("muon_cscsegs_gbl_dir_phi", muon_cscsegs_gbl_dir_phi);

  recoMuons->SetBranchAddress("muon_cscsegs_dxdz"   , muon_cscsegs_dxdz   );
  recoMuons->SetBranchAddress("muon_cscsegs_dydz"   , muon_cscsegs_dydz   );
  recoMuons->SetBranchAddress("muon_cscsegs_dxdzErr", muon_cscsegs_dxdzErr);
  recoMuons->SetBranchAddress("muon_cscsegs_dydzErr", muon_cscsegs_dydzErr);

  recoMuons->SetBranchAddress("muon_cscsegs_endcap" , muon_cscsegs_endcap );
  recoMuons->SetBranchAddress("muon_cscsegs_station", muon_cscsegs_station);
  recoMuons->SetBranchAddress("muon_cscsegs_ring"   , muon_cscsegs_ring   );
  recoMuons->SetBranchAddress("muon_cscsegs_chamber", muon_cscsegs_chamber);
  recoMuons->SetBranchAddress("muon_cscsegs_nhits"  , muon_cscsegs_nhits  );
  
  recoMuons->SetBranchAddress("muon_cscsegs_islctable", muon_cscsegs_islctable);
  recoMuons->SetBranchAddress("muon_cscsegs_ismatched", muon_cscsegs_ismatched);
  recoMuons->SetBranchAddress("muon_cscsegs_lctId"    , muon_cscsegs_lctId    );
  recoMuons->SetBranchAddress("muon_cscsegs_nmatched" , muon_cscsegs_nmatched );
  //---------------------------------------------------------------------

  int SizeTrk;
  std::vector<int>*    EndcapTrk = new vector<int>();  
  std::vector<int>*    SectorTrk = new vector<int>(); 
  std::vector<int>*    BxTrk     = new vector<int>();  

  std::vector<int>*    me1ID = new vector<int>(); 
  std::vector<int>*    me2ID = new vector<int>(); 
  std::vector<int>*    me3ID = new vector<int>(); 
  std::vector<int>*    me4ID = new vector<int>(); 
  std::vector<int>*    mb1ID = new vector<int>();     

  std::vector<int>*    OutputLinkTrk  = new vector<int>();  

  std::vector<int>*  ModeTrk = new vector<int>();  
  std::vector<float>* EtaTrk = new vector<float>();   
  std::vector<float>* PhiTrk = new vector<float>();
  std::vector<float>*  PtTrk = new vector<float>();

  std::vector<int>*      ChargeTrk = new vector<int>();  
  std::vector<int>* ChargeValidTrk = new vector<int>();  
  std::vector<int>*     QualityTrk = new vector<int>();  
  std::vector<int>*        ForRTrk = new vector<int>();  
  std::vector<int>*       Phi23Trk = new vector<int>();  
  std::vector<int>*       Phi12Trk = new vector<int>();    
  std::vector<int>*     PhiSignTrk = new vector<int>();    
                    
  std::vector<int>* EtaBitTrk = new vector<int>();    
  std::vector<int>* PhiBitTrk = new vector<int>();    
  std::vector<int>*  PtBitTrk = new vector<int>();    
 
  std::vector<int>* NumLCTsTrk = new vector<int>();  
  int       trLctEndcap[MAX_CSCTF_TRK][MAX_LCTS_PER_TRK]; 
  int       trLctSector[MAX_CSCTF_TRK][MAX_LCTS_PER_TRK]; 
  int    trLctSubSector[MAX_CSCTF_TRK][MAX_LCTS_PER_TRK]; 
  int           trLctBx[MAX_CSCTF_TRK][MAX_LCTS_PER_TRK]; 
  int          trLctBx0[MAX_CSCTF_TRK][MAX_LCTS_PER_TRK]; 
  int      trLctStation[MAX_CSCTF_TRK][MAX_LCTS_PER_TRK]; 
  int         trLctRing[MAX_CSCTF_TRK][MAX_LCTS_PER_TRK]; 
  int      trLctChamber[MAX_CSCTF_TRK][MAX_LCTS_PER_TRK]; 
  int trLctTriggerCSCID[MAX_CSCTF_TRK][MAX_LCTS_PER_TRK]; 
  int         trLctFpga[MAX_CSCTF_TRK][MAX_LCTS_PER_TRK];	  
  // note: the SPs return them in bits 
  int  trLctlocalPhi[MAX_CSCTF_TRK][MAX_LCTS_PER_TRK]; 
  int trLctglobalPhi[MAX_CSCTF_TRK][MAX_LCTS_PER_TRK];   
  int trLctglobalEta[MAX_CSCTF_TRK][MAX_LCTS_PER_TRK]; 
  int  trLctstripNum[MAX_CSCTF_TRK][MAX_LCTS_PER_TRK];   
  int trLctwireGroup[MAX_CSCTF_TRK][MAX_LCTS_PER_TRK];     
 
  csctfTTree->SetBranchAddress("SizeTrk"       , &SizeTrk       );
  csctfTTree->SetBranchAddress("EndcapTrk"     , &EndcapTrk     );
  csctfTTree->SetBranchAddress("SectorTrk"     , &SectorTrk     );
  csctfTTree->SetBranchAddress("BxTrk"         , &BxTrk         );
  csctfTTree->SetBranchAddress("me1ID"         , &me1ID         );
  csctfTTree->SetBranchAddress("me2ID"         , &me2ID         );
  csctfTTree->SetBranchAddress("me3ID"         , &me3ID         );
  csctfTTree->SetBranchAddress("me4ID"         , &me4ID         );
  csctfTTree->SetBranchAddress("mb1ID"         , &mb1ID         );
  csctfTTree->SetBranchAddress("OutputLinkTrk" , &OutputLinkTrk );
  csctfTTree->SetBranchAddress("ModeTrk"       , &ModeTrk       );
  csctfTTree->SetBranchAddress("EtaTrk"        , &EtaTrk        );
  csctfTTree->SetBranchAddress("PhiTrk_02PI"   , &PhiTrk        );
  csctfTTree->SetBranchAddress("PtTrk"         , &PtTrk         );
  csctfTTree->SetBranchAddress("ChargeTrk"     , &ChargeTrk     );
  csctfTTree->SetBranchAddress("ChargeValidTrk", &ChargeValidTrk);
  csctfTTree->SetBranchAddress("QualityTrk"    , &QualityTrk    );
  csctfTTree->SetBranchAddress("ForRTrk"       , &ForRTrk       );
  csctfTTree->SetBranchAddress("Phi23Trk"      , &Phi23Trk      );
  csctfTTree->SetBranchAddress("Phi12Trk"      , &Phi12Trk      );
  csctfTTree->SetBranchAddress("PhiSignTrk"    , &PhiSignTrk    );
  csctfTTree->SetBranchAddress("EtaBitTrk"     , &EtaBitTrk     );
  csctfTTree->SetBranchAddress("PhiBitTrk"     , &PhiBitTrk     );
  csctfTTree->SetBranchAddress("PtBitTrk"      , &PtBitTrk      );


  csctfTTree->SetBranchAddress("NumLCTsTrk"       , &NumLCTsTrk       );
  csctfTTree->SetBranchAddress("trLctEndcap"      ,  trLctEndcap      );
  csctfTTree->SetBranchAddress("trLctSector"      ,  trLctSector      );
  csctfTTree->SetBranchAddress("trLctSubSector"   ,  trLctSubSector   );
  csctfTTree->SetBranchAddress("trLctBx"          ,  trLctBx          );
  csctfTTree->SetBranchAddress("trLctBx0"         ,  trLctBx0         );
  csctfTTree->SetBranchAddress("trLctStation"     ,  trLctStation     );
  csctfTTree->SetBranchAddress("trLctRing"        ,  trLctRing        );
  csctfTTree->SetBranchAddress("trLctChamber"     ,  trLctChamber     );
  csctfTTree->SetBranchAddress("trLctTriggerCSCID",  trLctTriggerCSCID);
  csctfTTree->SetBranchAddress("trLctFpga"        ,  trLctFpga        );
  csctfTTree->SetBranchAddress("trLctlocalPhi"    ,  trLctlocalPhi    );
  csctfTTree->SetBranchAddress("trLctglobalPhi"   ,  trLctglobalPhi   );
  csctfTTree->SetBranchAddress("trLctglobalEta"   ,  trLctglobalEta   );
  csctfTTree->SetBranchAddress("trLctstripNum"    ,  trLctstripNum    );
  csctfTTree->SetBranchAddress("trLctwireGroup"   ,  trLctwireGroup   );

  // --------------------------------------------------------------------------- 
  int SizeLCTs;
  std::vector<int>*       lctEndcap = new vector<int>; 
  std::vector<int>*       lctSector = new vector<int>; 
  std::vector<int>*    lctSubSector = new vector<int>; 
  std::vector<int>*           lctBx = new vector<int>; 
  std::vector<int>*          lctBx0 = new vector<int>; 
  std::vector<int>*      lctStation = new vector<int>; 
  std::vector<int>*         lctRing = new vector<int>; 
  std::vector<int>*      lctChamber = new vector<int>; 
  std::vector<int>* lctTriggerCSCID = new vector<int>; 
  std::vector<int>*         lctFpga = new vector<int>;     
  
  // note: the SPs return them in bits 
  std::vector<int>*  lctlocalPhi = new vector<int>; 
  std::vector<int>* lctglobalPhi = new vector<int>;   
  std::vector<int>* lctglobalEta = new vector<int>; 
  std::vector<int>*  lctstripNum = new vector<int>;   
  std::vector<int>* lctwireGroup = new vector<int>;   
  
  csctfTTree->SetBranchAddress("SizeLCTs"       , &SizeLCTs       );
  csctfTTree->SetBranchAddress("lctEndcap"      , &lctEndcap      );
  csctfTTree->SetBranchAddress("lctSector"      , &lctSector      );
  csctfTTree->SetBranchAddress("lctSubSector"   , &lctSubSector   );
  csctfTTree->SetBranchAddress("lctBx"          , &lctBx          );
  csctfTTree->SetBranchAddress("lctBx0"         , &lctBx0         );
  csctfTTree->SetBranchAddress("lctStation"     , &lctStation     );
  csctfTTree->SetBranchAddress("lctRing"        , &lctRing        );
  csctfTTree->SetBranchAddress("lctChamber"     , &lctChamber     );
  csctfTTree->SetBranchAddress("lctTriggerCSCID", &lctTriggerCSCID);
  csctfTTree->SetBranchAddress("lctFpga"        , &lctFpga        );
  csctfTTree->SetBranchAddress("lctlocalPhi"    , &lctlocalPhi    );
  csctfTTree->SetBranchAddress("lctglobalPhi"   , &lctglobalPhi   );
  csctfTTree->SetBranchAddress("lctglobalEta"   , &lctglobalEta   );
  csctfTTree->SetBranchAddress("lctstripNum"    , &lctstripNum    );
  csctfTTree->SetBranchAddress("lctwireGroup"   , &lctwireGroup   );

  // ------------------------------------------------------------------
  // counters for diagnosys
  // ------------------------------------------------------------------
  int nMuons              = 0; // total # muons

  // # muons with no trigger info
  int nMuonsNoTrigger     = 0; // total # muons
  int nMuonsNoTriggerTime = 0; // # muons with no trigger info b/c of lcts timing
  int nMuonsNoTriggerCuts = 0; // # muons with no trigger info b/c of lcts extrapolation cuts
  int nMuonsNoTriggerEdge = 0; // # muons with no trigger info b/c close to the edge (Alex's snippet)
  int nMuonsNoTriggerElse = 0; // # muons with no trigger info b/c of other firmware nuances -> to skim and send to Alex

  // # muons with trigger info but not coming from the LCTs belonging to the global muon
  int nMuonsYesTriggerTime = 0; // lcts out-of-time
  int nMuonsYesTriggerCuts = 0; // the LCTs do not pass extrapolation cuts
  int nMuonsYesTriggerEdge = 0; // LCTs close to the edge (Alex's snippet)


  // ------------------------------------------------------------------
  // Loop over the events
  // ------------------------------------------------------------------
  //for (int iEvt=0; iEvt < recoMuons->GetEntries(); iEvt++) {
  //for (int iEvt=0; iEvt < 30000; iEvt++) {
  //for (int iEvt=400000; iEvt < 500000; iEvt++) {
  for (int iEvt=433764; iEvt < 433765; iEvt++) {

    recoMuons ->GetEntry(iEvt);
    csctfTTree->GetEntry(iEvt);
  
    //if ( ( iEvt % 10000) == 0 )
    if (printLevel > 0 ) {
      printf(" --- Event # %6d \n", iEvt+1);
      printf(" --- muonSize # %d \n", muonSize);
    } else {
      if ( ( iEvt % 10000) == 0 ) printf(" --- Event # %6d \n", iEvt+1);
    }

    for (int iReco=0; iReco < muonSize; iReco++) { 

      // this helps the id the muon on fireworks
      if (printLevel > 0 ) {
        std::cout << "muon="   << iReco+1 << std::endl;
        std::cout << "trkPt="  << trkPt->at(iReco) << std::endl;
        std::cout << "gmrEta=" << etaReco->at(iReco) << std::endl;
        std::cout << "gmrPhi=" << phiReco->at(iReco) << std::endl;
      }

      // -----------------------------------------------------------
      // global muon definition
      bool isGblMuon = true;
      if (!isGlobalMuon->at(iReco))            isGblMuon = false;
      if (!isGlobalMuonPromptTight->at(iReco)) isGblMuon = false;
      if (!isStandAloneMuon->at(iReco))        isGblMuon = false;
      if (!isTrackerMuon->at(iReco))           isGblMuon = false;
      if (gmrChi2Norm->at(iReco) > 10)         isGblMuon = false;
      if (fabs(gmrD0->at(iReco)) >  2)         isGblMuon = false;
      if ((rchEta->at(iReco)==-999))           isGblMuon = false;

      // standalone ONLY muon definition
      bool isStandAloneMuonOnly = true;
      if (isGlobalMuon->at(iReco))      isStandAloneMuonOnly = false;
      if (isTrackerMuon->at(iReco))     isStandAloneMuonOnly = false;
      if (!isStandAloneMuon->at(iReco)) isStandAloneMuonOnly = false;
      if ((rchEta->at(iReco)==-999))    isStandAloneMuonOnly = false;

      // tracker ONLY muon definition
      bool isTrackerMuonOnly = true;
      if (isGlobalMuon->at(iReco))           isTrackerMuonOnly=false;
      if (isStandAloneMuon->at(iReco))       isTrackerMuonOnly=false;
      if (muons_eta_me11->at(iReco) == -999) isTrackerMuonOnly=false;
      if (fabs(trkEta->at(iReco)) < 0.9)     isTrackerMuonOnly=false;
      if (trkNSegs->at(iReco) < 1)           isTrackerMuonOnly=false;
      // quarkonia tracker muons selections
      //if (trkValHits->at(iReco) < 12)        isTrackerMuonOnly=false;
      if (!isTMLastStationAngTight->at(iReco)) isTrackerMuonOnly=false;
      if (trkChi2Norm->at(iReco) > 5.0)        isTrackerMuonOnly=false;
      if (fabs(trkD0->at(iReco)) > 5.0)        isTrackerMuonOnly=false;
      if (fabs(trkDz->at(iReco)) > 20)         isTrackerMuonOnly=false;
    
      // Std + trk Mu != gbl
      bool isStdAndTrkButNotGbl = true;
      if (isGlobalMuon->at(iReco))       isStdAndTrkButNotGbl=false;
      if (!isStandAloneMuon->at(iReco))  isStdAndTrkButNotGbl=false;
      if (!isTrackerMuon->at(iReco))     isStdAndTrkButNotGbl=false;
      if ((rchEta->at(iReco)==-999))     isStdAndTrkButNotGbl=false;
    
      int counter=0;
      if (isGblMuon)            {counter++;}
      if (isStandAloneMuonOnly) {counter++;}
      if (isTrackerMuonOnly)    {counter++;}
      if (isStdAndTrkButNotGbl) {counter++;}
    
      // sanity check
      if (counter>1) {
        cout << "Error: the same muon fills two category!"
             << " Review their definitions\n";
        continue;
      }
      // --------------------------------------------------------------

      // --------------------------------------------------------------
      // pick your selection
      // --------------------------------------------------------------
      // here you select global, standalone only, tracker only
      if (!isGblMuon) continue;

      if (printLevel>0) {
        std::cout << "!!! is a Global Muon !!!\n";
        std::cout << "\n-------------------------------------------\n";
      }

      // --------------------------------------------------------------
      // routine to identify if the global muon has two segments
      // --------------------------------------------------------------
      bool hasTwoSegsMatched=false;

      int  lctEtaBitSeg[MAX_SEGS_STD];
      int  lctPhiBitSeg[MAX_SEGS_STD];
      int  lctEndcapSeg[MAX_SEGS_STD];
      int lctStationSeg[MAX_SEGS_STD];
      int  lctSectorSeg[MAX_SEGS_STD];
      int      lctBxSeg[MAX_SEGS_STD];

      int counterSegs=0;

      std::vector<csctftrack> csctftrack_vec;

      for (int iSeg=0;iSeg<muonNsegs->at(iReco);iSeg++) {
         
        if (muon_cscsegs_ismatched[iReco][iSeg]) {
          int id = muon_cscsegs_lctId[iReco][iSeg];
          lctEtaBitSeg[counterSegs]  = lctglobalEta->at(id); 
          lctPhiBitSeg[counterSegs]  = lctglobalPhi->at(id); 
          lctEndcapSeg[counterSegs]  = lctEndcap->at(id); 
          lctStationSeg[counterSegs] = lctStation->at(id); 
          lctSectorSeg[counterSegs]  = lctSector->at(id); 
          lctBxSeg[counterSegs]      = lctBx->at(id);

          if (printLevel>0) 
            std::cout << "id=" << id 
                      << ", lctEtaBitSeg="  << lctEtaBitSeg[counterSegs]
                      << ", lctPhiBitSeg="  << lctPhiBitSeg[counterSegs]
                      << ", lctStationSeg=" << lctStation->at(id)
                      << ", lctRingSeg="    << lctRing->at(id)
                      << ", lctSectorSeg="  << lctSector->at(id)
                      << ", lctEndcapSeg="  << lctEndcap->at(id)
                      << ", lctBxSeg="      << lctBx->at(id)
                      << std::endl;
           
          counterSegs++;

          // ----------------------------------------------------------------
          // identify if this information was already in a previous track
          bool IsAlreadyFilled = false;
          std::vector<csctftrack>::iterator csctftrack_it;

          for ( csctftrack_it = csctftrack_vec.begin(); csctftrack_it!= csctftrack_vec.end(); csctftrack_it++){
            
            // if endcap and sector matches
            if( csctftrack_it->endcap == lctEndcap->at(id) && 
                csctftrack_it->sector == lctSector->at(id)  ) {
              
              // does it match an already present LCT?
              bool isMatchingPreviousLCT = false;
              
              for (int iLCT=0; iLCT <4; iLCT++) {
                if(csctftrack_it->lctEta[iLCT] == -999) continue;
                if(csctftrack_it->lctPhi[iLCT] == -999) continue;
                
                int diffEta = abs (csctftrack_it->lctEta[iLCT] - lctglobalEta->at(id));
                int diffPhi = abs (csctftrack_it->lctPhi[iLCT] - lctglobalPhi->at(id));
                
                if (diffEta <= 6 && diffPhi <= 1024) isMatchingPreviousLCT = true;
              }
              
              if (!isMatchingPreviousLCT) continue;
              
              if (printLevel>0) {
                std::cout << "The candidate track is already there...\n";
                std::cout << "Adding a component\n";
              }
              
              int lctPos = lctStation->at(id)-1;
              
              csctftrack_it->lctEta[lctPos] = lctglobalEta->at(id);
              csctftrack_it->lctPhi[lctPos] = lctglobalPhi->at(id);
              
              IsAlreadyFilled = true;
            }
          }
          
          if (!IsAlreadyFilled) {
            if (printLevel > 0) 
              std::cout << "track not present: forming one" << std::endl;
           
            int lctPos = lctStation->at(id)-1;
            
            csctftrack track;
            track.endcap = lctEndcap->at(id);
            track.sector = lctSector->at(id);
            track.lctEta[lctPos] = lctglobalEta->at(id);
            track.lctPhi[lctPos] = lctglobalPhi->at(id);

            csctftrack_vec.push_back( track );
 
          }

          
        }
      }
      

      if (counterSegs>1) hasTwoSegsMatched=true;
       
      if (printLevel>0) 
        std::cout << " -> hasTwoSegsMatched?=" << hasTwoSegsMatched
                  << std::endl;
       
      // only Global muons with at least two segments
      if (!hasTwoSegsMatched) continue;


      // muons counter
      nMuons+=1;

      bool isTriggered = false;
      int    modeWinner=+999;      


      // loop over segments
      for (int iSeg=0;iSeg<muonNsegs->at(iReco);iSeg++) {
         
        if (isTriggered) continue;

        // only matched segments
        if (!muon_cscsegs_ismatched[iReco][iSeg]) continue;
        int segLCTid = muon_cscsegs_lctId[iReco][iSeg];

        if (printLevel>0) {
          // very verbose debugging
          std::cout  << "*********************************************************************\n";
          std::cout  << "lctEndcap->at(segLCTid)   =" << lctEndcap->at(segLCTid)   << std::endl
                     << "lctSector->at(segLCTid)   =" << lctSector->at(segLCTid)   << std::endl
                     << "lctStation->at(segLCTid)  =" << lctStation->at(segLCTid)  << std::endl
                     << "lctRing->at(segLCTid)     =" << lctRing->at(segLCTid)     << std::endl
                     << "lctChamber->at(segLCTid)  =" << lctChamber->at(segLCTid)  << std::endl
                     << "lctglobalPhi->at(segLCTid)=" << lctglobalPhi->at(segLCTid)<< std::endl
                     << "lctglobalEta->at(segLCTid)=" << lctglobalEta->at(segLCTid)<< std::endl;
        }

        //loop over the CSCTF raw information
        for (int iRaw=0; iRaw<SizeTrk; iRaw++) {

          //trigCounter++;
              
          //loop over  the LCT belonging to the track
          int nLcts = NumLCTsTrk->at(iRaw);

          for (int iLCT=0; iLCT<nLcts; iLCT++) {

         //    if (printLevel>0) {
//               std::cout  << "iLCT=" << iLCT << endl;
//               std::cout  << "trLctEndcap[" << iRaw << "][" << iLCT << "]   =" << trLctEndcap[iRaw][iLCT]   << std::endl
//                          << "trLctSector[" << iRaw << "][" << iLCT << "]   =" << trLctSector[iRaw][iLCT]   << std::endl
//                          << "trLctStation[" << iRaw << "][" << iLCT << "]  =" << trLctStation[iRaw][iLCT]  << std::endl
//                          << "trLctRing[" << iRaw << "][" << iLCT << "]     =" << trLctRing[iRaw][iLCT]     << std::endl
//                          << "trLctChamber[" << iRaw << "][" << iLCT << "]  =" << trLctChamber[iRaw][iLCT]  << std::endl
//                          << "trLctglobalPhi[" << iRaw << "][" << iLCT << "]=" << trLctglobalPhi[iRaw][iLCT] << std::endl
//                          << "trLctglobalEta[" << iRaw << "][" << iLCT << "]=" << trLctglobalEta[iRaw][iLCT]<< std::endl;
//               std::cout << "-----------------------------------------------------------\n";
//             }
            
            if (trLctEndcap[iRaw][iLCT]    != lctEndcap->at(segLCTid))    continue;
            if (trLctSector[iRaw][iLCT]    != lctSector->at(segLCTid))    continue;
            if (trLctStation[iRaw][iLCT]   != lctStation->at(segLCTid))   continue;
            if (trLctRing[iRaw][iLCT]      != lctRing->at(segLCTid))      continue;
            if (trLctChamber[iRaw][iLCT]   != lctChamber->at(segLCTid))   continue;
            if (trLctglobalPhi[iRaw][iLCT] != lctglobalPhi->at(segLCTid)) continue;
            if (trLctglobalEta[iRaw][iLCT] != lctglobalEta->at(segLCTid)) continue;
      
            
            if (printLevel>0) std::cout << " ----> IS TRIGGERD!\n";      
     

            isTriggered = true;
            modeWinner=ModeTrk->at(iRaw);      

          }
        }
      }
      

      // analysis of problematic triggers
      if (!isTriggered) {
        
        nMuonsNoTrigger+=1;
        
        if (printLevel>0) {
          std::cout << "\n\nHouston we have a problem" << std::endl;
          std::cout << "How many triggers?: " << SizeTrk << std::endl;           
          std::cout << "RUN=" << Run   << std::endl;
          std::cout << "EVT=" << Event << std::endl;
          std::cout << "LUMI=" << Lumi << std::endl;
        }
        
        // --------------------------------------------------------------
        bool atLeastTwoLCTinWindows=false;
        bool BxTimed=false;

        for (int i=0;i<counterSegs;i++) {
          
          for (int j=i;j<counterSegs;j++) {
            if ( abs( lctEtaBitSeg[i]-lctEtaBitSeg[j] ) <=    6 &&
                 abs( lctPhiBitSeg[i]-lctPhiBitSeg[j] ) <= 1024 &&
                 lctEndcapSeg[i] == lctEndcapSeg[j]             &&
                 lctStationSeg[i] != lctStationSeg[j]           &&
                 lctSectorSeg[i] == lctSectorSeg[j]              ) {

              atLeastTwoLCTinWindows=true;

              if (abs(lctBxSeg[i]-lctBxSeg[j]) < 3) BxTimed=true;
            }
          }
        }
        
          
        if (printLevel>0) 
          std::cout << " * atLeastTwoLCTinWindows?=" 
                    << atLeastTwoLCTinWindows << std::endl;
        

        bool stop = false;

        if ( !atLeastTwoLCTinWindows ) {
          if (SizeTrk==0) nMuonsNoTriggerCuts +=1;
          if (SizeTrk!=0) nMuonsYesTriggerCuts +=1;
        }
        else if ( !BxTimed ) {
          if (SizeTrk==0) nMuonsNoTriggerTime +=1;
          if (SizeTrk!=0) nMuonsYesTriggerTime +=1;
        }
        else{
          
          bool closeToEdge = false;
          std::vector<csctftrack>::iterator csctftrack_it;
          for ( csctftrack_it = csctftrack_vec.begin(); csctftrack_it!= csctftrack_vec.end(); csctftrack_it++){
            
            if (printLevel > 0) {
              std::cout << "**************************************************\n";
              std::cout << "csctftrack_it->endcap=" << csctftrack_it->endcap << std::endl;
              std::cout << "csctftrack_it->sector=" << csctftrack_it->sector << std::endl;
              
              for (int i=0;i<4;i++) {
                std::cout << "csctftrack_it->lctEta["<< i <<"]=" << csctftrack_it->lctEta[i] << std::endl;
                std::cout << "csctftrack_it->lctPhi["<< i <<"]=" << csctftrack_it->lctPhi[i] << std::endl;
              }
              
              std::cout << "csctftrack_it->mode=" << csctftrack_it->mode() << std::endl;
              std::cout << "csctftrack_it->phi="  << csctftrack_it->phi() << std::endl;
              std::cout << "IsCloseToEdge="  << IsCloseToEdge((*csctftrack_it)) << std::endl;
            }
            
            if ( closeToEdge ) continue;
            else closeToEdge = IsCloseToEdge((*csctftrack_it));
            
            
            
          }
          
          if (closeToEdge) {
            if (SizeTrk==0) nMuonsNoTriggerEdge  +=1;
            if (SizeTrk!=0) nMuonsYesTriggerEdge +=1;
          }
          else             {
            nMuonsNoTriggerElse +=1;   
            stop =true;
          }
        }
        
        if (/*printLevel>0 &&*/ stop) {
          
          // All CSCTF tracks information 
          for (int iRaw=0; iRaw<SizeTrk; iRaw++) 
            std::cout << " EndcapTrk->at(iRaw)=" << EndcapTrk->at(iRaw)
                      << " SectorTrk->at(iRaw)=" << SectorTrk->at(iRaw)
                      << " me1ID=" << me1ID->at(iRaw)
                      << " me2ID=" << me2ID->at(iRaw)
                      << " me3ID=" << me3ID->at(iRaw)
                      << " me4ID=" << me4ID->at(iRaw)
                      << " mb1ID=" << mb1ID->at(iRaw)
                      << " ModeTrk=" << ModeTrk->at(iRaw)
                      << std::endl;
          
          // All LCTs information
          for (int iLCT=0; iLCT<SizeLCTs; iLCT++) 


            std::cout << "    lctEndcap = " <<    lctEndcap->at(iLCT)
                      << "    lctSector = " <<    lctSector->at(iLCT)
                      << "   lctStation = " <<   lctStation->at(iLCT)
                      << "      lctRing = " <<      lctRing->at(iLCT)
                      << "   lctChamber = " <<   lctChamber->at(iLCT)
                      << "        lctBx = " <<        lctBx->at(iLCT)
                      << " lctglobalPhi = " << lctglobalPhi->at(iLCT)
                      << " lctglobalEta = " << lctglobalEta->at(iLCT) 
                      << std::endl;
          
          // stop to look at the event printout
          //int pippo;
          //cin >> pippo;
        }
      }
      
    }
  }// end Loop evts
  
  if (printLevel > 0)
    printf("\n----------------------------------------"
 	   "----------------------------------------\n");

  // ---------------------------------------------------------------------------
  // printout
  // ---------------------------------------------------------------------------
  std::cout << "nMuons              = " << nMuons              << std::endl;
  std::cout << "nMuonsNoTrigger     = " << nMuonsNoTrigger     << std::endl;
  std::cout << "nMuonsNoTriggerTime = " << nMuonsNoTriggerTime << std::endl;
  std::cout << "nMuonsNoTriggerCuts = " << nMuonsNoTriggerCuts << std::endl;
  std::cout << "nMuonsNoTriggerEdge = " << nMuonsNoTriggerEdge << std::endl;
  std::cout << "nMuonsNoTriggerElse = " << nMuonsNoTriggerElse << std::endl;
  
  std::cout << std::endl;

  std::cout << "nMuonsYesTriggerTime = " << nMuonsYesTriggerTime << std::endl;
  std::cout << "nMuonsYesTriggerCuts = " << nMuonsYesTriggerCuts << std::endl;
  std::cout << "nMuonsYesTriggerEdge = " << nMuonsYesTriggerEdge << std::endl;

}

double DR(double diffeta,
          double diffphi) {

  double diffetaSquare =diffeta*diffeta;  
  double diffphiSquare =diffphi*diffphi;  

  return sqrt(diffetaSquare+diffphiSquare);
}


bool IsCloseToEdge(csctftrack track) {

  bool retVal = false;

  int PHICUTL=128;
  int PHICUTH=(4095 - PHICUTL);
  
  int mode = track.mode();
  int phi = track.phi();

  if(
     (mode == 5 || mode == 8 || mode == 9 || mode == 10) &&
     (phi < PHICUTL || phi > PHICUTH)
     )
    retVal = true;

  return retVal;
}
