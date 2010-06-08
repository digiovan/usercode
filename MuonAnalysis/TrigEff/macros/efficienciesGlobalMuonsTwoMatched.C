//.L /afs/cern.ch/user/d/digiovan/code/CMSSW_3_5_4/tmp/slc5_ia32_gcc434/src/PhysicsTools/RooStatsCms/src/PhysicsToolsRooStatsCms/libPhysicsToolsRooStatsCms.so
//gSystem->CompileMacro("scripts/Collisions2010/trigEff/efficienciesGlobalMuonsTwoMatched.C","Ok")

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

using namespace std;

#include "/afs/cern.ch/user/d/digiovan/scripts/Init/modifiedStyle.C"

double DR(double diffeta,
          double diffphi);

double computeError(double num, double den);

TChain* recoMuons;
TChain* csctfTTree;

// ------------------------------------------------------------------------------
// which kind of muon?
// ------------------------------------------------------------------------------
TH1F* hTypeMu;
TH1F* hTriggerCounts;

// ------------------------------------------------------------------------------
// Phi/Eta Resolutions
// ------------------------------------------------------------------------------
TH1F* hDeltaPhi;
TH1F* hDeltaEta;
TH1F* hDR;
TH1F* hMode;

// ------------------------------------------------------------------------------
// Histograms for efficiencies
// ------------------------------------------------------------------------------
TH1F* hNTriggeredMuonsVsPt;  
TH1F* hNMuonsVsPt;  

TH1F* hNTriggeredMuonsVsEta;  
TH1F* hNMuonsVsEta;  

TH1F* hNTriggeredMuonsVsPhi;  
TH1F* hNMuonsVsPhi;  

// ------------------------------------------------------------------------------
// Efficiencies Graphs
// ------------------------------------------------------------------------------
TGraphAsymmErrors* gEffvsPt;
TGraphAsymmErrors* gEffvsEta;
TGraphAsymmErrors* gEffvsPhi;

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

void efficienciesGlobalMuonsTwoMatched(int  printLevel =     0,
                                           bool isMC   = !true,
                                           bool isSave = !true) {
    
  gROOT->Clear();
  gStyle->SetOptStat(111111);  
  
  //--------------------------------------------------------------------------
  // Construct the chain of input files
  //--------------------------------------------------------------------------
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
  
  int trkSegChamberId[MAX_MUONS][MAX_TRK_SEGS]; 
  int trkSegRing[MAX_MUONS][MAX_TRK_SEGS];    
  int trkSegStation[MAX_MUONS][MAX_TRK_SEGS]; 
  int trkSegEndcap[MAX_MUONS][MAX_TRK_SEGS];  
  int trkSegTriggerSector[MAX_MUONS][MAX_TRK_SEGS];
  int trkSegTriggerCscId[MAX_MUONS][MAX_TRK_SEGS]; 
  float trkSegXfromMatch[MAX_MUONS][MAX_TRK_SEGS];
  float trkSegYfromMatch[MAX_MUONS][MAX_TRK_SEGS];
  float trkSegPhifromMatch[MAX_MUONS][MAX_TRK_SEGS];

  int trkSegIsArb[MAX_MUONS][MAX_TRK_SEGS];
  float trkSegX[MAX_MUONS][MAX_TRK_SEGS];
  float trkSegY[MAX_MUONS][MAX_TRK_SEGS];
  float trkSegR[MAX_MUONS][MAX_TRK_SEGS];
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
  float muon_cscsegs_loc_x[MAX_MUONS][MAX_SEGS_STD];
  float muon_cscsegs_loc_y[MAX_MUONS][MAX_SEGS_STD];
  float muon_cscsegs_loc_eta[MAX_MUONS][MAX_SEGS_STD];
  float muon_cscsegs_loc_phi[MAX_MUONS][MAX_SEGS_STD];
  float muon_cscsegs_loc_dir_eta[MAX_MUONS][MAX_SEGS_STD];
  float muon_cscsegs_loc_dir_phi[MAX_MUONS][MAX_SEGS_STD];

  // segment position information, global
  float muon_cscsegs_gbl_x[MAX_MUONS][MAX_SEGS_STD];
  float muon_cscsegs_gbl_y[MAX_MUONS][MAX_SEGS_STD];
  float muon_cscsegs_gbl_eta[MAX_MUONS][MAX_SEGS_STD];
  float muon_cscsegs_gbl_phi[MAX_MUONS][MAX_SEGS_STD];
  float muon_cscsegs_gbl_dir_eta[MAX_MUONS][MAX_SEGS_STD];
  float muon_cscsegs_gbl_dir_phi[MAX_MUONS][MAX_SEGS_STD];

  // more on segment direction
  float muon_cscsegs_dxdz[MAX_MUONS][MAX_SEGS_STD];
  float muon_cscsegs_dydz[MAX_MUONS][MAX_SEGS_STD];
  float muon_cscsegs_dxdzErr[MAX_MUONS][MAX_SEGS_STD];
  float muon_cscsegs_dydzErr[MAX_MUONS][MAX_SEGS_STD];

  // general segment information
  int muon_cscsegs_endcap[MAX_MUONS][MAX_SEGS_STD];
  int muon_cscsegs_station[MAX_MUONS][MAX_SEGS_STD];
  int muon_cscsegs_ring[MAX_MUONS][MAX_SEGS_STD];
  int muon_cscsegs_chamber[MAX_MUONS][MAX_SEGS_STD];
  int muon_cscsegs_nhits[MAX_MUONS][MAX_SEGS_STD];

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

  recoMuons->SetBranchAddress("muonSize"           , &muonSize        );
  recoMuons->SetBranchAddress("isGlobalMuon"       , &isGlobalMuon    );
  recoMuons->SetBranchAddress("isStandAloneMuon"   , &isStandAloneMuon);
  recoMuons->SetBranchAddress("isTrackerMuon"      , &isTrackerMuon   );
  recoMuons->SetBranchAddress("isTMLastStationAngTight", &isTMLastStationAngTight   );
  recoMuons->SetBranchAddress("isGlobalMuonPromptTight", &isGlobalMuonPromptTight);
  
  recoMuons->SetBranchAddress("gmrPt" , &ptReco );
  recoMuons->SetBranchAddress("gmrEta", &etaReco);
  recoMuons->SetBranchAddress("gmrPhi", &phiReco);
  recoMuons->SetBranchAddress("gmrChi2Norm", &gmrChi2Norm);
  recoMuons->SetBranchAddress("gmrD0",  &gmrD0);
  recoMuons->SetBranchAddress("gmrDz",  &gmrDz);

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

  std::vector<int>*    ModeTrk = new vector<int>();  
  std::vector<float>*   EtaTrk = new vector<float>(); ;   
  std::vector<float>*   PhiTrk = new vector<float>(); ;   
  std::vector<float>*    PtTrk = new vector<float>(); ;   

  std::vector<int>*    ChargeTrk = new vector<int>();  
  std::vector<int>*    ChargeValidTrk = new vector<int>();  
  std::vector<int>*    QualityTrk = new vector<int>();  
  std::vector<int>*    ForRTrk = new vector<int>();  
  std::vector<int>*    Phi23Trk = new vector<int>();  
  std::vector<int>*    Phi12Trk = new vector<int>();    
  std::vector<int>*    PhiSignTrk = new vector<int>();    

  std::vector<int>*    EtaBitTrk = new vector<int>();    
  std::vector<int>*    PhiBitTrk = new vector<int>();    
  std::vector<int>*    PtBitTrk = new vector<int>();    
 
  std::vector<int>* NumLCTsTrk = new vector<int>();  
  int trLctEndcap[MAX_CSCTF_TRK][MAX_LCTS_PER_TRK]; 
  int trLctSector[MAX_CSCTF_TRK][MAX_LCTS_PER_TRK]; 
  int trLctSubSector[MAX_CSCTF_TRK][MAX_LCTS_PER_TRK]; 
  int trLctBx[MAX_CSCTF_TRK][MAX_LCTS_PER_TRK]; 
  int trLctBx0[MAX_CSCTF_TRK][MAX_LCTS_PER_TRK]; 
  int trLctStation[MAX_CSCTF_TRK][MAX_LCTS_PER_TRK]; 
  int trLctRing[MAX_CSCTF_TRK][MAX_LCTS_PER_TRK]; 
  int trLctChamber[MAX_CSCTF_TRK][MAX_LCTS_PER_TRK]; 
  int trLctTriggerCSCID[MAX_CSCTF_TRK][MAX_LCTS_PER_TRK]; 
  int trLctFpga[MAX_CSCTF_TRK][MAX_LCTS_PER_TRK];	  
  // note: the SPs return them in bits 
  int trLctlocalPhi[MAX_CSCTF_TRK][MAX_LCTS_PER_TRK]; 
  int trLctglobalPhi[MAX_CSCTF_TRK][MAX_LCTS_PER_TRK];   
  int trLctglobalEta[MAX_CSCTF_TRK][MAX_LCTS_PER_TRK]; 
  int trLctstripNum[MAX_CSCTF_TRK][MAX_LCTS_PER_TRK];   
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
  std::vector<int>* lctEndcap = new vector<int>; 
  std::vector<int>* lctSector = new vector<int>; 
  std::vector<int>* lctSubSector = new vector<int>; 
  std::vector<int>* lctBx = new vector<int>; 
  std::vector<int>* lctBx0 = new vector<int>; 
  std::vector<int>* lctStation = new vector<int>; 
  std::vector<int>* lctRing = new vector<int>; 
  std::vector<int>* lctChamber = new vector<int>; 
  std::vector<int>* lctTriggerCSCID = new vector<int>; 
  std::vector<int>* lctFpga = new vector<int>;     
  
  // note: the SPs return them in bits 
  std::vector<int>* lctlocalPhi = new vector<int>; 
  std::vector<int>* lctglobalPhi = new vector<int>;   
  std::vector<int>* lctglobalEta = new vector<int>; 
  std::vector<int>* lctstripNum = new vector<int>;   
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

  // --------------------------------------------------------------------------- 
  // Define Histograms
  // ---------------------------------------------------------------------------
  hTypeMu = new TH1F("hTypeMu","", 4,    0,   4); 
  hTriggerCounts = new TH1F("hTriggerCounts","",10,0,10); 
  
  hDeltaPhi = new TH1F("hDeltaPhi", "", 100,    Pi,   Pi);
  hDeltaEta = new TH1F("hDeltaEta", "", 100,  -1.0,  1.0);
  hDR       = new TH1F("hDR"    , "",   500,     0,  0.5);
  hMode     = new TH1F("hMode",     "",  15,     0,   15);


  // --------------------------------------------------------------------------- 
  // Efficiency Variables Pt
  // --------------------------------------------------------------------------- 
  // variable bin size Pt 
  Double_t scalePt[9] = {0, 2, 3, 4, 8, 15, 25, 35, 39.5};
  hNTriggeredMuonsVsPt = new TH1F("hNTriggeredMuonsVsPt","",8,scalePt);
  hNMuonsVsPt = new TH1F("hNMuonsVsPt","",8,scalePt);


  // --------------------------------------------------------------------------- 
  // Efficiency Variables Eta
  // --------------------------------------------------------------------------- 
  // variable bin size Eta 
  const int etaNBins = 8;
  Double_t scaleEta[etaNBins+1];
  
  double etaMin = 0.9;
  scaleEta[0]=etaMin;

  for (int iEta=1;iEta<etaNBins+1;iEta++)
    scaleEta[iEta] = etaMin + (1.6*iEta/etaNBins); 
  

  hNTriggeredMuonsVsEta = new TH1F("hNTriggeredMuonsVsEta","",etaNBins,scaleEta);
  hNMuonsVsEta = new TH1F("hNMuonsVsEta","",etaNBins,scaleEta);


  // --------------------------------------------------------------------------- 
  // Efficiency Variables Phi
  // --------------------------------------------------------------------------- 
  // variable bin size Phi 
  const int phiNBins = 8;
  Double_t scalePhi[phiNBins+1];
  
  double phiMin = -Pi;
  scalePhi[0]=phiMin;

  for (int iPhi=1;iPhi<phiNBins+1;iPhi++)
    scalePhi[iPhi] = phiMin + (2*Pi*iPhi/phiNBins); 
 

  hNTriggeredMuonsVsPhi = new TH1F("hNTriggeredMuonsVsPhi","",phiNBins,scalePhi);
  hNMuonsVsPhi = new TH1F("hNMuonsVsPhi","",phiNBins,scalePhi);

  
  // --------------------------------------------------------------------------- 
  // All global muons
  TString SaveExt("GblMu");
  TString TitleExt("Gbl #mu");
  
  // counters
  int nMuons              = 0; // total # muons
  int nMuonsNoTrigger     = 0; // # muons with no trigger info
  int nMuonsNoTriggerTime = 0; // # muons with no trigger info b/c of lct timing
  int nMuonsNoTriggerCuts = 0; // # muons with no trigger info b/c of extrapolation cuts
  int nMuonsNoTriggerElse = 0; // # muons with no trigger info b/c of firmware nuances

  // ------------------------------------------------------------------
  // Loop over the events
  // ------------------------------------------------------------------
  //for (int iEvt=0; iEvt < recoMuons->GetEntries(); iEvt++) {
  for (int iEvt=0; iEvt < 30000; iEvt++) {
    //for (int iEvt=0; iEvt < 100; iEvt++) {
    
    recoMuons ->GetEntry(iEvt);
    csctfTTree->GetEntry(iEvt);
  
    if ( ( iEvt % 10000) == 0 ) printf(" --- Event # %6d \n", iEvt+1);

    for (int iReco=0; iReco < muonSize; iReco++) { 

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
      if (isGblMuon)            {hTypeMu->Fill(0); counter++;}
      if (isStandAloneMuonOnly) {hTypeMu->Fill(1); counter++;}
      if (isTrackerMuonOnly)    {hTypeMu->Fill(2); counter++;}
      if (isStdAndTrkButNotGbl) {hTypeMu->Fill(3); counter++;}
    
      // sanity check
      if (counter>1) {
        cout << "Error: the same muon fills two category!"
             << " Review their definitions\n";
        continue;
      }
      // -----------------------------------------------------------

      // --------------------------------------------------------------
      // pick your selection
      // --------------------------------------------------------------
      // here you select global, standalone only, tracker only
      if (!isGblMuon) continue;

      if (printLevel>0) 
        std::cout << "\n-------------------------------------------\n";
       
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

      if (printLevel>0) std::cout << "List of Matched Segments:\n";
            
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
        }
         
      }
      

      if (counterSegs>1) hasTwoSegsMatched=true;
       
      if (printLevel>0) 
        std::cout << " -> hasTwoSegsMatched?=" << hasTwoSegsMatched
                  << std::endl;

       
      // only Global muons with at least two segments matched
      if (!hasTwoSegsMatched) continue;

      
      //fill the denominator    
      nMuons+=1;
      hNMuonsVsPt  -> Fill ( ptReco ->at(iReco) );
      hNMuonsVsEta -> Fill ( etaReco->at(iReco) );
      hNMuonsVsPhi -> Fill ( phiReco->at(iReco) );

    
      bool isTriggered = false;
      int  modeWinner=+999;      


      // loop over segments
      for (int iSeg=0;iSeg<muonNsegs->at(iReco);iSeg++) {
         
        if (isTriggered) continue;

        // look only at matched segments
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

          //loop over  the LCT belonging to the track
          int nLcts = NumLCTsTrk->at(iRaw);

          for (int iLCT=0; iLCT<nLcts; iLCT++) {

            if (printLevel>0) {
              std::cout << "-----------------------------------------------------------\n";
              std::cout  << "iLCT=" << iLCT << endl;
              std::cout  << "trLctEndcap[iRaw][iLCT]   =" << trLctEndcap[iRaw][iLCT]    << std::endl
                         << "trLctSector[iRaw][iLCT]   =" << trLctSector[iRaw][iLCT]    << std::endl
                         << "trLctStation[iRaw][iLCT]  =" << trLctStation[iRaw][iLCT]   << std::endl
                         << "trLctRing[iRaw][iLCT]     =" << trLctRing[iRaw][iLCT]      << std::endl
                         << "trLctChamber[iRaw][iLCT]  =" << trLctChamber[iRaw][iLCT]   << std::endl
                         << "trLctglobalPhi[iRaw][iLCT]=" << trLctglobalPhi[iRaw][iLCT] << std::endl
                         << "trLctglobalEta[iRaw][iLCT]=" << trLctglobalEta[iRaw][iLCT] << std::endl;

            }
            
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
      

//       // old method using DR
//       // identify the closest trigger if any
//       double dRWinner=+999;      
//       int  modeWinner=+999;      

//       int trigCounter=0;
    
//       for (int iRaw=0; iRaw<SizeTrk; iRaw++) {
      
//         if (EndcapTrk->at(iRaw) == -999) continue;
//         if (EndcapTrk->at(iRaw) * etaReco->at(iReco) < 0) continue;
  
//         float diffeta = rchEta->at(iReco) - (EndcapTrk->at(iRaw) * EtaTrk->at(iRaw));
//         float diffphi = rchPhi->at(iReco) - PhiTrk->at(iRaw);
     
//         if (diffphi < -Pi) diffphi += 2*Pi;
//         if (diffphi > +Pi) diffphi -= 2*Pi;
     
//         if (ModeTrk->at(iRaw) == 11 && NumLCTsTrk->at(iRaw) == 1) {
        
//           float lctPhi = fmod(trLctglobalPhi[iRaw][0]*PhiStep +
//                               ((trLctSector[iRaw][0]-1)*Pi/3) + //sector 1 starts at 15 degrees
//                               (Pi/12) , 2*Pi);
       
//           float lctEta = trLctEndcap[iRaw][0]*(trLctglobalEta[iRaw][0]*EtaStep + 0.9);
           
//           float diffetaLCT = rchEta->at(iReco) - lctEta;
//           float diffphiLCT = rchPhi->at(iReco) - lctPhi;
       
//           if (diffphiLCT < -Pi) diffphiLCT += 2*Pi;
//           if (diffphiLCT > +Pi) diffphiLCT -= 2*Pi;
       
//           hDeltaEta->Fill(diffetaLCT);
//           hDeltaPhi->Fill(diffphiLCT);
//           hDR->Fill(DR(diffetaLCT,diffphiLCT));
       
//           if (DR(diffetaLCT,diffphiLCT)<dRWinner) {
//             dRWinner=DR(diffetaLCT,diffphiLCT);
//             modeWinner=ModeTrk->at(iRaw);
//           }
       
//         }
//         else {
//           hDeltaEta->Fill(diffeta);
//           hDeltaPhi->Fill(diffphi);
//           hDR->Fill(DR(diffeta,diffphi));
       
//           if (DR(diffeta,diffphi)<dRWinner) {
//             dRWinner=DR(diffeta,diffphi);
//             modeWinner=ModeTrk->at(iRaw);
//           }
//         }
     
//         trigCounter++;
//       }

      if (!isTriggered) {
        
        nMuonsNoTrigger+=1;
        
        if (printLevel>0) {
          std::cout << "\n\nHouston we have a problem" << std::endl;
          std::cout << "SizeTrk: " << SizeTrk << std::endl;           
          std::cout << "RUN=" << Run   << std::endl;
          std::cout << "EVT=" << Event << std::endl;
        }
        
        // --------------------------------------------------------------
        bool atLeastTwoLCTinWindows=false;
        bool BxTimed=false;
        
        for (int i=0;i<counterSegs;i++) {
          
          //if (atLeastTwoLCTinWindows) continue;
          
          for (int j=i;j<counterSegs;j++) {
            if ( abs( lctEtaBitSeg[i]-lctEtaBitSeg[j] ) <    6 &&
                 abs( lctPhiBitSeg[i]-lctPhiBitSeg[j] ) < 1024 &&
                 lctEndcapSeg[i] == lctEndcapSeg[j]            &&
                 lctStationSeg[i] != lctStationSeg[j]          &&
                 lctSectorSeg[i] == lctSectorSeg[j]             ) {

              atLeastTwoLCTinWindows=true;

              if (abs(lctBxSeg[i]-lctBxSeg[j]) < 3) BxTimed=true;
            }
          }
        }
        
          
        if (printLevel>0) 
          std::cout << " * atLeastTwoLCTinWindows?=" 
                    << atLeastTwoLCTinWindows << std::endl;
        
        if (!atLeastTwoLCTinWindows) nMuonsNoTriggerCuts +=1;
        else if (!BxTimed)           nMuonsNoTriggerTime +=1;
        else                         nMuonsNoTriggerElse +=1;   
        // --------------------------------------------------------------
        if (printLevel>0) {
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
        
        
          // stop to look at the event printout
          int pippo;
          cin >> pippo;
        }

      }

      // Fill the numerator       
      //if (dRWinner != +999) {
      if (isTriggered) {

        hMode->Fill(modeWinner);
      
        hNTriggeredMuonsVsPt  -> Fill ( ptReco ->at(iReco) );
        hNTriggeredMuonsVsEta -> Fill ( etaReco->at(iReco) );
        hNTriggeredMuonsVsPhi -> Fill ( phiReco->at(iReco) );
         
      }// if (isTriggered) {

    } // end loop on reco muons
  }// end Loop evts
   
  if (printLevel > 0)
    printf("\n----------------------------------------"
 	   "----------------------------------------\n");
  
  TString isMCTitle("");
  isMC ? isMCTitle+=" (MC)" : isMCTitle+="";

   
  // --------------------------------------------------------------------------- 
  // hTypeMu
  // ---------------------------------------------------------------------------
  SetStyleh1(hTypeMu, 629, 1, 2, "#mu Type");

  hTypeMu->GetXaxis()->SetBinLabel(1,"GL");
  hTypeMu->GetXaxis()->SetBinLabel(2,"STA");
  hTypeMu->GetXaxis()->SetBinLabel(3,"TRK");
  hTypeMu->GetXaxis()->SetBinLabel(4,"TRK+STA != GBL");

  TCanvas* TypeMu = new TCanvas("TypeMu", "", 0, 0, 400, 400);

  hTypeMu->SetMinimum(0);
  hTypeMu->Draw("histo text");

  PrintIt(TypeMu, "Reco'd Muon Type"+isMCTitle);

  if (isSave) {
    TypeMu -> SaveAs(png+"TypeMuCollisions-TwoSegMatched.png");
    TypeMu -> SaveAs(eps+"TypeMuCollisions-TwoSegMatched.eps");
    TypeMu -> SaveAs(rootPlot+"TypeMuCollisions-TwoSegMatched.root");
  }


  // --------------------------------------------------------------------------- 
  // hMode
  // --------------------------------------------------------------------------- 
  SetStyleh1(hMode, 629, 1, 2, "Mode");

  hMode->SetMinimum(0);

  TCanvas* Mode = new TCanvas("Mode", "", 420, 0, 400, 400);

  hMode->Draw("histo text");

  PrintIt(Mode, "Mode"+isMCTitle);

  if (isSave) {
    Mode -> SaveAs(png+"GBL-ModeTriggerMu-TwoSegMatched.png");
    Mode -> SaveAs(eps+"GBL-ModeTriggerMu-TwoSegMatched.eps");
    Mode -> SaveAs(rootPlot+"GBL-ModeTriggerMu-TwoSegMatched.root");
  }

//   // --------------------------------------------------------------------------- 
//   // hTriggerCounts
//   // ---------------------------------------------------------------------------
//   SetStyleh1(hTriggerCounts, 629, 1, 2, "how many CSCTF triggers?");

//   hTriggerCounts->SetMinimum(0);

//   TCanvas* TriggerCounts = new TCanvas("TriggerCounts", "", 840, 0, 400, 400);

//   hTriggerCounts->Draw("histo text");

//   PrintIt(TriggerCounts, "CSCTF triggers per GBL #mu");

//   if (isSave) {
//     TriggerCounts -> SaveAs(png+"GBL-TriggerCounts-TwoSegMatched.png");
//     TriggerCounts -> SaveAs(eps+"GBL-TriggerCounts-TwoSegMatched.eps");
//     TriggerCounts -> SaveAs(rootPlot+"GBL-TriggerCounts-TwoSegMatched.root");
//   }

  // ---------------------------------------------------------------------------
  // ClopperPearsonBinomialInterval
  // ---------------------------------------------------------------------------
  ClopperPearsonBinomialInterval cp;
  //  alpha = 1 - CL
  const double alpha = (1-0.682);
  cp.init(alpha);


  // ---------------------------------------------------------------------------
  // Efficiency Vs Pt
  // ---------------------------------------------------------------------------
  double pt[8];
  double ptErr[8];
   
  double eff[8];
  double * eefflCP = new double[8];
  double * eeffhCP = new double[8];

  for(int iPt = 0; iPt < 8; ++iPt) {
  
    double lowEdge = hNMuonsVsPt->GetXaxis()->GetBinLowEdge(iPt+1);
    double upEdge  = hNMuonsVsPt->GetXaxis()->GetBinUpEdge(iPt+1);

    pt[iPt]    = (upEdge-lowEdge)/2 + lowEdge;
    ptErr[iPt] = (upEdge-lowEdge)/2;
    
    if (printLevel>0)
      cout << "pt="<<pt[iPt]<<", ptErr="<<ptErr[iPt];
    
    double num = hNTriggeredMuonsVsPt -> GetBinContent(iPt+1);
    double den = hNMuonsVsPt -> GetBinContent(iPt+1);

    // compute the efficiency
    if (den !=0 ) eff[iPt] = num/den;
    else          eff[iPt] = 0;
     
    // compute the error
    cp.calculate(num,den);
    eefflCP[iPt] = eff[iPt] - cp.lower();
    eeffhCP[iPt] = cp.upper() - eff[iPt];
    
    if (printLevel>0)
      std::cout << ", num=" << num
                << ", den=" << den
                << ", eff=" << eff[iPt]
                << ", eefflCP="<<eefflCP[iPt]
                << ", eeffhCP="<<eeffhCP[iPt]
                << std::endl;
  }

  // --- DRAW THE EFFICIENCY VS PT --- 
  gEffvsPt = new TGraphAsymmErrors(8, pt, eff, ptErr, ptErr, eefflCP, eeffhCP);
  gEffvsPt -> SetTitle("");
  gEffvsPt -> SetMinimum(0.);
  gEffvsPt -> SetMaximum(1.02);

  gEffvsPt -> SetLineColor(629);
  gEffvsPt -> SetLineWidth(2);
  gEffvsPt -> SetMarkerStyle(23);
  gEffvsPt -> SetMarkerSize(0.8);

  gEffvsPt -> GetXaxis()-> SetTitle("pt (GeV/c)");
  gEffvsPt -> GetYaxis()-> SetTitle("Efficiency");
  gEffvsPt -> GetYaxis()-> SetTitleOffset(1.35);

  gEffvsPt -> GetXaxis()-> SetNdivisions(509);
  gEffvsPt -> GetYaxis()-> SetNdivisions(514);

  //gEffvsPt -> GetXaxis()->SetRangeUser(0,20);
  TCanvas* EffvsPt = new TCanvas("EffvsPt", "", 0, 500, 600, 600);

  EffvsPt->SetGridx();
  EffvsPt->SetGridy();

  gEffvsPt -> Draw("AP");

  PrintIt(EffvsPt, "CSCTF Efficiency Curve, "+TitleExt+isMCTitle);
   
  if (isSave) {
    EffvsPt -> SaveAs(png+"GBL-EffvsPt-TwoSegMatched.png");
    EffvsPt -> SaveAs(eps+"GBL-EffvsPt-TwoSegMatched.eps");
    EffvsPt -> SaveAs(rootPlot+"GBL-EffvsPt-TwoSegMatched.root");
  }

  // ---------------------------------------------------------------------------
  // Efficiency Vs Eta
  // ---------------------------------------------------------------------------
  double eta[etaNBins];
  double etaErr[etaNBins];
   
  double effEta[etaNBins];
  double * eefflCPEta = new double[etaNBins];
  double * eeffhCPEta = new double[etaNBins];

  for(int iEta = 0; iEta < etaNBins; ++iEta) {
  
    double lowEdge = hNMuonsVsEta->GetXaxis()->GetBinLowEdge(iEta+1);
    double upEdge  = hNMuonsVsEta->GetXaxis()->GetBinUpEdge(iEta+1);

    eta[iEta]    = (upEdge-lowEdge)/2 + lowEdge;
    etaErr[iEta] = (upEdge-lowEdge)/2;
    
    if (printLevel>0)
      cout << "eta="<<eta[iEta]<<", etaErr="<<etaErr[iEta];
    
    double num = hNTriggeredMuonsVsEta -> GetBinContent(iEta+1);
    double den = hNMuonsVsEta -> GetBinContent(iEta+1);

    // compute the efficiency
    if (den !=0 ) effEta[iEta] = num/den;
    else          effEta[iEta] = 0;
     
    // compute the error
    cp.calculate(num,den);
    eefflCPEta[iEta] = effEta[iEta] - cp.lower();
    eeffhCPEta[iEta] = cp.upper() - effEta[iEta];
    
    if (printLevel>0)
      std::cout << ", num=" << num
                << ", den=" << den
                << ", eff=" << effEta[iEta]
                << ", eefflCP="<<eefflCPEta[iEta]
                << ", eeffhCP="<<eeffhCPEta[iEta]
                << std::endl;
  }

  gEffvsEta = new TGraphAsymmErrors(etaNBins, eta, effEta, etaErr, etaErr, eefflCPEta, eeffhCPEta);
  gEffvsEta -> SetTitle("");
  gEffvsEta -> SetMinimum(0.);
  gEffvsEta -> SetMaximum(1.02);

  gEffvsEta -> SetLineColor(629);
  gEffvsEta -> SetLineWidth(2);
  gEffvsEta -> SetMarkerStyle(23);
  gEffvsEta -> SetMarkerSize(0.8);

  gEffvsEta -> GetXaxis()-> SetTitle("|#eta|");
  gEffvsEta -> GetYaxis()-> SetTitle("Efficiency");
  gEffvsEta -> GetYaxis()-> SetTitleOffset(1.35);

  gEffvsEta -> GetXaxis()-> SetNdivisions(509);
  gEffvsEta -> GetYaxis()-> SetNdivisions(514);

  TCanvas* EffvsEta = new TCanvas("EffvsEta", "", 420, 500, 600, 600);

  EffvsEta->SetGridx();
  EffvsEta->SetGridy();

  gEffvsEta -> Draw("AP");

  PrintIt(EffvsEta, "CSCTF Efficiency Curve, "+TitleExt+isMCTitle);

  if (isSave) {
    EffvsEta -> SaveAs(png+"GBL-EffvsEta-TwoSegMatched.png");
    EffvsEta -> SaveAs(eps+"GBL-EffvsEta-TwoSegMatched.eps");
    EffvsEta -> SaveAs(rootPlot+"GBL-EffvsEta-TwoSegMatched.root");
  }


  // ---------------------------------------------------------------------------
  // Efficiency Vs Phi
  // ---------------------------------------------------------------------------
  double phi[phiNBins];
  double phiErr[phiNBins];
   
  double effPhi[phiNBins];
  double * eefflCPPhi = new double[phiNBins];
  double * eeffhCPPhi = new double[phiNBins];

  for(int iPhi = 0; iPhi < phiNBins; ++iPhi) {
  
    double lowEdge = hNMuonsVsPhi->GetXaxis()->GetBinLowEdge(iPhi+1);
    double upEdge  = hNMuonsVsPhi->GetXaxis()->GetBinUpEdge(iPhi+1);

    phi[iPhi]    = (upEdge-lowEdge)/2 + lowEdge;
    phiErr[iPhi] = (upEdge-lowEdge)/2;
    
    if (printLevel>0)
      cout << "phi="<<phi[iPhi]<<", phiErr="<<phiErr[iPhi];
    
    double num = hNTriggeredMuonsVsPhi -> GetBinContent(iPhi+1);
    double den = hNMuonsVsPhi -> GetBinContent(iPhi+1);

    // compute the efficiency
    if (den !=0 ) effPhi[iPhi] = num/den;
    else          effPhi[iPhi] = 0;
     
    // compute the error
    cp.calculate(num,den);
    eefflCPPhi[iPhi] = effPhi[iPhi] - cp.lower();
    eeffhCPPhi[iPhi] = cp.upper() - effPhi[iPhi];
    
    if (printLevel>0)
      std::cout << ", num=" << num
                << ", den=" << den
                << ", eff=" << effPhi[iPhi]
                << ", eefflCP="<<eefflCPPhi[iPhi]
                << ", eeffhCP="<<eeffhCPPhi[iPhi]
                << std::endl;
  }

  gEffvsPhi = new TGraphAsymmErrors(phiNBins, phi, effPhi, phiErr, phiErr, eefflCPPhi, eeffhCPPhi);
  gEffvsPhi -> SetTitle("");
  gEffvsPhi -> SetMinimum(0.);
  gEffvsPhi -> SetMaximum(1.02);

  gEffvsPhi -> SetLineColor(629);
  gEffvsPhi -> SetLineWidth(2);
  gEffvsPhi -> SetMarkerStyle(23);
  gEffvsPhi -> SetMarkerSize(0.8);

  gEffvsPhi -> GetXaxis()-> SetTitle("|#phi|");
  gEffvsPhi -> GetYaxis()-> SetTitle("Efficiency");
  gEffvsPhi -> GetYaxis()-> SetTitleOffset(1.35);

  gEffvsPhi -> GetXaxis()-> SetNdivisions(509);
  gEffvsPhi -> GetYaxis()-> SetNdivisions(514);

  TCanvas* EffvsPhi = new TCanvas("EffvsPhi", "", 840, 500, 600, 600);

  EffvsPhi->SetGridx();
  EffvsPhi->SetGridy();

  gEffvsPhi -> Draw("AP");

  PrintIt(EffvsPhi, "CSCTF Efficiency Curve, "+TitleExt+isMCTitle);

  if (isSave) {
    EffvsPhi -> SaveAs(png+"GBL-EffvsPhi-TwoSegMatched.png");
    EffvsPhi -> SaveAs(eps+"GBL-EffvsPhi-TwoSegMatched.eps");
    EffvsPhi -> SaveAs(rootPlot+"GBL-EffvsPhi-TwoSegMatched.root");
  }

  // ---------------------------------------------------------------------------
  // printout
  // ---------------------------------------------------------------------------
  std::cout << "nMuons              = " << nMuons              << std::endl;
  std::cout << "nMuonsNoTrigger     = " << nMuonsNoTrigger     << std::endl;
  std::cout << "nMuonsNoTriggerTime = " << nMuonsNoTriggerTime << std::endl;
  std::cout << "nMuonsNoTriggerCuts = " << nMuonsNoTriggerCuts << std::endl;
  std::cout << "nMuonsNoTriggerElse = " << nMuonsNoTriggerElse << std::endl;

}

double DR(double diffeta,
          double diffphi) {

  double diffetaSquare =diffeta*diffeta;  
  double diffphiSquare =diffphi*diffphi;  

  return sqrt(diffetaSquare+diffphiSquare);
}

double computeError(double num, double den) {
  
  if (den==0) return 0;
  double efficiency = num/den;
  
  return sqrt( ((1-efficiency) * efficiency)/den );
  
}


