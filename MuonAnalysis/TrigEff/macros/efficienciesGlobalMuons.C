//.L /afs/cern.ch/user/d/digiovan/code/CMSSW_3_5_4/tmp/slc5_ia32_gcc434/src/PhysicsTools/RooStatsCms/src/PhysicsToolsRooStatsCms/libPhysicsToolsRooStatsCms.so
//gSystem->CompileMacro("scripts/Collisions2010/trigEff/efficienciesGlobalMuons.C","Ok")

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


//TGraphErrors* gEffvsPt;   // DR cut
//TGraphErrors* gEffvsEta;  // DR cut
//TGraphErrors* gEffvsPhi;  // DR cut

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
const int MAX_CSC_RECHIT =  35;
const int MAX_TRK_SEGS = 100;
const int MAX_CSCTF_TRK = 36;
const int MAX_LCTS_PER_TRK = 4;


void efficienciesGlobalMuons(int  printLevel =     0,
                                 bool isMC   = !true,
                                 bool isSave = !true) {
    
  gROOT->Clear();
  gStyle->SetOptStat(111111);  
  
  //--------------------------------------------------------------------------
  // Construct the chain of input files
  //--------------------------------------------------------------------------
  cout << "Loading chain reco... \n";
    
  // get the Nutple  
  //#include "scripts/Collisions2010/trigEff/chains/collisionsChain_MinimumBiasOnlyCSCActivity_Reco.C"
#include "scripts/Collisions2010/trigEff/chains/collisionsChain_MinBiasMC_Reco.C"
#include "scripts/Collisions2010/trigEff/chains/collisionsChain_MinimumBias_Reco.C"
  
  if (isMC) recoMuons  = collisionsChainRecoMC;
  else      recoMuons  = collisionsChainReco;
  cout << "recoMuons->GetEntries(): " << recoMuons->GetEntries() << endl;

  // RAW
  //#include "scripts/Collisions2010/trigEff/chains/collisionsChain_MinimumBiasOnlyCSCActivity_Raw.C"
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

  vector<int>*    isGlobalMuon    = new vector<int>();
  vector<int>*    isStandAloneMuon= new vector<int>();
  vector<int>*    isTrackerMuon   = new vector<int>();
  vector<int>*    isTMLastStationAngTight = new vector<int>();

  vector<float>*  ptReco  = new vector<float>();
  vector<float>*  etaReco = new vector<float>();
  vector<float>*  phiReco = new vector<float>();
  vector<float>*  gmrChi2Norm = new vector<float>();
  vector<float>*        gmrDz = new vector<float>();
  vector<float>*        gmrD0 = new vector<float>();
  
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
  //float trkSegZfromMatch[MAX_MUONS][MAX_TRK_SEGS];
  //float trkSegRfromMatch[MAX_MUONS][MAX_TRK_SEGS];
  float trkSegPhifromMatch[MAX_MUONS][MAX_TRK_SEGS];
  //float trkSegEtafromMatch[MAX_MUONS][MAX_TRK_SEGS];

  int trkSegIsArb[MAX_MUONS][MAX_TRK_SEGS];
  float trkSegX[MAX_MUONS][MAX_TRK_SEGS];
  float trkSegY[MAX_MUONS][MAX_TRK_SEGS];
  //float trkSegZ[MAX_MUONS][MAX_TRK_SEGS];
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
    
   
  recoMuons->SetBranchAddress("Run"  , &Run  );
  recoMuons->SetBranchAddress("Event", &Event);
  recoMuons->SetBranchAddress("Bx"   , &Bx   );
  recoMuons->SetBranchAddress("Lumi" , &Lumi );

  recoMuons->SetBranchAddress("muonSize"           , &muonSize        );
  recoMuons->SetBranchAddress("isGlobalMuon"       , &isGlobalMuon    );
  recoMuons->SetBranchAddress("isStandAloneMuon"   , &isStandAloneMuon);
  recoMuons->SetBranchAddress("isTrackerMuon"      , &isTrackerMuon   );
  recoMuons->SetBranchAddress("isTMLastStationAngTight", &isTMLastStationAngTight   );

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

  // l1 ntuple
  //int  l1Size;
  //
  //vector<float>*  ptL1;
  //vector<float>*  etaL1;
  //vector<float>*  phiL1;
  //
  //NtupleL1->SetBranchAddress("l1Size", &l1Size);
  //NtupleL1->SetBranchAddress("l1Pt" ,  &ptL1 );
  //NtupleL1->SetBranchAddress("l1Eta",  &etaL1);
  //NtupleL1->SetBranchAddress("l1Phi",  &phiL1);
  // --------------------------------------------------------------------------- 
  // Define Histograms
  // ---------------------------------------------------------------------------
  hTypeMu = new TH1F("hTypeMu","", 4,    0,   4); 
  hTriggerCounts = new TH1F("hTriggerCounts","",10,0,10); 
  
  hDeltaPhi = new TH1F("hDeltaPhi", "", 100, -2*Pi, 2*Pi);
  hDeltaEta = new TH1F("hDeltaEta", "", 100,  -5.0,  5.0);
  hDR       = new TH1F("hDR"    , "",    50,     0,  1.0);
  hMode     = new TH1F("hMode",     "",  15,     0,   15);

  // --------------------------------------------------------------------------- 
  // Efficiency Variables Pt
  // --------------------------------------------------------------------------- 
  double num[8]      = {0,0,0,0,0,0,0,0};
  double den[8]      = {0,0,0,0,0,0,0,0};
  double effic[8]    = {0,0,0,0,0,0,0,0};
  double errorEff[8] = {0,0,0,0,0,0,0,0};
  
  double pt[8]      = { 1.25,   3,   4,   5,   8,  15,   25,  35};
  double ptError[8] = { 1.25, 0.5, 0.5, 0.5, 2.5, 4.5,  5.5, 4.5};


  // --------------------------------------------------------------------------- 
  // Efficiency Variables Eta
  // --------------------------------------------------------------------------- 
  double numEta[7]      = {0,0,0,0,0,0,0}; 
  double denEta[7]      = {0,0,0,0,0,0,0}; 
  double efficEta[7]    = {0,0,0,0,0,0,0}; 
  double errorEffEta[7] = {0,0,0,0,0,0,0}; 

  double etaMin[7] = {0.9,1.3,1.5,1.7,1.9,2.1,2.3}; 
  double etaMax[7] = {1.3,1.5,1.7,1.9,2.1,2.3,2.5}; 

  double eta[7];
  double etaError[7];

  for (int iEta=0; iEta<7; iEta++) {
    etaError[iEta]=(etaMax[iEta]-etaMin[iEta])/2;
    eta[iEta]=etaMin[iEta]+etaError[iEta];
  }


  // --------------------------------------------------------------------------- 
  // Efficiency Variables Phi
  // --------------------------------------------------------------------------- 
  double numPhi[8]      = {0,0,0,0,0,0,0,0}; 
  double denPhi[8]      = {0,0,0,0,0,0,0,0}; 
  double efficPhi[8]    = {0,0,0,0,0,0,0,0}; 
  double errorEffPhi[8] = {0,0,0,0,0,0,0,0}; 

  double phiMin[8] = {    -Pi, -3*Pi/4, -Pi/2, -Pi/4,    0, +Pi/4,   +Pi/2, +3*Pi/4};
  double phiMax[8] = {-3*Pi/4,   -Pi/2, -Pi/4,     0,+Pi/4, +Pi/2, +3*Pi/4,    Pi  }; 

  double phi[8];
  double phiError[8];

  for (int iPhi=0; iPhi<8; iPhi++) {
    phiError[iPhi]=(phiMax[iPhi]-phiMin[iPhi])/2;
    phi[iPhi]=phiMin[iPhi]+phiError[iPhi];
  }
  
  // All global muons
  TString SaveExt("GblMu");
  TString TitleExt("Gbl #mu");
  
  // Only Gbl mu with key station
  //TString SaveExt("GblWithKeySta");
  //TString TitleExt("Gbl #mu w/ key station");

  // Only Gbl mu with NO key station
  //TString SaveExt("GblWithNoKeySta");
  //TString TitleExt("Gbl #mu w/o key station");

  // All global muons
  //TString SaveExt("StaOnlyMu");
  //TString TitleExt("STA Only #mu");
  
  // Only Gbl mu with key station
  //TString SaveExt("StaOnlyWithKeySta");
  //TString TitleExt("STA only #mu w/ key station");

  //TString SaveExt("StaOnlyWithKeyStaEtaGt2.1");
  //TString TitleExt("STA only #mu w/ key station |#eta| > 2.1");

  // Only Gbl mu with NO key station
  //TString SaveExt("StaOnlyWithNoKeySta");
  //TString TitleExt("STA Only #mu w/o key station");

  // ------------------------------------------------------------------
  // Loop over the events
  // ------------------------------------------------------------------
  for (int iEvt=0; iEvt < recoMuons->GetEntries(); iEvt++) {
  //for (int iEvt=0; iEvt < 30000; iEvt++) {
    
    recoMuons ->GetEntry(iEvt);
    csctfTTree->GetEntry(iEvt);
    
    if ( ( iEvt % 10000) == 0 ) printf(" --- Event # %6d \n", iEvt+1);

    //for (int iReco=0; iReco < gmrSize; iReco++) { 
    for (int iReco=0; iReco < muonSize; iReco++) { 
    
      // --------------------------------------------------------------
      // Analyze the event ONLY if CSC raw info available
      // --------------------------------------------------------------
      // here you select if csc hit, if csc hit in key station or not
      // All global muons
      //if ( rchEta->at(iReco) == -999 ) continue;
      
      // Only mu with key station
      //if ( rchCSCtype->at(iReco) > 15) continue;
     
      // Only mu withOUT key station
      //if ( rchCSCtype->at(iReco) <= 15) continue;     

      // -----------------------------------------------------------
      // global muon definition
      bool isGblMuon = true;
      if (!isGlobalMuon->at(iReco))     isGblMuon = false;
      if (!isStandAloneMuon->at(iReco)) isGblMuon = false;
      if (!isTrackerMuon->at(iReco))    isGblMuon = false;
      if (gmrChi2Norm->at(iReco) > 10)  isGblMuon = false;
      if (fabs(gmrD0->at(iReco)) >  2)  isGblMuon = false;
      if ((rchEta->at(iReco)==-999))    isGblMuon = false;

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
     
      // bool onlySta = true;
      // if ( isGlobalMuon  -> at(iReco) == 1) onlySta = false;
      // if ( isTrackerMuon -> at(iReco) == 1) onlySta = false;
      // if ( isStandAloneMuon -> at(iReco) == 0) onlySta = false;
      // //// noticed that some standalone have eta values above 2.5 up to ~5!
      // //if ( fabs(stdEta -> at(iReco)) > 2.5) onlySta = false;
      // //
      // //if ( fabs(stdEta -> at(iReco)) < 2.1) onlySta = false;
      // if (!onlySta) continue;

      //fill the denominator
      // pt
      if (                          ptReco->at(iReco)<= 2.5) den[0]++;
      if (ptReco->at(iReco)> 2.5 && ptReco->at(iReco)<= 3.5) den[1]++;
      if (ptReco->at(iReco)> 3.5 && ptReco->at(iReco)<= 4.5) den[2]++;
      if (ptReco->at(iReco)> 4.5 && ptReco->at(iReco)<= 5.5) den[3]++;
      if (ptReco->at(iReco)> 5.5 && ptReco->at(iReco)<=10.5) den[4]++;
      if (ptReco->at(iReco)>10.5 && ptReco->at(iReco)<=19.5) den[5]++;
      if (ptReco->at(iReco)>19.5 && ptReco->at(iReco)<=30.5) den[6]++;
      if (ptReco->at(iReco)>30.5 && ptReco->at(iReco)<=39.5) den[7]++;

      //eta
      for (int iEta=0; iEta<7; iEta++) {
        if (fabs(etaReco->at(iReco)) >= etaMin[iEta] &&  
            fabs(etaReco->at(iReco)) <  etaMax[iEta]  )
          denEta[iEta]++;  
      }


      //phi
      for (int iPhi=0; iPhi<8; iPhi++) {
        if (phiReco->at(iReco) >= phiMin[iPhi] &&  
            phiReco->at(iReco) <  phiMax[iPhi]  )
          denPhi[iPhi]++;  
      }

      //cout << "\n\nGlobal muon" << endl;
      //cout << "etaReco->at(iReco)=" << etaReco->at(iReco) << endl;
      //cout << "rchEta->at(iReco)="  << rchEta->at(iReco)  << endl;
      //cout << "phiReco->at(iReco)=" << phiReco->at(iReco) << endl;
      //cout << "rchPhi->at(iReco)="  << rchPhi->at(iReco)  << endl;

      double dRWinner=+999;      
      int  modeWinner=+999;      

      int trigCounter=0;
      
      for (int iRaw=0; iRaw<SizeTrk; iRaw++) {
        
        if (EndcapTrk->at(iRaw) == -999) continue;
        if (EndcapTrk->at(iRaw) * etaReco->at(iReco) < 0) continue;
    
        float diffeta = rchEta->at(iReco) - (EndcapTrk->at(iRaw) * EtaTrk->at(iRaw));
        float diffphi = rchPhi->at(iReco) - PhiTrk->at(iRaw);
       
        if (diffphi < -Pi) diffphi += 2*Pi;
        if (diffphi > +Pi) diffphi -= 2*Pi;
       
        if (ModeTrk->at(iRaw) == 11 && NumLCTsTrk->at(iRaw) == 1) {
          
          float lctPhi = fmod(trLctglobalPhi[iRaw][0]*PhiStep +
                              ((trLctSector[iRaw][0]-1)*Pi/3) + //sector 1 starts at 15 degrees
                              (Pi/12) , 2*Pi);
         
          float lctEta = trLctEndcap[iRaw][0]*(trLctglobalEta[iRaw][0]*EtaStep + 0.9);
             
          float diffetaLCT = rchEta->at(iReco) - lctEta;
          float diffphiLCT = rchPhi->at(iReco) - lctPhi;
         
          if (diffphiLCT < -Pi) diffphiLCT += 2*Pi;
          if (diffphiLCT > +Pi) diffphiLCT -= 2*Pi;
         
          hDeltaEta->Fill(diffetaLCT);
          hDeltaPhi->Fill(diffphiLCT);
          hDR->Fill(DR(diffetaLCT,diffphiLCT));
         
          if (DR(diffetaLCT,diffphiLCT)<0.5 && DR(diffetaLCT,diffphiLCT)<dRWinner) {
            dRWinner=DR(diffetaLCT,diffphiLCT);
            modeWinner=ModeTrk->at(iRaw);
          }
         
        }
        else {
          hDeltaEta->Fill(diffeta);
          hDeltaPhi->Fill(diffphi);
          hDR->Fill(DR(diffeta,diffphi));
         
          if (DR(diffeta,diffphi)<0.5 && DR(diffeta,diffphi)<dRWinner) {
            dRWinner=DR(diffeta,diffphi);
            modeWinner=ModeTrk->at(iRaw);
          }
        }
       
        trigCounter++;
      }
      
      hTriggerCounts->Fill(trigCounter);
         
      if (dRWinner != +999) {

        hMode->Fill(modeWinner);
        
        if (                          ptReco->at(iReco)<= 2.5) num[0]++;
        if (ptReco->at(iReco)> 2.5 && ptReco->at(iReco)<= 3.5) num[1]++;
        if (ptReco->at(iReco)> 3.5 && ptReco->at(iReco)<= 4.5) num[2]++;
        if (ptReco->at(iReco)> 4.5 && ptReco->at(iReco)<= 5.5) num[3]++;
        if (ptReco->at(iReco)> 5.5 && ptReco->at(iReco)<=10.5) num[4]++;
        if (ptReco->at(iReco)>10.5 && ptReco->at(iReco)<=19.5) num[5]++;
        if (ptReco->at(iReco)>19.5 && ptReco->at(iReco)<=30.5) num[6]++;
        if (ptReco->at(iReco)>30.5 && ptReco->at(iReco)<=39.5) num[7]++;
        
        //eta
        for (int iEta=0; iEta<7; iEta++) {
          if (fabs(etaReco->at(iReco)) >= etaMin[iEta] &&  
              fabs(etaReco->at(iReco)) <  etaMax[iEta]  )
            numEta[iEta]++;  
        }

        //phi
        for (int iPhi=0; iPhi<8; iPhi++) {
          if (phiReco->at(iReco) >= phiMin[iPhi] &&  
              phiReco->at(iReco) <  phiMax[iPhi]  )
            numPhi[iPhi]++;  
        }
      }// if (dRWinner != +999) {
    }
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
    TypeMu -> SaveAs(png+"TypeMuCollisions.png");
    TypeMu -> SaveAs(eps+"TypeMuCollisions.eps");
    TypeMu -> SaveAs(rootPlot+"TypeMuCollisions.root");
  }

  // --------------------------------------------------------------------------- 
  // hDeltaPhi
  // ---------------------------------------------------------------------------  
  SetStyleh1(hDeltaPhi, 629, 1, 2, "#phi_{key}-#phi_{CSCTF}");
  
  TCanvas* DeltaPhi 
    = new TCanvas("DeltaPhi", "", 840, 0, 400, 400);

  hDeltaPhi->SetMinimum(0);
  hDeltaPhi->Draw("histo");
  
  PrintIt(DeltaPhi, "#Delta #phi = #phi_{key}-#phi_{CSCTF}"+isMCTitle);
  
  if (isSave) {
    DeltaPhi -> SaveAs(png+"GBL-DeltaPhiKeyMenCSCTF.png");
    DeltaPhi -> SaveAs(eps+"GBL-DeltaPhiKeyMenCSCTF.eps");
    DeltaPhi -> SaveAs(rootPlot+"GBL-DeltaPhiKeyMenCSCTF.root");
  }
  

  // --------------------------------------------------------------------------- 
  // hDeltaEta
  // ---------------------------------------------------------------------------  
  SetStyleh1(hDeltaEta, 629, 1, 2, "#eta_{key}-#eta_{CSCTF}");
    
  TCanvas* DeltaEta 
    = new TCanvas("DeltaEta", "", 1260, 0, 400, 400);

  hDeltaEta->SetMinimum(0);
  hDeltaEta->Draw("histo");
  
  PrintIt(DeltaEta, "#Delta #eta = #eta_{key}-#eta_{CSCTF}"+isMCTitle);

  if (isSave) {
    DeltaEta -> SaveAs(png+"GBL-DeltaEtaKeyMenCSCTF.png");
    DeltaEta -> SaveAs(eps+"GBL-DeltaEtaKeyMenCSCTF.eps");
    DeltaEta -> SaveAs(rootPlot+"GBL-DeltaEtaKeyMenCSCTF.root");
  }

  // --------------------------------------------------------------------------- 
  // hDR
  // ---------------------------------------------------------------------------  
  SetStyleh1(hDR, 629, 1, 2, "#Delta R_{key-CSCTF}");
  
  TCanvas* DRC = new TCanvas("DRC", "", 0, 500, 400, 400);

  hDR->SetMinimum(0);  
  hDR->Draw("histo");
  
  PrintIt(DRC, "#Delta R_{key-CSCTF}"+isMCTitle);
  
  if (isSave) {
    DRC -> SaveAs(png+"GBL-DRKeyMenCSCTF.png");
    DRC -> SaveAs(eps+"GBL-DRKeyMenCSCTF.eps");
    DRC -> SaveAs(rootPlot+"GBL-DRKeyMenCSCTF.root");
  }

  // --------------------------------------------------------------------------- 
  // hMode
  // ---------------------------------------------------------------------------  
  SetStyleh1(hMode, 629, 1, 2, "Mode");
  
  hMode->SetMinimum(0);

  TCanvas* Mode 
    = new TCanvas("Mode", "", 420, 500, 400, 400);

  hMode->Draw("histo text");

  PrintIt(Mode, "Mode"+isMCTitle);

  if (isSave) {
    Mode -> SaveAs(png+"GBL-ModeTriggerMu.png");
    Mode -> SaveAs(eps+"GBL-ModeTriggerMu.eps");
    Mode -> SaveAs(rootPlot+"GBL-ModeTriggerMu.root");
  }

  // --------------------------------------------------------------------------- 
  // hTriggerCounts
  // ---------------------------------------------------------------------------  
  SetStyleh1(hTriggerCounts, 629, 1, 2, "how many CSCTF triggers?");
  
  hTriggerCounts->SetMinimum(0);

  TCanvas* TriggerCounts = new TCanvas("TriggerCounts", "", 840, 500, 400, 400);

  hTriggerCounts->Draw("histo text");

  PrintIt(TriggerCounts, "CSCTF triggers per GBL #mu");

  if (isSave) {
    TriggerCounts -> SaveAs(png+"GBL-TriggerCounts.png");
    TriggerCounts -> SaveAs(eps+"GBL-TriggerCounts.eps");
    TriggerCounts -> SaveAs(rootPlot+"GBL-TriggerCounts.root");
  }

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
  double * eefflCP = new double[8];
  double * eeffhCP = new double[8];

  for (int iPt=0; iPt<8; iPt++) {
    cout << "num[" << iPt <<"]: " << num[iPt] << endl;
    cout << "den[" << iPt <<"]: " << den[iPt] << endl;
    cout << "errorEff[" << iPt <<"]: " << computeError(num[iPt], den[iPt]) << endl;
    errorEff[iPt] = computeError(num[iPt], den[iPt]);
    
    if (den[iPt]!=0) effic[iPt] = num[iPt]/den[iPt];
    else             effic[iPt] = 0;

    cp.calculate(num[iPt],den[iPt]);
    eefflCP[iPt] = effic[iPt] - cp.lower();
    eeffhCP[iPt] = cp.upper() - effic[iPt];

  }
  
  gEffvsPt = new TGraphAsymmErrors(8, pt, effic, ptError, ptError, eefflCP, eeffhCP);
  //gEffvsPt = new TGraphErrors(8, pt, effic, ptError, errorEff);
  gEffvsPt -> SetTitle("");
  gEffvsPt -> SetMinimum(0.8);
  gEffvsPt -> SetMaximum(1.1);
  
  gEffvsPt -> SetLineColor(629);
  gEffvsPt -> SetLineWidth(2);
  gEffvsPt -> SetMarkerStyle(23);
  gEffvsPt -> SetMarkerSize(0.8);
  
  gEffvsPt -> GetXaxis()-> SetTitle("pt (GeV/c)");
  gEffvsPt -> GetYaxis()-> SetTitle("Efficiency [%]");
  gEffvsPt -> GetYaxis()-> SetTitleOffset(1.35);
  
  gEffvsPt -> GetXaxis()-> SetNdivisions(509);
  gEffvsPt -> GetYaxis()-> SetNdivisions(514);

  gEffvsPt -> GetXaxis()->SetRangeUser(0,20);
  TCanvas* EffvsPt = new TCanvas("EffvsPt", "", 0, 500, 600, 600);
  
  EffvsPt->SetGridx();
  EffvsPt->SetGridy();
  
  gEffvsPt -> Draw("AP");
  
  PrintIt(EffvsPt, "CSCTF Efficiency Curve, "+TitleExt+isMCTitle);
  
  if (isSave) {
    EffvsPt -> SaveAs(png+"GBL-EffvsPt"+SaveExt+".png");
    EffvsPt -> SaveAs(eps+"GBL-EffvsPt"+SaveExt+".eps");
    EffvsPt -> SaveAs(rootPlot+"GBL-EffvsPt"+SaveExt+".root");
  }
  
  
  // ---------------------------------------------------------------------------
  // Efficiency Vs Eta
  // ---------------------------------------------------------------------------
  double * eefflCPEta = new double[7];
  double * eeffhCPEta = new double[7];

  for (int iEta=0; iEta<7; iEta++) {
    //     cout << "numEta[" << iEta <<"]: " << numEta[iEta] << endl;
    //     cout << "denEta[" << iEta <<"]: " << denEta[iEta] << endl;    
    //     cout << "errorEffEta[" << iEta <<"]: " << computeError(numEta[iEta], denEta[iEta]) << endl;
    errorEffEta[iEta] = computeError(numEta[iEta], denEta[iEta]);
  
    if (denEta[iEta]!=0) efficEta[iEta] = numEta[iEta]/denEta[iEta];
    else             efficEta[iEta] = 0;

    cp.calculate(numEta[iEta],denEta[iEta]);
    eefflCPEta[iEta] = efficEta[iEta] - cp.lower();
    eeffhCPEta[iEta] = cp.upper() - efficEta[iEta];
  }
  
  gEffvsEta = new TGraphAsymmErrors(7, eta, efficEta, etaError, etaError, eefflCPEta, eeffhCPEta);
  //gEffvsEta = new TGraphErrors(8, eta, efficEta, etaError, errorEffEta);
  gEffvsEta -> SetTitle("");
  gEffvsEta -> SetMinimum(0.8);
  gEffvsEta -> SetMaximum(1.1);
  
  gEffvsEta -> SetLineColor(629);
  gEffvsEta -> SetLineWidth(2);
  gEffvsEta -> SetMarkerStyle(23);
  gEffvsEta -> SetMarkerSize(0.8);
  
  gEffvsEta -> GetXaxis()-> SetTitle("|#eta|");
  gEffvsEta -> GetYaxis()-> SetTitle("Efficiency [%]");
  gEffvsEta -> GetYaxis()-> SetTitleOffset(1.35);

  gEffvsEta -> GetXaxis()-> SetNdivisions(509);
  gEffvsEta -> GetYaxis()-> SetNdivisions(514);
  
  TCanvas* EffvsEta = new TCanvas("EffvsEta", "", 420, 500, 600, 600);
  
  EffvsEta->SetGridx();
  EffvsEta->SetGridy();
  
  gEffvsEta -> Draw("AP");
  
  PrintIt(EffvsEta, "CSCTF Efficiency Curve, "+TitleExt+isMCTitle);
  
  if (isSave) {
    EffvsEta -> SaveAs(png+"GBL-EffvsEta"+SaveExt+".png");
    EffvsEta -> SaveAs(eps+"GBL-EffvsEta"+SaveExt+".eps");
    EffvsEta -> SaveAs(rootPlot+"GBL-EffvsEta"+SaveExt+".root");
  }
  

  // ---------------------------------------------------------------------------
  // Efficiency Vs Phi
  // ---------------------------------------------------------------------------
  double * eefflCPPhi = new double[8];
  double * eeffhCPPhi = new double[8];

  for (int iPhi=0; iPhi<8; iPhi++) {
    //     cout << "numPhi[" << iPhi <<"]: " << numPhi[iPhi] << endl;
    //     cout << "denPhi[" << iPhi <<"]: " << denPhi[iPhi] << endl;    
    //     cout << "errorEffPhi[" << iPhi <<"]: " << computeError(numPhi[iPhi], denPhi[iPhi]) << endl;
    errorEffPhi[iPhi] = computeError(numPhi[iPhi], denPhi[iPhi]);
  
    if (denPhi[iPhi]!=0) efficPhi[iPhi] = numPhi[iPhi]/denPhi[iPhi];
    else             efficPhi[iPhi] = 0;

    cp.calculate(numPhi[iPhi],denPhi[iPhi]);
    eefflCPPhi[iPhi] = efficPhi[iPhi] - cp.lower();
    eeffhCPPhi[iPhi] = cp.upper() - efficPhi[iPhi];

  }

  gEffvsPhi = new TGraphAsymmErrors(8, phi, efficPhi, phiError, phiError, eefflCPPhi, eeffhCPPhi);
  //gEffvsPhi = new TGraphErrors(8, phi, efficPhi, phiError, errorEffPhi);
  gEffvsPhi -> SetTitle("");
  gEffvsPhi -> SetMinimum(0.8);
  gEffvsPhi -> SetMaximum(1.1);

  gEffvsPhi -> SetLineColor(629);
  gEffvsPhi -> SetLineWidth(2);
  gEffvsPhi -> SetMarkerStyle(23);
  gEffvsPhi -> SetMarkerSize(0.8);

  gEffvsPhi -> GetXaxis()-> SetTitle("#phi (rad)");
  gEffvsPhi -> GetYaxis()-> SetTitle("Efficiency [%]");
  gEffvsPhi -> GetYaxis()-> SetTitleOffset(1.35);

  gEffvsPhi -> GetXaxis()-> SetNdivisions(509);
  gEffvsPhi -> GetYaxis()-> SetNdivisions(514);

  TCanvas* EffvsPhi = new TCanvas("EffvsPhi", "", 840, 500, 600, 600);

  EffvsPhi->SetGridx();
  EffvsPhi->SetGridy();

  gEffvsPhi -> Draw("AP");

  PrintIt(EffvsPhi, "CSCTF Efficiency Curve, "+TitleExt+isMCTitle);

  if (isSave) {
    EffvsPhi -> SaveAs(png+"GBL-EffvsPhi"+SaveExt+".png");
    EffvsPhi -> SaveAs(eps+"GBL-EffvsPhi"+SaveExt+".eps");
    EffvsPhi -> SaveAs(rootPlot+"GBL-EffvsPhi"+SaveExt+".root");
  }

} // end 


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


