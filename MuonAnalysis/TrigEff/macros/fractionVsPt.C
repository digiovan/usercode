//gSystem->CompileMacro("scripts/Collisions2010/trigEff/fractionVsPt.C","Ok")

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
using namespace std;

#include "/afs/cern.ch/user/d/digiovan/scripts/Init/modifiedStyle.C"

// ------------------------------------------------------------------------------
// useful functions
void FillFraction(TH1F* histoOrig,
                  TH1F* histoNorm,
                  double fra[],
                  double fraErr[]);

double computeError(double num, double den);
// ------------------------------------------------------------------------------

TChain* recoMuons;
TChain* csctfTTree;

// ------------------------------------------------------------------------------
// which kind of muon?
// ------------------------------------------------------------------------------
TGraphErrors* gFractionGblVsPt;
TGraphErrors* gFractionTrkOnlyVsPt;
TGraphErrors* gFractionStaOnlyVsPt;
TGraphErrors* gFractionTrkStaNoGblVsPt;

TH1F* hPt;
TH1F* hPtGbl;
TH1F* hPtStaOnly;
TH1F* hPtTrkOnly;
TH1F* hPtTrkStaNoGbl;

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

void fractionVsPt(int  printLevel =     0,
                      bool isMC   =  true,
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
  TString eps = "eps/Collisions2010/trigEff/Fractions/";
  TString png = "png/Collisions2010/trigEff/Fractions/";
  TString rootPlot = "rootPlot/Collisions2010/trigEff/Fractions/";
  
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
 
  csctfTTree->SetBranchAddress("SizeTrk"       , &SizeTrk      );
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

  // ------------------------------------------------------------------ 
  // Define Histograms
  // ------------------------------------------------------------------
  const int nBins=20;

  hPt            = new TH1F("hPt"           , "", nBins, 0, nBins);
  hPtGbl         = new TH1F("hPtGbl"        , "", nBins, 0, nBins);
  hPtStaOnly     = new TH1F("hPtStaOnly"    , "", nBins, 0, nBins);
  hPtTrkOnly     = new TH1F("hPtTrkOnly"    , "", nBins, 0, nBins);
  hPtTrkStaNoGbl = new TH1F("hPtTrkStaNoGbl", "", nBins, 0, nBins);

  // --------------------------------------------------------------------------- 
  // Useful value to compute the fractions
  // --------------------------------------------------------------------------- 
  double fraGbl[nBins], fraTrkOnly[nBins], fraStaOnly[nBins];
  double fraErrGbl[nBins], fraErrTrkOnly[nBins], fraErrStaOnly[nBins];

  double pt[nBins],ptError[nBins];

  //initialization  
  for (int i=0;i<nBins;i++){
    fraGbl[i]=0;
    fraTrkOnly[i]=0;
    fraStaOnly[i]=0;
    fraErrGbl[i]=0;
    fraErrTrkOnly[i]=0;
    fraErrStaOnly[i]=0;

    pt[i]=0.5+i;
    ptError[i]=0.5;
  }

  // ------------------------------------------------------------------
  // Loop over the events
  // ------------------------------------------------------------------
  for (int iEvt=0; iEvt < recoMuons->GetEntries(); iEvt++) {
    //for (int iEvt=0; iEvt < 10000; iEvt++) {
    recoMuons ->GetEntry(iEvt);
    csctfTTree->GetEntry(iEvt);
    
    if ( ( iEvt % 5000) == 0 ) printf(" --- Event # %6d \n", iEvt+1);
    
    for (int iReco=0; iReco < muonSize; iReco++) { 
            
      // -----------------------------------------------------------
      // global muon definition
      bool isGblMuon = true;
      if (!isGlobalMuon->at(iReco))     isGblMuon = false;
      if (!isStandAloneMuon->at(iReco)) isGblMuon = false;
      if (!isTrackerMuon->at(iReco))    isGblMuon = false;
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
      if (isGblMuon)            {/*hTypeMu->Fill(0);*/ counter++;}
      if (isStandAloneMuonOnly) {/*hTypeMu->Fill(1);*/ counter++;}
      if (isTrackerMuonOnly)    {/*hTypeMu->Fill(2);*/ counter++;}
      if (isStdAndTrkButNotGbl) {/*hTypeMu->Fill(3);*/ counter++;}
      
      // sanity check
      if (counter>1) {
        cout << "Error: the same muon fills two category!"
             << " Review their definitions\n";
        continue;
      }
      // -----------------------------------------------------------

    
      // -----------------------------------------------------------
      // for this efficiency studies use only tracker muons
      if (isGblMuon) {
        hPtGbl          ->Fill(trkPt->at(iReco));
        
        hPt->Fill(trkPt->at(iReco));
      }
      if (isStandAloneMuonOnly) {
        hPtStaOnly          ->Fill(stdPt->at(iReco));
  
        hPt->Fill(stdPt->at(iReco));
      }
      if (isTrackerMuonOnly) {
        hPtTrkOnly          ->Fill(trkPt->at(iReco));
  
        hPt->Fill(trkPt->at(iReco));
      }
      // -----------------------------------------------------------

    }
  
  }// end Loop evts
  
  FillFraction(hPtGbl,hPt,fraGbl,fraErrGbl);
  FillFraction(hPtTrkOnly,hPt,fraTrkOnly,fraErrTrkOnly);
  FillFraction(hPtStaOnly,hPt,fraStaOnly,fraErrStaOnly);
  
  if (printLevel > 0)
       printf("\n----------------------------------------"
       "----------------------------------------\n");

  TString isMCTitle("");
  isMC ? isMCTitle+=" (MC)" : isMCTitle+="";
  
//   // --------------------------------------------------------------------------- 
//   // hTypeMu
//   // ---------------------------------------------------------------------------  
//   SetStyleh1(hTypeMu, 629, 1, 2, "#mu Type");
  
//   hTypeMu->GetXaxis()->SetBinLabel(1,"GL");
//   hTypeMu->GetXaxis()->SetBinLabel(2,"STA");
//   hTypeMu->GetXaxis()->SetBinLabel(3,"TRK");
//   hTypeMu->GetXaxis()->SetBinLabel(4,"TRK+STA != GBL");
  
//   hTypeMu->SetMinimum(0);

//   TCanvas* TypeMu = new TCanvas("TypeMu", "", 0, 0, 400, 400);

//   hTypeMu->Draw("histo text");

//   PrintIt(TypeMu, "Reco'd Muon Type");

//   if (isSave) {
//     TypeMu -> SaveAs(png+"TypeMuCollisions.png");
//     TypeMu -> SaveAs(eps+"TypeMuCollisions.eps");
//     TypeMu -> SaveAs(rootPlot+"TypeMuCollisions.root");
//   }


  // --------------------------------------------------------------------------- 
  // hPt
  // ---------------------------------------------------------------------------  
  SetStyleh1(hPt, 1, 1, 3, "P_{T} (GeV/c)");
  
  // hPt->SetMinimum(0);
  
  TCanvas* Pt = new TCanvas("Pt", "", 420, 0, 400, 400);
  Pt->SetLogy();

  hPt->Draw("histo");

  PrintItLog(Pt, "P_{T}^{#mu}"+isMCTitle);
  
  if (isSave) {
    Pt -> SaveAs(png+"PtAllMusCollisions.png");
    Pt -> SaveAs(eps+"PtAllMusCollisions.eps");
    Pt -> SaveAs(rootPlot+"PtAllMusCollisions.root");
  }

  // --------------------------------------------------------------------------- 
  // hPtCompare
  // ---------------------------------------------------------------------------  
  SetStyleh1(hPtGbl, 629, 3001, 2, "P_{T} (GeV/c)");
  SetStyleh1(hPtTrkOnly, 597, 3001, 2, "P_{T} (GeV/c)");
  SetStyleh1(hPtStaOnly, 414, 3001, 2, "P_{T} (GeV/c)");
  
  TCanvas* PtCompare = new TCanvas("PtCompare", "", 840, 0, 400, 400);
  
  hPt       ->Draw("histo");
  hPtTrkOnly->Draw("histo same");
  hPtGbl    ->Draw("histo same");
  hPtStaOnly->Draw("histo same");

  tl = SetLegend(0.6, 0.67, 0.82, 0.87);
  tl->AddEntry(hPt       , " Total"   , "l");
  tl->AddEntry(hPtGbl    , " Gbl"     , "l");
  tl->AddEntry(hPtTrkOnly, " Trk Only", "l");
  tl->AddEntry(hPtStaOnly, " Sta Only", "l");
  tl -> Draw("same");

  PtCompare->SetLogy();
  PrintItLog(PtCompare, "P_{T}^{#mu} Composition"+isMCTitle);
  
  if (isSave) {
    PtCompare -> SaveAs(png+"PtMuCompositionCollisions.png");
    PtCompare -> SaveAs(eps+"PtMuCompositionCollisions.eps");
    PtCompare -> SaveAs(rootPlot+"PtMuCompositionCollisions.root");
  }


  // --------------------------------------------------------------------------- 
  // hFractionCompare
  // ---------------------------------------------------------------------------  
  gFractionGblVsPt = new TGraphErrors(nBins, pt, fraGbl, 
                                       ptError, fraErrGbl);
  gFractionGblVsPt -> SetTitle("");
  gFractionGblVsPt -> SetMinimum(0.0);
  gFractionGblVsPt -> SetMaximum(1.4);

  gFractionGblVsPt -> SetLineColor(629);
  gFractionGblVsPt -> SetLineWidth(2);
  gFractionGblVsPt -> SetMarkerStyle(23);
  gFractionGblVsPt -> SetMarkerSize(0.8);

  gFractionGblVsPt -> GetXaxis()-> SetTitle("pt (GeV/c)");
  gFractionGblVsPt -> GetYaxis()-> SetTitle("Percentage [%]");
  gFractionGblVsPt -> GetYaxis()-> SetTitleOffset(1.35);

  gFractionGblVsPt -> GetXaxis()-> SetNdivisions(509);
  gFractionGblVsPt -> GetYaxis()-> SetNdivisions(514);


  gFractionTrkOnlyVsPt = new TGraphErrors(nBins, pt, fraTrkOnly, 
                                          ptError, fraErrTrkOnly);
  gFractionTrkOnlyVsPt -> SetTitle("");
  gFractionTrkOnlyVsPt -> SetMinimum(0.0);
  gFractionTrkOnlyVsPt -> SetMaximum(1.4);

  gFractionTrkOnlyVsPt -> SetLineColor(597);
  gFractionTrkOnlyVsPt -> SetLineWidth(2);
  gFractionTrkOnlyVsPt -> SetMarkerStyle(23);
  gFractionTrkOnlyVsPt -> SetMarkerSize(0.8);

  gFractionTrkOnlyVsPt -> GetXaxis()-> SetTitle("pt (GeV/c)");
  gFractionTrkOnlyVsPt -> GetYaxis()-> SetTitle("Percentage [%]");
  gFractionTrkOnlyVsPt -> GetYaxis()-> SetTitleOffset(1.35);

  gFractionTrkOnlyVsPt -> GetXaxis()-> SetNdivisions(509);
  gFractionTrkOnlyVsPt -> GetYaxis()-> SetNdivisions(514);


  gFractionStaOnlyVsPt = new TGraphErrors(nBins, pt, fraStaOnly, 
                                          ptError, fraErrStaOnly);
  gFractionStaOnlyVsPt -> SetTitle("");
  gFractionStaOnlyVsPt -> SetMinimum(0.0);
  gFractionStaOnlyVsPt -> SetMaximum(1.4);

  gFractionStaOnlyVsPt -> SetLineColor(414);
  gFractionStaOnlyVsPt -> SetLineWidth(2);
  gFractionStaOnlyVsPt -> SetMarkerStyle(23);
  gFractionStaOnlyVsPt -> SetMarkerSize(0.8);

  gFractionStaOnlyVsPt -> GetXaxis()-> SetTitle("pt (GeV/c)");
  gFractionStaOnlyVsPt -> GetYaxis()-> SetTitle("Percentage [%]");
  gFractionStaOnlyVsPt -> GetYaxis()-> SetTitleOffset(1.35);

  gFractionStaOnlyVsPt -> GetXaxis()-> SetNdivisions(509);
  gFractionStaOnlyVsPt -> GetYaxis()-> SetNdivisions(514);

  
  TCanvas* FractionCompare = new TCanvas("FractionCompare", "", 
                                         1260, 0, 400, 400);
  
  FractionCompare->SetGridx();
  FractionCompare->SetGridy();

  gFractionGblVsPt -> Draw("AP");
  gFractionTrkOnlyVsPt -> Draw("P same");
  gFractionStaOnlyVsPt -> Draw("P same");

  tl = SetLegend(0.6, 0.67, 0.82, 0.87);
  tl->AddEntry(gFractionGblVsPt    , " Gbl"     , "l");
  tl->AddEntry(gFractionTrkOnlyVsPt, " Trk Only", "l");
  tl->AddEntry(gFractionStaOnlyVsPt, " Sta Only", "l");
tl -> Draw("same");

  PrintIt(FractionCompare, "Percentage Vs Pt"+isMCTitle);
  
  if (isSave) {
    FractionCompare -> SaveAs(png+"PercentageVsPtCollisions.png");
    FractionCompare -> SaveAs(eps+"PercentageVsPtCollisions.eps");
    FractionCompare -> SaveAs(rootPlot+"PercentageVsPtCollisions.root");
  }
  
}// end 


void FillFraction(TH1F* histoOrig,
                  TH1F* histoNorm,
                  double fra[],
                  double fraErr[]) {

  if (histoOrig->GetNbinsX()!=histoNorm->GetNbinsX()) {
    std::cout<<"Error: the two histograms must have the same number of bins\n";
    return;
  }
  
  //cout<<"Nbins="<<histoOrig->GetNbinsX()<<endl;
  
  for (int iBinX=1;iBinX<histoOrig->GetNbinsX()+1; iBinX++) {
    double percentage=0;
    double error=0;
    
    if (histoNorm->GetBinContent(iBinX)!=0) {
      percentage 
        = histoOrig->GetBinContent(iBinX)/histoNorm->GetBinContent(iBinX);
      error=computeError(histoOrig->GetBinContent(iBinX),
                         histoNorm->GetBinContent(iBinX));
    }
    
    //fill the array   
    fra[iBinX-1]=percentage;
    fraErr[iBinX-1]=error;
    
    //cout<<"fra="<<fra[iBinX-1]<<endl;
    //cout<<"fraErr="<<fraErr[iBinX-1]<<endl;
  }
}

double computeError(double num, double den) {
  
  if (den==0) return 0;
  double efficiency = num/den;
  
  return sqrt( ((1-efficiency) * efficiency)/den );
  
}


