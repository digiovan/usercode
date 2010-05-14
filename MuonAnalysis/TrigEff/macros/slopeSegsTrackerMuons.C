//gSystem->CompileMacro("scripts/Collisions2010/trigEff/slopeSegsTrackerMuons.C","Ok")

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
// TTree
// ------------------------------------------------------------------------------
TTree* recoMuons;
TTree* csctfTTree;

// ------------------------------------------------------------------------------
// which kind of muon?
// ------------------------------------------------------------------------------
TH1F* hTrkSegsStationNoTrig;
TH1F* hNSegsNoTriggers;

TH1F* hTrkSegsDxDz;
TH1F* hTrkSegsDyDz;

TH1F* hTrkSegsDxDzNoTrig;
TH1F* hTrkSegsDyDzNoTrig;

TH2F* hSlopeCorr;
TH2F* hSlopeCorrNoTrig;

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

void slopeSegsTrackerMuons(int  printLevel =     0,
                               bool isMC   = !true,
                               bool isSave =  !true) {
  
  gROOT->Clear();
  gStyle->SetOptStat(111111);  
  
  //--------------------------------------------------------------------------
  // Construct the chain of input files
  //--------------------------------------------------------------------------
  cout << "Loading chain reco... \n";
    
  // get the Nutple  
  //#include "scripts/Collisions2010/trigEff/chains/collisionsChain_MinimumBiasOnlyCSCActivity_Reco.C"
  #include "scripts/Collisions2010/trigEff/chains/collisionsChain_MinimumBias_Reco.C"
  #include "scripts/Collisions2010/trigEff/chains/collisionsChain_MinBiasMC_Reco.C"
  if (isMC) recoMuons  = collisionsChainRecoMC;
  else      recoMuons  = collisionsChainReco;
  cout << "recoMuons->GetEntries(): " << recoMuons->GetEntries() << endl;

  // RAW
  cout << "Loading chain raw... \n";

  //#include "scripts/Collisions2010/trigEff/chains/collisionsChain_MinimumBiasOnlyCSCActivity_Raw.C"
  #include "scripts/Collisions2010/trigEff/chains/collisionsChain_MinimumBias_Raw.C"
  #include "scripts/Collisions2010/trigEff/chains/collisionsChain_MinBiasMC_Raw.C"
    
  // get the Nutple  
  if (isMC) csctfTTree  = collisionsChainRawMC;
  else      csctfTTree  = collisionsChainRaw;

  cout << "csctfTTree->GetEntries(): " << csctfTTree->GetEntries() << endl;

  // ===========================================================================
  TString eps = "eps/Collisions2010/trigEff/TrkMuOnly/";
  TString png = "png/Collisions2010/trigEff/TrkMuOnly/";
  TString rootPlot = "rootPlot/Collisions2010/trigEff/TrkMuOnly/";

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
 
  float trkSegDxDz[MAX_MUONS][MAX_TRK_SEGS];
  float trkSegDyDz[MAX_MUONS][MAX_TRK_SEGS];
  
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
    
  recoMuons->SetBranchAddress("trkSegDxDz" , &trkSegDxDz );
  recoMuons->SetBranchAddress("trkSegDyDz" , &trkSegDyDz );
  
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
  //csctfTTree->SetBranchAddress("PhiTrk"        , &PhiTrk        );
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
  csctfTTree->SetBranchAddress("trLctEndcap"      , trLctEndcap      );
  csctfTTree->SetBranchAddress("trLctSector"      , trLctSector      );
  csctfTTree->SetBranchAddress("trLctSubSector"   , trLctSubSector   );
  csctfTTree->SetBranchAddress("trLctBx"          , trLctBx          );
  csctfTTree->SetBranchAddress("trLctBx0"         , trLctBx0         );
  csctfTTree->SetBranchAddress("trLctStation"     , trLctStation     );
  csctfTTree->SetBranchAddress("trLctRing"        , trLctRing        );
  csctfTTree->SetBranchAddress("trLctChamber"     , trLctChamber     );
  csctfTTree->SetBranchAddress("trLctTriggerCSCID", trLctTriggerCSCID);
  csctfTTree->SetBranchAddress("trLctFpga"        , trLctFpga        );
  csctfTTree->SetBranchAddress("trLctlocalPhi"    , trLctlocalPhi    );
  csctfTTree->SetBranchAddress("trLctglobalPhi"   , trLctglobalPhi   );
  csctfTTree->SetBranchAddress("trLctglobalEta"   , trLctglobalEta   );
  csctfTTree->SetBranchAddress("trLctstripNum"    , trLctstripNum    );
  csctfTTree->SetBranchAddress("trLctwireGroup"   , trLctwireGroup   );

  
  // ------------------------------------------------------------------ 
  // Define Histograms
  // ------------------------------------------------------------------
  hTrkSegsStationNoTrig = new TH1F("hTrkSegsStationNoTrig","", 10, 0, 11); 
  hNSegsNoTriggers = new TH1F("hNSegsNoTriggers","", 20,    0,  20); 

  hTrkSegsDxDz = new TH1F("hTrkSegsDxDz","", 100, -2, 2); 
  hTrkSegsDyDz = new TH1F("hTrkSegsDyDz","", 100, -2, 2); 

  hTrkSegsDxDzNoTrig = new TH1F("hTrkSegsDxDzNoTrig","", 100, -2, 2); 
  hTrkSegsDyDzNoTrig = new TH1F("hTrkSegsDyDzNoTrig","", 100, -2, 2); 

  hSlopeCorr = new TH2F("hSlopeCorr","",100,-2,2,100,-2,2);  
  hSlopeCorrNoTrig = new TH2F("hSlopeCorrNoTrig","",100,-2,2,100,-2,2);  


  // ------------------------------------------------------------------
  // Loop over the events
  // ------------------------------------------------------------------
  for (int iEvt=0; iEvt < recoMuons->GetEntries(); iEvt++) {
  //for (int iEvt=0; iEvt < 30000; iEvt++) {
    
    recoMuons ->GetEntry(iEvt);
    csctfTTree->GetEntry(iEvt);
    //l1extraMuons->GetEntry(iEvt);
    
    if ( ( iEvt % 5000) == 0 ) printf(" --- Event # %6d \n", iEvt+1);
    //printf(" --- Event # %6d \n", iEvt+1);

    //cout << muonSize << endl;
    for (int iReco=0; iReco < muonSize; iReco++) { 
    
      // -----------------------------------------------------------
      // global muon definition
      bool isGblMuon = true;
      if (!isGlobalMuon->at(iReco))  isGblMuon = false;
      if ((rchEta->at(iReco)==-999)) isGblMuon = false;

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
      if (counter>1) cout << "Error: the same muon fills two category!"
                          << " Review their definitions\n";
      // -----------------------------------------------------------


      // -----------------------------------------------------------
      // for this efficiency studies use only tracker muons
      if (!isTrackerMuonOnly) continue;
      
      int trigCounter=0;

      for (int iRaw=0; iRaw<SizeTrk; iRaw++) {
        
        if (EndcapTrk->at(iRaw) == -999) continue;
        if (EndcapTrk->at(iRaw) * muons_eta_me11->at(iReco) < 0) continue;
          
        trigCounter++;
      }
      
      int StaId=-999;
      float dxdz=-999;
      float dydz=-999;

      for (int iSeg=0;iSeg<trkNSegs->at(iReco);iSeg++) {
        int isArb = trkSegIsArb[iReco][iSeg];
        if (isArb) {
          
          if(trigCounter==0) { 
            //hTrkSegsDxDzNoTrig->Fill(trkSegDxDz[iReco][iSeg]);
            //hTrkSegsDyDzNoTrig->Fill(trkSegDyDz[iReco][iSeg]);
            //
            //hSlopeCorrNoTrig->Fill(trkSegDxDz[iReco][iSeg],
            //                       trkSegDyDz[iReco][iSeg]);

            int StaIdTmp=-999;
            
            if (trkSegStation[iReco][iSeg]==1) StaIdTmp=trkSegRing[iReco][iSeg];
            else
              StaIdTmp=3+(trkSegStation[iReco][iSeg]-2)+trkSegRing[iReco][iSeg];
            
            if(StaId<StaIdTmp) {
              StaId=StaIdTmp;
              dxdz=trkSegDxDz[iReco][iSeg];
              dydz=trkSegDyDz[iReco][iSeg];
            }
            //if (StaId<4) {
            //  hTrkSegsDxDzNoTrig->Fill(trkSegDxDz[iReco][iSeg]);
            //  hTrkSegsDyDzNoTrig->Fill(trkSegDyDz[iReco][iSeg]);
              
            //  hSlopeCorrNoTrig->Fill(trkSegDxDz[iReco][iSeg],
            //                         trkSegDyDz[iReco][iSeg]);
            //}
            
          }
          else {
            hTrkSegsDxDz->Fill(trkSegDxDz[iReco][iSeg]);
            hTrkSegsDyDz->Fill(trkSegDyDz[iReco][iSeg]);

            hSlopeCorr->Fill(trkSegDxDz[iReco][iSeg],
                             trkSegDyDz[iReco][iSeg]);

          }
        }
      }
      
      if(trigCounter==0) {
        hNSegsNoTriggers->Fill(trkNSegs->at(iReco));
        hTrkSegsStationNoTrig->Fill(StaId);
        if (StaId<5) {
          hTrkSegsDxDzNoTrig->Fill(dxdz);
          hTrkSegsDyDzNoTrig->Fill(dydz);
          hSlopeCorrNoTrig->Fill(dxdz,dydz);
        }
      }
    }
  }// end Loop evts
  
  
  if (printLevel > 0)
    printf("\n----------------------------------------"
           "----------------------------------------\n");
  
  
  TString isMCTitle("");
  isMC ? isMCTitle+=" (MC)" : isMCTitle+="";

  // --------------------------------------------------------------------------- 
  // hTrkSegsDxDz
  // ---------------------------------------------------------------------------  
  SetStyleh1(hTrkSegsDxDz, 1, 0, 2, "dX/dZ");
  hTrkSegsDxDz->SetLineStyle(7);
  hTrkSegsDxDz->SetMinimum(0);

  TCanvas* TrkSegsDxDz = new TCanvas("TrkSegsDxDz", "", 0, 0, 400, 400);

  hTrkSegsDxDz->Draw("histo");

  PrintIt(TrkSegsDxDz, "Reco'd Segs dX/dZ"+isMCTitle);

  if (isSave) {
    TrkSegsDxDz -> SaveAs(png+"TrkSegDxDzTrkMu.png");
    TrkSegsDxDz -> SaveAs(eps+"TrkSegDxDzTrkMu.eps");
    TrkSegsDxDz -> SaveAs(rootPlot+"TrkSegDxDzTrkMu.root");
  }

  // --------------------------------------------------------------------------- 
  // hTrkSegsDyDz
  // ---------------------------------------------------------------------------  
  SetStyleh1(hTrkSegsDyDz, 1, 0, 2, "dY/dZ");
  hTrkSegsDyDz->SetLineStyle(7);
  hTrkSegsDyDz->SetMinimum(0);

  TCanvas* TrkSegsDyDz = new TCanvas("TrkSegsDyDz", "", 420, 0, 400, 400);

  hTrkSegsDyDz->Draw("histo");

  PrintIt(TrkSegsDyDz, "Reco'd Segs dY/dZ"+isMCTitle);

  if (isSave) {
    TrkSegsDyDz -> SaveAs(png+"TrkSegDyDzTrkMu.png");
    TrkSegsDyDz -> SaveAs(eps+"TrkSegDyDzTrkMu.eps");
    TrkSegsDyDz -> SaveAs(rootPlot+"TrkSegDyDzTrkMu.root");
  }


  // --------------------------------------------------------------------------- 
  // hTrkSegsDxDzNoTrig
  // ---------------------------------------------------------------------------  
  SetStyleh1(hTrkSegsDxDzNoTrig, 629, 1, 2, "dX/dZ");
  hTrkSegsDxDzNoTrig->SetMinimum(0);

  TCanvas* TrkSegsDxDzNoTrig = new TCanvas("TrkSegsDxDzNoTrig", "", 
                                           0, 500, 400, 400);

  hTrkSegsDxDzNoTrig->Draw("histo");

  PrintIt(TrkSegsDxDzNoTrig, "Reco'd Segs dX/dZ"+isMCTitle);

  if (isSave) {
    TrkSegsDxDzNoTrig -> SaveAs(png+"TrkSegDxDzNoTrigTrkMu.png");
    TrkSegsDxDzNoTrig -> SaveAs(eps+"TrkSegDxDzNoTrigTrkMu.eps");
    TrkSegsDxDzNoTrig -> SaveAs(rootPlot+"TrkSegDxDzNoTrigTrkMu.root");
  }

  // --------------------------------------------------------------------------- 
  // hTrkSegsDyDzNoTrig
  // ---------------------------------------------------------------------------  
  SetStyleh1(hTrkSegsDyDzNoTrig, 414, 1, 2, "dY/dZ");
  hTrkSegsDyDzNoTrig->SetMinimum(0);

  TCanvas* TrkSegsDyDzNoTrig = new TCanvas("TrkSegsDyDzNoTrig", "", 
                                           420, 500, 400, 400);

  hTrkSegsDyDzNoTrig->Draw("histo");

  PrintIt(TrkSegsDyDzNoTrig, "Reco'd Segs dY/dZ"+isMCTitle);

  if (isSave) {
    TrkSegsDyDzNoTrig -> SaveAs(png+"TrkSegDyDzNoTrigTrkMu.png");
    TrkSegsDyDzNoTrig -> SaveAs(eps+"TrkSegDyDzNoTrigTrkMu.eps");
    TrkSegsDyDzNoTrig -> SaveAs(rootPlot+"TrkSegDyDzNoTrigTrkMu.root");
  }

  
  // --------------------------------------------------------------------------- 
  // hSlopeCorr
  // ---------------------------------------------------------------------------  
  SetStyleh2(hSlopeCorr, 1, 0.1, "dX/dZ", "dY/dZ");
  
  TCanvas* SlopeCorr = new TCanvas("SlopeCorr", "", 840, 0, 400, 400);
  hSlopeCorr ->Draw();
  
  PrintIt(SlopeCorr, "dY/dZ vs dX/dZ"+isMCTitle);
  
  if (isSave) {
    SlopeCorr -> SaveAs(png+"SlopeCorr.png");
    SlopeCorr -> SaveAs(eps+"SlopeCorr.eps");
    SlopeCorr -> SaveAs(rootPlot+"SlopeCorr.root");
  }


  // --------------------------------------------------------------------------- 
  // hSlopeCorrNoTrig
  // ---------------------------------------------------------------------------  
  SetStyleh2(hSlopeCorrNoTrig, 1, 0.1, "dX/dZ", "dY/dZ");
  
  TCanvas* SlopeCorrNoTrig = new TCanvas("SlopeCorrNoTrig", "", 
                                         840, 500, 400, 400);
  hSlopeCorrNoTrig ->Draw();
  
  PrintIt(SlopeCorrNoTrig, "dY/dZ vs dX/dZ"+isMCTitle);
  
  if (isSave) {
    SlopeCorrNoTrig -> SaveAs(png+"SlopeCorrNoTrig.png");
    SlopeCorrNoTrig -> SaveAs(eps+"SlopeCorrNoTrig.eps");
    SlopeCorrNoTrig -> SaveAs(rootPlot+"SlopeCorrNoTrig.root");
  }

  
  // --------------------------------------------------------------------------- 
  // hTrkSegsDxDzCompare
  // ---------------------------------------------------------------------------  
  TCanvas* TrkSegsDxDzCompare = new TCanvas("TrkSegsDxDzCompare", "", 
                                            0, 1000, 400, 400);

  compareHistos(hTrkSegsDxDz,hTrkSegsDxDzNoTrig,TrkSegsDxDzCompare,true);
  
  PrintIt(TrkSegsDxDzCompare, "Reco'd Segs dX/dZ"+isMCTitle);
  
  if (isSave) {
    TrkSegsDxDzCompare -> SaveAs(png+"TrkSegDxDzCompareTrkMu.png");
    TrkSegsDxDzCompare -> SaveAs(eps+"TrkSegDxDzCompareTrkMu.eps");
    TrkSegsDxDzCompare -> SaveAs(rootPlot+"TrkSegDxDzCompareTrkMu.root");
  }


  // --------------------------------------------------------------------------- 
  // hTrkSegsDyDzCompare
  // ---------------------------------------------------------------------------  
  TCanvas* TrkSegsDyDzCompare = new TCanvas("TrkSegsDyDzCompare", "", 
                                            420, 1000, 400, 400);

  compareHistos(hTrkSegsDyDz,hTrkSegsDyDzNoTrig,TrkSegsDyDzCompare,true);

  PrintIt(TrkSegsDyDzCompare, "Reco'd Segs dY/dZ"+isMCTitle);
  
  if (isSave) {
    TrkSegsDyDzCompare -> SaveAs(png+"TrkSegDyDzCompareTrkMu.png");
    TrkSegsDyDzCompare -> SaveAs(eps+"TrkSegDyDzCompareTrkMu.eps");
    TrkSegsDyDzCompare -> SaveAs(rootPlot+"TrkSegDyDzCompareTrkMu.root");
  }


  // --------------------------------------------------------------------------- 
  // hTrkSegsStationNoTrig
  // ---------------------------------------------------------------------------  
  SetStyleh1(hTrkSegsStationNoTrig, 414, 1, 2, "Station/Ring");
  TCanvas* TrkSegsStationNoTrig = new TCanvas("TrkSegsStationNoTrig", "", 
                                            840, 1000, 400, 400);

  hTrkSegsStationNoTrig->GetXaxis()->SetBinLabel( 1,"ME1/1");
  hTrkSegsStationNoTrig->GetXaxis()->SetBinLabel( 2,"ME1/2");
  hTrkSegsStationNoTrig->GetXaxis()->SetBinLabel( 3,"ME1/3");
  hTrkSegsStationNoTrig->GetXaxis()->SetBinLabel( 4,"ME1/1a");
  hTrkSegsStationNoTrig->GetXaxis()->SetBinLabel( 5,"ME2/1");
  hTrkSegsStationNoTrig->GetXaxis()->SetBinLabel( 6,"ME2/2");
  hTrkSegsStationNoTrig->GetXaxis()->SetBinLabel( 7,"ME3/1");
  hTrkSegsStationNoTrig->GetXaxis()->SetBinLabel( 8,"ME3/2");
  hTrkSegsStationNoTrig->GetXaxis()->SetBinLabel( 9,"ME4/1");
  hTrkSegsStationNoTrig->GetXaxis()->SetBinLabel(10,"ME4/2");

  hTrkSegsStationNoTrig->Draw("histo");

  PrintIt(TrkSegsStationNoTrig, "Station/Ring Segs, no CSCTF trigger info"+isMCTitle);
  
  if (isSave) {
    TrkSegsStationNoTrig -> SaveAs(png+"TrkSegsStationNoTrigTrkMu.png");
    TrkSegsStationNoTrig -> SaveAs(eps+"TrkSegsStationNoTrigTrkMu.eps");
    TrkSegsStationNoTrig -> SaveAs(rootPlot+"TrkSegsStationNoTrigTrkMu.root");
  }
  
  // --------------------------------------------------------------------------- 
  // hNSegsNoTriggers
  // ---------------------------------------------------------------------------  
  SetStyleh1(hNSegsNoTriggers, 629, 1, 2, "# of segs");
  
  TCanvas* NSegsNoTriggers = new TCanvas("NSegsNoTriggers", "", 1260, 1000, 400, 400);

  hNSegsNoTriggers->Draw();

  PrintIt(NSegsNoTriggers, "# segs (trk only #mu) no CSCTF trig info"+isMCTitle);

  if (isSave) {
    NSegsNoTriggers -> SaveAs(png+"NSegsNoTriggersTrkMuCollisions.png");
    NSegsNoTriggers -> SaveAs(eps+"NSegsNoTriggersTrkMuCollisions.eps");
    NSegsNoTriggers -> SaveAs(rootPlot+"NSegsNoTriggersTrkMuCollisions.root");
  }

} // end 



