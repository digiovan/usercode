// ------------------------------------------------------------------------------
// TTree
// ------------------------------------------------------------------------------
TTree* Ntuple;

// ------------------------------------------------------------------------------
// histograms
// ------------------------------------------------------------------------------
TH1F* hL1ABXN; //L1A BXN from GT

// ------------------------------------------------------------------------------
// variables
// ------------------------------------------------------------------------------
double Pi  = TMath::Pi();

const int MAXCSCTFTR = 60;
const int MAXCSCTFLCTSTR = 4;
const int MAXCSCTFSPS = 12;

void timeBasic(TString runNumber  = "123893",
	       int     printLevel =      0 ,
	       bool    isSave     =      !true) {


  
  gROOT->Clear();
  gInterpreter->LoadMacro("/afs/cern.ch/user/d/digiovan/scripts/Init/"
			  "modifiedStyle.C");

  
  
  //--------------------------------------------------------------------------
  // Construct the chain of input files
  //--------------------------------------------------------------------------
  gInterpreter->LoadMacro("./scripts/cosmicsChain.C");
    
  TChain* Chain  = cosmicsChain();
  
  // get the Nutple
  Ntuple  = (TTree*) Chain->Get("Ntuple");
  cout << "Ntuple->GetEntries(): " << Ntuple->GetEntries() << endl;
  
  //----------------------------------------------------------------------------
  // Access the needed variables
  //----------------------------------------------------------------------------
  int    run, event, Bx;

  int    csctf_trSize;

  Ntuple -> SetBranchAddress("run", &run );
  Ntuple -> SetBranchAddress("event", &event );
  Ntuple -> SetBranchAddress("Bx", &Bx );
  Ntuple -> SetBranchAddress("csctf_trSize", &csctf_trSize);

  // --------------------------------------------------------------------------- 
  // Histograms
  // ---------------------------------------------------------------------------  
  hL1ABXN    = new TH1F("hL1ABXN", "", 3600, 0 , 3600);

 // --------------------------------------------------------------------------- 
  // Loop over the events
  // ---------------------------------------------------------------------------  
  for (int iEvt=0; iEvt < Ntuple->GetEntries(); iEvt++) {
   
    Ntuple->GetEntry(iEvt);

    if (csctf_trSize
     hDeltaEta -> Fill(deltaEta);    
	
    
  }// end Loop evts
  if (printLevel > 0)
    printf("\n----------------------------------------"
	   "----------------------------------------\n");
    
  // --------------------------------------------------------------------------- 
  // Delta Eta
  // ---------------------------------------------------------------------------  
  SetStyleh1(hDeltaEta   , 622, 1001, 1, "#phi");
  TCanvas* DeltaEta 
    = new TCanvas("DeltaEta", "Delta eta",   0, 0, 400, 400);
  
  hDeltaEta -> Draw();
  
  PrintIt(DeltaEta, "#Delta #eta = #eta_{RecHit} - #eta_{Raw}");

 
  // --------------------------------------------------------------------------- 
  // Delta Phi
  // ---------------------------------------------------------------------------  
  SetStyleh1(hDeltaPhi   , 622, 1001, 1, "#phi");
  TCanvas* DeltaPhi 
    = new TCanvas("DeltaPhi", "Delta phi", 420, 0, 400, 400);
  
  hDeltaPhi -> Draw();
  PrintIt(DeltaPhi, "#Delta #phi = #phi_{RecHit} - #phi_{Raw}");


  // --------------------------------------------------------------------------- 
  // Pt
  // ---------------------------------------------------------------------------  
  SetStyleh2(hPt, 1, 0.1, "P_{T}^{raw} [Gev]", "P_{T}^{reco} [Gev]");
  TCanvas* Pt 
    = new TCanvas("Pt", "Pt", 840, 0, 400, 400);
  
  hPt -> Draw("COLZ");
  SetMargin(Pt);

  titleCanv = "P_{T}^{raw} vs P_{T}^{reco}";
  if (phiMin ==0  && phiMax == Pi  ) titleCanv += (" (Top Sectors)"   );
  if (phiMin ==Pi && phiMax == 2*Pi) titleCanv += (" (Bottom Sectors)");
  
  PrintIt(Pt, titleCanv.Data());

  if (isSave) {
    Pt->SaveAs("eps/"+runNumber+"/PtRawVsPtReco-"+runNumber+"-"+phiCut+".eps");
    Pt->SaveAs("png/"+runNumber+"/PtRawVsPtReco-"+runNumber+"-"+phiCut+".png");
  }


  // --------------------------------------------------------------------------- 
  // Pt Profile X
  // --------------------------------------------------------------------------- 
  hPt -> ProfileX();
  
  TCanvas* Pt_pfx 
    = new TCanvas("Pt_pfx", "Pt ProfileX", 1260, 0, 400, 400);
  
  // Additional makeup...
  hPt_pfx -> SetMarkerSize (0.8);
  hPt_pfx -> SetMarkerColor(629);
  hPt_pfx -> SetLineWidth  (1);
  
  hPt_pfx -> GetYaxis()    -> SetTitle("P_{T}^{reco} [Gev]");
  hPt_pfx -> GetYaxis()    -> SetTitleOffset(1);

  hPt_pfx -> Draw("");

  titleCanv = "Profile in P_{T}^{raw}";
  if (phiMin ==0  && phiMax == Pi  ) titleCanv += (" (Top Sectors)"   );
  if (phiMin ==Pi && phiMax == 2*Pi) titleCanv += (" (Bottom Sectors)");

  PrintIt(Pt_pfx, titleCanv.Data());

  if (isSave) {
    Pt_pfx->SaveAs("eps/"+runNumber+"/PtRawVsPtReco_pfx-"
		   +runNumber+"-"+phiCut+".eps");
    Pt_pfx->SaveAs("png/"+runNumber+"/PtRawVsPtReco_pfx-"
		   +runNumber+"-"+phiCut+".png");
  }


  // --------------------------------------------------------------------------- 
  // Pt Profile Y
  // --------------------------------------------------------------------------- 
  hPt -> ProfileY();
  
  TCanvas* Pt_pfy 
    = new TCanvas("Pt_pfy", "Pt ProfileY", 0, 500, 400, 400);
  
  // Additional makeup...
  hPt_pfy -> SetMarkerSize (0.8);
  hPt_pfy -> SetMarkerColor(629);
  hPt_pfy -> SetLineWidth  (1);
  hPt_pfy -> GetXaxis() -> SetTitleOffset(1);

  hPt_pfy -> GetYaxis()    -> SetTitle("P_{T}^{raw} [Gev]");
  hPt_pfy -> GetYaxis()    -> SetTitleOffset(1);
  
  hPt_pfy -> Draw("");

  titleCanv = "Profile in P_{T}^{reco'd}";
  if (phiMin ==0  && phiMax == Pi  ) titleCanv += (" (Top Sectors)"   );
  if (phiMin ==Pi && phiMax == 2*Pi) titleCanv += (" (Bottom Sectors)");

  PrintIt(Pt_pfy, titleCanv.Data());

  if (isSave) {
    Pt_pfy->SaveAs("eps/"+runNumber+"/PtRawVsPtReco_pfy-"
		   +runNumber+"-"+phiCut+".eps");
    Pt_pfy->SaveAs("png/"+runNumber+"/PtRawVsPtReco_pfy-"
		   +runNumber+"-"+phiCut+".png");
  }
  
  
  // --------------------------------------------------------------------------- 
  // Q/P
  // ---------------------------------------------------------------------------  
  SetStyleh2(hQoverP, 1, 0.1, "q/P_{T}^{raw} [1/Gev]", 
	                      "q/P_{T}^{reco} [1/Gev]");
  TCanvas* QoverP 
    = new TCanvas("QoverP", "QoverP", 420, 500, 400, 400);
  
  hQoverP -> Draw("COLZ");
  SetMargin(QoverP);

  titleCanv = "q/P_{T}^{raw} vs q/P_{T}^{reco}";
  if (phiMin ==0  && phiMax == Pi  ) titleCanv += (" (Top Sectors)"   );
  if (phiMin ==Pi && phiMax == 2*Pi) titleCanv += (" (Bottom Sectors)");

  PrintIt(QoverP, titleCanv.Data());

  if (isSave) {
    QoverP->SaveAs("eps/"+runNumber+"/qPtRawVsqPtReco-"
		   +runNumber+"-"+phiCut+".eps");
    QoverP->SaveAs("png/"+runNumber+"/qPtRawVsqPtReco-"
		   +runNumber+"-"+phiCut+".png");
  }


  // --------------------------------------------------------------------------- 
  // Q/P Profile X
  // --------------------------------------------------------------------------- 
  hQoverP -> ProfileX();
  
  TCanvas* QoverP_pfx 
    = new TCanvas("QoverP_pfx", "q/P_{T} ProfileX", 840, 500, 400, 400);
  
  // Additional makeup...
  hQoverP_pfx -> SetMarkerSize (0.8);
  hQoverP_pfx -> SetMarkerColor(629);
  hQoverP_pfx -> SetLineWidth  (1);

  hQoverP_pfx -> GetYaxis()    -> SetTitle("q/P_{T}^{reco} [1/Gev]");
  hQoverP_pfx -> GetYaxis()    -> SetTitleOffset(1);

  hQoverP_pfx -> Draw("");

  titleCanv = "Profile in q/P_{T}^{raw}";
  if (phiMin ==0  && phiMax == Pi  ) titleCanv += (" (Top Sectors)"   );
  if (phiMin ==Pi && phiMax == 2*Pi) titleCanv += (" (Bottom Sectors)");

  PrintIt(QoverP_pfx, titleCanv.Data());

  if (isSave) {
    QoverP_pfx->SaveAs("eps/"+runNumber+"/qPtRawVsqPtReco_pfx-"
		       +runNumber+"-"+phiCut+".eps");
    QoverP_pfx->SaveAs("png/"+runNumber+"/qPtRawVsqPtReco_pfx-"
		       +runNumber+"-"+phiCut+".png");
  }


  // --------------------------------------------------------------------------- 
  // Q/P Profile Y
  // --------------------------------------------------------------------------- 
  hQoverP -> ProfileY();
  
  TCanvas* QoverP_pfy 
    = new TCanvas("QoverP_pfy", "q/P_{T} ProfileY", 1260, 500, 400, 400);
  
  // Additional makeup...
  hQoverP_pfy -> SetMarkerSize (0.8);
  hQoverP_pfy -> SetMarkerColor(629);
  hQoverP_pfy -> SetLineWidth  (1);
  hQoverP_pfy -> GetXaxis()    -> SetTitleOffset(1);

  hQoverP_pfy -> GetYaxis()    -> SetTitle("q/P_{T}^{raw} [1/Gev]");
  hQoverP_pfy -> GetYaxis()    -> SetTitleOffset(1);

  hQoverP_pfy -> Draw("");

  titleCanv = "Profile in q/P_{T}^{reco'd}";
  if (phiMin ==0  && phiMax == Pi  ) titleCanv += (" (Top Sectors)"   );
  if (phiMin ==Pi && phiMax == 2*Pi) titleCanv += (" (Bottom Sectors)");

  PrintIt(QoverP_pfy, titleCanv.Data());

  if (isSave) {
    QoverP_pfy->SaveAs("eps/"+runNumber+"/qPtRawVsqPtReco_pfy-"
		       +runNumber+"-"+phiCut+".eps");
    QoverP_pfy->SaveAs("png/"+runNumber+"/qPtRawVsqPtReco_pfy-"
		       +runNumber+"-"+phiCut+".png");
  }


  // --------------------------------------------------------------------------- 
  // ChrgValid
  // ---------------------------------------------------------------------------  
  SetStyleh1(hChargeValid   , 622, 1001, 1, "charge");
  TCanvas* ChargeValid 
    = new TCanvas("ChargeValid", "is charge valid?", 0, 1000, 400, 400);
  
  hChargeValid -> Draw();

  titleCanv = "Is charge valid?";
  if (phiMin ==0  && phiMax == Pi  ) titleCanv += (" (Top Sectors)"   );
  if (phiMin ==Pi && phiMax == 2*Pi) titleCanv += (" (Bottom Sectors)");

  PrintIt(ChargeValid, titleCanv.Data());

  if (isSave) {
    ChargeValid->SaveAs("eps/"+runNumber+"/chargeValid-"
			+runNumber+"-"+phiCut+".eps");
    ChargeValid->SaveAs("png/"+runNumber+"/chargeValid-"
			+runNumber+"-"+phiCut+".png");
  }


  // --------------------------------------------------------------------------- 
  // ChrgRaw Vs ChrgReco
  // ---------------------------------------------------------------------------  
  SetStyleh2(hChrgRawVsChrgReco, 1, 0.1, "q^{raw}", "q^{reco}");
  TCanvas* ChrgRawVsChrgReco 
    = new TCanvas("ChrgRawVsChrgReco", 
		  "ChrgRawVsChrgReco", 420, 1000, 400, 400);
  
  hChrgRawVsChrgReco -> Draw("COLZ");

  SetMargin(ChrgRawVsChrgReco);

  titleCanv = "q^{raw} vs q^{reco}";
  if (phiMin ==0  && phiMax == Pi  ) titleCanv += (" (Top Sectors)"   );
  if (phiMin ==Pi && phiMax == 2*Pi) titleCanv += (" (Bottom Sectors)");

  PrintIt(ChrgRawVsChrgReco , titleCanv.Data());

  if (isSave) { 
    ChrgRawVsChrgReco->SaveAs("eps/"+runNumber+"/qRawVsqReco-"
			      +runNumber+"-"+phiCut+".eps");
    ChrgRawVsChrgReco->SaveAs("png/"+runNumber+"/qRawVsqReco-"
			      +runNumber+"-"+phiCut+".png");
  }

  cout << "PtRawVsPtReco has "
       <<  hPt -> GetEntries() 
       << " entries\n";
       
} // end 
