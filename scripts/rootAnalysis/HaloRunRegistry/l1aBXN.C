// ------------------------------------------------------------------------------
// TTree
// ------------------------------------------------------------------------------
TTree* Ntuple;

// ------------------------------------------------------------------------------
// histograms
// ------------------------------------------------------------------------------
TH1F* hL1ABXN; // L1A BXN from GT
TH1F* hL1ABXN_Halo; // L1A BXN from GT for Halo
TH1F* hL1ABXN_Coinc; // L1A BXN from GT for Coinc
TH1F* hL1ABXN_Singles; // L1A BXN from GT for Singles
TH1F* hMode; // Mode

// ------------------------------------------------------------------------------
// variables
// ------------------------------------------------------------------------------
double Pi  = TMath::Pi();

const int MAXCSCTFTR = 60;
const int MAXCSCTFLCTSTR = 4;
const int MAXCSCTFSPS = 12;

void l1aBXN(TString runNumber  = "123893",
	    int     printLevel =       0 ,
	    bool    isSave     =    !true) {


  
  gROOT->Clear();
  gInterpreter->LoadMacro("/afs/cern.ch/user/d/digiovan/scripts/Init/"
			  "modifiedStyle.C");
  
  //--------------------------------------------------------------------------
  // Construct the chain of input files
  //--------------------------------------------------------------------------
  TString chainMacro("scripts/HaloRunRegistry/chains/");
  chainMacro += "collsionsChain";
  chainMacro += runNumber;
  chainMacro += ".C";

  cout << "Loading chain: " << chainMacro.Data() << endl;
  gInterpreter->ExecuteMacro(chainMacro.Data());
  
  // get the Nutple  
  TChain* Ntuple  = collisionsChain;
  cout << "Ntuple->GetEntries(): " << Ntuple->GetEntries() << endl;
  
  //----------------------------------------------------------------------------
  // Access the needed variables
  //----------------------------------------------------------------------------
  int    Run, Event, Bx;
  
  int    csctf_trSize;
  int    csctf_trMode[MAXCSCTFTR];

  Ntuple -> SetBranchAddress("Run", &Run );
  Ntuple -> SetBranchAddress("Event", &Event );
  Ntuple -> SetBranchAddress("Bx", &Bx );

  Ntuple -> SetBranchAddress("csctf_trSize", &csctf_trSize);
  Ntuple -> SetBranchAddress("csctf_trMode", &csctf_trMode);

  // --------------------------------------------------------------------------- 
  // Histograms
  // ---------------------------------------------------------------------------  
  hL1ABXN = new TH1F("hL1ABXN", "", 3600, 0, 3600);
  hL1ABXN_Halo = new TH1F("hL1ABXN_Halo", "", 3600, 0, 3600);
  hL1ABXN_Coinc = new TH1F("hL1ABXN_Coinc", "", 3600, 0, 3600);
  hL1ABXN_Singles = new TH1F("hL1ABXN_Singles", "", 3600, 0, 3600);
  hMode = new TH1F("hMode", "", 16, 0, 16);

  // --------------------------------------------------------------------------- 
  // Loop over the events
  // ---------------------------------------------------------------------------  
  for (int iEvt=0; iEvt < Ntuple->GetEntries(); iEvt++) {
   
    Ntuple->GetEntry(iEvt);
     
    for (int iTrk=0; iTrk<csctf_trSize; iTrk++) {
      hL1ABXN -> Fill(Bx);

      switch (csctf_trMode[iTrk]) {
      case 15:
	hL1ABXN_Halo->Fill(Bx);
	break;
      case 11:
	hL1ABXN_Singles->Fill(Bx);
	break;
      default:
	hL1ABXN_Coinc->Fill(Bx);
	break;
      }
      
      if(Bx>48 && Bx<53)
	hMode -> Fill(csctf_trMode[iTrk]);
    }
    
  }// end Loop evts
  if (printLevel > 0)
    printf("\n----------------------------------------"
	   "----------------------------------------\n");
   
  // --------------------------------------------------------------------------- 
  // hL1ABXN
  // ---------------------------------------------------------------------------  
  SetStyleh1(hL1ABXN, 629, 3001, 1, "L1A (BX)");
  TCanvas* L1ABXN = new TCanvas("L1ABXN", "", 0, 0, 400, 400);

  hL1ABXN -> Draw();

  PrintIt(L1ABXN, "L1A BXN Distribution");

  if (isSave) {
    L1ABXN -> SaveAs("png/L1ABXN.png");
    L1ABXN -> SaveAs("eps/L1ABXN.eps");
    L1ABXN -> SaveAs("rootPlot/L1ABXN.root");
  }

  // --------------------------------------------------------------------------- 
  // hL1ABXN_Halo
  // ---------------------------------------------------------------------------  
  SetStyleh1(hL1ABXN_Halo, 629, 3001, 1, "L1A (BX)");
  TCanvas* L1ABXN_Halo = new TCanvas("L1ABXN_Halo", "", 420, 0, 400, 400);

  hL1ABXN_Halo -> Draw();

  PrintIt(L1ABXN_Halo, "L1A BXN Distribution (HALO)");

  if (isSave) {
    L1ABXN_Halo -> SaveAs("png/L1ABXN_Halo.png");
    L1ABXN_Halo -> SaveAs("eps/L1ABXN_Halo.eps");
    L1ABXN_Halo -> SaveAs("rootPlot/L1ABXN_Halo.root");
  }

  // --------------------------------------------------------------------------- 
  // hL1ABXN_Coinc
  // ---------------------------------------------------------------------------  
  SetStyleh1(hL1ABXN_Coinc, 629, 3001, 1, "L1A (BX)");
  TCanvas* L1ABXN_Coinc = new TCanvas("L1ABXN_Coinc", "", 840, 0, 400, 400);

  hL1ABXN_Coinc -> Draw();

  PrintIt(L1ABXN_Coinc, "L1A BXN Distribution (COINC)");

  if (isSave) {
    L1ABXN_Coinc -> SaveAs("png/L1ABXN_Coinc.png");
    L1ABXN_Coinc -> SaveAs("eps/L1ABXN_Coinc.eps");
    L1ABXN_Coinc -> SaveAs("rootPlot/L1ABXN_Coinc.root");
  }

  // --------------------------------------------------------------------------- 
  // hL1ABXN_Singles
  // ---------------------------------------------------------------------------  
  SetStyleh1(hL1ABXN_Singles, 629, 3001, 1, "L1A (BX)");
  TCanvas* L1ABXN_Singles = new TCanvas("L1ABXN_Singles", "", 1260, 0, 400, 400);

  hL1ABXN_Singles -> Draw();

  PrintIt(L1ABXN_Singles, "L1A BXN Distribution (SINGLES)");

  if (isSave) {
    L1ABXN_Singles -> SaveAs("png/L1ABXN_Singles.png");
    L1ABXN_Singles -> SaveAs("eps/L1ABXN_Singles.eps");
    L1ABXN_Singles -> SaveAs("rootPlot/L1ABXN_Singles.root");
  }

  // --------------------------------------------------------------------------- 
  // hMode
  // ---------------------------------------------------------------------------  
  SetStyleh1(hMode, 629, 3001, 1, "L1A (BX)");
  TCanvas* Mode = new TCanvas("Mode", "", 0, 500, 400, 400);

  hMode -> Draw();

  PrintIt(Mode, "Mode Distribution");
       
} // end 
