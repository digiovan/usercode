// ------------------------------------------------------------------------------
// TTree
// ------------------------------------------------------------------------------
TTree* Ntuple;

// ------------------------------------------------------------------------------
// histograms
// ------------------------------------------------------------------------------
TH1F* hL1ABXN; // L1A BXN from GT
TH1F* hL1ABXN_Halo; // L1A BXN from GT for Halo

TH1F* hNumLCTsHalo; // Num of LCTs forming the Halo (expected 2)
TH1F* hHaloCombo; // Halo Extrapolation Combo
TH1F* hHaloComboBkg; // Halo Extrapolation Combo Bkg

TH1F* hDPhiME12; // Phi_ME1-Phi_ME2
TH1F* hDPhiME13; // Phi_ME1-Phi_ME3
TH1F* hDPhiME112; // Phi_ME1/1-Phi_ME2
TH1F* hDPhiME113; // Phi_ME1/1-Phi_ME3

TH1F* hDPhiME12Bkg; // Phi_ME1-Phi_ME2 Bkg
TH1F* hDPhiME13Bkg; // Phi_ME1-Phi_ME3 Bkg
TH1F* hDPhiME112Bkg; // Phi_ME1/1-Phi_ME2 Bkg
TH1F* hDPhiME113Bkg; // Phi_ME1/1-Phi_ME3 Bkg

TH1F* hDEtaME12; // Eta_ME1-Eta_ME2
TH1F* hDEtaME13; // Eta_ME1-Eta_ME3
TH1F* hDEtaME112; // Eta_ME1/1-Eta_ME2
TH1F* hDEtaME113; // Eta_ME1/1-Eta_ME3

TH1F* hDEtaME12Bkg; // Eta_ME1-Eta_ME2 Bkg
TH1F* hDEtaME13Bkg; // Eta_ME1-Eta_ME3 Bkg
TH1F* hDEtaME112Bkg; // Eta_ME1/1-Eta_ME2 Bkg
TH1F* hDEtaME113Bkg; // Eta_ME1/1-Eta_ME3 Bkg

// ------------------------------------------------------------------------------
// legend
// ------------------------------------------------------------------------------
TLegend *tl;

// ------------------------------------------------------------------------------
// variables
// ------------------------------------------------------------------------------
double Pi  = TMath::Pi();
double PhiStep = (62*Pi/180)/4096;
double EtaStep = 1.6/128; //(2.5-0.9)/(2^7)

const int MAXCSCTFTR = 60;
const int MAXCSCTFLCTSTR = 4;
const int MAXCSCTFSPS = 12;

// ------------------------------------------------------------------------------
// macro
// ------------------------------------------------------------------------------
void bkgSeparation(TString runNumber  = "123893",
		   int     printLevel =        0,
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
  int Run, Event, Bx;
  
  int csctf_trSize;
  int csctf_trMode[MAXCSCTFTR];
  int csctf_trNumLCTs[MAXCSCTFTR]; 

  int csctf_lctStation[MAXCSCTFTR][MAXCSCTFLCTSTR];
  int csctf_lctRing[MAXCSCTFTR][MAXCSCTFLCTSTR];
  
  // note: the SPs return them in bits
  int csctf_lctglobalPhi[MAXCSCTFTR][MAXCSCTFLCTSTR];
  int csctf_lctglobalEta[MAXCSCTFTR][MAXCSCTFLCTSTR];

  Ntuple -> SetBranchAddress("Run", &Run );
  Ntuple -> SetBranchAddress("Event", &Event );
  Ntuple -> SetBranchAddress("Bx", &Bx );

  Ntuple -> SetBranchAddress("csctf_trSize", &csctf_trSize);
  Ntuple -> SetBranchAddress("csctf_trMode", &csctf_trMode);
  Ntuple -> SetBranchAddress("csctf_trNumLCTs", &csctf_trNumLCTs);

  Ntuple -> SetBranchAddress("csctf_lctStation", &csctf_lctStation);
  Ntuple -> SetBranchAddress("csctf_lctRing", &csctf_lctRing);
  Ntuple -> SetBranchAddress("csctf_lctglobalPhi", &csctf_lctglobalPhi);
  Ntuple -> SetBranchAddress("csctf_lctglobalEta", &csctf_lctglobalEta);

  // --------------------------------------------------------------------------- 
  // Histograms
  // ---------------------------------------------------------------------------  
  hL1ABXN = new TH1F("hL1ABXN", "", 3600, 0, 3600);
  hL1ABXN_Halo = new TH1F("hL1ABXN_Halo", "", 3600, 0, 3600);

  hNumLCTsHalo = new TH1F("hNumLCTsHalo", "", 4, 0, 4);
  hHaloCombo = new TH1F("hHaloCombo", "", 4, 1, 5);
  hHaloComboBkg = new TH1F("hHaloComboBkg", "", 4, 1, 5);

  hDPhiME12 = new TH1F("hDPhiME12", "", 601, -300*PhiStep, 300*PhiStep);
  hDPhiME12Bkg = new TH1F("hDPhiME12Bkg", "", 601, -300*PhiStep, 300*PhiStep);
  hDPhiME13 = new TH1F("hDPhiME13", "", 601, -300*PhiStep, 300*PhiStep);
  hDPhiME13Bkg = new TH1F("hDPhiME13Bkg", "", 601, -300*PhiStep, 300*PhiStep);
  hDPhiME112 = new TH1F("hDPhiME112", "", 601, -300*PhiStep, 300*PhiStep);
  hDPhiME112Bkg = new TH1F("hDPhiME112Bkg", "", 601, -300*PhiStep, 300*PhiStep);
  hDPhiME113 = new TH1F("hDPhiME113", "", 601, -300*PhiStep, 300*PhiStep);
  hDPhiME113Bkg = new TH1F("hDPhiME113Bkg", "", 601, -300*PhiStep, 300*PhiStep);

  hDEtaME12 = new TH1F("hDEtaME12", "", 44, 0, 44*EtaStep);
  hDEtaME12Bkg = new TH1F("hDEtaME12Bkg", "", 44, 0, 44*EtaStep);
  hDEtaME13 = new TH1F("hDEtaME13", "", 44, 0, 44*EtaStep);
  hDEtaME13Bkg = new TH1F("hDEtaME13Bkg", "", 44, 0, 44*EtaStep);
  hDEtaME112 = new TH1F("hDEtaME112", "", 44, 0, 44*EtaStep);
  hDEtaME112Bkg = new TH1F("hDEtaME112Bkg", "", 44, 0, 44*EtaStep);
  hDEtaME113 = new TH1F("hDEtaME113", "", 44, 0, 44*EtaStep);
  hDEtaME113Bkg = new TH1F("hDEtaME113Bkg", "", 44, 0, 44*EtaStep);


  // --------------------------------------------------------------------------- 
  // Loop over the events
  // ---------------------------------------------------------------------------  
  for (int iEvt=0; iEvt < Ntuple->GetEntries(); iEvt++) {
   
    Ntuple->GetEntry(iEvt);
     
    for (int iTrk=0; iTrk<csctf_trSize; iTrk++) {
      hL1ABXN -> Fill(Bx);

      if (csctf_trMode[iTrk] == 15) { 
	hL1ABXN_Halo->Fill(Bx);

	// haloCombo: 1 -> ME1-ME2
        //            2 -> ME1-ME3     
	//            3 -> ME1/1-ME2     
	//            4 -> ME1/1-ME3     
	int haloCombo = -999; 
      
	// protection: csctf does not preserve info about the lcts belonging
	// to the track. This is done at the software level, so any bug could
	// cause this to happen...
	if (csctf_trNumLCTs[iTrk] != 2) {
	  cout << "Weird: the halo candidate does not have LCTs... ";
	  cout << "Bx: " << Bx << endl;
	  continue;
	}

	if (printLevel) {
	  
	  // loop over lcts
	  for (int iLCT=0; iLCT < csctf_trNumLCTs[iTrk]; iLCT++) {
	    cout << "csctf_lctStation[" << iTrk << "][" << iLCT << "]: " 
		 << csctf_lctStation[iTrk][iLCT];
	    cout << "csctf_lctRing[" << iTrk << "][" << iLCT << "]: " 
		 << csctf_lctRing[iTrk][iLCT] << endl << endl;
	  }// end loop over lcts
	} 


	if (csctf_lctStation[iTrk][1] == 2)
	  if (csctf_lctRing[iTrk][0] == 1) haloCombo = 3;
	  else 		                   haloCombo = 1;		     

	if (csctf_lctStation[iTrk][1] == 3)
	  if (csctf_lctRing[iTrk][0] == 1) haloCombo = 4;
	  else 		                   haloCombo = 2;		     


	double dphiHalo = (csctf_lctglobalPhi[iTrk][0]-csctf_lctglobalPhi[iTrk][1])*PhiStep; 
	double detaHalo = fabs(csctf_lctglobalEta[iTrk][0]-csctf_lctglobalEta[iTrk][1])*EtaStep; 

	if(Bx>48 && Bx<54) {
	  hNumLCTsHalo -> Fill(csctf_trNumLCTs[iTrk]);
	  hHaloCombo -> Fill(haloCombo);

	  if (haloCombo == 1) hDPhiME12 -> Fill(dphiHalo);
	  if (haloCombo == 2) hDPhiME13 -> Fill(dphiHalo);
	  if (haloCombo == 3) hDPhiME112 -> Fill(dphiHalo);
	  if (haloCombo == 4) hDPhiME113 -> Fill(dphiHalo);

	  if (haloCombo == 1) hDEtaME12 -> Fill(detaHalo);
	  if (haloCombo == 2) hDEtaME13 -> Fill(detaHalo);
	  if (haloCombo == 3) hDEtaME112 -> Fill(detaHalo);
	  if (haloCombo == 4) hDEtaME113 -> Fill(detaHalo);

	}// signal bx
	else{
	  hHaloComboBkg -> Fill(haloCombo);

	  if (haloCombo == 1) hDPhiME12Bkg -> Fill(dphiHalo);
	  if (haloCombo == 2) hDPhiME13Bkg -> Fill(dphiHalo);
	  if (haloCombo == 3) hDPhiME112Bkg -> Fill(dphiHalo);
	  if (haloCombo == 4) hDPhiME113Bkg -> Fill(dphiHalo);

	  if (haloCombo == 1) hDEtaME12Bkg -> Fill(detaHalo);
	  if (haloCombo == 2) hDEtaME13Bkg -> Fill(detaHalo);
	  if (haloCombo == 3) hDEtaME112Bkg -> Fill(detaHalo);
	  if (haloCombo == 4) hDEtaME113Bkg -> Fill(detaHalo);
	  
	}// bkg


      }//only halo
    
    }//trk loop
    
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
  // hNumLCTsHalo
  // ---------------------------------------------------------------------------  
  SetStyleh1(hNumLCTsHalo, 629, 3001, 1, "# LCTs");
  TCanvas* NumLCTsHalo = new TCanvas("NumLCTsHalo", "", 840, 0, 400, 400);

  hNumLCTsHalo -> Draw();

  PrintIt(NumLCTsHalo, "Number of LCTs forming the Halo candidate");

  if (isSave) {
    NumLCTsHalo -> SaveAs("png/NumLCTsHalo.png");
    NumLCTsHalo -> SaveAs("eps/NumLCTsHalo.eps");
    NumLCTsHalo -> SaveAs("rootPlot/NumLCTsHalo.root");
  }

  // --------------------------------------------------------------------------- 
  // hHaloCombo
  // ---------------------------------------------------------------------------  
  SetStyleh1(hHaloCombo, 629, 0, 3, "LCTs Combo");
  SetStyleh1(hHaloComboBkg, 597, 0, 3, "LCTs Combo");

  hHaloCombo -> GetXaxis() -> SetBinLabel(1,"ME1-ME2");
  hHaloCombo -> GetXaxis() -> SetBinLabel(2,"ME1-ME3");
  hHaloCombo -> GetXaxis() -> SetBinLabel(3,"ME1/1-ME2");
  hHaloCombo -> GetXaxis() -> SetBinLabel(4,"ME1/1-ME3");

  hHaloComboBkg -> SetLineStyle(2);

  hHaloComboBkg -> GetXaxis() -> SetBinLabel(1,"ME1-ME2");
  hHaloComboBkg -> GetXaxis() -> SetBinLabel(2,"ME1-ME3");
  hHaloComboBkg -> GetXaxis() -> SetBinLabel(3,"ME1/1-ME2");
  hHaloComboBkg -> GetXaxis() -> SetBinLabel(4,"ME1/1-ME3");

  TCanvas* HaloCombo = new TCanvas("HaloCombo", "", 1260, 0, 400, 400);

  if (hHaloCombo->GetEntries() == 0 ||
      hHaloComboBkg->GetEntries() == 0)
    cout << "hHaloCombo->GetEntries(): " << hHaloCombo->GetEntries() << " -- "
	 << "hHaloComboBkg->GetEntries(): " << hHaloComboBkg->GetEntries() << endl;
  
  if (hHaloCombo->GetMaximum() > hHaloComboBkg->GetMaximum()) {
    if (hHaloCombo->GetEntries())    hHaloCombo -> DrawNormalized("", 1);
    if (hHaloComboBkg->GetEntries()) hHaloComboBkg -> DrawNormalized("same", 1);
  }
  else {
    if (hHaloComboBkg->GetEntries()) hHaloComboBkg -> DrawNormalized("", 1);
    if (hHaloCombo->GetEntries())    hHaloCombo -> DrawNormalized("same", 1);
  }

  PrintIt(HaloCombo, "Halo Candidate Combos");

  TString collisions = (" COLLISIONS:");
  collisions += hHaloCombo->GetEntries();
  collisions += " entries";

  TString cosmics = (" COSMICS:");
  cosmics += hHaloComboBkg->GetEntries();
  cosmics += " entries";

  tl = SetLegend(0.50, 0.73, 0.64, 0.87);
  tl->AddEntry(hHaloCombo, collisions, "l");
  tl->AddEntry(hHaloComboBkg, cosmics, "l");

  tl -> Draw("same");

  if (isSave) {
    HaloCombo -> SaveAs("png/HaloCombo.png");
    HaloCombo -> SaveAs("eps/HaloCombo.eps");
    HaloCombo -> SaveAs("rootPlot/HaloCombo.root");
  }


  // --------------------------------------------------------------------------- 
  // hDPhiME12
  // ---------------------------------------------------------------------------  
  SetStyleh1(hDPhiME12, 629, 3001, 3, "#Delta #phi_{ME1-ME2}");
  SetStyleh1(hDPhiME12Bkg, 597, 3001, 3, "#Delta #phi_{ME1-ME2}");

  hDPhiME12Bkg -> SetLineStyle(2);

  TCanvas* DPhiME12 = new TCanvas("DPhiME12", "", 0, 500, 400, 400);

  Double_t normColl =  hDPhiME12->GetEntries();
  if (normColl) hDPhiME12->Scale(1/normColl);
  else cout << "hDPhiME12->GetEntries(): " << normColl << endl; 

  Double_t normCosm =  hDPhiME12Bkg->GetEntries();
  if (normCosm) hDPhiME12Bkg->Scale(1/normCosm);
  else cout << "hDPhiME12Bkg->GetEntries(): " << normCosm << endl; 

  if (hDPhiME12->GetMaximum() > hDPhiME12Bkg->GetMaximum()) {
    if (hDPhiME12->GetEntries())    hDPhiME12 -> Draw();
    if (hDPhiME12Bkg->GetEntries()) hDPhiME12Bkg -> Draw("same");
  }
  else {
    if (hDPhiME12Bkg->GetEntries()) hDPhiME12Bkg -> Draw();
    if (hDPhiME12->GetEntries())    hDPhiME12 -> Draw("same");
  }

  PrintIt(DPhiME12, "Halo #Delta #phi (ME1-ME2 extrapolation)");

  TString collisions = (" COLLISIONS:");
  collisions += hDPhiME12->GetEntries();
  collisions += " entries";

  TString cosmics = (" COSMICS:");
  cosmics += hDPhiME12Bkg->GetEntries();
  cosmics += " entries";

  tl = SetLegend(0.50, 0.73, 0.64, 0.87);
  tl->AddEntry(hDPhiME12, collisions, "l");
  tl->AddEntry(hDPhiME12Bkg, cosmics, "l");

  tl -> Draw("same");

  if (isSave) {
    DPhiME12 -> SaveAs("png/DPhiME12.png");
    DPhiME12 -> SaveAs("eps/DPhiME12.eps");
    DPhiME12 -> SaveAs("rootPlot/DPhiME12.root");
  }


  // --------------------------------------------------------------------------- 
  // hDPhiME13
  // ---------------------------------------------------------------------------  
  SetStyleh1(hDPhiME13, 629, 3001, 3, "#Delta #phi_{ME1-ME3}");
  SetStyleh1(hDPhiME13Bkg, 597, 3001, 3, "#Delta #phi_{ME1-ME3}");

  hDPhiME13Bkg -> SetLineStyle(2);

  TCanvas* DPhiME13 = new TCanvas("DPhiME13", "", 420, 500, 400, 400);

  Double_t normColl =  hDPhiME13->GetEntries();
  if (normColl) hDPhiME13->Scale(1/normColl);
  else cout << "hDPhiME13->GetEntries(): " << normColl << endl; 

  Double_t normCosm =  hDPhiME13Bkg->GetEntries();
  if (normCosm) hDPhiME13Bkg->Scale(1/normCosm);
  else cout << "hDPhiME13Bkg->GetEntries(): " << normCosm << endl; 

  if (hDPhiME13->GetMaximum() > hDPhiME13Bkg->GetMaximum()) {
    if (hDPhiME13->GetEntries())    hDPhiME13 -> Draw();
    if (hDPhiME13Bkg->GetEntries()) hDPhiME13Bkg -> Draw("same");
  }
  else {
    if (hDPhiME13Bkg->GetEntries()) hDPhiME13Bkg -> Draw();
    if (hDPhiME13->GetEntries())    hDPhiME13 -> Draw("same");
  }

  PrintIt(DPhiME13, "Halo #Delta #phi (ME1 #rightarrow ME3 extrapolation)");

  TString collisions = (" COLLISIONS:");
  collisions += hDPhiME13->GetEntries();
  collisions += " entries";

  TString cosmics = (" COSMICS:");
  cosmics += hDPhiME13Bkg->GetEntries();
  cosmics += " entries";

  tl = SetLegend(0.50, 0.73, 0.64, 0.87);
  tl->AddEntry(hDPhiME13, collisions, "l");
  tl->AddEntry(hDPhiME13Bkg, cosmics, "l");

  tl -> Draw("same");

  if (isSave) {
    DPhiME13 -> SaveAs("png/DPhiME13.png");
    DPhiME13 -> SaveAs("eps/DPhiME13.eps");
    DPhiME13 -> SaveAs("rootPlot/DPhiME13.root");
  }


  // --------------------------------------------------------------------------- 
  // hDPhiME112
  // ---------------------------------------------------------------------------  
  SetStyleh1(hDPhiME112, 629, 3001, 3, "#Delta #phi_{ME1/1-ME2}");
  SetStyleh1(hDPhiME112Bkg, 597, 3001, 3, "#Delta #phi_{ME1/1-ME2}");

  hDPhiME112Bkg -> SetLineStyle(2);

  TCanvas* DPhiME112 = new TCanvas("DPhiME112", "", 840, 500, 400, 400);

  Double_t normColl =  hDPhiME112->GetEntries();
  if (normColl) hDPhiME112->Scale(1/normColl);
  else cout << "hDPhiME112->GetEntries(): " << normColl << endl; 

  Double_t normCosm =  hDPhiME112Bkg->GetEntries();
  if (normCosm) hDPhiME112Bkg->Scale(1/normCosm);
  else cout << "hDPhiME112Bkg->GetEntries(): " << normCosm << endl; 

  if (hDPhiME112->GetMaximum() > hDPhiME112Bkg->GetMaximum()) {
    if (hDPhiME112->GetEntries())    hDPhiME112 -> Draw();
    if (hDPhiME112Bkg->GetEntries()) hDPhiME112Bkg -> Draw("same");
  }
  else {
    if (hDPhiME112Bkg->GetEntries()) hDPhiME112Bkg -> Draw();
    if (hDPhiME112->GetEntries())    hDPhiME112 -> Draw("same");
  }

  PrintIt(DPhiME112, "Halo #Delta #phi (ME1/1 #rightarrow ME2 extrapolation)");

  TString collisions = (" COLLISIONS:");
  collisions += hDPhiME112->GetEntries();
  collisions += " entries";

  TString cosmics = (" COSMICS:");
  cosmics += hDPhiME112Bkg->GetEntries();
  cosmics += " entries";

  tl = SetLegend(0.50, 0.73, 0.64, 0.87);
  tl->AddEntry(hDPhiME112, collisions, "l");
  tl->AddEntry(hDPhiME112Bkg, cosmics, "l");

  tl -> Draw("same");

  if (isSave) {
    DPhiME112 -> SaveAs("png/DPhiME112.png");
    DPhiME112 -> SaveAs("eps/DPhiME112.eps");
    DPhiME112 -> SaveAs("rootPlot/DPhiME112.root");
  }


  // --------------------------------------------------------------------------- 
  // hDPhiME113
  // ---------------------------------------------------------------------------  
  SetStyleh1(hDPhiME113, 629, 3001, 3, "#Delta #phi_{ME1/1-ME3}");
  SetStyleh1(hDPhiME113Bkg, 597, 3001, 3, "#Delta #phi_{ME1/1-ME3}");

  hDPhiME113Bkg -> SetLineStyle(2);

  TCanvas* DPhiME113 = new TCanvas("DPhiME113", "", 1260, 500, 400, 400);

  Double_t normColl =  hDPhiME113->GetEntries();
  if (normColl) hDPhiME113->Scale(1/normColl);
  else cout << "hDPhiME113->GetEntries(): " << normColl << endl; 

  Double_t normCosm =  hDPhiME113Bkg->GetEntries();
  if (normCosm) hDPhiME113Bkg->Scale(1/normCosm);
  else cout << "hDPhiME113Bkg->GetEntries(): " << normCosm << endl; 

  if (hDPhiME113->GetMaximum() > hDPhiME113Bkg->GetMaximum()) {
    if (hDPhiME113->GetEntries())    hDPhiME113 -> Draw();
    if (hDPhiME113Bkg->GetEntries()) hDPhiME113Bkg -> Draw("same");
  }
  else {
    if (hDPhiME113Bkg->GetEntries()) hDPhiME113Bkg -> Draw();
    if (hDPhiME113->GetEntries())    hDPhiME113 -> Draw("same");
  }

  PrintIt(DPhiME113, "Halo #Delta #phi (ME1/1 #rightarrow ME3 extrapolation)");

  TString collisions = (" COLLISIONS:");
  collisions += hDPhiME113->GetEntries();
  collisions += " entries";

  TString cosmics = (" COSMICS:");
  cosmics += hDPhiME113Bkg->GetEntries();
  cosmics += " entries";

  tl = SetLegend(0.50, 0.73, 0.64, 0.87);
  tl->AddEntry(hDPhiME113, collisions, "l");
  tl->AddEntry(hDPhiME113Bkg, cosmics, "l");

  tl -> Draw("same");

  if (isSave) {
    DPhiME113 -> SaveAs("png/DPhiME113.png");
    DPhiME113 -> SaveAs("eps/DPhiME113.eps");
    DPhiME113 -> SaveAs("rootPlot/DPhiME113.root");
  }



  // --------------------------------------------------------------------------- 
  // hDEtaME12
  // ---------------------------------------------------------------------------  
  SetStyleh1(hDEtaME12, 629, 3001, 3, "#Delta #eta_{ME1-ME2}");
  SetStyleh1(hDEtaME12Bkg, 597, 3001, 3, "#Delta #eta_{ME1-ME2}");

  hDEtaME12Bkg -> SetLineStyle(2);

  TCanvas* DEtaME12 = new TCanvas("DEtaME12", "", 0, 1000, 400, 400);

  Double_t normColl =  hDEtaME12->GetEntries();
  if (normColl) hDEtaME12->Scale(1/normColl);
  else cout << "hDEtaME12->GetEntries(): " << normColl << endl; 

  Double_t normCosm =  hDEtaME12Bkg->GetEntries();
  if (normCosm) hDEtaME12Bkg->Scale(1/normCosm);
  else cout << "hDEtaME12Bkg->GetEntries(): " << normCosm << endl; 

  if (hDEtaME12->GetMaximum() > hDEtaME12Bkg->GetMaximum()) {
    if (hDEtaME12->GetEntries())    hDEtaME12 -> Draw();
    if (hDEtaME12Bkg->GetEntries()) hDEtaME12Bkg -> Draw("same");
  }
  else {
    if (hDEtaME12Bkg->GetEntries()) hDEtaME12Bkg -> Draw();
    if (hDEtaME12->GetEntries())    hDEtaME12 -> Draw("same");
  }

  PrintIt(DEtaME12, "Halo #Delta #eta (ME1-ME2 extrapolation)");

  TString collisions = (" COLLISIONS:");
  collisions += hDEtaME12->GetEntries();
  collisions += " entries";

  TString cosmics = (" COSMICS:");
  cosmics += hDEtaME12Bkg->GetEntries();
  cosmics += " entries";

  tl = SetLegend(0.50, 0.73, 0.64, 0.87);
  tl->AddEntry(hDEtaME12, collisions, "l");
  tl->AddEntry(hDEtaME12Bkg, cosmics, "l");

  tl -> Draw("same");

  if (isSave) {
    DEtaME12 -> SaveAs("png/DEtaME12.png");
    DEtaME12 -> SaveAs("eps/DEtaME12.eps");
    DEtaME12 -> SaveAs("rootPlot/DEtaME12.root");
  }


  // --------------------------------------------------------------------------- 
  // hDEtaME13
  // ---------------------------------------------------------------------------  
  SetStyleh1(hDEtaME13, 629, 3001, 3, "#Delta #eta_{ME1-ME3}");
  SetStyleh1(hDEtaME13Bkg, 597, 3001, 3, "#Delta #eta_{ME1-ME3}");

  hDEtaME13Bkg -> SetLineStyle(2);

  TCanvas* DEtaME13 = new TCanvas("DEtaME13", "", 420, 1000, 400, 400);

  Double_t normColl =  hDEtaME13->GetEntries();
  if (normColl) hDEtaME13->Scale(1/normColl);
  else cout << "hDEtaME13->GetEntries(): " << normColl << endl; 

  Double_t normCosm =  hDEtaME13Bkg->GetEntries();
  if (normCosm) hDEtaME13Bkg->Scale(1/normCosm);
  else cout << "hDEtaME13Bkg->GetEntries(): " << normCosm << endl; 

  if (hDEtaME13->GetMaximum() > hDEtaME13Bkg->GetMaximum()) {
    if (hDEtaME13->GetEntries())    hDEtaME13 -> Draw();
    if (hDEtaME13Bkg->GetEntries()) hDEtaME13Bkg -> Draw("same");
  }
  else {
    if (hDEtaME13Bkg->GetEntries()) hDEtaME13Bkg -> Draw();
    if (hDEtaME13->GetEntries())    hDEtaME13 -> Draw("same");
  }

  PrintIt(DEtaME13, "Halo #Delta #eta (ME1-ME3 extrapolation)");

  TString collisions = (" COLLISIONS:");
  collisions += hDEtaME13->GetEntries();
  collisions += " entries";

  TString cosmics = (" COSMICS:");
  cosmics += hDEtaME13Bkg->GetEntries();
  cosmics += " entries";

  tl = SetLegend(0.50, 0.73, 0.64, 0.87);
  tl->AddEntry(hDEtaME13, collisions, "l");
  tl->AddEntry(hDEtaME13Bkg, cosmics, "l");

  tl -> Draw("same");

  if (isSave) {
    DEtaME13 -> SaveAs("png/DEtaME13.png");
    DEtaME13 -> SaveAs("eps/DEtaME13.eps");
    DEtaME13 -> SaveAs("rootPlot/DEtaME13.root");
  }


  // --------------------------------------------------------------------------- 
  // hDEtaME112
  // ---------------------------------------------------------------------------  
  SetStyleh1(hDEtaME112, 629, 3001, 3, "#Delta #eta_{ME1/1-ME2}");
  SetStyleh1(hDEtaME112Bkg, 597, 3001, 3, "#Delta #eta_{ME1/1-ME2}");

  hDEtaME112Bkg -> SetLineStyle(2);

  TCanvas* DEtaME112 = new TCanvas("DEtaME112", "", 840, 1000, 400, 400);

  Double_t normColl =  hDEtaME112->GetEntries();
  if (normColl) hDEtaME112->Scale(1/normColl);
  else cout << "hDEtaME112->GetEntries(): " << normColl << endl; 

  Double_t normCosm =  hDEtaME112Bkg->GetEntries();
  if (normCosm) hDEtaME112Bkg->Scale(1/normCosm);
  else cout << "hDEtaME112Bkg->GetEntries(): " << normCosm << endl; 

  if (hDEtaME112->GetMaximum() > hDEtaME112Bkg->GetMaximum()) {
    if (hDEtaME112->GetEntries())    hDEtaME112 -> Draw();
    if (hDEtaME112Bkg->GetEntries()) hDEtaME112Bkg -> Draw("same");
  }
  else {
    if (hDEtaME112Bkg->GetEntries()) hDEtaME112Bkg -> Draw();
    if (hDEtaME112->GetEntries())    hDEtaME112 -> Draw("same");
  }

  PrintIt(DEtaME112, "Halo #Delta #eta (ME1/1-ME2 extrapolation)");

  TString collisions = (" COLLISIONS:");
  collisions += hDEtaME112->GetEntries();
  collisions += " entries";

  TString cosmics = (" COSMICS:");
  cosmics += hDEtaME112Bkg->GetEntries();
  cosmics += " entries";

  tl = SetLegend(0.50, 0.73, 0.64, 0.87);
  tl->AddEntry(hDEtaME112, collisions, "l");
  tl->AddEntry(hDEtaME112Bkg, cosmics, "l");

  tl -> Draw("same");

  if (isSave) {
    DEtaME112 -> SaveAs("png/DEtaME112.png");
    DEtaME112 -> SaveAs("eps/DEtaME112.eps");
    DEtaME112 -> SaveAs("rootPlot/DEtaME112.root");
  }


  // --------------------------------------------------------------------------- 
  // hDEtaME113
  // ---------------------------------------------------------------------------  
  SetStyleh1(hDEtaME113, 629, 3001, 3, "#Delta #eta_{ME1/1-ME3}");
  SetStyleh1(hDEtaME113Bkg, 597, 3001, 3, "#Delta #eta_{ME1/1-ME3}");

  hDEtaME113Bkg -> SetLineStyle(2);

  TCanvas* DEtaME113 = new TCanvas("DEtaME113", "", 1260, 1000, 400, 400);

  Double_t normColl =  hDEtaME113->GetEntries();
  if (normColl) hDEtaME113->Scale(1/normColl);
  else cout << "hDEtaME113->GetEntries(): " << normColl << endl; 

  Double_t normCosm =  hDEtaME113Bkg->GetEntries();
  if (normCosm) hDEtaME113Bkg->Scale(1/normCosm);
  else cout << "hDEtaME113Bkg->GetEntries(): " << normCosm << endl; 

  if (hDEtaME113->GetMaximum() > hDEtaME113Bkg->GetMaximum()) {
    if (hDEtaME113->GetEntries())    hDEtaME113 -> Draw();
    if (hDEtaME113Bkg->GetEntries()) hDEtaME113Bkg -> Draw("same");
  }
  else {
    if (hDEtaME113Bkg->GetEntries()) hDEtaME113Bkg -> Draw();
    if (hDEtaME113->GetEntries())    hDEtaME113 -> Draw("same");
  }

  PrintIt(DEtaME113, "Halo #Delta #eta (ME1/1-ME3 extrapolation)");

  TString collisions = (" COLLISIONS:");
  collisions += hDEtaME113->GetEntries();
  collisions += " entries";

  TString cosmics = (" COSMICS:");
  cosmics += hDEtaME113Bkg->GetEntries();
  cosmics += " entries";

  tl = SetLegend(0.50, 0.73, 0.64, 0.87);
  tl->AddEntry(hDEtaME113, collisions, "l");
  tl->AddEntry(hDEtaME113Bkg, cosmics, "l");

  tl -> Draw("same");

  if (isSave) {
    DEtaME113 -> SaveAs("png/DEtaME113.png");
    DEtaME113 -> SaveAs("eps/DEtaME113.eps");
    DEtaME113 -> SaveAs("rootPlot/DEtaME113.root");
  }


       
} // end 
