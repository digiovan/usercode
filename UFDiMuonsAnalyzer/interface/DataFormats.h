// event info
typedef struct {
  int run;
  int lumi;
  int event;
  int bx;
  int orbit;
} _EventInfo;


// vertex info
typedef struct{
  int nVertices;
  int isValid[20];
  float x[20];	
  float y[20];	
  float z[20];	
  float xErr[20];	
  float yErr[20];	
  float zErr[20];	
  float chi2[20];
  int ndf[20];
  float normChi2[20];
} _VertexInfo;


// basic track info
typedef struct {
  int charge;
  float pt;
  float ptErr;
  float eta;
  float phi;
} _TrackInfo;


// muon info
typedef struct {

  int isTracker;
  int isStandAlone;
  int isGlobal;

  int charge;
  float pt;
  float ptErr;
  float eta;
  float phi;

  float trkPt;
  float trkPtErr;
  float trketa;
  float trkPhi;
  
  float normChiSquare;
  float d0;
  float dz;
  
  int numPixelLayers;   //number of pixel layers with valid hits
  int numTrackerLayers; //number of tracker layers with valid hits 
  int numStripLayers;   //number of strip layers with valid hits

  float validFracTracker; //valid fraction of tracker hits

  int numValidMuonHits;
  int numValidPixelHits;
  int numValidTrackerHits;
  int numValidStripHits;
  int numSegmentMatches;
  int numOfMatchedStations;

  float trackIsoSumPt;
  float hcalIso;
  float ecalIso;
  float relCombIso;

  // PF information
  int isPFMuon;

  float pfPt;
  float pfEta;
  float pfPhi;
  
  float sumChargedHadronPtR03; // sum-pt of charged Hadron 
  float sumChargedParticlePtR03; // sum-pt of charged Particles(inludes e/mu) 
  float sumNeutralHadronEtR03;  // sum pt of neutral hadrons
  float sumPhotonEtR03;  // sum pt of PF photons
  float sumPUPtR03;  // sum pt of charged Particles not from PV  (for Pu corrections)

  float sumChargedHadronPtR04; 
  float sumChargedParticlePtR04;
  float sumNeutralHadronEtR04;  
  float sumPhotonEtR04;
  float sumPUPtR04;

  int isHltMatched[2];
  
} _MuonInfo;


//MET
typedef struct {
  float metPx;
  float metPy;
  float metPt;
  float metPhi;
  float sumEt;
} _MetInfo;


// pf Jets
typedef struct {
  int   nPFjets;
  float pfJetPx[10];
  float pfJetPy[10];
  float pfJetPz[10];
  float pfJetPt[10];
  float pfJetEta[10];
  float pfJetPhi[10];
  float pfJetM[10];
  int   pfJetCharge[10];
  /////// Energy Fractions //////
  //Charged Hadron
  float pfJetCHF[10];
  //NeutralHadron
  float pfJetNHF[10];
  //Charged EM
  float pfJetCEF[10];
  //Neutral EM
  float pfJetNEF[10];
  //Mu
  float pfJetMuF[10];
  // HF Hadron Fraction
  float pfJetHFHF[10];
  // HF EM Fraction
  float pfJetHFEF[10];
  /////// Multiplicities //////
  // Total Charged
  int pfJetCM[10];
  //Charged Hadron
  int pfJetCHM[10];
  //NeutralHadron
  int pfJetNHM[10];
  //Charged EM
  int pfJetCEM[10];
  //Neutral EM
  int pfJetNEM[10];
  //Mu
  int pfJetMuM[10];
  // HF Hadron Fraction
  int pfJetHFHM[10];
  // HF EM Fraction
  int pfJetHFEM[10];
  // Jet Correction Factor--Above momentum is already corrected!!
  // This factor will return momentum to uncorrected value!!
  float pfJetJECFactor[10];
  // Gen Jet Values
  bool pfJetGenMatched[10];
  float pfJetGenPx[10];
  float pfJetGenPy[10];
  float pfJetGenPz[10];
  float pfJetGenPt[10];
  float pfJetGenEta[10];
  float pfJetGenPhi[10];
  float pfJetGenM[10];
  ///// Gen Jet Energy Fractions ///////
  // EM Fraction
  float pfJetGenEMF[10];
  // Had Fraction
  float pfJetGenHadF[10];
  // Invisible Fraction
  float pfJetGenInvF[10];
  // Auxiliary Fraction (Undecayed Sigmas, etc.)
  float pfJetGenAuxF[10];
} _PFJetInfo;


// generator level jets
typedef struct {
  int nGenJets;
  float genJetPx[10];
  float genJetPy[10];
  float genJetPz[10];
  float genJetPt[10];
  float genJetEta[10];
  float genJetPhi[10];
  float genJetM[10];
  int   genJetCharge[10];
} _GenJetInfo;
