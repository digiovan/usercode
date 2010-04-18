#ifndef TrigEff_h
#define TrigEff_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

#include <string>
#include <TFile.h>
#include <TTree.h>

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
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
//
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
// L1GTUtils
#include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtUtils.h"

// Geometry
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/CSCRecHit/interface/CSCRecHit2D.h"
#include "FWCore/Framework/interface/ESHandle.h"

// for track propagation
#include "DataFormats/GeometrySurface/interface/Plane.h"

// Transient tracks (for extrapolations)
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TrackTransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"

//CSCTF
#include "DataFormats/L1CSCTrackFinder/interface/L1CSCTrackCollection.h"
#include "L1Trigger/CSCTrackFinder/interface/CSCTFPtLUT.h"
#include "L1Trigger/CSCTrackFinder/interface/CSCSectorReceiverLUT.h"
#include "CondFormats/L1TObjects/interface/L1MuTriggerScales.h"
#include "CondFormats/DataRecord/interface/L1MuTriggerScalesRcd.h"
#include "CondFormats/L1TObjects/interface/L1MuTriggerPtScale.h"
#include "CondFormats/DataRecord/interface/L1MuTriggerPtScaleRcd.h"
#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigiCollection.h"

#include <vector>
#include "TMatrixF.h"
#include "TMatrixD.h"

//---------------------------------------------------------------------------
#include "TMath.h"
#define PI TMath::Pi()
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
#define MAX_MUONS 100 
#define MAX_CSC_RECHIT 35 // 4 stations x 6 layers + some overlap between chambers
#define MAX_TRK_SEGS 10   // max # segments to a given tracker muon
#define MAX_CSCTF_TRK 36  // max # of CSCTF tracks per BX
#define MAX_LCTS_PER_TRK 4  // max # of LCTS which form a CSCTF track
//---------------------------------------------------------------------------

using namespace trigger;
using namespace reco;
using namespace std;
using namespace edm;

class TrigEff : public edm::EDAnalyzer {
private:
    edm::InputTag genTag;
    edm::InputTag L1extraTag;
    edm::InputTag Level1Tag;
    edm::InputTag Level2Tag;
    edm::InputTag Level3Tag;
    edm::InputTag tracksTag;
    edm::InputTag muonsTag;
    edm::InputTag saMuonsTag;
    edm::InputTag TriggerEventTag;
    edm::InputTag HLTriggerTag;
    double eta_cut, pt_cut, matching_cone;
    std::string outputFile;
    std::string level2module;
    std::string level3module;

    edm::InputTag csctfTag;

    // to Remove
    edm::InputTag gpTag;

    //---------------------------------------------------------------------
    // Default global muon collection
    //---------------------------------------------------------------------
    TTree *recoMuons;
    void fillMuons(const edm::Handle<reco::MuonCollection> muons);
 
    //---------------------------------------------------------------------
    // general information 
    //--------------------------------------------------------------------- 
    int Run, Event, Lumi, Bx, Orbit;

    int muonSize; 

    std::vector<int>*   isGlobalMuon;
    std::vector<int>*   isTrackerMuon;
    std::vector<int>*   isStandAloneMuon;
    std::vector<int>*   isMuonAllArbitrated;

    std::vector<int>*   isEnergyValid;
    std::vector<float>* caloCompatibility; 
    std::vector<float>* em; 
    std::vector<float>* emS9;
    std::vector<float>* emS25;
    std::vector<float>* emMax;
    std::vector<float>* had;
    std::vector<float>* hadS9;
    std::vector<float>* hadMax;
    std::vector<float>* ho;
    std::vector<float>* hoS9;
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    // Global Muon Block
    //---------------------------------------------------------------------
    std::vector<float>* gmrPt;
    std::vector<float>* gmrEta;
    std::vector<float>* gmrPhi;
    std::vector<float>* gmrP;
    std::vector<float>* gmrPx;
    std::vector<float>* gmrPy;
    std::vector<float>* gmrPz;
    std::vector<float>* gmrTheta;
    std::vector<float>* gmrVx;
    std::vector<float>* gmrVy;
    std::vector<float>* gmrVz;
    std::vector<float>* gmrCharge;
    std::vector<float>* gmrNDoF;
    std::vector<float>* gmrChi2;
    std::vector<float>* gmrChi2Norm;
    std::vector<float>* gmrDXY;
    std::vector<float>* gmrDTheta;
    std::vector<float>* gmrDPt;
    std::vector<float>* gmrDEta;
    std::vector<float>* gmrDPhi;
    std::vector<float>* gmrDDXY;
    std::vector<float>* gmrIso03nTracks;
    std::vector<float>* gmrIso03sumPt;
    std::vector<float>* gmrDz;
    std::vector<float>* gmrD0;
    std::vector<float>* gmrDsz;
    std::vector<float>* gmrDDz;
    std::vector<float>* gmrDD0;
    std::vector<float>* gmrDDsz;
    std::vector<float>* gmrInnerX;
    std::vector<float>* gmrInnerY;
    std::vector<float>* gmrInnerZ;
    std::vector<float>* gmrOuterX;
    std::vector<float>* gmrOuterY;
    std::vector<float>* gmrOuterZ;
    std::vector<int>*   gmrValHits;   
 
    //---------------------------------------------------------------------
    // Standalone Muon Block
    //---------------------------------------------------------------------
    std::vector<float>* stdEnergy;
    std::vector<float>* stdPt;
    std::vector<float>* stdEta;
    std::vector<float>* stdPhi;
    std::vector<float>* stdPx;
    std::vector<float>* stdPy;
    std::vector<float>* stdPz;
    std::vector<float>* stdVx;
    std::vector<float>* stdVy;
    std::vector<float>* stdVz;
    std::vector<float>* stdCharge;
    std::vector<float>* stdDPt;
    std::vector<float>* stdDEta;
    std::vector<float>* stdDPhi;
    std::vector<float>* stdDz;
    std::vector<float>* stdD0;    
    std::vector<float>* stdNDoF;
    std::vector<float>* stdChi2;
    std::vector<float>* stdChi2Norm;
    std::vector<float>* stdDXY;
    std::vector<float>* stdTheta;
    std::vector<float>* stdDTheta;
    std::vector<float>* stdDDz;
    std::vector<float>* stdDD0;
    std::vector<int>*   stdValHits;
    std::vector<float>* stdInnerX;
    std::vector<float>* stdInnerY;
    std::vector<float>* stdInnerZ;
    std::vector<float>* stdOuterX;
    std::vector<float>* stdOuterY;
    std::vector<float>* stdOuterZ;
   
    //---------------------------------------------------------------------
    // Tracker Muon Block
    //---------------------------------------------------------------------
    std::vector<float>* trkEnergy;
    std::vector<float>* trkPt;
    std::vector<float>* trkEta;
    std::vector<float>* trkPhi;
    std::vector<float>* trkPx;
    std::vector<float>* trkPy;
    std::vector<float>* trkPz;
    std::vector<float>* trkVx;
    std::vector<float>* trkVy;
    std::vector<float>* trkVz;
    std::vector<float>* trkCharge;
    std::vector<float>* trkDPt;
    std::vector<float>* trkDEta;
    std::vector<float>* trkDPhi;
    std::vector<float>* trkDz;
    std::vector<float>* trkD0;    
    std::vector<float>* trkNDoF;
    std::vector<float>* trkChi2;
    std::vector<float>* trkChi2Norm;
    std::vector<float>* trkDXY;
    std::vector<float>* trkTheta;
    std::vector<float>* trkDTheta;
    std::vector<float>* trkDDz;
    std::vector<float>* trkDD0;
    std::vector<int>*   trkValHits;    

    // CSC segment for the tracker muon
    //chamber
    std::vector<int>* trkNchambers; 
    std::vector<int>* trkNofMatches;
    std::vector<int>* trkIsMatchValid;

    //segment
    std::vector<int>* trkNSegs;
    TMatrixF trkSegChamberId; 
    TMatrixF trkSegRing;    
    TMatrixF trkSegStation; 
    TMatrixF trkSegEndcap;  
    TMatrixF trkSegTriggerSector;
    TMatrixF trkSegTriggerCscId; 
    TMatrixF trkSegXfromMatch;
    TMatrixF trkSegYfromMatch;
    TMatrixF trkSegPhifromMatch;
 
    TMatrixF trkSegIsArb;
    TMatrixF trkSegX;
    TMatrixF trkSegY;
    TMatrixF trkSegR;
    TMatrixF trkSegPhi;

    //---------------------------------------------------------------------
    // RECHIT information: only for standalone/global muons!
    //---------------------------------------------------------------------
    std::vector<int>*   rchCSCtype;
    std::vector<float>* rchEta;
    std::vector<float>* rchPhi;
    std::vector<float>* rchPhi_02PI;
    std::vector<int>*   rchStation;
    std::vector<int>*   rchChamber;
    std::vector<int>*   rchRing;
    std::vector<int>*   rchLayer;
    
    std::vector<int>* rchMuonSize;//# rechits per muon
    TMatrixF     rchEtaMatrix;
    TMatrixF     rchPhiMatrix;
    TMatrixF     rchPhi02PIMatrix;
    TMatrixF     rchStationMatrix;
    TMatrixF     rchChamberMatrix;
    TMatrixF     rchRingMatrix;
    TMatrixF     rchLayerMatrix;
    TMatrixF     rchTypeMatrix;

    //--------------------------------------------------------------------- 
    // old format: keep it until you are sure the TMatrixD works
    Int_t     nMu_nCscRchHits; // nMu x MAX_CSC_RECHIT
    Double_t *rchEtaList;
    Double_t *rchPhiList;
    Double_t *rchPhiList_02PI;

    void resizeRchHits(int nMu_nMaxCscRchHits);
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    // Propagation block
    //---------------------------------------------------------------------
    // propagation to ME1/1
    std::vector<float>*    muons_x_mep11;
    std::vector<float>*    muons_y_mep11;
    std::vector<float>*    muons_z_mep11;
    std::vector<float>*  muons_phi_mep11;
    std::vector<float>*  muons_eta_mep11;
                                 
    std::vector<float>*    muons_x_mem11;
    std::vector<float>*    muons_y_mem11;
    std::vector<float>*    muons_z_mem11;
    std::vector<float>*  muons_phi_mem11;
    std::vector<float>*  muons_eta_mem11;

    // propagation to ME1
    std::vector<float>*    muons_x_mep1;
    std::vector<float>*    muons_y_mep1;
    std::vector<float>*    muons_z_mep1;
    std::vector<float>*  muons_phi_mep1;
    std::vector<float>*  muons_eta_mep1;
                                 
    std::vector<float>*    muons_x_mem1;
    std::vector<float>*    muons_y_mem1;
    std::vector<float>*    muons_z_mem1;
    std::vector<float>*  muons_phi_mem1;
    std::vector<float>*  muons_eta_mem1;

    // propagation to ME2
    std::vector<float>*    muons_x_mep2;
    std::vector<float>*    muons_y_mep2;
    std::vector<float>*    muons_z_mep2;
    std::vector<float>*  muons_phi_mep2;
    std::vector<float>*  muons_eta_mep2;
      
    std::vector<float>*    muons_x_mem2;
    std::vector<float>*    muons_y_mem2;
    std::vector<float>*    muons_z_mem2;
    std::vector<float>*  muons_phi_mem2;
    std::vector<float>*  muons_eta_mem2;

    // propagation to ME3
    std::vector<float>*    muons_x_mep3;
    std::vector<float>*    muons_y_mep3;
    std::vector<float>*    muons_z_mep3;
    std::vector<float>*  muons_phi_mep3;
    std::vector<float>*  muons_eta_mep3;
                                 
    std::vector<float>*    muons_x_mem3;
    std::vector<float>*    muons_y_mem3;
    std::vector<float>*    muons_z_mem3;
    std::vector<float>*  muons_phi_mem3;
    std::vector<float>*  muons_eta_mem3;
    //---------------------------------------------------------------------


    //---------------------------------------------------------------------
    // l1 extra muon collection
    //---------------------------------------------------------------------
    TTree *l1extraMuons;
    void fillExtra(l1extra::L1MuonParticleCollection::const_iterator l1muon);

    int l1Size;
    
    std::vector<float>* l1Eta;
    std::vector<float>* l1Pt;
    std::vector<float>* l1Phi;

    std::vector<int>* isIsolated;
    std::vector<int>* isMip;
    std::vector<int>* isForward;
    std::vector<int>* isRPC;     
    std::vector<int>* detectorType;     
    std::vector<int>* rank;     
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    // Service
    //---------------------------------------------------------------------
    edm::ESHandle<CSCGeometry> cscGeom;
    //
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    // TRACK information 
    //---------------------------------------------------------------------
    TTree *csctfTTree; 
    void fillCSCTF(const edm::Handle<L1CSCTrackCollection> tracks,
                   const L1MuTriggerScales  *ts, 
                   const L1MuTriggerPtScale *tpts, 
                   CSCSectorReceiverLUT* srLUTs_[5][2]); 
    
    int SizeTrk;
    std::vector<int>*    EndcapTrk; 
    std::vector<int>*    SectorTrk; 
    std::vector<int>*    BxTrk; 
    
    std::vector<int>*    me1ID; 
    std::vector<int>*    me2ID; 
    std::vector<int>*    me3ID; 
    std::vector<int>*    me4ID; 
    std::vector<int>*    mb1ID;     
    
    std::vector<int>*    OutputLinkTrk; 
    
    std::vector<int>*    ModeTrk; 
    std::vector<float>*  EtaTrk;   
    std::vector<float>*  PhiTrk;   
    std::vector<float>*  PhiTrk_02PI; 
    std::vector<float>*  PtTrk;   
    
    std::vector<int>*    ChargeTrk; 
    std::vector<int>*    ChargeValidTrk; 
    std::vector<int>*    QualityTrk; 
    std::vector<int>*    ForRTrk; 
    std::vector<int>*    Phi23Trk; 
    std::vector<int>*    Phi12Trk;   
    std::vector<int>*    PhiSignTrk;   
    
    std::vector<int>*    EtaBitTrk;   
    std::vector<int>*    PhiBitTrk;   
    std::vector<int>*    PtBitTrk;   

    CSCSectorReceiverLUT* srLUTs_[5][2];
    const L1MuTriggerScales  *ts;
    const L1MuTriggerPtScale *tpts;
    unsigned long long m_scalesCacheID ;
    unsigned long long m_ptScaleCacheID ;

    // LCT (STUBS FORMING THE TRACK)  
    // it contains the number of LCT forming a track 
    std::vector<int>* NumLCTsTrk; 
    
    TMatrixD trLctEndcap; 
    TMatrixD trLctSector; 
    TMatrixD trLctSubSector; 
    TMatrixD trLctBx; 
    TMatrixD trLctBx0; 
       
    TMatrixD trLctStation; 
    TMatrixD trLctRing; 
    TMatrixD trLctChamber; 
    TMatrixD trLctTriggerCSCID; 
    TMatrixD trLctFpga;	  

     // note: the SPs return them in bits 
    TMatrixD trLctlocalPhi; 
    TMatrixD trLctglobalPhi;   
    TMatrixD trLctglobalEta; 

    TMatrixD trLctstripNum;   
    TMatrixD trLctwireGroup;     

    //---------------------------------------------------------------------

    //snippet from HLTEventAnalyzerRAW
    /// module config parameters
    std::string   processName_;
    std::string   triggerName_;
    edm::InputTag triggerResultsTag_;
    edm::InputTag triggerEventWithRefsTag_;
    
    /// additional class data memebers
    edm::Handle<edm::TriggerResults>           triggerResultsHandle_;
    edm::Handle<trigger::TriggerEventWithRefs> triggerEventWithRefsHandle_;
    HLTConfigProvider hltConfig_;
    
    /// payload extracted from TriggerEventWithRefs
    trigger::Vids        muonIds_;
    trigger::VRmuon      muonRefs_;
    trigger::VRmuon      muonRefsl2_;
    trigger::VRmuon      muonRefsl3_;
    trigger::Vids        l1muonIds_;
    trigger::VRl1muon    l1muonRefs_;
    ////////////////////////////

    // L1GTUtils
    L1GtUtils m_l1GtUtils;
    std::string m_nameAlgTechTrig;
    //

    // ttree
    TFile *file;
    TTree *trig1, *trig2, *trig3, *trig123, *trig321, *mcTruth;

    int l1trigger, l2trigger, l3trigger, trueMuonPair;
    double etaTrg1, phiTrg1, ptTrg1;
    double etaTrg2, phiTrg2, ptTrg2;
    double etaTrg3, phiTrg3, ptTrg3;
    double etaRec,  phiRec,  ptRec;
    double dR1, dR2, dR3;
    double minv, minv1, minv2, minv3, pt3, eta3;

    double minvLow, minvHigh;

    int printLevel;
    char processName[200], triggerName[200];
    int wasRun;
    int wasAccepted;
    int whichError;

    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob(void){}

    //virtual void analyzeTrigger(const edm::Event&, const std::string&);
    //virtual void beginRun(edm::Run const&, edm::EventSetup const&);
    //virtual void analyzeL1GtUtils(const edm::Event&, const edm::EventSetup&);

    class TrigInfo {
      public:
        l1extra::L1MuonParticleCollection::const_iterator l1cand;
        VRmuon::const_iterator l2cand;
        VRmuon::const_iterator l3cand;
        MuonRef tag;
        double minv12;
        double minv3;
        double pt3;
        double eta3;
    };

    // The Magnetic field
    edm::ESHandle<MagneticField> theBField;
    // The GlobalTrackingGeometry
    edm::ESHandle<GlobalTrackingGeometry> theTrackingGeometry;
    // Extrapolator to cylinder
    edm::ESHandle<Propagator> propagatorAlong;
    edm::ESHandle<Propagator> propagatorOpposite;
    
    FreeTrajectoryState freeTrajStateMuon(reco::TrackRef track);

public:
    explicit TrigEff(const edm::ParameterSet&);
    ~TrigEff(void);

    void muonsInit();
    void muonsDel (); 

    void l1extraInit();
    void l1extraDel (); 

    TrajectoryStateOnSurface  surfExtrapTrkSam (reco::TrackRef track, double z);

    void csctfInit(); 
    void csctfDel ();  

};

#endif
