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

#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"

#include <vector>
#include "TMatrixF.h"
#include "TMatrixD.h"

#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"

#include "SegmentLCTMatchBox.h"

//---------------------------------------------------------------------------
#include "TMath.h"
#define PI TMath::Pi()
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
#define MAX_MUONS 100 
#define MAX_CSC_RECHIT 48 // 4 stations x 6 layers + (some overlap between chambers: added other 24 hits to be safe)
#define MAX_TRK_SEGS 100   // max # segments to a given tracker muon
#define MAX_CSCTF_TRK 36  // max # of CSCTF tracks per BX
#define MAX_LCTS_PER_TRK 4  // max # of LCTS which form a CSCTF track
#define MAX_SEGS_STD 16 // MAX number of segments which can be associated to StandAlone component of the GBL muon
//---------------------------------------------------------------------------


using namespace trigger;
using namespace reco;
using namespace std;
using namespace edm;

class TrigEff : public edm::EDAnalyzer {
private:
    edm::InputTag L1extraTag;
    edm::InputTag muonsTag;
    std::string outputFile;

    edm::InputTag csctfTag;
    edm::InputTag csctfLctsTag;

    edm::InputTag cscSegTag;

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
    std::vector<int>*   isTMLastStationAngTight;//for trk muons only
    std::vector<int>*   isGlobalMuonPromptTight;//for gbl muons only

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
    //std::vector<float>* gmrEnergy;
    //std::vector<float>* gmrDEnergy;
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
    //std::vector<float>* stdEnergy;
    //std::vector<float>* stdDEnergy;
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
    //std::vector<float>* trkEnergy;
    //std::vector<float>* trkDEnergy;
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
    int trkSegChamberId[MAX_MUONS][MAX_TRK_SEGS]; 
    int trkSegRing[MAX_MUONS][MAX_TRK_SEGS];    
    int trkSegStation[MAX_MUONS][MAX_TRK_SEGS]; 
    int trkSegEndcap[MAX_MUONS][MAX_TRK_SEGS];  
    int trkSegTriggerSector[MAX_MUONS][MAX_TRK_SEGS];
    int trkSegTriggerCscId[MAX_MUONS][MAX_TRK_SEGS]; 
    float trkSegXfromMatch[MAX_MUONS][MAX_TRK_SEGS];
    float trkSegYfromMatch[MAX_MUONS][MAX_TRK_SEGS];
    float trkSegZfromMatch[MAX_MUONS][MAX_TRK_SEGS];
    float trkSegRfromMatch[MAX_MUONS][MAX_TRK_SEGS];
    float trkSegPhifromMatch[MAX_MUONS][MAX_TRK_SEGS];
    float trkSegEtafromMatch[MAX_MUONS][MAX_TRK_SEGS];
      
    int trkSegIsArb[MAX_MUONS][MAX_TRK_SEGS];
    float trkSegX[MAX_MUONS][MAX_TRK_SEGS];
    float trkSegY[MAX_MUONS][MAX_TRK_SEGS];
    float trkSegZ[MAX_MUONS][MAX_TRK_SEGS];
    float trkSegR[MAX_MUONS][MAX_TRK_SEGS];
    float trkSegPhi[MAX_MUONS][MAX_TRK_SEGS];
    float trkSegEta[MAX_MUONS][MAX_TRK_SEGS];

    float trkSegDxDz[MAX_MUONS][MAX_TRK_SEGS];
    float trkSegDyDz[MAX_MUONS][MAX_TRK_SEGS];
    float trkSegDxDzErr[MAX_MUONS][MAX_TRK_SEGS];
    float trkSegDyDzErr[MAX_MUONS][MAX_TRK_SEGS];

    //---------------------------------------------------------------------
    // RECHIT information: only for standalone/global muons!
    //---------------------------------------------------------------------
    std::vector<int>*   rchCSCtype;
    std::vector<float>* rchEtaLocal;
    std::vector<float>* rchPhiLocal;
    std::vector<float>* rchEta;
    std::vector<float>* rchPhi;
    std::vector<float>* rchPhi_02PI;
    std::vector<int>*   rchStation;
    std::vector<int>*   rchChamber;
    std::vector<int>*   rchRing;
    std::vector<int>*   rchLayer;
    
    std::vector<int>* rchMuonSize;//# rechits per muon
    float rchEtaMatrixLocal[MAX_MUONS][MAX_CSC_RECHIT];
    float rchPhiMatrixLocal[MAX_MUONS][MAX_CSC_RECHIT];
    float rchEtaMatrix[MAX_MUONS][MAX_CSC_RECHIT];
    float rchPhiMatrix[MAX_MUONS][MAX_CSC_RECHIT];
    float rchPhi02PIMatrix[MAX_MUONS][MAX_CSC_RECHIT];
    int   rchStationMatrix[MAX_MUONS][MAX_CSC_RECHIT];
    int   rchChamberMatrix[MAX_MUONS][MAX_CSC_RECHIT];
    int   rchRingMatrix[MAX_MUONS][MAX_CSC_RECHIT];
    int   rchLayerMatrix[MAX_MUONS][MAX_CSC_RECHIT];
    int   rchTypeMatrix[MAX_MUONS][MAX_CSC_RECHIT];

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
    std::vector<float>*    muons_x_me11;
    std::vector<float>*    muons_y_me11;
    std::vector<float>*    muons_z_me11;
    std::vector<float>*  muons_phi_me11;
    std::vector<float>*  muons_eta_me11;
                                 
    // propagation to ME1
    std::vector<float>*    muons_x_me1;
    std::vector<float>*    muons_y_me1;
    std::vector<float>*    muons_z_me1;
    std::vector<float>*  muons_phi_me1;
    std::vector<float>*  muons_eta_me1;
                                 
    // propagation to ME2
    std::vector<float>*    muons_x_me2;
    std::vector<float>*    muons_y_me2;
    std::vector<float>*    muons_z_me2;
    std::vector<float>*  muons_phi_me2;
    std::vector<float>*  muons_eta_me2;
      
    // propagation to ME3
    std::vector<float>*    muons_x_me3;
    std::vector<float>*    muons_y_me3;
    std::vector<float>*    muons_z_me3;
    std::vector<float>*  muons_phi_me3;
    std::vector<float>*  muons_eta_me3;
                                 
    //---------------------------------------------------------------------
    // segment information
    //---------------------------------------------------------------------

    int segsSize; 
    std::vector<float>* cscsegs_loc_x;
    std::vector<float>* cscsegs_loc_y;
    std::vector<float>* cscsegs_loc_z;

    std::vector<float>* cscsegs_loc_theta;
    std::vector<float>* cscsegs_loc_eta;
    std::vector<float>* cscsegs_loc_phi;

    std::vector<float>* cscsegs_loc_dir_theta;
    std::vector<float>* cscsegs_loc_dir_eta;
    std::vector<float>* cscsegs_loc_dir_phi;

    std::vector<float>* cscsegs_gbl_x;
    std::vector<float>* cscsegs_gbl_y;
    std::vector<float>* cscsegs_gbl_z;

    std::vector<float>* cscsegs_gbl_theta;
    std::vector<float>* cscsegs_gbl_eta;
    std::vector<float>* cscsegs_gbl_phi;

    std::vector<float>* cscsegs_gbl_dir_theta;
    std::vector<float>* cscsegs_gbl_dir_eta;
    std::vector<float>* cscsegs_gbl_dir_phi;

    std::vector<int>* cscsegs_endcap;
    std::vector<int>* cscsegs_station;
    std::vector<int>* cscsegs_ring;
    std::vector<int>* cscsegs_chamber;

    void  fillSegments(edm::Handle<CSCSegmentCollection> cscSegments, 
                       edm::ESHandle<CSCGeometry> cscGeom);


    //--------------------------------------------------------------------------
    // Record information about segments belonging to the STD muon component
    //--------------------------------------------------------------------------
    // how many segments are associated to the muon candidate
    std::vector<int>* muonNsegs; 

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


    // fill the variables of the segments belonging to the muon
    void  fillSegmentsMuons ( const edm::Handle<reco::MuonCollection> muons,
                              edm::Handle<CSCSegmentCollection> cscSegments, 
                              edm::ESHandle<CSCGeometry> cscGeom,
                              const edm::Handle<CSCCorrelatedLCTDigiCollection> CSCTFlcts);

    // useful methods...
    // 1) new matching code from Ivan to see if the segment is LCTAble
    bool isLCTAble ( const CSCSegment &segment);
    // 2) given a segment, is there any LCT in the collection which matches it?
    bool isMatched ( const CSCSegment &segment, 
                     edm::Handle<CSCCorrelatedLCTDigiCollection> );
  
    // 3) inline function for the half strip and wiregroup determination
    static int wireGroup ( const CSCRecHit2D& hit ); 
    static int halfStrip ( const CSCRecHit2D& hit );
    /* int halfStrip ( int channel, float position ){  */
/*       int retVal = 2*channel; */
/*       if (position < 0) retVal -=1; */
/*       return retVal; */
/*     } */

/*     int wireGroup ( const CSCRecHit2D &hit ){ */
/*       return ( hit.wgroups()[0] -1 ); */
/*     } */


    // class container for the segment and number of hits matching the
    // global muon
    class Segment {
    public:      
      const CSCSegment cscsegcand;
      int nMatchedHits;
      Segment(const CSCSegment &seg,int nMHits):cscsegcand(seg),nMatchedHits(nMHits){};
    };

    typedef std::vector< Segment * > SegmentVector;
    
    // return a vector containing the segments associated to the muon
    TrigEff::SegmentVector* SegmentsInMuon(const reco::Muon* muon, const CSCSegmentCollection* segments );

    // importing directly Ivan's object...
    SegmentLCTMatchBox _matchBox;  
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
    
    // ALL LCTS
    void fillAllLCTs(const edm::Handle<CSCCorrelatedLCTDigiCollection> corrlcts,
                     CSCSectorReceiverLUT* srLUTs_[5][2]); 
   
    int SizeLCTs;
    std::vector<int>* lctEndcap; 
    std::vector<int>* lctSector; 
    std::vector<int>* lctSubSector; 
    std::vector<int>* lctBx; 
    std::vector<int>* lctBx0; 
    std::vector<int>* lctStation; 
    std::vector<int>* lctRing; 
    std::vector<int>* lctChamber; 
    std::vector<int>* lctTriggerCSCID; 
    std::vector<int>* lctFpga;     

    // note: the SPs return them in bits 
    std::vector<int>* lctlocalPhi; 
    std::vector<int>* lctglobalPhi;   
    std::vector<int>* lctglobalEta; 
    std::vector<int>* lctstripNum;   
    std::vector<int>* lctwireGroup;   


    //---------------------------------------------------------------------
    // L1GTUtils
    L1GtUtils m_l1GtUtils;
    std::string m_nameAlgTechTrig;
    
    // ttree
    TFile *file;
    
    int printLevel;
    
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob(void);//{}

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

    MuonServiceProxy* theService;


};

#endif
