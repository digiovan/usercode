import FWCore.ParameterSet.Config as cms

##from RecoMuon.TrackingTools.MuonServiceProxy_cff import *
#
# yay! the analyzer!
#
DiMuons = cms.EDAnalyzer('UFDiMuonsAnalyzer',
                         #ServiceParameters = cms.PSet(MuonServiceProxy),
                         getFilename = cms.untracked.string("tmpName.root"),
                         muonColl = cms.InputTag("muons"),
                         
                         isVerbose = cms.untracked.bool(False),
                         isMonteCarlo = cms.bool(False),
                         nMuons = cms.int32(2),
                         
                         # muon kinematic cuts
                         isGlobal    =cms.int32(0),
                         isStandAlone=cms.int32(0),
                         isTracker   =cms.int32(0),
                         ptMin       = cms.double(20),
                         etaMax      = cms.double(2.4),
                         normChiSquareMax= cms.double(999),

                         # number of hits cuts
                         minMuonHits	  = cms.int32(1),
                         minPixelHits	  = cms.int32(1),
                         minStripHits	  = cms.int32(-999),
                         minTrackerHits	  = cms.int32(-999),
                         minSegmentMatches= cms.int32(2),
                         minNumOfMatchedStations = cms.int32(-999),

                         minPixelLayers  = cms.int32(-999),
                         minTrackerLayers= cms.int32(8),
                         minStripLayers  = cms.int32(-999),
                         validFracTrackerMin= cms.int32(-999),   
    
                         # if beamSpot is not specified
                         # assumes you are using "offlineBeamSpot"
                         beamSpot = cms.untracked.InputTag("offlineBeamSpot"),
                         d0Max = cms.double(0.2), 
                         
                         #track isolation
                         trackIsoMaxSumPt = cms.double(9999),
                         relCombIsoMax    = cms.double(0.1),
                         
                         #triggerName = cms.string("@"),#wild card: all triggers
                         # HLT trigger info
                         checkTrigger = cms.bool(True),

                         selectLowestSingleMuTrigger = cms.untracked.bool(False),
                         processName = cms.string("HLT"),
                         triggerNames = cms.vstring("HLT_Mu15"),
                         #triggerName = cms.string("I want to crash"),
                         triggerResults = cms.InputTag("TriggerResults","","HLT"),
                         triggerEvent   = cms.InputTag("hltTriggerSummaryAOD","","HLT"),
                         
                         metTag         = cms.InputTag("null"),
                         pfJetsTag      = cms.InputTag("null"),
                         genJetsTag     = cms.InputTag("null")
                         
                         )
