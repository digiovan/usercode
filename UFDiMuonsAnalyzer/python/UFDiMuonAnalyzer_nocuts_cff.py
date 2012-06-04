import FWCore.ParameterSet.Config as cms

DiMuons = cms.EDAnalyzer('UFDiMuonsAnalyzer',
                         getFilename = cms.untracked.string("tmpName.root"),
                         muonColl = cms.InputTag("muons"),
                         
                         isVerbose = cms.untracked.bool(False),
                         isMonteCarlo = cms.bool(False),
                         nMuons = cms.int32(0),
                         
                         # muon kinematic cuts
                         isGlobal    =cms.int32(0),
                         isStandAlone=cms.int32(0),
                         isTracker   =cms.int32(0),
                         ptMin       = cms.double(-999),
                         etaMax      = cms.double(999),
                         normChiSquareMax= cms.double(999),

                         # number of hits cuts
                         minMuonHits	  = cms.int32(-999),
                         minPixelHits	  = cms.int32(-999),
                         minStripHits	  = cms.int32(-999),
                         minTrackerHits	  = cms.int32(-999),
                         minSegmentMatches= cms.int32(-999),
                         minNumOfMatchedStations = cms.int32(-999),

                         minPixelLayers  = cms.int32(-999),
                         minTrackerLayers= cms.int32(-999),
                         minStripLayers  = cms.int32(-999),
                         validFracTrackerMin= cms.int32(-999),   
    
                         # if beamSpot is not specified
                         # assumes you are using "offlineBeamSpot"
                         beamSpot = cms.untracked.InputTag("offlineBeamSpot"),
                         d0Max = cms.double(999), 
                         
                         #track isolation
                         trackIsoMaxSumPt = cms.double(9999),
                         relCombIsoMax    = cms.double(9999),
                         
                         #triggerName = cms.string("@"),#wild card: all triggers
                         # HLT trigger info
                         checkTrigger = cms.bool(False),

                         selectLowestSingleMuTrigger = cms.untracked.bool(False),
                         processName = cms.string("HLT"),
                         triggerNames = cms.vstring("HLT_Mu30"),
                         #triggerName = cms.string("I want to crash"),
                         triggerResults = cms.InputTag("TriggerResults","","HLT"),
                         triggerEvent   = cms.InputTag("hltTriggerSummaryAOD","","HLT"),
                         
                         metTag         = cms.InputTag("null"),
                         pfJetsTag      = cms.InputTag("null"),
                         genJetsTag     = cms.InputTag("null")
                         
                         )
