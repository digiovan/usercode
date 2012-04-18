import FWCore.ParameterSet.Config as cms

process = cms.Process("UFDiMuonAnalyzer")

# no pat routine as it is already in the skims
thisIsData = True

if thisIsData:
    print 'Running over data sample'
else:
    print 'Running over MC sample'



process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.GeometryExtended_cff")
process.load('Configuration.EventContent.EventContent_cff')

# global tag
if thisIsData:
    print 'Loading Global Tag For Data'
    process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
    process.GlobalTag.globaltag = "GR_R_52_V7::All"
else:
    print 'Loading Global Tag For MC'
    process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
    process.GlobalTag.globaltag = "START52_V5::All"


# ------------ PoolSource -------------
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring('file:/data/uftrig01b/digiovan/root/CMSSW_5_2_3_patch1/SkimsSingleMuRun2012A-PromptReco-v1_LowPU/mergedFiles/Run2012A-PromptReco-v1_run190949_00001.root')
                            )
# -------- PoolSource END -------------


#===============================================================================
process.load("ZPtToMuMu.UFDiMuonsAnalyzer.UFDiMuonAnalyzer_nocuts_cff")
process.dimuons = process.DiMuons.clone()
process.dimuons.getFilename = cms.untracked.string("Run190949.root")
process.dimuons.isMonteCarlo = cms.bool(False) 

## process.dimuons.nMuons = cms.int32(0)
## process.dimuons.ptMin= cms.double(-999)
## process.dimuons.etaMaxLoose = cms.double(999)
## process.dimuons.etaMaxTight = cms.double(999)
## process.dimuons.normChiSquareMax= cms.double(999)
## # number of hits cuts
## process.dimuons.minMuonHits	  = cms.int32(-999)
## process.dimuons.minPixelHits	  = cms.int32(-999)
## process.dimuons.minStripHits	  = cms.int32(-999)
## process.dimuons.minTrackerHits	  = cms.int32(-999)
## process.dimuons.minSegmentMatches = cms.int32(-999)
## 
## process.dimuons.d0Max = cms.double(999) 
## # isolation
## process.dimuons.trackIsoMaxSumPt = cms.double(9999)
## process.dimuons.relCombIsoMax    = cms.double(9999)
## 
## process.dimuons.checkTrigger   = cms.bool(False)

#process.dimuons.selectLowestSingleMuTrigger = cms.untracked.bool(True)
process.dimuons.processName    = cms.string("HLT")
process.dimuons.triggerNames   = cms.vstring("HLT_Mu15_eta2p1")
process.dimuons.triggerResults = cms.InputTag("TriggerResults","","HLT")
process.dimuons.triggerEvent   = cms.InputTag("hltTriggerSummaryAOD","","HLT")
 
## process.dimuons.metTag         = cms.InputTag("pfMet")
## process.dimuons.pfJetsTag      = cms.InputTag("cleanPatJets")
## process.dimuons.genJetsTag     = cms.InputTag("null")

#===============================================================================

process.p = cms.Path( process.dimuons )

process.outpath = cms.EndPath()
#===============================================================================

