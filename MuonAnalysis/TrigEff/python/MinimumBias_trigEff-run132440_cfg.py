import FWCore.ParameterSet.Config as cms

process = cms.Process("TrigEff")

# import of standard configurations
process.load('Configuration/StandardSequences/Services_cff')
process.load('Configuration/StandardSequences/GeometryIdeal_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/SimL1Emulator_cff')
process.load("Configuration.StandardSequences.RawToDigi_Data_cff")
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/EventContent/EventContent_cff')
process.load('Configuration.StandardSequences.ReconstructionCosmics_cff')

process.load("RecoMuon.TrackingTools.MuonServiceProxy_cff")
process.load("RecoMuon.TrackingTools.MuonTrackLoader_cff")

##process.load('Configuration/StandardSequences/Services_cff')
##process.load('Configuration/StandardSequences/GeometryIdeal_cff')
##process.load('Configuration/StandardSequences/MagneticField_38T_cff')
###process.load('Configuration/StandardSequences/SimL1Emulator_cff')
###process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
###process.load('Configuration/StandardSequences/EndOfProcess_cff')
##process.load('Configuration/EventContent/EventContent_cff')
###process.load('Configuration.StandardSequences.ReconstructionCosmics_cff')
##
# global tag
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
#process.GlobalTag.globaltag = 'GR09_R_35X_V3::All'
process.GlobalTag.globaltag = 'GR10_P_V4::All'

# message logger
process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cout.placeholder = cms.untracked.bool(False)
#process.MessageLogger.cout.threshold = cms.untracked.string('WARNING')
#process.MessageLogger.cout.threshold = cms.untracked.string('INFO')
process.MessageLogger.cerr.threshold = cms.untracked.string('ERROR')
process.MessageLogger.debugModules = cms.untracked.vstring('*')

#Tracer: uncomment helps debugging
#process.Tracer = cms.Service("Tracer")

# how many events?
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.load("Configuration.StandardSequences.VtxSmearedBetafuncEarlyCollision_cff")
process.load("Configuration.StandardSequences.Generator_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.L1Emulator_cff")

# uncommenting does not work...
#process.load("Configuration.StandardSequences.L1TriggerDefaultMenu_cff")

#run HLT (the load option only will crash the script)
from HLTrigger.Configuration.HLT_8E29_cff import *



# ------------ Trig Eff ----------------
# example for J/Psi
process.te = cms.EDAnalyzer("TrigEff",
                            process.MuonServiceProxy,
genTag     = cms.InputTag("genParticles"),
L1extraTag = cms.InputTag("l1extraParticles"),
Level1Tag  = cms.InputTag("null"),#using l1extra
Level2Tag  = cms.InputTag("hltSingleMu9L2Filtered7"),
Level3Tag  = cms.InputTag("hltSingleMu9L3Filtered9"),
tracksTag  = cms.InputTag("generalTracks"),
saMuonsTag = cms.InputTag("standAloneMuons","UpdatedAtVtx"),
muonsTag   = cms.InputTag("muons"),
#csctfTag = cms.InputTag("null"),
#csctfLctsTag = cms.InputTag("null"),
csctfTag     = cms.InputTag("csctfDigis"),
csctfLctsTag = cms.InputTag("csctfDigis"),
TriggerEventTag = cms.InputTag("null"),
HLTriggerTag = cms.InputTag("null"),
eta_cut= cms.double(2.4),
pt_cut= cms.double(0.0),
matching_cone= cms.double(0.2),
minvLow  = cms.double(0),
minvHigh = cms.double(20),
outputFile = cms.string("MinimumBias_trigEff-run132440_1.root"),
gpTag = cms.InputTag("hltTriggerSummaryAOD","","HLT"),
## mimic HLTEventAnalyzerRaw
processName = cms.string("HLT"),
#triggerName = cms.string("@"), # all paths
triggerName = cms.string("HLT_Mu0_L1MuOpen"),
triggerResults = cms.InputTag("TriggerResults","","HLT"),
triggerEventWithRefs = cms.InputTag("hltTriggerSummaryRAW","","HLT"),
# mimic L1GTUtils
AlgorithmName = cms.string('L1_SingleMu7'),
#
level2module  = cms.string("hltMu0L1MuOpenL2Filtered0"),
level3module  = cms.string("hltMu0L1MuOpenL3Filtered0"),
# printLevel
printLevel = cms.untracked.int32(1)
)
# --------- Trig Eff END ----------------



# skimming taken from DPGAnalysis/Skims/python/MinBiasPDSkim_cfg.py?revision=1.15
#===============================================================================
# CSC Activity
from DPGAnalysis.Skims.CSCSkim_cfi import cscSkim
process.cscSkimLower = cscSkim.clone()
#set to minimum activity (1 segment and 1 hit in one of the chambers)
process.cscSkimLower.minimumSegments = 1
process.cscSkimLower.minimumHitChambers = 1
process.cscSkim = cms.Path(process.cscSkimLower)
#===============================================================================

#===============================================================================
# GOOD COLLISIONS
process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskAlgoTrigConfig_cff')
process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff')
process.load('HLTrigger/HLTfilters/hltLevel1GTSeed_cfi')

process.L1T1coll=process.hltLevel1GTSeed.clone()
process.L1T1coll.L1TechTriggerSeeding = cms.bool(True)
process.L1T1coll.L1SeedsLogicalExpression = cms.string('0 AND (40 OR 41) AND NOT (36 OR 37 OR 38 OR 39) AND NOT ((42 AND NOT 43) OR (43 AND NOT 42))')
#===============================================================================

#===============================================================================
# Good Vertex Filter
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
vertexCollection = cms.InputTag('offlinePrimaryVertices'),
minimumNDOF = cms.uint32(4) ,
maxAbsZ = cms.double(15),
maxd0 = cms.double(2)
)
#===============================================================================

#===============================================================================
# No Scraping
process.scraping = cms.EDFilter("FilterOutScraping",
applyfilter = cms.untracked.bool(True),
debugOn = cms.untracked.bool(True),
numtrack = cms.untracked.uint32(10),
thresh = cms.untracked.double(0.2))
#===============================================================================


#===============================================================================
process.load( "HLTrigger.HLTcore.triggerSummaryAnalyzerAOD_cfi" )

# triggerSummaryAnalyzerAOD,
# singleMuNoIso & hltEnd,
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

#from PhysicsTools.HepMCCandAlgos.genParticleCandidatesFast_cfi import *
#process.p = cms.Path(process.genParticleCandidates*process.te)

### L1 Extra needed for the analysis
##process.load('L1Trigger.Configuration.L1Extra_cff')
##
### path
##process.p = cms.Path(process.L1Extra*process.te)

process.p = cms.Path(
    process.L1T1coll
    +process.primaryVertexFilter
    +process.scraping
    +process.cscSkimLower
    #+process.gtDigis
    #+process.gtEvmDigis
    +process.csctfDigis
    +process.te)

# ------------ PoolSource -------------
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()
process.source = cms.Source ("PoolSource",
fileNames = readFiles,
secondaryFileNames = secFiles
)
# -------- PoolSource END -------------

## readFiles.extend([
## 
## 
##     '/store/data/Commissioning10/MinimumBias/RECO/Apr1ReReco-v2/0139/3ADE63D6-923E-DF11-B92A-001A92971BD8.root',
##     '/store/data/Commissioning10/MinimumBias/RECO/Apr1ReReco-v2/0133/6E74370E-6D3E-DF11-890B-0018F3D09684.root',
##     '/store/data/Commissioning10/MinimumBias/RECO/Apr1ReReco-v2/0132/C851858F-613E-DF11-A403-00261894387E.root',
##     '/store/data/Commissioning10/MinimumBias/RECO/Apr1ReReco-v2/0132/C0938B90-613E-DF11-8F45-00261894386D.root',
##     '/store/data/Commissioning10/MinimumBias/RECO/Apr1ReReco-v2/0132/9A9D3396-613E-DF11-8643-00261894398B.root',
##     '/store/data/Commissioning10/MinimumBias/RECO/Apr1ReReco-v2/0132/02CA6895-613E-DF11-A4F0-003048678B92.root',
##     '/store/data/Commissioning10/MinimumBias/RECO/Apr1ReReco-v2/0131/E024013C-5D3E-DF11-8F34-0018F3D09682.root',
##     '/store/data/Commissioning10/MinimumBias/RECO/Apr1ReReco-v2/0131/DE499234-603E-DF11-94A6-003048678B06.root',
##     '/store/data/Commissioning10/MinimumBias/RECO/Apr1ReReco-v2/0131/BAFFF7F1-5A3E-DF11-9F13-00248C0BE018.root',
##     '/store/data/Commissioning10/MinimumBias/RECO/Apr1ReReco-v2/0131/A8F6AE35-5D3E-DF11-8601-002354EF3BE0.root',
##     '/store/data/Commissioning10/MinimumBias/RECO/Apr1ReReco-v2/0131/8C884021-5B3E-DF11-9D5E-003048678ADA.root',
##     '/store/data/Commissioning10/MinimumBias/RECO/Apr1ReReco-v2/0131/7AA74636-603E-DF11-A5DE-003048678F8A.root',
##     '/store/data/Commissioning10/MinimumBias/RECO/Apr1ReReco-v2/0131/62DE6421-5B3E-DF11-B53B-003048678FA6.root',
##     '/store/data/Commissioning10/MinimumBias/RECO/Apr1ReReco-v2/0131/58C598F2-5A3E-DF11-9B89-00304867BEE4.root',
##     '/store/data/Commissioning10/MinimumBias/RECO/Apr1ReReco-v2/0131/40EA1340-5D3E-DF11-B49E-0026189438D7.root',
##     '/store/data/Commissioning10/MinimumBias/RECO/Apr1ReReco-v2/0130/FAED0412-5A3E-DF11-B5B7-002618943870.root',
##     '/store/data/Commissioning10/MinimumBias/RECO/Apr1ReReco-v2/0130/EC8E3A0D-5A3E-DF11-A4CC-002618943900.root',
##     '/store/data/Commissioning10/MinimumBias/RECO/Apr1ReReco-v2/0130/D2DB4C0A-593E-DF11-97AE-00304867902C.root',
##     '/store/data/Commissioning10/MinimumBias/RECO/Apr1ReReco-v2/0130/C82F420B-5A3E-DF11-B9FD-002618943948.root',
##     '/store/data/Commissioning10/MinimumBias/RECO/Apr1ReReco-v2/0130/BED9F431-593E-DF11-925F-003048D15E24.root',
##     '/store/data/Commissioning10/MinimumBias/RECO/Apr1ReReco-v2/0130/98FEF80A-5A3E-DF11-B1AF-00248C0BE018.root',
##     '/store/data/Commissioning10/MinimumBias/RECO/Apr1ReReco-v2/0130/7CF5670C-5A3E-DF11-8C52-003048678F1C.root',
##     '/store/data/Commissioning10/MinimumBias/RECO/Apr1ReReco-v2/0130/5AAA5A6D-553E-DF11-AEE5-002354EF3BE6.root',
##     '/store/data/Commissioning10/MinimumBias/RECO/Apr1ReReco-v2/0130/50BFA301-593E-DF11-9B2E-00304867918E.root',
##     '/store/data/Commissioning10/MinimumBias/RECO/Apr1ReReco-v2/0130/50178414-593E-DF11-851A-0030486790B0.root',
##     '/store/data/Commissioning10/MinimumBias/RECO/Apr1ReReco-v2/0130/14C3AF58-563E-DF11-B19B-00261894396F.root',
##     '/store/data/Commissioning10/MinimumBias/RECO/Apr1ReReco-v2/0130/1092AD71-973E-DF11-8F4E-0018F3D096D8.root',
##     '/store/data/Commissioning10/MinimumBias/RECO/Apr1ReReco-v2/0129/58274A60-523E-DF11-A311-0030486792F0.root',
##     '/store/data/Commissioning10/MinimumBias/RECO/Apr1ReReco-v2/0129/322D2491-533E-DF11-9711-003048679000.root' 
##  ]
## 
## )
## 
## secFiles.extend( [
##  ]
## 
## )

readFiles.extend([
    'rfio:/castor/cern.ch/user/d/digiovan/Collisions2010/MinimumBiasSkim/132440/root/MinimumBiasSkim-run132440_1.root'
 ]
)

secFiles.extend([])

