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

# global tag
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
#process.GlobalTag.globaltag = 'GR10_P_V4::All'
process.GlobalTag.globaltag = 'GR_R_35X_V8B::All'

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
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

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
muonsTag   = cms.InputTag("muons"),
csctfTag     = cms.InputTag("csctfDigis"),
csctfLctsTag = cms.InputTag("csctfDigis"),
outputFile = cms.string("MinimumBiasTrigEffNtuple-run132440to134987_1.root"),
cscSegTag = cms.InputTag("cscSegments"), 
# mimic L1GTUtils
AlgorithmName = cms.string('L1_SingleMu7'),
# printLevel
printLevel = cms.untracked.int32(1)
)
# --------- Trig Eff END ----------------


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

process.p = cms.Path(process.csctfDigis
+process.te)

# ------------ PoolSource -------------
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()
process.source = cms.Source ("PoolSource",
fileNames = readFiles,
secondaryFileNames = secFiles#,
#eventsToProcess =cms.untracked.VEventRange('132601:24017629')
)
# -------- PoolSource END -------------

readFiles.extend([

    #'rfio:/castor/cern.ch/user/d/digiovan/Skims/MinimumBiasMay6thPDSkim/MinBiasSkim_167_1.root'
    #'rfio:/castor/cern.ch/user/d/digiovan/Skims/MinimumBiasMay6thPDSkim_2attempt/MinBiasSkim_110_1.root'
    'rfio:/castor/cern.ch/user/d/digiovan/Skims/MinimumBiasMay6thPDSkim_2attempt/MinBiasSkim_205_1.root'
 ]

)

secFiles.extend( [
 ]

)

