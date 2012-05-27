import FWCore.ParameterSet.Config as cms

process = cms.Process("UFMuonEle")

thisIsData = False

if thisIsData:
    print 'Running over data sample'
else:
    print 'Running over MC sample'



process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.GeometryExtended_cff")
process.load('Configuration.EventContent.EventContent_cff')


# global tag
if thisIsData:
    print 'Loading Global Tag For Data'
    process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
    process.GlobalTag.globaltag = "GR_R_42_V19::All"
else:
    print 'Loading Global Tag For MC'
    process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
    process.GlobalTag.globaltag = "START52_V5::All"#new JEC correction w.r.t. V11


# ------------ PoolSource -------------
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring())
#process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring('file:TTJets_1_1_ReS.root'))
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
process.source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)

# list of files
# from 
readFiles.extend([
    'file:TTJets_1_1_mkk.root'
])
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
#process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange()
# -------- PoolSource END -------------


#===============================================================================
#process.load("ZPtToMuMu.UFDiMuonsAnalyzer.UFDiMuonAnalyzer_nocuts_cff")
process.load("ZPtToMuMu.UFMuonElectronAnalyzer.UFMuonElectronAnalyzer_cff")
process.muonele = process.MuonElectron.clone()
process.muonele.isVerbose = cms.untracked.bool(True)
process.muonele.isMonteCarlo = cms.bool(True)
process.muonele.getFilename = cms.untracked.string("test.root")

#===============================================================================

process.p = cms.Path( process.muonele )

process.outpath = cms.EndPath()
#===============================================================================

