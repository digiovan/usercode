import FWCore.ParameterSet.Config as cms

process = cms.Process("UFDiMuonAnalyzer")

from PhysicsTools.PatAlgos.patTemplate_cfg import *
from PhysicsTools.PatAlgos.tools.coreTools import *
from PhysicsTools.PatAlgos.tools.pfTools import *

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
    process.GlobalTag.globaltag = "START52_V5::All"


# ------------ PoolSource -------------
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring('file:DYJetsToLL_2_1_RQm.root')
                            )
#process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
# -------- PoolSource END -------------


#===============================================================================
process.load("UserCode.UFDiMuonsAnalyzer.UFDiMuonAnalyzer_nocuts_cff")
process.dimuons = process.DiMuons.clone()
process.dimuons.getFilename = cms.untracked.string("DYJetsToLL.root")
process.dimuons.isMonteCarlo = cms.bool(True) 

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
process.dimuons.processName    = cms.string("HLT")
process.dimuons.triggerNames   = cms.vstring("HLT_Mu30_eta2p1","HLT_Mu40_eta2p1")
process.dimuons.triggerResults = cms.InputTag("TriggerResults","","HLT")
process.dimuons.triggerEvent   = cms.InputTag("hltTriggerSummaryAOD","","HLT")
## 
## process.dimuons.metTag         = cms.InputTag("pfMet")
## process.dimuons.pfJetsTag      = cms.InputTag("cleanPatJets")
## process.dimuons.genJetsTag     = cms.InputTag("null")

## ADDING PAT
removeMCMatching(process, ['All'])

if thisIsData:
    print "\nData Jet Corrections"
    switchJetCollection(process,cms.InputTag('ak5PFJets'),
                        doJTA        = True,
                        doBTagging   = True,
                        jetCorrLabel = ('AK5PF', cms.vstring([ 'L1Offset', 'L2Relative', 'L3Absolute', 'L2L3Residual'])),
                        doType1MET   = True,
                        #   genJetCollection=cms.InputTag("ak5GenJets"),
                        doJetID      = True
                        )

else:
    print "\nMC Jet Corrections"
    switchJetCollection(process,cms.InputTag('ak5PFJets'),
                        doJTA        = True,
                        doBTagging   = True,
                        jetCorrLabel = ('AK5PF', cms.vstring(['L2Relative', 'L3Absolute' ])),
                        doType1MET   = True,
                        genJetCollection=cms.InputTag("ak5GenJets"),
                        doJetID      = True
                        )

# Clean the Jets from the seleted leptons, and apply loose(???) btag cuts and loose id cuts!
process.cleanPatJets = cms.EDProducer("PATJetCleaner",
          src = cms.InputTag("patJets"),
          preselection = cms.string('( (neutralEmEnergy/energy < 0.99) &&  (neutralHadronEnergy/energy < 0.99) && numberOfDaughters>1) '),
          checkOverlaps = cms.PSet(
             muons = cms.PSet(
               #src       = cms.InputTag("userDataMuons"),
               src       = cms.InputTag("muons"),
               algorithm = cms.string("byDeltaR"),
               preselection        = cms.string(""),
               deltaR              = cms.double(0.5),
               checkRecoComponents = cms.bool(False),
               pairCut             = cms.string(""),
               requireNoOverlaps   = cms.bool(True),
             ),
         ),
###applying 2 loose cutson the b-tagged jets
         finalCut = cms.string('')
)


#===============================================================================

process.p = cms.Path( process.patDefaultSequence*
                      process.cleanPatJets*
                      process.dimuons
                    )

process.outpath = cms.EndPath()
#===============================================================================

