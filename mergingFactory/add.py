#!/bin/python
import FWCore.ParameterSet.Config as cms
import sys,os

#
# first argment is cmsRun
# second argument is the script that's called
# third argument is the file name for the list to be loaded
#
filelist = "test"
if len(sys.argv) != 3: 
	print "not enough command line arguments!"
	print "number of arguments: %d"%len(sys.argv)
	for arg in sys.argv: print arg
else: 
	filelist = sys.argv[2]
	
#
# run the main program
#
process = cms.Process("ADDFILES1")
process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 500 
#
# load necessary processes
#
process.load(filelist)

## process.options = cms.untracked.PSet(
##         fileMode = cms.untracked.string('NOMERGE')
## )


process.source.inputCommands = cms.untracked.vstring( "keep *" )

process.out1 = cms.OutputModule("PoolOutputModule",
                                fileName = cms.untracked.string("/data/uftrig01b/digiovan/root/CMSSW_5_2_3_patch1/SkimsSingleMuRun2012A-PromptReco-v1_AOD/mergedFiles/"+filelist+".root"),
                                )
                                
#process.out1.fileName = "test.root"
process.outpath = cms.EndPath(process.out1)
