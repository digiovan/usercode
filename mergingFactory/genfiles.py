#!/bin/python
import ROOT
from ROOT import *
def prestuff(thefile):
	thefile.write("import FWCore.ParameterSet.Config as cms\n")
	thefile.write("source = cms.Source(\"PoolSource\",\n")
	thefile.write("fileNames = cms.untracked.vstring(\n")


inFile = open("thefiles.list", "r")

ifile = 1
iout = 1

# cms523_patch1
outFile = open("Run2012A-PromptReco-v1_%05d.py"%int(iout),"w")
print "nohup cmsRun add.py Run2012A-PromptReco-v1_%05d >& logfile_Run2012A-PromptReco-v1_%05d < /dev/null &" % (int(iout), int(iout))

prestuff(outFile)
for line in inFile.readlines():
	file = line.rstrip("\n")
		
	outFile.write("'file:"+file+"\',\n")
	ifile += 1
	#if ifile == 10:
        if ifile == 151:
		ifile = 1
		iout += 1
		outFile.write("))\n")
		outFile.close()

                # CMSSW_5_2_3_patch1
                outFile = open("Run2012A-PromptReco-v1_%05d.py"%int(iout),"w")
                print "nohup cmsRun add.py Run2012A-PromptReco-v1_%05d >& logfile_Run2012A-PromptReco-v1_%05d < /dev/null &" % (int(iout), int(iout))

                prestuff(outFile)

outFile.write("))\n")
outFile.close()


