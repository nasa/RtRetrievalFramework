#!/usr/bin/env python

import glob
import sys
import shutil

tcDirBase = '/home/mcduffie/oco_l2/auto_testing/testcases'
specSubDir = 'in/l1b/spec/OCO/'

if (len(sys.argv) < 4):
    print "usage:\n", sys.argv[0], "<glob_list_file> <source_testset_dir> <dest_testset_dir>\n"
    print "Copies L1B spectra generared by FM testcases into their retrieval"
    print "testcase directories"
    sys.exit(1)

tcNameGlobFile = sys.argv[1]

testsetSourceDir = tcDirBase + '/' + sys.argv[2]
testsetDestDir = tcDirBase + '/' + sys.argv[3]

globListFileObj = open(tcNameGlobFile, 'r')
globList = globListFileObj.readlines()
globListFileObj.close()

for globTxt in globList:
    globTxt = globTxt.strip()
    
    if globTxt == None or len(globTxt) == 0:
        next

    sourceGlob = testsetSourceDir + '/' + globTxt
    destGlob   = testsetDestDir + '/' + globTxt
    
    sourceDir = glob.glob( sourceGlob )
    destDirs = glob.glob( destGlob )

    if len(sourceDir) > 1:
        raise Exception, "Multiple source directories matched from glob: " + sourceGlob
    elif len(sourceDir) == 0:
        raise Exception, "Could not find source directories for glob: " + sourceGlob

    for destTCDir in destDirs:
        print sourceDir[0] + ' => ' + destTCDir
        
        shutil.rmtree(destTCDir + '/' + specSubDir)
        shutil.copytree(sourceDir[0] + '/' + specSubDir, destTCDir + '/' + specSubDir)
