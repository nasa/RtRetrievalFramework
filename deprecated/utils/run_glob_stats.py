#!/usr/bin/env python

import sys
import re
import os
import glob
from optparse import OptionParser

sys.path.append('/home/mcduffie/oco_l2/auto_testing/generator')

import L2_Input

def get_out_info_data(tcDirList):

    maxVal = float(sys.maxint)
    minVal = float(-sys.maxint - 1)

    rmsMin = [maxVal for x in range(3)]
    chiSqMin = [maxVal for x in range(3)]
    dSigmaSqMin = maxVal

    rmsMax = [minVal for x in range(3)]
    chiSqMax = [minVal for x in range(3)]
    dSigmaSqMax = minVal

    foundOne = False
    for tcDir in tcDirList:
        infoFile = tcDir + '/out/out_info.dat'
        
        if os.path.exists(infoFile):
            foundOne = True
            
            infoObj = L2_Input.Input_File(infoFile)
            infoMatrix = infoObj.rootNode.Get_Matrix_Data()
            numRows = len(infoMatrix)

            dataRow = infoMatrix[numRows-1]

            rmsData      = [float(x) for x in dataRow[3:6]]
            chiSqData    = [float(x) for x in dataRow[6:9]]
            dSigmaSqData = float(dataRow[9])

            for rmsIdx in range(len(rmsData)):
                if rmsData[rmsIdx] < rmsMin[rmsIdx]:
                    rmsMin[rmsIdx] = rmsData[rmsIdx]
                if rmsData[rmsIdx] > rmsMax[rmsIdx]:
                    rmsMax[rmsIdx] = rmsData[rmsIdx]

            for chiSqIdx in range(len(chiSqData)):
                if chiSqData[chiSqIdx] < chiSqMin[chiSqIdx]:
                    chiSqMin[chiSqIdx] = chiSqData[chiSqIdx]
                if chiSqData[chiSqIdx] > chiSqMax[chiSqIdx]:
                    chiSqMax[chiSqIdx] = chiSqData[chiSqIdx]

            if dSigmaSqData < dSigmaSqMin:
                dSigmaSqMin = dSigmaSqData
            if dSigmaSqData > dSigmaSqMax:
                dSigmaSqMax = dSigmaSqData

    if not foundOne:
        return ['NO_OUT_INFO' for x in range(14)]

    outputRanges = []
    for rmsIdx in range(len(rmsMin)):
        outputRanges.append('%.2e' % rmsMin[rmsIdx])
        outputRanges.append('%.2e' % rmsMax[rmsIdx])

    for chiSqIdx in range(len(chiSqMin)):
        outputRanges.append('%.2e' % chiSqMin[chiSqIdx])
        outputRanges.append('%.2e' % chiSqMax[chiSqIdx])

    outputRanges.append('%.2e' % dSigmaSqMin)
    outputRanges.append('%.2e' % dSigmaSqMax)

    return outputRanges

def get_log_info_data(tcDirList, logDir):

    maxVal = float(sys.maxint)
    minVal = float(-sys.maxint - 1)

    xTargMin = maxVal
    xTargMax = minVal

    radRuntimeMin = maxVal
    radRuntimeMax = minVal

    numRadCallsMin = maxVal
    numRadCallsMax = minVal

    foundOne = False
    for tcDir in tcDirList:
        tcName = (re.findall('testcase_[^/]+', tcDir))[0]

        logFile = ''
        tcLogFileSearch = glob.glob(tcDir + '/stdout*')
        logDirFileSearch = glob.glob(logDir + '/stdout*' + tcName + '*')

        if len(tcLogFileSearch)  > 0:
            logFile = tcLogFileSearch[0]
        elif len(logDirFileSearch)  > 0:
            logFile = logDirFileSearch[0]
        else:
            continue
        
        foundOne = True

        logFileObj = open(logFile, 'r')
        logFile = logFileObj.readlines()
        logFileObj.close()

        xTargVal = 0
        numRadCallsVal = 0
        for logLine in logFile:
            logLine = logLine.strip()
            if logLine.find('Xtarg') >= 0:
                lineParts = logLine.split()
                xTargVal = float(lineParts[4])
            if logLine.find('calculate_radiance') >= 0:
                numRadCallsVal = numRadCallsVal + 1

        errFile = None
        tcErrFileSearch = glob.glob(tcDir + '/stderr*')
        errDirFileSearch = glob.glob(logDir + '/stderr*' + tcName + '*')

        if len(tcErrFileSearch)  > 0:
            errFile = tcErrFileSearch[0]
        elif len(errDirFileSearch)  > 0:
            errFile = errDirFileSearch[0]

        radRuntimeVal = 0
        if errFile != None:
            errFileObj = open(errFile, 'r')
            errFile = errFileObj.readlines()
            errFileObj.close()

            for errLine in errFile:
                errLine = errLine.strip()
                if errLine.find('Total time spent in radiant') >= 0:
                    lineParts = errLine.split()
                    radRuntimeVal = radRuntimeVal + float(lineParts[6])

        if xTargVal < xTargMin:
            xTargMin = xTargVal
        if xTargVal > xTargMax:
            xTargMax = xTargVal

        if radRuntimeVal < radRuntimeMin:
            radRuntimeMin = radRuntimeVal
        if radRuntimeVal > radRuntimeMax:
            radRuntimeMax = radRuntimeVal

        if numRadCallsVal < numRadCallsMin:
            numRadCallsMin = numRadCallsVal
        if numRadCallsVal > numRadCallsMax:
            numRadCallsMax = numRadCallsVal

    if not foundOne:
        return ['NO_LOG_FILE' for x in range(6)]

    logResults = []

    logResults.append('%.2e' % xTargMin)
    logResults.append('%.2e' % xTargMax)

    logResults.append('%.2e' % radRuntimeMin)
    logResults.append('%.2e' % radRuntimeMax)

    logResults.append('%d' % numRadCallsMin)
    logResults.append('%d' % numRadCallsMax)

    return logResults

#############
# Main level

# Set up command line arguments
parser = OptionParser("usage: %prog [options] glob_file list_file [list_file...]")

parser.add_option( "-c", "--count", dest="count",
                   action="store_true",
                   help="count the occurance of a blob in each list file",
                   default=False)

parser.add_option( "-o", "--out_info", dest="out_info",
                   action="store_true",
                   help="display ranges for parameters from out_info.dat",
                   default=False)

parser.add_option( "-l", "--log_info", dest="log_info",
                   action="store_true",
                   help="display ranges for parameters from log files",
                   default=False)

parser.add_option( "-d", "--log_file_dir", dest="log_dir",
                   metavar="DIR",
                   help="directory where log files are located if not found in the testcase directory",
                   default="./")

# Parse command line arguments
(options, args) = parser.parse_args()

if len(args) < 2:
    parser.error("incorrect number of arguments")

# Get filenames used for gathering stats
globListFilename = args[0]
listFilenames    = args[1:]

# Load the glob file into memory
globListObj = open(globListFilename, 'r')
globList = globListObj.readlines()
globListObj.close()

# Create format string for glob text
globListMaxLen = max([len(x) for x in globList])
globFormat = "%" + str(globListMaxLen) + "s\t"

# Create format for list filenames
listFilesMaxLen = max([len(x) for x in listFilenames])
listStrFormat = "%" + str(listFilesMaxLen) + 's\t'
listNumFormat = "%" + str(listFilesMaxLen) + 'd\t'

print globFormat % 'glob',
print 'total\t',

# If counting then list names of files being counted
if( options.count ):
    for listName in listFilenames:
        print listStrFormat % listName,

outInfoStrFormat = '%11s\t%11s\t%11s\t%11s\t%11s\t%11s\t%11s\t%11s\t%11s\t%11s\t%11s\t%11s\t%14s\t%14s\t'
if options.out_info:
    print outInfoStrFormat % ('rms_1_min', 'rms_1_max',
                              'rms_2_min', 'rms_2_max',
                              'rms_3_min', 'rms_3_max',
                              'chi2_1_min', 'chi2_1_max',
                              'chi2_2_min', 'chi2_2_max',
                              'chi2_3_min', 'chi2_3_max',
                              'd_sigma_sq_min', 'd_sigma_sq_max'),

logInfoStrFormat = '%11s\t%11s\t%11s\t%11s\t%13s\t%13s\t'
if options.log_info:
    print logInfoStrFormat % ('xtarg_min', 'xtarg_max',
                              'rad_time_min', 'rad_time_max',
                              'rad_calls_min', 'rad_calls_max'
                              ),

print ''


for globValue in globList:
    globValue = globValue.strip()
    globRe = globValue.replace('.', '\.').replace('*', '.*')

    listCounts = []
    allMatches = []
    for listFile in listFilenames:
        listFileObj = open(listFile, 'r')
        listContents = listFileObj.readlines()
        listFileObj.close()

        listMatches = []
        for listLine in listContents:
            listLine = listLine.strip()
            reResult = re.findall(globRe, listLine)

            if len(reResult) > 0:
                listMatches.append(listLine)
                allMatches.append(listLine)

        listCounts.append(len(listMatches))

    print globFormat % globValue,
    print '%5d\t' % len(allMatches), 
    
    if options.count:
        for listIdx in range(len(listFilenames)):
            print listNumFormat % listCounts[listIdx],

    if options.out_info:
        outInfoData = get_out_info_data(allMatches)
        print outInfoStrFormat % tuple(outInfoData),

    if options.log_info:
        logInfoData = get_log_info_data(allMatches, options.log_dir)
        print logInfoStrFormat % tuple(logInfoData),
            
    print ''

