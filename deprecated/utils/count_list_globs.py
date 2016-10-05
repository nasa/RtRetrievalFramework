#!/usr/bin/env python

import sys
import re

if len(sys.argv) < 3:
    print 'usage: ' + sys.argv[0] + ' <glob_file> [<list_file>...]'
    sys.exit(0)

globListFilename = sys.argv[1]
listFilenames = sys.argv[2:]

listFilesMaxLen = max([len(x) for x in listFilenames])

globListObj = open(globListFilename, 'r')
globList = globListObj.readlines()
globListObj.close()

globListMaxLen = max([len(x) for x in globList])

globFormat = "%" + str(globListMaxLen) + "s\t"

listStrFormat = "%" + str(listFilesMaxLen) + 's\t'
listNumFormat = "%" + str(listFilesMaxLen) + 'd\t'

print globFormat % 'glob',
for listName in listFilenames:
    print listStrFormat % listName,
print ''

for glob in globList:
    glob = glob.strip()
    globRe = glob.replace('.', '\.').replace('*', '.*')

    listCounts = []
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

        listCounts.append(len(listMatches))
            
    print globFormat % glob,
    for listIdx in range(len(listFilenames)):
        print listNumFormat % listCounts[listIdx],
    print ''
