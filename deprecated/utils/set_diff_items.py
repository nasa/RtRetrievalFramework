#!/usr/bin/env python

import sys

if (len(sys.argv) < 2):
    print "usage:\n\t", sys.argv[0], " <complete_set_file> <subset_file>\n"
    print "Displays the items of complete_set_file that are not in subset_file\n"
    sys.exit(1)


allItemsFile = sys.argv[1]
subItemsFile = sys.argv[2]

allItemsFObj = open(allItemsFile, 'r')
subItemsFObj = open(subItemsFile, 'r')

allItemsList = allItemsFObj.readlines()
subItemsList = subItemsFObj.readlines()

for currItem in allItemsList:
    if currItem not in subItemsList:
        print currItem,
