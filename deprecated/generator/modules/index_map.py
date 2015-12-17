import os
import logging

def Process_File(source, destination, fileKeywords, moduleSections, valuesDict, mapDict):
   logger = logging.getLogger(os.path.basename(__file__))

   if type(source) is str:
      srcFileObj = open(source, 'r')
   elif hasattr(source, 'read'):
      srcFileObj = source
   else:
      raise Exception('Unrecognized source object: %s' % source)

   if type(destination) is str:
      dstFileObj = open(destination, 'w')
   elif hasattr(destination, 'write'):
      dstFileObj = destination
   else: 
      raise Exception('Unrecognized destination object: %s' % destination)     

   file_idx = 0
   for src_line in srcFileObj.readlines():
       print >>dstFileObj, '%s\t%d' % (src_line.strip(), file_idx)
       file_idx += 1

   if type(source) is str:
      srcFileObj.close()
   if type(destination) is str:
      dstFileObj.close()
