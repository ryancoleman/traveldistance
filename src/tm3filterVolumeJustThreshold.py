#!/usr/bin/env python2.5

import tm3, string, sys

def deleteLeaves(tmData, attributeName, threshold):
  column = tmData.attributeTitles.index(attributeName)
  for node in tmData.bottomUpTraversal(): 
    if node in tmData.tree: #may have already been merged
      if tmData.isLeaf(node.getId()):
        if node.attributes[column] < threshold:
          tmData.deleteLeaf(node.getId())
    
#main is only run for testing import of tm3 from commandline        
if -1 != string.find(sys.argv[0], "tm3filterVolumeJustThreshold.py"):
  try:
    threshold = float(sys.argv[1])
    if 2 == len(sys.argv):
      print "no input files specified"
    else:
      for filename in sys.argv[2:]:
        tmData = tm3.tmTreeFromFile(filename)
        deleteLeaves(tmData, "Volume", threshold)
        tmData.write(filename+".v"+str(threshold)+".jt.tm3") 
  except (TypeError, IndexError):
    print "tm3filterVolumeJustThreshold.py threshold file [more files]"

