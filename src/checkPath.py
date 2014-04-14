#!/usr/bin/env python

import string
import sys
import comparePaths
import paths
import tstTopology
import tstdata
import tstdebug

def checkPathTst(tstFileName, pathFileName):
  path = comparePaths.readCGOPath(pathFileName)
  tstD = tstdata.tstData(tstFileName)
  loopTris, loopPoints = tstTopology.getLoopsAndVerticesTst(tstD)
  throughTris, throughPts = paths.checkPathPoints(
      path, loopPoints, tstD.dict['POINT_XYZ'])
  if throughTris:  # it worked... make some more debugging output
    tstdebug.debugTrianglesNotOrig(
        throughTris, tstD.dict['POINT_XYZ'],
        tstFileName + ".path.through.loop.py", ptList=throughPts)
    print "path goes through a hole"
  else:
    print "path does NOT go through a hole"

#main
if -1 != string.find(sys.argv[0], "checkPath"):
  if len(sys.argv) >= 3:
    tst = sys.argv[1]
    path = sys.argv[2]
    checkPathTst(tst, path)
  else:
    print "Usage: checkPath.py tstFile pathFile"
    print "Checks whether path goes through any hole in the protein"
