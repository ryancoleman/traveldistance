#!/usr/bin/env python2.5
#ryan g. coleman ryangc@mail.med.upenn.edu

import sys,string,os

import tstCreate
import tstConvexHull
import oldTravelDist
import tstdata
import cavity
from tstTravelDepth import runTravelDepthCompletely

#this is main
if -1 != string.find(sys.argv[0], "tstCreateNoChunnel.py"):
  if 6 < len(sys.argv):
    pdbFileName = sys.argv[1]
    gridSize = sys.argv[2]
    whichSurface = sys.argv[3]
    probe = float(sys.argv[4])
    radscale = sys.argv[5]
    path = sys.argv[6]
    runTravelDepthCompletely(pdbFileName, gridSize, whichSurface, \
                             probe, radscale, path, check2trisPerEdge=False, \
                             doTravelDepth=False)
  elif 5 < len(sys.argv):
    pdbFileName = sys.argv[1]
    gridSize = sys.argv[2]
    whichSurface = sys.argv[3]
    probe = float(sys.argv[4])
    radscale = sys.argv[5]
    runTravelDepthCompletely(pdbFileName, gridSize, whichSurface, \
                             probe, radscale, check2trisPerEdge=False, \
                             doTravelDepth=False)
  elif 4 < len(sys.argv):
    pdbFileName = sys.argv[1]
    gridSize = sys.argv[2]
    whichSurface = sys.argv[3]
    probe = float(sys.argv[4])
    runTravelDepthCompletely(pdbFileName, gridSize, whichSurface, probe, \
                             check2trisPerEdge=False, \
                             doTravelDepth=False)
  elif 3 < len(sys.argv):
    pdbFileName = sys.argv[1]
    gridSize = sys.argv[2]
    whichSurface = sys.argv[3]
    runTravelDepthCompletely(pdbFileName, gridSize, whichSurface, \
                             check2trisPerEdge=False, doTravelDepth=False)
  elif 2 < len(sys.argv):
    pdbFileName = sys.argv[1]
    gridSize = sys.argv[2]
    runTravelDepthCompletely(pdbFileName, gridSize, check2trisPerEdge=False, \
                             doTravelDepth=False)
  elif 1 < len(sys.argv):
    pdbFileName = sys.argv[1]
    runTravelDepthCompletely(pdbFileName, check2trisPerEdge=False, \
                             doTravelDepth=False)
  else:
    print "Usage: tstCreateNoChunnel.py file.pdb [gridsize] [tri|mesh] [probe] [radius scale] [pathToExecs]"
