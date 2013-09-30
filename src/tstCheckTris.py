#!/usr/bin/env python2.5
#script to read .tst format and make surfaces in pymol
#ryan coleman

#from pymol.cgo import *
#from pymol import cmd
import math,  string
#for running commandline instead of through pymol
import sys, os
from tstdata import tstData

  

#this function is the opening to pymol to open a tst file, read in the data
def tstOpen(filename):
  tstData(filename)

def tstCheckTris(tstFileName):
  tstD = tstData(tstFileName) #read the file into the data structure
  okay = True
  pointTris = tstD.dict['POINT_TRIANGLE']
  for pointTriRec in pointTris:
    if okay and pointTriRec[1] >= 13:
      okay = False
  if okay:
    print tstFileName + " checks out okay"
    return False
  else:
    print tstFileName + " FAILS check, 13 or more triangles share one point!"
    return tstFileName

#this is where main will go
if len(sys.argv) > 1: #else do nothing, read in as module
  badlist = []
  for tstFile in sys.argv[1:]:
    checked = tstCheckTris(tstFile)
    if checked:
      badlist.append(checked)
  if len(badlist) > 0:
    print str(len(badlist)) + " files did  not check out"
    for bad in badlist:
      print bad
else:
  print "Usage: tstCheckTris tstFile [list of additional tst files]"
