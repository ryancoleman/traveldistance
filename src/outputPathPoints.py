#!/usr/bin/env python2.5
#ryan g coleman ryangc@mail.med.upenn.edu

import pdb
from geometry import dist
import sys, string
import tstdebug
import comparePaths
from mesh import meshFromSpheres

def readPathWritePoints(original, txtfile):
  '''outputs a txtfile of distance vs radius'''
  origPath = comparePaths.readCGOPathWithRadius(original) #0 is radius
  #make 3 radius instead of 0
  newPath = []
  for pathUnit in origPath:
    radius = pathUnit.pop(0)
    pathUnit.append(radius)
    newPath.append(pathUnit)
  meshFromSpheres(newPath, 1.0, txtfile)

if -1 != string.find(sys.argv[0], "outputPathPoints"):
  if len(sys.argv) > 2:
    orig = sys.argv[1]
    txtfile = sys.argv[2]
    readPathWritePoints(orig, txtfile)