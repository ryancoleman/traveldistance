#!/usr/bin/env python2.5
#ryan g coleman ryangc@mail.med.upenn.edu

import pdb
from geometry import dist
import sys, string
import tstdebug
import comparePaths
from paths import outputRadiiTxt

#tstdebug.debugSetGridSpheres(pathIn, stepSize, root
def readOutputRadiiTxt(original, txtfile):
  '''outputs a txtfile of distance vs radius'''
  origPath = comparePaths.readCGOPathWithRadius(original) #0 is radius
  outputRadiiTxt(origPath, txtfile)  

if -1 != string.find(sys.argv[0], "outputPathTxt"):
  if len(sys.argv) > 2:
    orig = sys.argv[1]
    txtfile = sys.argv[2]
    readOutputRadiiTxt(orig, txtfile)
