#!/usr/bin/env python

#output a pdb file of just the residues between the membrane barrier

import string
import sys
import comparePaths
import os
import pdb
import geometry

def getJustMembranePdb(inputFileName, outputFileName):
  pdbBarriers = pdb.pdbData(inputFileName)
  #get the barriers read in and defined
  barrierAtomList = [[], []]
  for index, resName in enumerate(pdbBarriers.resNames):
    if resName == "DUM":
      if pdbBarriers.atoms[index][0] == "O":
        barrierAtomList[0].append(pdbBarriers.coords[index])
      elif pdbBarriers.atoms[index][0] == "N":
        barrierAtomList[1].append(pdbBarriers.coords[index])
  barrierZ = [barrierAtomList[0][0][2], barrierAtomList[1][0][2]]
  barrierZ.sort()
  barrierSep = geometry.distL2(barrierAtomList[0][0], barrierAtomList[1][0])
  zCoord = barrierZ[1]
  goodResChain = []
  for index, thisResNum in enumerate(pdbBarriers.resNums):
    chain = pdbBarriers.chains[index]
    resChain = str(thisResNum) + str(chain)
    if resChain not in goodResChain:
      #otherwise don't need to check, already in
      zTest = pdbBarriers.coords[index][2]
      if abs(zTest) <= zCoord:
        goodResChain.append(resChain)
  newPdb = pdbBarriers.getListResiduesChains(goodResChain)
  newPdb.write(outputFileName)

if -1 != string.find(sys.argv[0], "getMembranePdb"):
  if len(sys.argv) >= 2:
    for prefix in sys.argv[1:]:
      getJustMembranePdb(prefix, "just_membrane_" + prefix)
  else:
    print "getMembranePdb.py file_planes.pdb [more pdbs]"
