#!/usr/bin/env python
#ryan g. coleman ryangc ATSYMBOL mail.med.upenn.edu

#analyzes the b-factor column of pdb file... for travel-in analysis

import analyzePdbB
import string
import sys
import pdb

def analyzePdbBCompareFiles(firstFilename):
  pdbs = {}  # each list is a pdb, which has a dictionary of residues
  filenameList = []
  if firstFilename:
    listFile = open(firstFilename, 'r')
    try:
      for line in listFile:
        if len(line) == 6:
          filenameList.append(
              line[0:4] + "." + line[4] + ".nocav.tst.mesh.travelin.pdb")
        elif len(line) == 5:  # analyze non-domains
          filenameList.append(line[0:4] + ".nocav.tst.mesh.travelin.pdb")
        elif len(line) > 8:
          filenameList.append(line[:-1])
        else:  # yeah just allow anything
          filenameList.append(line[:-1])
    except StopIteration:
      pass
    listFile.close()
  if filenameList:
    for filename in filenameList:
      pdbs[filename] = []
      pdbD = pdb.pdbData(filename)
      tempRes = {}
      for index, resName in enumerate(pdbD.resNames):
        if pdbD.radii[index] > 0.:
          if resName not in tempRes:
            tempRes[resName] = {}
          if string.strip(pdbD.atoms[index]) not in tempRes[resName]:
            tempRes[resName][string.strip(pdbD.atoms[index])] = []
          tempRes[resName][string.strip(pdbD.atoms[index])].append(
              pdbD.factors[index][1])
      pdbs[filename] = tempRes
      #residues now contains all the b-factor (travelin) data
  #need to change after this line to call new methods
  listFilename = firstFilename
  analyzePdbB.makeMeanPerProteinReport(
      pdbs, outName=listFilename + ".meansbypdb.txt")

#this is main
if -1 != string.find(sys.argv[0], "analyzePdbB_meanPdbs.py"):
  if len(sys.argv) > 1:
    for filename1 in sys.argv[1:]:
      analyzePdbBCompareFiles(filename1)
