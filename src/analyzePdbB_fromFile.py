#!/usr/bin/env python2.5
#ryan g. coleman ryangc@mail.med.upenn.edu

#analyzes the b-factor column of pdb file... for travel-in analysis

import sys, string, pdb
import analyzePdbB 

def analyzePdbBfromFile(listFilename, buriedThresh=2.):
  residues = {} # a dict of dicts where the sub-dicts are keyed on atom name
  buriedResidues, surfaceResidues = {},{}
  filenameList = []
  if listFilename:
    listFile = open(listFilename, 'r')
    try:
      for line in listFile:
        if len(line) == 6:
          #filenameList.append(line[0:4] + "." + line[4] + ".nocav.tst.travelin.pdb")
          filenameList.append(line[0:4] + "." + line[4] + ".nocav.tst.mesh.travelin.pdb")
    except StopIteration:
      pass
    listFile.close()
  if filenameList:
    for filename in filenameList:
      pdbD = pdb.pdbData(filename)
      for index,resName in enumerate(pdbD.resNames):
        if pdbD.radii[index] > 0.:
          atomName = string.strip(pdbD.atoms[index])
          bfactors = pdbD.getFactorsByResidueChain(pdbD.resNums[index], \
                                           pdbD.chains[index])
          if buriedThresh > min(bfactors):
            if resName not in surfaceResidues:
              surfaceResidues[resName] = {} #init sub-dict
            if atomName not in surfaceResidues[resName]:
              surfaceResidues[resName][atomName] = []
            surfaceResidues[resName][atomName].append(pdbD.factors[index][1])
          else:
            if resName not in buriedResidues:
              buriedResidues[resName] = {} #init sub-dict
            if atomName not in buriedResidues[resName]:
              buriedResidues[resName][atomName] = []
            buriedResidues[resName][atomName].append(pdbD.factors[index][1])
          #add it to this either way
          if resName not in residues:
            residues[resName] = {} #init sub-dict
          if atomName not in residues[resName]:
            residues[resName][atomName] = []
          residues[resName][atomName].append(pdbD.factors[index][1])

    #residues now contains all the b-factor (travelin) data
  analyzePdbB.makeResidueReport(residues, listFilename + ".residue.bfactor", \
                    maxY=6, maxYBeta=6) #hardcoded for now...
  analyzePdbB.makeAtomReport(residues, listFilename + ".atom.bfactor")
  analyzePdbB.makeAtomReport(surfaceResidues, listFilename + ".atom.bfactor.surface")
  analyzePdbB.makeAtomReport(buriedResidues, listFilename + ".atom.bfactor.buried")
  analyzePdbB.makeHistogramReport(residues, listFilename+ ".histogram.bfactor")

#this is main
if -1 != string.find(sys.argv[0], "analyzePdbB_fromFile.py"): 
  if len(sys.argv) > 1:
    for filename in sys.argv[1:]:
      analyzePdbBfromFile(filename) #passing in a filename now

