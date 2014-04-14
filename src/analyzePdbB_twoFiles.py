#!/usr/bin/env python
#ryan g. coleman ryangc@mail.med.upenn.edu

#analyzes the b-factor column of pdb file... for travel-in analysis

import analyzePdbB
import string
import sys
import pdb

def analyzePdbBCompareFiles(firstFilename, compareFilename):
  residues = ({}, {})  # a dict of dicts where the sub-dicts keys are atom names
  pdbs = ([], [])  # each list is a pdb, which has a dictionary of residues
  for lindex, listFilename in enumerate((firstFilename, compareFilename)):
    filenameList = []
    if listFilename:
      listFile = open(listFilename, 'r')
      try:
        for line in listFile:
          if len(line) == 6:
            filenameList.append(
                line[0:4] + "." + line[4] + ".nocav.tst.mesh.travelin.pdb")
          if len(line) == 5:  # analyze non-domains
            filenameList.append(line[0:4] + ".nocav.tst.mesh.travelin.pdb")
      except StopIteration:
        pass
      listFile.close()
    if filenameList:
      for filename in filenameList:
        pdbD = pdb.pdbData(filename)
        tempRes = {}
        for index, resName in enumerate(pdbD.resNames):
          if pdbD.radii[index] > 0.:
            if resName not in tempRes:
              tempRes[resName] = {}
            if resName not in residues[lindex]:
              residues[lindex][resName] = {}  # init sub-dict
            if string.strip(pdbD.atoms[index]) not in residues[lindex][resName]:
              residues[lindex][resName][string.strip(pdbD.atoms[index])] = []
            if string.strip(pdbD.atoms[index]) not in tempRes[resName]:
              tempRes[resName][string.strip(pdbD.atoms[index])] = []
            #if pdbD.factors[index][1] > 30:
            #  #print pdbD.factors[index][1], filename #debugging
            #print pdbD.atoms[index],
            residues[lindex][resName][string.strip(pdbD.atoms[index])].append(
                pdbD.factors[index][1])
            tempRes[resName][string.strip(pdbD.atoms[index])].append(
                pdbD.factors[index][1])
        pdbs[lindex].append(tempRes)
      #residues now contains all the b-factor (travelin) data
  #need to change after this line to call new methods
  listFilename = firstFilename + "." + compareFilename
  analyzePdbB.makeCompareResidueReport(
      residues, listFilename + ".residue.bfactor")
  #do alternate pvalue tests and by file comparisons with pdbs lists
  corrAll, corrBeta = analyzePdbB.makeCompareResidueReportAlternate(
      pdbs, listFilename + ".residue.alt.bfactor")
  #now do a p-value comparison that is corrected for the overall depth diffs
  analyzePdbB.makeCompareResidueReportAlternate(
      pdbs, listFilename + ".residue.alt.corr.bfactor",
      correctionAll=corrAll, correctionBeta=corrBeta)
  analyzePdbB.makeMeanPerProteinReport(
      pdbs, outName=listFilename + ".meansbypdb.txt")

#this is main
if -1 != string.find(sys.argv[0], "analyzePdbB_twoFiles.py"):
  if len(sys.argv) > 2:
    filename1 = sys.argv[1]
    filename2 = sys.argv[2]
    analyzePdbBCompareFiles(filename1, filename2)  # passing in a filename now
