#!/usr/bin/env python
#ryan g. coleman ryangc ATSYMBOL mail.med.upenn.edu

#finds nmr model closest to the average model, outputs it as best, doesn't clean

#grab system, string,  regular expression, and operating system modules
import sys
import string
import re
import os
import math
import pdb  # for chain sorting ease

def getClosestToAverage(fileName, skipAlign=True):
  '''takes a pdb file of an nmr structure, superposes all structures, finds
  euclidean average, then finds closest original model to the average, then
  picks it out as the one to use'''
  if not skipAlign:
    import superimpose  # superposition function, calls fortran code
  pdbD = pdb.pdbData(fileName)
  modelNums = pdbD.getModelNumbers()
  firstModel = False
  alignedModels = []
  for modelNum in modelNums:
    newPdb = pdbD.getOneModel(modelNum)
    outputFileName = fileName[:-4] + ".model." + str(modelNum) + ".pdb"
    newPdb.write(outputFileName)
    if not firstModel:
      firstModel = outputFileName  # don't align first
    else:
      if not skipAlign:
        alignFileName = "align." + outputFileName
        superimpose.superposition(firstModel, outputFileName, alignFileName)
        alignedModels.append(alignFileName)
      else:
        alignedModels.append(outputFileName)
  firstModelPdbD = pdb.pdbData(firstModel)
  averageFileName = fileName[:-4] + ".average.pdb"
  averagePdbD = firstModelPdbD.getAverageCoords(alignedModels)
  averagePdbD.write(averageFileName)
  bestModelName = firstModel
  bestRMSD = averagePdbD.calcRMSDfile(firstModel)
  for alignModel in alignedModels:
    otherRMSD = averagePdbD.calcRMSDfile(alignModel)
    if otherRMSD < bestRMSD:
      bestRMSD = otherRMSD
      bestModelName = alignModel
  #print bestRMSD, alignModel
  bestFileName = fileName[:-4] + ".best.pdb"
  modelNumber = bestModelName.replace(".pdb", "").replace(
      "align."+fileName[:-4]+".model.", "").replace(fileName[:-4]+".model.", "")
  pdbD = pdb.pdbData(fileName)
  newPdb = pdbD.getOneModel(int(modelNumber))
  newPdb.write(bestFileName)  # writes all other file info not just coords

if -1 != string.find(sys.argv[0], "getClosestToAverageNMRpdb.py"):
  fileName = sys.argv[1]
  getClosestToAverage(fileName, skipAlign=True)
