#!/usr/bin/env python

import tm3
import string
import sys
import statistics
import pdb

columnNames = [
    "Volume", "Apolar Surface Area", "Positive Surface Area",
    "Negative Surface Area", "Fraction Apolar Surface Area",
    "Fraction Positive Surface Area", "Fraction Negative Surface Area",
    "maxTD", "Residue Name List", "Residue Number List",
    "Surface Area", "threshold", "meanTD", "height", "mean height",
    "mouths", "Area of Biggest Mouth", "Diameter of Biggest Mouth",
    "mean Curvature", "mean absolute Curvature"]
params = [{
    "volMin": 300., "volMax": 800., "apolarMin": 150.,
    "polarMin": 100., "apolarFractionMax": 0.8, "polarFractionMax": 0.5}, {
    "volMin": 200., "volMax": 1000., "apolarMin": 100.,
    "polarMin": 50., "apolarFractionMax": 0.9, "polarFractionMax": 0.7}, {
    "volMin": 200., "volMax": 1000., "apolarMin": 50.,
    "polarMin": 50., "apolarFractionMax": 0.9, "polarFractionMax": 0.9}]

def pickGoodPockets(tmData, goodNodeList, columnNums, params):
  '''picks the good pockets according to params'''
  for aNode in tmData.bottomUpTraversal():
    thisVolume = aNode.attributes[columnNums[0]]
    thisApolarSA = aNode.attributes[columnNums[1]]
    thisPolarSA = aNode.attributes[columnNums[2]] + \
        aNode.attributes[columnNums[3]]
    thisApolarFraction = aNode.attributes[columnNums[4]]
    thisPolarFraction = aNode.attributes[columnNums[5]] + \
        aNode.attributes[columnNums[6]]
    #all logic is embedded here
    if thisVolume >= params["volMin"] and thisVolume <= params["volMax"] and \
        thisApolarSA >= params["apolarMin"] and \
        thisPolarSA >= params["polarMin"] and \
        thisApolarFraction <= params["apolarFractionMax"] and \
        thisPolarFraction <= params["polarFractionMax"]:  # meets all criteria
      if tmData.areNotOffspring(aNode, goodNodeList):
        goodNodeList.append(aNode)
      else:
        #had a deeper child already meet criteria
        pass

def pickPocketsPdb(pdbData, tmData, userParams):
  '''finds pockets that meet certain criteria'''
  if userParams is not None:
    params.append(userParams)
  columnNums = tmData.titlesToColumns(columnNames)
  for paramCount, paramSet in enumerate(params):
    goodNodeList = []
    pickGoodPockets(tmData, goodNodeList, columnNums, paramSet)
    if len(goodNodeList) == 0:
      print "WARNING: no good pockets found with parameter set", paramCount
    else:
      print "parameter set", paramCount, "found", len(goodNodeList), "pockets"
    for goodCount, goodNode in enumerate(goodNodeList):
      resNumStr = goodNode.attributes[columnNums[9]]
      newPdb = pdbData.getListResiduesChains(resNumStr)
      #open the output file
      outFile = open(
          str(paramCount) + "." + str(goodCount) + ".xtal-lig.pdb", 'w')
      outFile.write("REMARK parameter set: " + str(paramCount) + "\n")
      for name, value in params[paramCount].iteritems():
        outFile.write("REMARK paramater " + name + ": " + str(value) + "\n")
      outFile.write("REMARK pocket number: " + str(goodCount) + "\n")
      #write a bunch of remarks about all the attributes
      for colCount in xrange(len(columnNames)):
        outFile.write(
            "REMARK " + str(columnNames[colCount]) + ": " +
            str(goodNode.attributes[columnNums[colCount]]) + "\n")
      newPdb.outputLines(outFile)  # write the actual PDB file
      outFile.close()

if -1 != string.find(sys.argv[0], "tm3pickerPdb.py"):
  pdbData = pdb.pdbData(sys.argv[1])
  tmData = tm3.tmTreeFromFile(sys.argv[2])
  paramsNew = None
  if len(sys.argv) > 8:  # user wants to declare some an extra set of parameters
    paramsNew = {}
    paramNames = [
        "volMin", "volMax", "apolarMin",
        "polarMin", "apolarFractionMax", "polarFractionMax"]
    for count in xrange(3, 9):
      paramsNew[paramNames[count - 3]] = float(sys.argv[count])
    print "user declare pocket picker parameters:", paramsNew
  pickPocketsPdb(pdbData, tmData, paramsNew)
