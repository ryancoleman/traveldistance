#!/usr/bin/env python

import tm3
import string
import sys
import statistics

paramsDefault = {
    "volMin": 300., "volMax": 800., "apolarMin": 150.,
    "polarMin": 100., "apolarFractionMax": 0.8,
    "polarFractionMax": 0.5}
paramsLooser = {
    "volMin": 200., "volMax": 1000., "apolarMin": 100.,
    "polarMin": 50., "apolarFractionMax": 0.9,
    "polarFractionMax": 0.7}

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

def pickPockets(tmData):
  '''finds pockets that meet certain criteria'''
  columnNames = ["Volume", "Apolar Surface Area", "Positive Surface Area",
                 "Negative Surface Area", "Fraction Apolar Surface Area",
                 "Fraction Positive Surface Area",
                 "Fraction Negative Surface Area", "maxTD",
                 "Residue Name List"]
  columnNums = tmData.titlesToColumns(columnNames)
  goodNodeList = []
  pickGoodPockets(tmData, goodNodeList, columnNums, paramsDefault)
  if len(goodNodeList) == 0:
    pickGoodPockets(tmData, goodNodeList, columnNums, paramsLooser)
  for goodNode in goodNodeList:
    print goodNode.attributes[columnNums[8]]
  if len(goodNodeList) == 0:
    print "no good pockets found"

if -1 != string.find(sys.argv[0], "tm3picker.py"):
  tmData = tm3.tmTreeFromFile(sys.argv[1])
  pickPockets(tmData)
