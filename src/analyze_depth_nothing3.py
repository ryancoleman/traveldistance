#!/usr/bin/env python2.5
#ryan g. coleman ryangc@mail.med.upenn.edu

#statistical analysis of depth script

import math,  string, sets, sys #standards for old functions
import glob #easy way to grab all ligand files
import random #permutations, random
import tstdata

def analyzeDepth(prefix): 
  pointXYZD = []
  if len(prefix) == 4:
    tstFiles = glob.glob("*" + prefix + "*nocav*" + "tst*")
  elif len(prefix) == 5:
    tstFiles = glob.glob("*" + prefix[:4] + "*" + prefix[4] + "*" + "tst*")
  elif len(prefix) > 5: #just do anything
    tstFiles = glob.glob(prefix)
  for tstFile in tstFiles:
    tstDataP = tstdata.tstData(tstFile,necessaryKeys=['DEPTH_TRAVEL_DIST'])
    pointTravelDist = tstDataP.dict['DEPTH_TRAVEL_DIST']
    for pointDepth in pointTravelDist:
      distT = pointDepth[1]
      pointXYZD.append(distT)  
    #no matter what, return the constructed array    
  return pointXYZD

def analyzeAll(listPrefixes):
  allData = []
  for prefix in listPrefixes:
    allData.append(analyzeDepth(prefix))
  #print allData
  #now we have the data for all the prefixes (structures),
  #first guess, look at false-bar and true-bar for each protein
  # and overall...
  #pickle the data.
  #pickleData = [listPrefixes, allData]
  #pickle.dump(pickleData, open("pickled.analysis.data",'w'))  
  allhists = []
  prefMean = {}
  for index,surfPoints in enumerate(allData):
    histogramBins = [0 for x in range(2000)] #more than that, in trouble
    runningTotal, numberPoints = 0.,0
    for surfPoint in surfPoints:
      histogramBins[int(surfPoint*10)] += 1
      runningTotal += surfPoint
      numberPoints += 1
    prefMean[listPrefixes[index]]= float(runningTotal)/float(numberPoints)
    allhists.append(histogramBins)
  for index,histogramBins in enumerate(allhists):
    print listPrefixes[index], prefMean[listPrefixes[index]],
    for histBin in histogramBins:
      print histBin, 
    print ""

#this is main
if len(sys.argv) > 1:  
  analyzeAll(sys.argv[1:])
