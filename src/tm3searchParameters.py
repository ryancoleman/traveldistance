#!/usr/bin/env python2.5

#reads in tm3 files, ignores tree structure, just looks for nodes/pockets
#that are similar

import tm3, string, sys

if -1 != string.find(sys.argv[0], "tm3searchParameters.py"):
  tmDataList = []
  try:
    alpha = float(sys.argv[1]) #the second parameter should be a float, 
                      #which is multiplied by the redudancy scores
    minVolume = float(sys.argv[2])
    maxVolume = float(sys.argv[3])
    startConn = float(sys.argv[4])
    endConn = float(sys.argv[5])
    stepConn = float(sys.argv[6])
    useNormalMeans = bool(sys.argv[7]) 
    calculateNormalMeans = not useNormalMeans #actually want to invert
  except (IndexError, ValueError):
    print "error with parameters, all are required..."
    print "tm3searchParameters.py alpha minVolume maxVolume startConnections"+ \
         " endConnections stepConnections useNormalMeansStds(boolean)" + \
         " [tm3 file list]"
    print "defaults are tm3searchParameters.py 1 25 2000 10000 20000 1000 True"
    exit(1)
  filenames = sys.argv[8:]
  filenames.sort() #so coloring is consistent
  for filename in filenames:
    tmDataList.append(tm3.tmTreeFromFile(filename))
  tm3.findSimilarNodes(tmDataList, ["Surface Area","Volume","height", \
          "mean Curvature", "mouths",\
          "longest dimension", "middle dimension", "short dimension", \
          "Area of Biggest Mouth","Diameter of Biggest Mouth","mean height"], \
          startConn, endConn, stepConn, \
          resCols=["Atom Name List"], sizeColName = "Volume", \
          sizeMin=minVolume,sizeMax=maxVolume, \
          calcColMeansStd=calculateNormalMeans, alpha=alpha)
