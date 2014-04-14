#!/usr/bin/env python

#reads in tm3 files, ignores tree structure, just looks for nodes/pockets
#that are similar

import tm3, string, sys

if -1 != string.find(sys.argv[0], "tm3searchAlphaExamineClusters.py"):
  tmDataList = []
  alpha = float(sys.argv[1]) #the second parameter should be a float,
                      #which is multiplied by the redudancy scores
  filenames = sys.argv[2:]
  filenames.sort() #so coloring is consistent
  for filename in filenames:
    tmDataList.append(tm3.tmTreeFromFile(filename))
  tm3.findSimilarNodes(tmDataList, ["Surface Area","Volume","height", \
          "mean Curvature", "mouths",\
          "longest dimension", "middle dimension", "short dimension", \
          "Area of Biggest Mouth","Diameter of Biggest Mouth","mean height"], \
          100000, 20000, 1000, \
          resCols=["Atom Name List"], sizeColName = "Volume", \
          sizeMin=25,sizeMax=2000, calcColMeansStd=False, alpha=alpha, \
          clusterOutput=True, printThings=False)
