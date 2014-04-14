#!/usr/bin/env python
#ryan g coleman ryangc@mail.med.upenn.edu ryan.g.coleman@gmail.com

import string
import sys

def readRadiiTxtFile(filename):
  '''reads a file that is a pair of floats separated by a comma, returns list'''
  distRadPairs = []
  try:
    openfile = open(filename, 'r')
    for line in openfile:
      tokens = string.split(string.strip(line), ',')
      distRadPairs.append((float(tokens[0]), float(tokens[1])))
  except StopIteration:
    pass  # it is okay
  return distRadPairs

def combinePairs(allDistRadPairs, names):
  '''combines dist, rad pairs into one master list'''
  allFirst = set()
  for distRadPairs in allDistRadPairs:
    for dist, rad in distRadPairs:
      allFirst.add(dist)
  distances = list(allFirst)
  distances.sort()  # now have all the unique distances in sorted order
  radii = []
  for index, distRadPairs in enumerate(allDistRadPairs):
    name = names[index]
    theseRadii = []
    for dist, rad in distRadPairs:
      found = distances.index(dist)
      while found > len(theseRadii):  # must interpolate to produce new data
        lastRadius = theseRadii[-1]
        lastDistIndex = len(theseRadii)-1
        lastDist = distances[lastDistIndex]
        thisDist = distances[lastDistIndex+1]  # next one
        difference = (thisDist-lastDist)/(dist-lastDist)  # between 0 and 1
        curRad = lastRadius + ((rad-lastRadius)*difference)
        theseRadii.append(curRad)
      theseRadii.append(rad)  # add the current one
    while len(theseRadii) < len(distances):
      theseRadii.append("-")
    theseRadii.insert(0, name)  # easier to do now due to indexing
    radii.append(theseRadii)
  distances.insert(0, 'distances')
  radii.insert(0, distances)
  return radii

def outputNewData(newTxtData):
  #made lists wrong, so messy to output
  columns = len(newTxtData)
  rows = len(newTxtData[0])
  for row in range(rows):
    for col in range(columns):
      print str(newTxtData[col][row]) + ",",
    print " "

#main function
if -1 != string.find(sys.argv[0], "combineRadiiFiles"):
  allDistRadPairs = []
  for filename in sys.argv[1:]:
    distRadPairs = readRadiiTxtFile(filename)
    allDistRadPairs.append(distRadPairs)
  newTxtData = combinePairs(allDistRadPairs, sys.argv[1:])
  outputNewData(newTxtData)
