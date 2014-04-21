#!/usr/bin/env python

import tm3
import string
import sys
import statistics

def summarizeOneFeature(tmDataList, columnName, intervals=50, outName="a.txt"):
  '''takes that column, makes a histogram for each structure'''
  outFile = open(outName, 'w')
  columnNum = tmDataList[0].titleToColumn(columnName)
  treeData = {}
  overallMax = 0.
  for tm3tree in tmDataList:
    data = tm3tree.getListColumn(columnNum)
    overallMax = max(overallMax, max(data))
    treeData[tm3tree] = data
  if intervals == "max":
    intervals = overallMax  # 1 per
  interval = overallMax/intervals  # number of intervals desired
  #print a header
  outFile.write("name\tcount\tmean\tstddev\t")
  currentOut = 0.
  while currentOut < overallMax:
    outFile.write(str(currentOut) + "\t")
    currentOut += interval
  outFile.write("\n")
  for tm3tree in tmDataList:
    tm3data = treeData[tm3tree]
    avgData = statistics.computeMean(tm3data)
    stddevData = statistics.computeStdDev(tm3data, avgData)
    histo, outMax = statistics.computeHistogram(tm3data, interval, overallMax)
    outFile.write(tm3tree.inputFileName + "\t")
    outFile.write(str(len(tm3data)) + "\t")
    outFile.write(str(avgData) + "\t")
    outFile.write(str(stddevData) + "\t")
    for histoCount in histo:
      outFile.write(str(histoCount) + "\t")
    outFile.write("\n")
  outFile.close()

if -1 != string.find(sys.argv[0], "tm3summarize.py"):
  tmDataList = []
  for filename in sys.argv[1:]:
    tmData = tm3.tmTreeFromFile(filename)
    tmDataList.append(tmData)
  summarizeOneFeature(tmDataList, "Volume", outName="volume.txt")
  summarizeOneFeature(tmDataList, "Surface Area", outName="surfarea.txt")
  summarizeOneFeature(
      tmDataList, "mouths", outName="mouths.txt", intervals="max")
