#!/usr/bin/env python

#ryan g. coleman ryan.g.coleman ATSYMBOL gmail.com ryangc ATSYMBOL mail.med.upenn.edu
#kim sharp lab http://crystal.med.upenn.edu

#analyzes and computes p-values for difference of means test, reads data
#from a file with column data and 2 files for the 2 data sets to switch to
#construct p-values.

import string
import sys
import statistics

def readSummaryFile(inputFileName):
  '''reads a space-delimited file into a dictionary keyed on column then row'''
  if inputFileName:
    names = []
    data = []
    inputFile = open(inputFileName, 'r')
    try:
      for line in inputFile:
        tokens = string.split(line)
        if len(data) == 0:
          countDicts = len(tokens) - 1
          data = [{} for count in range(countDicts)]
          names = tokens[1:]
        else:
          name = tokens[0]
          others = tokens[1:]
          for dataIndex, dataDict in enumerate(data):
            dataDict[name] = float(others[dataIndex])
    except StopIteration:
      pass  # EOF
    newData = {}
    for count in range(len(names)):
      newData[names[count]] = data[count]
    return newData

def readList(inputFileName):
  '''assumes each line has a single token, returns list of such'''
  if inputFileName:
    data = []
    inputFile = open(inputFileName, 'r')
    try:
      for line in inputFile:
        tokens = string.split(line)
        if len(tokens) >= 1:
          data.append(tokens[0])  # just add first... assume all there is
    except StopIteration:
      pass  # EOF
    return data

def getMean(summaryDataColumn, listPdb):
  '''for each entry in the second argument, finds the matching value in the
  first argument, adds these up and divides by the length'''
  summary = 0.
  for pdb in listPdb:
    summary += summaryDataColumn[pdb]
  return summary / len(listPdb)

def getDiffMean(summaryDataColumn, listPdbs):
  '''get the diff of the first list - second by using the first argument'''
  return getMean(
      summaryDataColumn, listPdbs[0]) - getMean(summaryDataColumn, listPdbs[1])

def analyzeColumn(
    summaryDataColumn, listPdbs, outputFileName, numTests=1000000):
  '''takes a column of summary data, computes average and permutes, outputs'''
  origDiffMean = getDiffMean(summaryDataColumn, listPdbs)
  pValCounts = [0., 0.]  # above, below
  for test in range(numTests):
    newLists = statistics.permuteLists(listPdbs)
    testMean = getDiffMean(summaryDataColumn, newLists)
    if testMean >= origDiffMean:
      pValCounts[0] += 1.
    if testMean <= origDiffMean:
      pValCounts[1] += 1.
  pVals = [pValCount/float(numTests) for pValCount in pValCounts]
  outputFile = open(outputFileName, 'w')
  outputFile.write("origDiffMean\tmean1\tmean2\tpValAbove\tpValBelow\n")
  outputFile.write(str(origDiffMean) + "\t")
  outputFile.write(str(getMean(summaryDataColumn, listPdbs[0])) + "\t")
  outputFile.write(str(getMean(summaryDataColumn, listPdbs[1])) + "\t")
  for pVal in pVals:
    outputFile.write(str(pVal) + "\t")
  outputFile.write("\n")
  outputFile.close()

if -1 != string.find(sys.argv[0], "analyze_packing.py"):
  try:
    summaryData = readSummaryFile(sys.argv[1])
    listPdbs = [readList(fileName) for fileName in sys.argv[2:4]]
    bothName = str(sys.argv[2]) + "." + str(sys.argv[3]) + "."
    for name in summaryData.keys():
      analyzeColumn(summaryData[name], listPdbs, bothName + name + ".txt")
  except IndexError:
    print "analyze_packing.py summary.txt file1.txt file2.txt"
    print "uses information in summary.txt to compare matching files in 1 and 2"
    sys.exit(1)
