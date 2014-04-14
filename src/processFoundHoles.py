#!/usr/bin/env python
#ryan g coleman, ryangc@mail.med.upenn.edu
#copyright 2006-7 ryan g coleman, kim sharp crystal.med.upenn.edu

#processes output from findholes routine... in *findholes.log files
#file format is whitespace delimited

import sys
import string
import statistics
import comparePaths
import paths  # get path stats methods
gnuplotAvailable, plotter = True, False
try:
  import Gnuplot
  plotter = Gnuplot.Gnuplot(debug=0)
except ImportError:
  gnuplotAvailable = False

colNames = ["stepLength", "pathLength", "pathMinRadius", "pathMaxInsideRadius"]
colNamesBonus = [
    "stepLength", "pathLength", "pathMinRadius", "pathMaxInsideRadius",
    "pRMSD", "coverage", "span", "wRMSD", "less1", "lessRad", "radiicomp"]

def getDataFromFile(fileName):
  '''processes a single file into a list of lists for later interpretation'''
  file = open(fileName, 'r')
  returnData = []
  try:
    for line in file:
      data = string.split(line)
      if data[0] == 'number':  # opening line
        pass #do nothing
      else: #process it
        #convert to floats
        floatData = []
        for datum in data:
          try:
            floatData.append(float(datum))
          except ValueError:  # plugs column
            floatData.append(str(datum))
        while len(floatData) < 15:  # pad
          floatData.append(0)
        returnData.append(floatData)
  except StopIteration:
    pass #EOF
  file.close()
  return returnData

def getOneColumn(data, column):
  '''returns the one column as a list, useful for histograms'''
  returnList = []
  for line in data:
    returnList.append(line[column])
  return returnList

def getFromColumnIfMin(data, column, ifCol, valueCol):
  '''returns the one column as a list,
     if each entry in ifCol is less than or equal to valueCol'''
  returnList = []
  for line in data:
    if line[ifCol] <= valueCol:
      returnList.append(line[column])
  return returnList

def getFromColumnIfMax(data, column, ifCol, valueCol):
  '''
  returns the one column as a list,
  if each entry in ifCol is less than or equal to valueCol
  '''
  returnList = []
  for line in data:
    if line[ifCol] >= valueCol:
      returnList.append(line[column])
  return returnList

def getRankColumn(data, column, row):
  '''returns the number of data in that column less than or equal to
  the row index'''
  #first get all the column data
  columnData = []
  for line in data:
    columnData.append(line[column])
  columnData.sort()
  piece = data[row][column]
  for index,columnDatum in enumerate(columnData):
    if piece <= columnDatum:
      return index + 1

def getMinColumn(data, column):
  '''returns the row entry and index of the minimum'''
  bestIndex, bestValue = 0, data[0][column]
  for rowIndex,line in enumerate(data):
    if line[column] < bestValue:
      bestIndex, bestValue = rowIndex,line[column]
  return bestValue, bestIndex

def getMaxColumn(data, column):
  '''returns the row entry and index of the maximum'''
  bestIndex, bestValue = 0, data[0][column]
  for rowIndex,line in enumerate(data):
    if line[column] > bestValue:
      bestIndex, bestValue = rowIndex,line[column]
  return bestValue, bestIndex

def sortByColumnReturnOther(data, column, other):
  '''sorts by the column column, returns the column other'''
  #make tuples
  tuples = []
  for line in data:
    tuples.append((line[column], line[other]))
  tuples.sort()
  #now just get the last one
  returnArray = []
  for oneTuple in tuples:
    returnArray.append(oneTuple[1])
  return returnArray

def processDataOne(data, filename, sortStat, sortMethod):
  if sortMethod == 'min':
    bestVal, bestIndex = getMinColumn(data, sortStat)  # index 8 is the prmsd
  else:
    bestVal, bestIndex = getMaxColumn(data, sortStat)  # index 8 is the prmsd
  bestRankStr = filename + ", "
  sortedEvals = []
  backwardsEvals = []
  for column in range(4,8):
    rank = getRankColumn(data, column, bestIndex)
    revRank = len(data) - rank + 1
    bestRankStr +=  str(rank) + ", " +  str(revRank) + ", "
    sorted = filename + ", "
    backwards = filename + ", "
    sortedOther = sortByColumnReturnOther(data, column, sortStat)
    for additional in sortedOther:
      sorted += str(additional) + ", "
    sortedOther.reverse()
    for additional in sortedOther:
      backwards += str(additional) + ", "
    sortedEvals.append(sorted)
    backwardsEvals.append(backwards)
  bestRankStr += str(bestVal)
  return bestRankStr, sortedEvals, backwardsEvals

def processData(dataList, nameList, listPaths, \
                outputFileName="processed.foundholes."):
  #first do the one big summary output file
  bestLists = []
  compCols = [8,9,10,11,12,13,14]
  comps = ['min', 'max', 'max', 'min', 'max', 'max', 'min']
  compThreshs = [5., .4, .8, 2.5, .5, .8, 2.5]
  for index,colIdx in enumerate(compCols):
    bestList = []
    for data in dataList:
      if comps[index] == 'min':
        bestVal, bestIndex = getMinColumn(data, colIdx)  # index 8 is the prmsd
      else:
        bestVal, bestIndex = getMaxColumn(data, colIdx)  # index 8 is the prmsd
      bestList.append(bestVal)
    bestLists.append(bestList)
    #fileOut = open(outputFileName + colNamesBonus[colIdx-4] + ".best.log", 'w')
    #for index,name in enumerate(nameList):
    #  fileOut.write(name + " " + str(bestList[index]) + "\n")
    #fileOut.close()
  fileOut = open(outputFileName + "overall.best.log", 'w')
  fileOut.write("name pRMSD coverage span wRMSD less1 lessRad radiicomp\n")
  for index,name in enumerate(nameList):
    fileOut.write(name + " ")
    for index2 in range(len(compCols)):
      fileOut.write(str(bestLists[index2][index]) + " ")
    fileOut.write("\n")
  fileOut.close()
  for colIdx, sortColNumber in enumerate(compCols):
    bestRankStrings, sortedEvals, backwardsEvals = [], [], []
    for index, data in enumerate(dataList):
      bestRankString, sortedEval, backEval = processDataOne(
          data, nameList[index], sortColNumber, comps[colIdx])
      bestRankStrings.append(bestRankString)
      sortedEvals.append(sortedEval)
      backwardsEvals.append(backEval)
    fileOut = open(
        outputFileName + colNamesBonus[sortColNumber-4] + ".best.rankings.log",
        'w')
    for bestRankStr in bestRankStrings:
      fileOut.write(bestRankStr + "\n")
    fileOut.close()
    for index,colName in enumerate(colNames):
      fileOut = open(
          outputFileName + "rankings." + colName + "." +
          colNamesBonus[sortColNumber-4] + ".log", 'w')
      for line in sortedEvals:
        fileOut.write(line[index] + "\n")
      fileOut.close()
      fileOut = open(
          outputFileName + "rankings.reverse." + colName + "." +
          colNamesBonus[sortColNumber-4] + ".log", 'w')
      for line in backwardsEvals:
        fileOut.write(line[index] + "\n")
      fileOut.close()
  drilledData = [
      False, False, False, False, [], [], [], [], [], [], [], [], [], []]
  for drilledPath in listPaths:
    drilledData[4].append(float(len(drilledPath)))
    drilledData[5].append(float(paths.pathLength(drilledPath)))
    drilledData[6].append(float(paths.pathMinRadius(drilledPath)))
    drilledData[7].append(float(paths.pathMaxInsideRadius(drilledPath)))
  for column in range(4, 14):
    columnName = colNamesBonus[column - 4]
    for colIdx, sortColNumber in enumerate(compCols):
      #print columnName, colNamesBonus[sortColNumber - 4]
      columnData = []
      selectColumnData = []
      bestColumnData = []
      for data in dataList:
        if comps[colIdx] == 'min':
          bestVal, bestIndex = getMinColumn(data, sortColNumber)  
          selectColumnData.extend(
              getFromColumnIfMin(data, column, sortColNumber,
              compThreshs[colIdx]))
          bestColumnData.extend(
              getFromColumnIfMin(data, column, sortColNumber, bestVal))
        else:
          bestVal, bestIndex = getMaxColumn(data, sortColNumber)
          selectColumnData.extend(getFromColumnIfMax(
              data, column, sortColNumber, compThreshs[colIdx]))
          bestColumnData.extend(
              getFromColumnIfMax(data, column, sortColNumber, bestVal))
        columnData.extend(getOneColumn(data, column))
      if len(drilledData[column]) > 0:
        maxHere = (max(max(columnData),max(drilledData[column])) + 1.)
      else:
        maxHere = max(columnData) + 1.
      interval = maxHere/40.
      histogram = statistics.computeHistogram(columnData,interval,maxData=maxHere)
      #selectHistogram = statistics.computeHistogram(selectColumnData,interval,maxData=histogram[1])
      bestHistogram = statistics.computeHistogram(bestColumnData,interval,maxData=histogram[1])
      if len(drilledData[column]) > 0:
        realHistogram = statistics.computeHistogram(drilledData[column],interval,maxData=histogram[1])
      else:
        realHistogram = False
      #need to scale select part of data to be same height as histogram
      maxHeight = max(histogram[0])
      #maxSelectHeight = max(selectHistogram[0])
      maxBestHeight = max(bestHistogram[0])
      if realHistogram:
        maxRealHeight = max(realHistogram[0])
        realScaledHistogram = [[],realHistogram[1]]
        for histPoint in realHistogram[0]:
          realScaledHistogram[0].append(histPoint*maxHeight/maxRealHeight)
      #selectScaledHistogram = [[], selectHistogram[1]]
      bestScaledHistogram = [[], bestHistogram[1]]
      #for histPoint in selectHistogram[0]:
        #selectScaledHistogram[0].append(histPoint*maxHeight/maxSelectHeight)
      for histPoint in bestHistogram[0]:
        bestScaledHistogram[0].append(histPoint*maxHeight/maxBestHeight)
      #print histogram, len(histogram[0])
      #print selectHistogram, len(selectHistogram[0])
      #make gnuplot version if possible
      if gnuplotAvailable:
        #plotter = Gnuplot.Gnuplot(debug=0)
        xVals = []
        for cutoff in range(2+int(histogram[1]/interval)):
          xVals.append(cutoff*interval)
        graphData = Gnuplot.Data(xVals, histogram[0], title="All")
        #if comps[colIdx] == 'min':
          #graphSelectData = Gnuplot.Data(xVals, selectScaledHistogram[0], title="<" + str(compThreshs[colIdx]) + " " + str(colNamesBonus[sortColNumber-4]) + " Scaled by " + str(maxHeight/maxSelectHeight))
        #else:
          #graphSelectData = Gnuplot.Data(xVals, selectScaledHistogram[0], title=">" + str(compThreshs[colIdx]) + " " + str(colNamesBonus[sortColNumber-4]) + " Scaled by " + str(maxHeight/maxSelectHeight))
        graphBestData = Gnuplot.Data(xVals, bestScaledHistogram[0], title="Best "+str(colNamesBonus[sortColNumber-4])+" Scaled by " + str(maxHeight/maxBestHeight))
        if realHistogram:
          graphRealData = Gnuplot.Data(xVals, realScaledHistogram[0], title="Drilled Holes Scaled by " + str(maxHeight/maxRealHeight))
        graphDataCum = Gnuplot.Data(xVals, statistics.computeCumulativeHistogram(histogram[0]))
        plotter('set terminal png')
        plotter('set output "' + outputFileName + columnName + '.' + colNamesBonus[sortColNumber-4] + '.png"')
        plotter('set data style linespoints')
        plotter('set xrange [' + str(min(xVals)-.01) + ':' + str(max(xVals)+.01) + ']')
        plotter.xlabel(columnName)
        plotter.ylabel('Count')
        plotter('set multiplot')
        plotter('set key right top')
        if realHistogram:
          #plotter.plot(graphData, graphSelectData, graphBestData, graphRealData)
          plotter.plot(graphData, graphBestData, graphRealData)
        else:
          #plotter.plot(graphData, graphSelectData, graphBestData)
          plotter.plot(graphData, graphBestData)
        plotter('unset multiplot')

#this is main
if len(sys.argv) > 1:
  listData = []
  drilledFileNames,listPaths = [], []
  for fileName in sys.argv[1:]:
    drilledName = string.replace(fileName, ".nocav.tst.findholes.log", ".py")
    drilledFileNames.append(drilledName)
    data = getDataFromFile(fileName)
    listData.append(data)
    drilledPath = comparePaths.readCGOPathWithRadius(drilledName)
    listPaths.append(drilledPath) #radius, x, y, z
  processData(listData, sys.argv[1:], listPaths)
else:
  print "call with list of findholes.log files"
