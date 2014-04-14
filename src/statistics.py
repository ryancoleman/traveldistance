#just a bunch of handy methods for basic stats

import math
import random

def computeMean(totalList):
  if len(totalList) > 0:
    return sum(totalList)/float(len(totalList))
  else:
    return 0.

def computeStdDev(totalList, mean):
  if len(totalList) > 0:
    sumSq = 0.
    for datum in totalList:
      sumSq += (datum - mean)**2.
    return (sumSq/float(len(totalList)))**0.5
  else:
    return 0.

def computeHistogram(totalList, interval, maxData=False):
  if not maxData:
    maxData = max(totalList)
  histogram = [0 for x in xrange(1 + int(maxData/interval))]
  for data in totalList:
    place = int(math.floor(data/interval))
    try:
      histogram[place] += 1
    except IndexError:
      pass  # if user wants max set to be less than real max
  return histogram, maxData

def computeNormalizedHistogram(totalList, interval, maxData=False):
  '''histogram where counts sum to 1'''
  if not maxData:
    maxData = max(totalList)
  histogram = [0. for x in xrange(2 + int(maxData/interval))]  # 2=beg&end bin
  for data in totalList:
    place = int(math.floor(data/interval))
    try:
      histogram[place] += 1./len(totalList)
    except IndexError:
      pass  # if user wants max set to be less than real max
  return histogram, maxData

def computeCumulativeHistogram(histogram):
  total = sum(histogram)
  cumulative = []
  runningTotal = 0.
  for value in histogram:
    runningTotal += value
    if 0. < total:
      cumulative.append(runningTotal/total)
    else:
      cumulative.append(0.)
  return cumulative

def pvalueDiffMeansLazy(list1, list2, numTests):
  '''calculate the difference in the means of the list for you, pass to
  pvalueDiffMeans'''
  diffMean = computeMean(list1) - computeMean(list2)
  pVal1, pVal2 = pvalueDiffMeans(list1, list2, diffMean, numTests)
  return diffMean, pVal1, pVal2

def pvalueDiffMeans(list1, list2, diffMean, numTests):
  '''does a pvalue test without actually making new lists over and over.
  does both one-sided tests, returns the above then below p-value.'''
  lens = (len(list1), len(list2))
  combinedList = list1 + list2
  testsAboveBelow = [1, 1]
  for test in xrange(numTests):
    sums = [0, 0]
    counts = [0, 0]
    for position in combinedList:
      if counts[0] == lens[0]:  # means 0th is full
        sums[1] += position
      elif counts[1] == lens[1]:  # means 1st is full
        sums[0] += position
      else:  # actually make choice
        chosen = random.randint(0, 1)
        sums[chosen] += position
        counts[chosen] += 1
    averages = [sums[0]/float(lens[0]), sums[1]/float(lens[1])]
    diff = averages[0] - averages[1]
    #print sums, averages, diff, diffMean
    if diff >= diffMean:
      testsAboveBelow[0] += 1
    if diff <= diffMean:
      testsAboveBelow[1] += 1
  return testsAboveBelow[0]/float(1.+numTests), \
      testsAboveBelow[1]/float(1.+numTests)

def pvalueDiffMeansUnpairedLazy(list1, list2, numTests):
  '''calculate the difference in the means of the list for you, pass to
  pvalueDiffMeans'''
  diffMean = computeMean(list1) - computeMean(list2)
  pVal1, pVal2 = pvalueDiffMeansUnpaired(list1, list2, diffMean, numTests)
  return diffMean, pVal1, pVal2

def pvalueDiffMeansUnpaired(list1, list2, diffMean, numTests):
  '''does a pvalue test without actually making new lists over and over.
  does both one-sided tests, returns the above then below p-value.
  for unpaired data.'''
  lens = (len(list1), len(list2))
  combinedList = list1 + list2
  sumList = sum(combinedList)
  testsAboveBelow = [1, 1]
  for test in xrange(numTests):
    sums = [0, 0]
    counts = [0, 0]
    listOne = random.sample(combinedList, len(list1))
    sums[0] = sum(listOne)
    sums[1] = sumList - sums[0]
    means = [sums[0]/float(lens[0]), sums[1]/float(lens[1])]
    diff = means[0] - means[1]
    #print sums, averages, diff, diffMean
    if diff >= diffMean:
      testsAboveBelow[0] += 1
    if diff <= diffMean:
      testsAboveBelow[1] += 1
  return testsAboveBelow[0]/float(1.+numTests), \
      testsAboveBelow[1]/float(1.+numTests)

def calculateIntervals(rowExtrema, interval):
  '''returns a list of min, gap1, gap2, ... gapn-2, gapn-1'''
  retList = [rowExtrema[0]]
  tempVal = retList[-1]
  while tempVal < rowExtrema[1]:
    tempVal += interval
    retList.append(tempVal)
  return retList

def compute2dHistogramRows(row1data, row2data, interval1, interval2):
  '''computes intersection of 2 rows as 2d histogram counts'''
  row1ext = [min(row1data), max(row1data)]  # extrema
  row2ext = [min(row2data), max(row2data)]
  row1span = row1ext[1] - row1ext[0]
  row2span = row2ext[1] - row2ext[0]
  histogramLine = [0 for x in xrange(1+int(row1span/interval1))]
  histogram = [histogramLine[:] for x in xrange(1+int(row2span/interval2))]
  for index in xrange(len(row1data)):   # initialized, now actually increment
    pair = row1data[index], row2data[index]
    place1 = int(math.floor((pair[0]-row1ext[0])/interval1))
    place2 = int(math.floor((pair[1]-row2ext[0])/interval2))
    histogram[place2][place1] += 1
  return histogram, calculateIntervals(row2ext, interval2), \
      calculateIntervals(row1ext, interval1)

def compute2dHistogramData(data, row1, row2, interval1, interval2):
  '''the data is a list of lists, this counts the intersections of the 2 rows'''
  row1data, row2data = [], []  # extract rows from data
  for line in data:
    row1data.append(line[row1])
    row2data.append(line[row2])
  return compute2dHistogramRows(row1data, row2data, interval1, interval2)

def compute2dHistogramRowsAlphabetic(row1data, row2data, interval1):
  '''row2 is alphabetic, not numeric data'''
  row1ext = [min(row1data), max(row1data)]  # extrema
  row1span = row1ext[1] - row1ext[0]
  row2set = set(row2data)
  row2heads = list(row2set)
  row2heads.sort()
  histogramLine = [0 for x in xrange(1+int(row1span/interval1))]
  histogram = [histogramLine[:] for x in xrange(len(row2heads))]
  for index in xrange(len(row1data)):   # initialized, now actually increment
    pair = row1data[index], row2data[index]
    place1 = int(math.floor((pair[0]-row1ext[0])/interval1))
    place2 = int(row2heads.index(row2data[index]))
    histogram[place2][place1] += 1
  return histogram, row2heads, calculateIntervals(row1ext, interval1)

def compute2dHistogramAlphabetic(data, row1, row2, interval1):
  '''row2 is alphabetic, not numeric data'''
  row1data, row2data = [], []  # extract rows from data
  for line in data:
    row1data.append(line[row1])
    row2data.append(line[row2])
  return compute2dHistogramRowsAlphabetic(row1data, row2data, interval1)

def write2dHistogram(outFile, histogram, intervalY, intervalX):
  '''takes an open file, writes the histogram to it'''
  outFile.write("-\t")
  for data in intervalX:
    outFile.write(str(data) + "\t")
  outFile.write("\n")
  for index, row in enumerate(histogram):
    outFile.write(str(intervalY[index]) + "\t")
    for data in row:
      outFile.write(str(data) + "\t")
    outFile.write("\n")

def permuteLists(pdbs):
  '''input a length 2 list where both sublists have same length.
  permute the lists by randomly switching pairs in same position.
  return new length 2 list of permuted sublists'''
  switch = []  # a 1, 0 list that defines which pairs are switched
  for position in range(len(pdbs[0])):
    switch.append(random.randint(0, 1))
  newPdbs = ([], [])
  for position in range(len(pdbs[0])):
    if switch[position] == 0:  # don't switch
      newPdbs[0].append(pdbs[0][position])
      newPdbs[1].append(pdbs[1][position])
    else:  # switch
      newPdbs[1].append(pdbs[0][position])
      newPdbs[0].append(pdbs[1][position])
  return newPdbs

def listOrderCorrectness(inList):
  '''sorts the inlist. computes the difference at each position. list should be
  integers. square the differences, divide by the length, return that'''
  sortList = inList[:]  # copy
  sortList.sort()
  squaredSum = 0.
  for index in xrange(len(inList)):
    diff = abs(sortList[index] - inList[index])
    squaredSum += diff**2.
  return (squaredSum / float(len(inList)))

def getAllSubLists(inList, minLength=1):
  '''calls getsublists with all possible values of length greater than X'''
  masterList = []
  lastList = None
  for count in xrange(len(inList)):
    if count >= minLength:
      thisList = getSubLists(inList, count + 1, smallerLists=lastList)
      masterList.extend(thisList)
      lastList = thisList
  return masterList

def getSubLists(inList, retLength, smallerLists=None):
  '''enumerates all possible ways of choosing retLength items from inList
  smallerLists can be the list of possibilities of length 1 less than retlen'''
  if retLength == len(inList):
    return [inList]
  elif retLength == len(inList) - 1:  # all ways of removing 1 is easy
    retLists = []
    for index in xrange(len(inList)):
      newList = inList[:]
      del newList[index]
      retLists.append(newList)
    return retLists
  elif retLength == 1:
    retLists = []
    for index in xrange(len(inList)):
      retLists.append([inList[index]])
    return retLists
  else:  # recurse
    retLists = []
    if smallerLists is None:
      smallLists = getSubLists(inList, retLength-1)  # recursive step
    else:
      smallLists = smallerLists
    for aList in smallLists:
      lastIndex = inList.index(aList[-1])
      for lastItem in inList[lastIndex + 1:]:
        newList = aList[:] + [lastItem]
        retLists.append(newList)
    return retLists

def cohenEffectSize(list1, list2):
  '''uses list sizes, list means and list standard deviations to compute
  effect size. definitions <0.2 trivial; 0.2-0.5 small; 0.5-0.8 medium;
  >0.8 large. with much credit to Paul Hawkins.'''
  meanX = computeMean(list1)
  meanY = computeMean(list2)
  sizeX = len(list1)
  sizeY = len(list2)
  stdDevX = computeStdDev(list1, meanX)
  stdDevY = computeStdDev(list2, meanY)
  cohenD = (meanX - meanY) * math.sqrt(sizeX + sizeY - 2) / \
      math.sqrt(
          ((sizeX - 1) * stdDevX ** 2.0 + (sizeY - 1) * stdDevY ** 2.0)/2.)
  return abs(cohenD)

def statisticalPower(list1, list2):
  '''assumes alpha=0.05, beta=0.2 and calculates cohend &
  then uses 2*stddev^2 (.84-1.96/2)^2/cohenD^2 to calculate
  necessary sample size. easiest formula i could find, happy to get a better
  one to use!'''
  meanX = computeMean(list1)
  meanY = computeMean(list2)
  stdDevX = computeStdDev(list1, meanX)
  stdDevY = computeStdDev(list2, meanY)
  meanStdDev = computeMean([stdDevX, stdDevY])
  samplePower = ((2*meanStdDev**2.0)*((.84+1.96/2.)**2.0))/((meanX-meanY)**2.0)
  return math.ceil(samplePower)
