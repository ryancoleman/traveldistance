#!/usr/bin/env python
#ryan g. coleman ryangc@mail.med.upenn.edu

#analyzes the b-factor column of pdb file... for travel-in analysis
#uses gnuplot for graphing... needs gnuplot 4 (apparently)

import pdb
import sys
import string
import math
import random  # for pvalue tests
gnuplotAvailable = True
try:
  import Gnuplot
except ImportError:
  gnuplotAvailable = False
import statistics

#useful lists
aminoAcid3Codes = [
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]

aminoAcidCodes = [
    'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P',
    'S', 'T', 'W', 'Y', 'V']

backboneAtomCodes = ['C', 'CA', 'N', 'O', 'OXT']

cbCode = 'CB'
caCode = 'CA'

carbonBetaCodes = {}

for aminoAcid in aminoAcid3Codes:
  carbonBetaCodes[aminoAcid] = cbCode
carbonBetaCodes["GLY"] = caCode

#hopefully non-arbitrary split for analyzing travel in distance
highCodes = ["CYS", "HIS", "ILE", "LEU", "MET", "PHE", "TRP", "TYR", "VAL"]

lowCodes = [
    "ALA", "ARG", "ASN", "ASP", "GLN", "GLU", "GLY", "LYS", "PRO", "SER", "THR"]

def makeResidueReport(
    residueData, outputFilename="residue.bfactor",
    maxY=False, maxYBeta=False, runGraphs=False):
  #residueNames = residueData.keys()
  residueNames = aminoAcid3Codes
  residueNames.sort()
  fileTemp = open(outputFilename + ".txt", 'w')
  fileTemp.write("ResidueName Mean StdDev Low High Count\n")
  averages, stddevs = {}, {}
  betaAverages, betaStddevs = {}, {}
  for residueName in residueNames:
    #assemble into one big list
    totalList = []
    if residueName in residueData:
      for data in residueData[residueName].values():
        totalList.extend(data)
    average = statistics.computeMean(totalList)
    averages[residueName] = average
    stddev = statistics.computeStdDev(totalList, average)
    stddevs[residueName] = stddev
    betaList = []
    if residueName in residueData:
      data = residueData[residueName]
      betaList.extend(data[carbonBetaCodes[residueName]])
    else:
      data = []
    if len(betaList) > 0:
      betaAvg = statistics.computeMean(betaList)
      betaAverages[residueName] = betaAvg
      betaStddevs[residueName] = statistics.computeStdDev(betaList, betaAvg)
    if len(totalList) > 0:
      fileTemp.write(
          residueName + " " + str(average) + " " + str(stddev) + " " +
          str(min(totalList)) + " " + str(max(totalList)) + " " +
          str(len(totalList)) + "\n")
    else:
      fileTemp.write(
          residueName + " " + str(average) + " " + str(stddev) + " " +
          str(0.) + " " + str(0.) + " " + str(0.) + "\n")
  fileTemp.close()
  if gnuplotAvailable and runGraphs:
    plotter = Gnuplot.Gnuplot(debug=0)
    yLabels = '('
    yData, yError, yMin, yMax = [], [], 10, 0
    yBetaData, yBetaError, yBetaMin, yBetaMax = [], [], 10, 0
    for index, code in enumerate(aminoAcid3Codes):
      yLabels += '"' + str(code) + '" ' + str(index)
      if index != len(aminoAcid3Codes) - 1:
        yLabels += ', '
      if code in averages:
        yData.append(averages[code])
        yError.append(stddevs[code])
        yMin = min(yMin, yData[-1]-yError[-1])
        yMax = max(yMax, yData[-1]+yError[-1])
        yBetaData.append(betaAverages[code])
        yBetaError.append(betaStddevs[code])
        yBetaMin = min(yBetaMin, yBetaData[-1]-yBetaError[-1])
        yBetaMax = max(yBetaMax, yBetaData[-1]+yBetaError[-1])
      else:  # none of that residue
        yData.append(0)
        yError.append(0)
        yBetaData.append(0)
        yBetaError.append(0)
    yLabels += ')'
    graphData = Gnuplot.Data(range(20), yData, yError)
    plotter('set terminal png')
    plotter('set output "' + outputFilename + '.png"')
    plotter('set data style yerrorbars')
    plotter('set boxwidth 0.9 absolute')
    plotter('set xtics ' + yLabels)
    if maxY is False:
      plotter('set yrange [' + str(yMin-0.2) + ':' + str(yMax+0.2) + ']')
    else:
      plotter('set yrange [0:' + str(maxY) + ']')
    plotter('set xrange [-1:20]')
    plotter.xlabel('Residue')
    plotter.ylabel('Mean Travel In Distance')
    plotter.plot(graphData)
    #do another graph with just carbon-betas
    plotter('set output "' + outputFilename + '.beta.png"')
    graphDataBeta = Gnuplot.Data(range(20), yBetaData, yBetaError)
    plotter.ylabel('Mean Travel In Distance of Carbon Beta')
    if maxYBeta is False:
      plotter(
          'set yrange [' + str(yBetaMin-0.2) + ':' + str(yBetaMax+0.2) + ']')
    else:
      plotter('set yrange [0:' + str(maxYBeta) + ']')
    plotter.plot(graphDataBeta)

def makeCompareResidueReportAlternate(
    pdbs, outputFilename="residue.bfactor", numTests=9999,
    correctionAll=0., correctionBeta=0.):
  '''different way to do p-vals, instead of permuting all data, permute the
  pairs of hyp/meso pdb files.'''
  residueNames = aminoAcid3Codes  # for now ignore what is in the files
  fileTemp2 = open(outputFilename + ".pvals.txt", 'w')
  fileTemp2.write("ResidueName DiffMeans MeanA MeanB pValAbove pValBelow\n")
  fileTemp3 = open(outputFilename + ".pvals.beta.txt", 'w')
  fileTemp3.write("ResidueName DiffMeans MeanA MeanB pValAbove pValBelow\n")
  #first find means
  means = [{}, {}]
  betaMeans = [{}, {}]
  overallList = [[], []]
  overallBetaList = [[], []]
  totalMeans, totalBetaMeans = [0., 0.], [0., 0.]
  for code in residueNames:
    betaKsLists = [[], []]
    for lindex in range(2):  # either a or b
      totalList, betaList = [], []
      for pdbResidues in pdbs[lindex]:
        if code in pdbResidues:
          for atomValues in pdbResidues[code].values():
            totalList.extend(atomValues)
      means[lindex][code] = statistics.computeMean(totalList)
      for pdbResidues in pdbs[lindex]:
        if code in pdbResidues:
          betaList.extend(pdbResidues[code][carbonBetaCodes[code]])
      betaKsLists[lindex] = betaList
      betaMeans[lindex][code] = statistics.computeMean(betaList)
      overallList[lindex].extend(totalList)
      overallBetaList[lindex].extend(betaList)
    #use betaKsLists to compute ks stuff
  for lindex in range(2):  # either a or b
    totalMeans[lindex] = statistics.computeMean(overallList[lindex])
    totalBetaMeans[lindex] = statistics.computeMean(overallBetaList[lindex])
  #print means, betaMeans
  pValueCounts = [{}, {}]  # first is above, second is below
  pValueBetaCounts = [{}, {}]
  for code in residueNames+["ALL"]:  # initialize counts, even for overall total
    for aboveBelow in range(2):
      pValueCounts[aboveBelow][code] = 1
      pValueBetaCounts[aboveBelow][code] = 1
  for test in xrange(numTests):
    testMeans = [{}, {}]
    testBetaMeans = [{}, {}]
    overallList = [[], []]
    overallBetaList = [[], []]
    totalTestMeans, totalTestBetaMeans = [0., 0.], [0., 0.]
    newPdbs = statistics.permuteLists(pdbs)
    for code in residueNames:
      for lindex in range(2):  # either a or b
        totalList, betaList = [], []
        for pdbResidues in newPdbs[lindex]:
          if code in pdbResidues:
            for atomValues in pdbResidues[code].values():
              totalList.extend(atomValues)
        testMeans[lindex][code] = statistics.computeMean(totalList)
        for pdbResidues in newPdbs[lindex]:
          if code in pdbResidues:
            betaList.extend(pdbResidues[code][carbonBetaCodes[code]])
        testBetaMeans[lindex][code] = statistics.computeMean(betaList)
        overallList[lindex].extend(totalList)
        overallBetaList[lindex].extend(betaList)
    for lindex in range(2):  # either a or b
      totalTestMeans[lindex] = statistics.computeMean(overallList[lindex])
      totalTestBetaMeans[lindex] = \
          statistics.computeMean(overallBetaList[lindex])
    for code in residueNames:  # calc pval for each residue
      testMeanDiff = testMeans[0][code] - testMeans[1][code]
      origMeanDiff = means[0][code] - means[1][code] - correctionAll
      if origMeanDiff <= testMeanDiff:
        pValueCounts[0][code] += 1
      if origMeanDiff >= testMeanDiff:
        pValueCounts[1][code] += 1
      testMeanDiff = testBetaMeans[0][code] - testBetaMeans[1][code]
      origMeanDiff = betaMeans[0][code] - betaMeans[1][code] - correctionBeta
      if origMeanDiff <= testMeanDiff:
        pValueBetaCounts[0][code] += 1
      if origMeanDiff >= testMeanDiff:
        pValueBetaCounts[1][code] += 1
    code = "ALL"  # fake residue name for overall
    testMeanDiff = totalTestMeans[0] - totalTestMeans[1]
    origMeanDiff = totalMeans[0] - totalMeans[1] - correctionAll
    if origMeanDiff <= testMeanDiff:
      pValueCounts[0][code] += 1
    if origMeanDiff >= testMeanDiff:
      pValueCounts[1][code] += 1
    testMeanDiff = totalTestBetaMeans[0] - totalTestBetaMeans[1] - \
        correctionBeta
    origMeanDiff = totalBetaMeans[0] - totalBetaMeans[1]
    if origMeanDiff <= testMeanDiff:
      pValueBetaCounts[0][code] += 1
    if origMeanDiff >= testMeanDiff:
      pValueBetaCounts[1][code] += 1
  for code in residueNames:  # output time
    fileTemp2.write(code + " " + str(means[0][code]-means[1][code]) + " ")
    fileTemp2.write(str(means[0][code]) + " " + str(means[1][code]) + " ")
    fileTemp2.write(str(pValueCounts[0][code]/float(1+numTests)) + " ")
    fileTemp2.write(str(pValueCounts[1][code]/float(1+numTests)) + " ")
    fileTemp2.write("\n")
    fileTemp3.write(
        code + " " + str(betaMeans[0][code]-betaMeans[1][code]) + " ")
    fileTemp3.write(
        str(betaMeans[0][code]) + " " + str(betaMeans[1][code]) + " ")
    fileTemp3.write(str(pValueBetaCounts[0][code]/float(1+numTests)) + " ")
    fileTemp3.write(str(pValueBetaCounts[1][code]/float(1+numTests)) + " ")
    fileTemp3.write("\n")
  code = "ALL"  # fake for overall
  fileTemp2.write("ALL " + str(totalMeans[0]-totalMeans[1]) + " ")
  fileTemp2.write(str(totalMeans[0]) + " " + str(totalMeans[1]) + " ")
  fileTemp2.write(str(pValueCounts[0][code]/float(1+numTests)) + " ")
  fileTemp2.write(str(pValueCounts[1][code]/float(1+numTests)) + " ")
  fileTemp2.write("\n")
  fileTemp3.write("ALL " + str(totalBetaMeans[0]-totalBetaMeans[1]) + " ")
  fileTemp3.write(str(totalBetaMeans[0]) + " " + str(totalBetaMeans[1]) + " ")
  fileTemp3.write(str(pValueBetaCounts[0][code]/float(1+numTests)) + " ")
  fileTemp3.write(str(pValueBetaCounts[1][code]/float(1+numTests)) + " ")
  fileTemp3.write("\n")
  fileTemp2.close()
  fileTemp3.close()
  return totalMeans[0]-totalMeans[1], totalBetaMeans[0]-totalBetaMeans[1]

def makeCompareResidueReportVersus(pdbs, outputFilename="residue.each.bfactor"):
  '''outputs the averages from each residue across each pdb file, for
  making comparison graphs'''
  pass

def makeCompareResidueReport(
    residueBoth, outputFilename="residue.bfactor", maxY=False, maxYBeta=False,
    numTests=9):
  ranges = [-0.3, 0.6]
  residueNames = []
  for residueName in residueBoth[0].keys() + residueBoth[1].keys():
    if residueName not in residueNames:
      residueNames.append(residueName)
  residueNames.sort()
  #residueNames = aminoAcid3Codes #for now ignore what is in the files
  fileTemp = open(outputFilename + ".txt", 'w')
  fileTemp.write("ResidueName AtomName Mean StdDev Low High Count\n")
  fileTemp2 = open(outputFilename + ".pvals.txt", 'w')
  fileTemp2.write("ResidueName DiffMeans MeanA MeanB pValAbove pValBelow\n")
  fileTemp3 = open(outputFilename + ".pvals.beta.txt", 'w')
  fileTemp3.write("ResidueName DiffMeans MeanA MeanB pValAbove pValBelow\n")
  averages, stddevs = ({}, {}), ({}, {})
  betaAverages, betaStddevs = ({}, {}), ({}, {})
  totalLists, betaLists = ({}, {}), ({}, {})
  for residueName in residueNames:
    totalList = [], []
    betaList = [], []
    for indexSet, residueData in enumerate(residueBoth):
      try:
        for data in residueData[residueName].values():
          totalList[indexSet].extend(data)
        totalLists[indexSet][residueName] = totalList[indexSet]
        average = statistics.computeMean(totalList[indexSet])
        averages[indexSet][residueName] = average
        #print average, residueName
        stddev = statistics.computeStdDev(totalList[indexSet], average)
        stddevs[indexSet][residueName] = stddev
        data = residueData[residueName]
        betaList[indexSet].extend(data[carbonBetaCodes[residueName]])
        betaLists[indexSet][residueName] = betaList[indexSet]
        if len(betaList[indexSet]) > 0:
          betaAvg = statistics.computeMean(betaList[indexSet])
          #print betaAvg, residueName
          betaAverages[indexSet][residueName] = betaAvg
          betaStddevs[indexSet][residueName] = statistics.computeStdDev(
              betaList[indexSet], betaAvg)
        fileTemp.write(
            residueName + " " + str(average) + " " + str(stddev) + " " +
            str(min(totalList)) + " " + str(max(totalList)) + " " +
            str(len(totalList)) + "\n")
      except (ZeroDivisionError, KeyError):
        pass  # probably don't really need this residue anyway
  fileTemp.close()
  for index, code in enumerate(aminoAcid3Codes):  # now do the pvalue tests
    meanA = averages[0][code]
    meanB = averages[1][code]
    listA = totalLists[0][code]
    listB = totalLists[1][code]
    pvals = statistics.pvalueDiffMeans(listA, listB, meanA-meanB, numTests)
    #fileTemp2.write("ResidueName DiffMeans MeanA MeanB pValAbove pValBelow\n")
    fileTemp2.write(code + " " + str(meanA-meanB) + " " + str(meanA) + " ")
    fileTemp2.write(str(meanB) + " " + str(pvals[0]) + " " + str(pvals[1]))
    fileTemp2.write("\n")
    meanA = betaAverages[0][code]
    meanB = betaAverages[1][code]
    listA = betaLists[0][code]
    listB = betaLists[1][code]
    pvals = statistics.pvalueDiffMeans(listA, listB, meanA-meanB, numTests)
    fileTemp3.write(code + " " + str(meanA-meanB) + " " + str(meanA) + " ")
    fileTemp3.write(str(meanB) + " " + str(pvals[0]) + " " + str(pvals[1]))
    fileTemp3.write("\n")
  fileTemp2.close()
  fileTemp3.close()
  if gnuplotAvailable:
    plotter = Gnuplot.Gnuplot(debug=0)
    yLabels = '('
    yData, yError, yMin, yMax = [], [], 10, 0
    yBetaData, yBetaError, yBetaMin, yBetaMax = [], [], 10, 0
    for index, code in enumerate(aminoAcid3Codes):
      yLabels += '"' + str(code) + '" ' + str(index)
      if index != len(aminoAcid3Codes) - 1:
        yLabels += ', '
      yData.append(averages[0][code] - averages[1][code])
      #yError.append(stddevs[0][code])
      #yMin = min(yMin, yData[-1] - yError[-1])
      #yMax = max(yMax, yData[-1] + yError[-1])
      yMin = min(yMin, yData[-1])
      yMax = max(yMax, yData[-1])
      #print betaAverages[0][code]
      #print betaAverages[1][code]
      betaAvgDiff = 0.
      try:
        betaAvg0 = betaAverages[0][code]
        betaAvg1 = betaAverages[1][code]
        betaAvgDiff = betaAvg0 - betaAvg1
      except KeyError:
        print code
        betaAvgDiff = 0.
      yBetaData.append(betaAvgDiff)
      #yBetaError.append(betaStddevs[0][code])
      yBetaMin = min(yBetaMin, yBetaData[-1])
      yBetaMax = max(yBetaMax, yBetaData[-1])
    yLabels += ')'
    graphData = Gnuplot.Data(range(20), yData)
    plotter('set terminal png')
    plotter('set output "' + outputFilename + '.png"')
    plotter('set data style points')
    plotter('set boxwidth 0.9 absolute')
    plotter('set xtics ' + yLabels)
    if ranges:
      plotter('set yrange [' + str(ranges[0]) + ':' + str(ranges[1]) + ']')
    elif maxY is False:
      plotter('set yrange [' + str(yMin-0.2) + ':' + str(yMax+0.2) + ']')
    else:
      plotter('set yrange [0:' + str(maxY) + ']')
    plotter('set xrange [-1:20]')
    plotter.xlabel('Residue')
    plotter.ylabel('Mean Travel In Distance')
    plotter.plot(graphData)
    #do another graph with just carbon-betas
    plotter('set output "' + outputFilename + '.beta.png"')
    graphDataBeta = Gnuplot.Data(range(20), yBetaData)
    plotter.ylabel('Mean Travel In Distance of Carbon Beta')
    if ranges:
      plotter('set yrange [' + str(ranges[0]) + ':' + str(ranges[1]) + ']')
    elif maxYBeta is False:
      plotter(
          'set yrange [' + str(yBetaMin-0.2) + ':' + str(yBetaMax+0.2) + ']')
    else:
      plotter('set yrange [0:' + str(maxYBeta) + ']')
    plotter.plot(graphDataBeta)

def makeAminoAcidHistogram(
    plotter,  backValues, sideValues, caValues, cbValues, outNameAmino,
    interval=0.5):
  '''makes one histogram for one amino acid for the backbone/sidechain values'''
  outTextFileName = outNameAmino + ".res.txt"
  maxBoth = max(backValues+sideValues+caValues+cbValues)
  #print outNameAmino,
  #print len(backValues), len(sideValues), len(caValues), len(cbValues), maxBoth
  histoBack, maxBothRet = statistics.computeHistogram(
      backValues, interval, maxBoth)
  histoSide, maxBothRet = statistics.computeHistogram(
      sideValues, interval, maxBoth)
  histoCa, maxBothRet = statistics.computeHistogram(
      caValues, interval, maxBoth)
  histoCb, maxBothRet = statistics.computeHistogram(cbValues, interval, maxBoth)
  xVals = []
  for cutoff in range(2+int(maxBoth/interval)):
    xVals.append(cutoff*interval)
  graphDataBack = Gnuplot.Data(xVals, histoBack, title="Backbone")
  graphDataSide = Gnuplot.Data(xVals, histoSide, title="Sidechain")
  graphDataCa = Gnuplot.Data(xVals, histoCa, title="C-alpha")
  graphDataCb = Gnuplot.Data(xVals, histoCb, title="C-beta")
  if outTextFileName:
    outTextFile = open(outTextFileName, 'w')
    outTextFile.write("xVal\tback\tside\tca\tcb\tbackN\tsideN\tcaN\tcbN\n")
    for index in range(max(len(xVals), len(histoBack))):
      outTextFile.write(str(xVals[index]) + "\t")
      outTextFile.write(str(histoBack[index]) + "\t")
      outTextFile.write(str(histoSide[index]) + "\t")
      outTextFile.write(str(histoCa[index]) + "\t")
      outTextFile.write(str(histoCb[index]) + "\t")
      if sum(histoBack) > 0.:
        outTextFile.write(str(float(histoBack[index])/sum(histoBack)) + "\t")
      else:
        outTextFile.write("0.\t")
      if sum(histoSide) > 0.:
        outTextFile.write(str(float(histoSide[index])/sum(histoSide)) + "\t")
      else:
        outTextFile.write("0.\t")
      if sum(histoCa) > 0.:
        outTextFile.write(str(float(histoCa[index])/sum(histoCa)) + "\t")
      else:
        outTextFile.write("0.\t")
      if sum(histoCb) > 0.:
        outTextFile.write(str(float(histoCb[index])/sum(histoCb)) + "\n")
      else:
        outTextFile.write("0.\n")
    outTextFile.close()
  plotter('set terminal png')
  plotter('set output "' + outNameAmino + '.png"')
  #plotter('set data style boxes')
  plotter('set data style linespoints')
  plotter('set key right top')
  plotter('set xrange [' + str(min(xVals)-1) + ':' + str(max(xVals)+1) + ']')
  plotter(
      'set yrange [' + str(0.) + ':' + str(1.05*max(histoBack+histoSide)) + ']')
  plotter.xlabel('Travel In Distance')
  plotter.ylabel('Atom Count')
  plotter.plot(graphDataBack, graphDataSide, graphDataCa, graphDataCb)

def makeAtomReport(residueData, outputFilename="atom.bfactor", runGraphs=True):
  residueNames = residueData.keys()
  residueNames.sort()
  fileTemp = open(outputFilename + '.txt', 'w')
  fileTemp.write("ResidueName AtomName Mean StdDev Low High Count\n")
  resAtomAverage = {}
  for residueName in residueNames:
    resAtomAverage[residueName] = {}
    atomNames = residueData[residueName].keys()
    atomNames.sort()
    for atomName in atomNames:
      totalList = residueData[residueName][atomName]
      average = statistics.computeMean(totalList)
      resAtomAverage[residueName][atomName] = average
      stddev = statistics.computeStdDev(totalList, average)
      fileTemp.write(
          residueName + " " + atomName + " " + str(average) + " " +
          str(stddev) + " " + str(min(totalList)) + " " +
          str(max(totalList)) + " " + str(len(totalList)) + "\n")
  fileTemp.close()
  if gnuplotAvailable and runGraphs:
    #first make backbone-sidechain report
    plotter = Gnuplot.Gnuplot(debug=0)
    yLabels = '('
    yDataBackbone, yDataSidechain = [], []
    yDataCa, yDataCb = [], []
    for index, code in enumerate(aminoAcid3Codes):
      yLabels += '"' + str(code) + '" ' + str(index)
      if index != len(aminoAcid3Codes) - 1:
        yLabels += ', '
      backValues, sideValues = [], []
      caValues, cbValues = [], []
      try:
        for key, values in residueData[code].iteritems():
          if string.strip(key) in backboneAtomCodes:
            backValues.extend(values)
          else:
            sideValues.extend(values)
          if string.strip(key) == caCode:
            caValues.extend(values)
          elif string.strip(key) == cbCode:
            cbValues.extend(values)
      except KeyError:  # sometimes one residue won't be represented
        pass  # but that is okay
      if len(backValues) == 0:
        yDataBackbone.append(0)
      else:
        yDataBackbone.append(sum(backValues)/float(len(backValues)))
      if len(sideValues) == 0:
        yDataSidechain.append(0)
      else:
        yDataSidechain.append(sum(sideValues)/float(len(sideValues)))
      if len(caValues) == 0:
        yDataCa.append(0)
      else:
        yDataCa.append(sum(caValues)/float(len(caValues)))
      if len(cbValues) == 0:
        yDataCb.append(0)
      else:
        yDataCb.append(sum(cbValues)/float(len(cbValues)))
      if len(backValues + sideValues + caValues + cbValues) > 0:
        makeAminoAcidHistogram(
            plotter, backValues, sideValues, caValues, cbValues,
            outputFilename + "." + str(code))
    yLabels += ')'
    graphDataBackbone = Gnuplot.Data(range(20), yDataBackbone, title="Backbone")
    graphDataSidechain = Gnuplot.Data(
        range(20), yDataSidechain, title="Sidechain")
    graphDataCa = Gnuplot.Data(range(20), yDataCa, title="C-alpha")
    graphDataCb = Gnuplot.Data(range(20), yDataCb, title="C-beta")
    plotter('set terminal png')
    plotter('set output "' + outputFilename + '.png"')
    plotter('set data style points')
    plotter('set key right top')
    plotter('set xtics ' + yLabels)
    plotter(
        'set yrange [' + str(min(yDataBackbone + yDataSidechain) - 0.5) +
        ':' + str(max(yDataBackbone+yDataSidechain)+0.5) + ']')
    plotter('set xrange [-1:20]')
    plotter.xlabel('Residue')
    plotter.ylabel('Mean Travel In Distance')
    plotter.plot(graphDataBackbone, graphDataSidechain)
    plotter('set output "' + outputFilename + '.ab.png"')
    if "buried" in outputFilename:
      plotter('set yrange [' + str(min(yDataCa + yDataCb) - 0.5) + ':6.]')
    else:
      plotter(
          'set yrange [' + str(min(yDataCa + yDataCb) - 0.5) + ':' +
          str(max(yDataCa + yDataCb) + 0.5) + ']')
    plotter.plot(graphDataCa, graphDataCb)

def makeHistogramReport(
    residueData, outputFilename="histogram.bfactor", interval=0.1,
    yMaxHisto=0.16):
  fileTemp = open(outputFilename + ".txt", 'w')
  fileTemp.write("LowEndInterval Count\n")
  totalList = []
  betaList = []
  for oneResName, oneResData in residueData.iteritems():
    #assemble into one big list
    if oneResName in aminoAcid3Codes:
      for data in oneResData.values():
        totalList.extend(data)
      betaList.extend(oneResData[carbonBetaCodes[oneResName]])
  resList = {}
  betaResList = {}
  for oneResKey in aminoAcid3Codes:
    resList[oneResKey] = []
    betaResList[oneResKey] = []
    if oneResKey in residueData:
      for data in residueData[oneResKey].values():
        resList[oneResKey].extend(data)
      betaResList[oneResKey].extend(
          residueData[oneResKey][carbonBetaCodes[oneResKey]])
  #now do histogram stuff
  histogram, maxData = statistics.computeNormalizedHistogram(
      totalList, interval, 8.)
  betaHistogram, betaMax = statistics.computeNormalizedHistogram(
      betaList, interval, 8.)
  for index, data in enumerate(histogram):
    fileTemp.write(str(index*interval) + " " + str(data) + "\n")
  fileTemp.close()
  #make gnuplot version if possible
  if gnuplotAvailable:
    plotter = Gnuplot.Gnuplot(debug=0)
    xVals = []
    for cutoff in range(2+int(maxData/interval)):
      xVals.append(cutoff*interval)
    graphData = Gnuplot.Data(xVals, histogram)
    graphDataCum = Gnuplot.Data(
        xVals, statistics.computeCumulativeHistogram(histogram))
    graphDataBeta = Gnuplot.Data(xVals, betaHistogram)
    graphDataBetaCum = Gnuplot.Data(
        xVals, statistics.computeCumulativeHistogram(betaHistogram))
    plotter('set terminal png')
    plotter('set output "' + outputFilename + '.png"')
    plotter('set data style linespoints')
    plotter('set xrange [' + str(min(xVals)-1) + ':' + str(max(xVals)+1) + ']')
    plotter('set yrange [' + str(0.) + ':' + str(yMaxHisto) + ']')
    plotter.xlabel('Travel In Distance')
    plotter.ylabel('Atom Count')
    plotter.plot(graphData)
    plotter('set output "' + outputFilename + '.beta.png"')
    plotter.ylabel('Carbon Beta Atom Count')
    plotter.plot(graphDataBeta)
    plotter('set output "' + outputFilename + '.cumulative.png"')
    plotter('set yrange []')  # automatic
    plotter.ylabel('Cumulative Atom Count')
    plotter.plot(graphDataCum)
    plotter.ylabel('Cumulative Beta Atom Count')
    plotter('set output "' + outputFilename + '.cumulative.beta.png"')
    plotter.plot(graphDataBetaCum)
    plotter.ylabel('Atom Count')
    #now do one for each residue
    plotter('set key right top')
    plotter('set data style lines')
    ylabels = ('Atom Count', 'Carbon Beta Atom Count')
    outputNames = (outputFilename, outputFilename + '.beta')
    histogramData = (resList, betaResList)
    for index in range(len(histogramData)):
      thisResList = histogramData[index]
      outputName = outputNames[index]
      ylabel = ylabels[index]
      plotter.ylabel(ylabel)
      histograms = {}
      maxOverRes = 0.
      for resName in aminoAcid3Codes:
        histogramRes, maxData = statistics.computeHistogram(
            thisResList[resName], interval, maxData)
        histograms[resName] = histogramRes
        maxOverRes = max(maxOverRes, max(histogramRes))
      plotter('set yrange [0:' + str(maxOverRes+1000) + ']')
      resGraphDatas = [], []
      lowGraphDatas = [], []
      highGraphDatas = [], []
      for resName in aminoAcid3Codes:
        plotter('set output "' + outputName + "." + resName + '.png"')
        plotter('set yrange [0:' + str(maxOverRes + 1000) + ']')
        histogramRes = histograms[resName]
        resGraphDatas[0].append(
            Gnuplot.Data(xVals, histogramRes, title=resName))
        plotter.plot(resGraphDatas[0][-1])
        plotter(
            'set output "' + outputName + "." + resName + '.cumulative.png"')
        plotter('set yrange [0:1]')
        cumData = Gnuplot.Data(
            xVals, statistics.computeCumulativeHistogram(histogramRes),
            title=resName)
        plotter.plot(cumData)
        resGraphDatas[1].append(cumData)
        if resName in highCodes:
          highGraphDatas[0].append(resGraphDatas[0][-1])
          highGraphDatas[1].append(cumData)
        elif resName in lowCodes:
          lowGraphDatas[0].append(resGraphDatas[0][-1])
          lowGraphDatas[1].append(cumData)
      #very stupid hack... plot() is dumb
      outNames = (
          'set output "' + outputName + ".residues" + '.png"',
          'set output "' + outputName + ".residues" + '.cumulative.png"')
      lowNames = (
          'set output "' + outputName + ".residues.low" + '.png"',
          'set output "' + outputName + ".residues.low" + '.cumulative.png"')
      highNames = (
          'set output "' + outputName + ".residues.high" + '.png"',
          'set output "' + outputName + ".residues.high" + '.cumulative.png"')
      ranges = (
          'set yrange [0:' + str(maxOverRes+1000) + ']', 'set yrange [0:1]')
      for count, resGraphData in enumerate(resGraphDatas):
        plotter(outNames[count])
        plotter('set key right bottom')
        plotter(ranges[count])
        plotter.plot(
            resGraphData[0], resGraphData[1], resGraphData[2],
            resGraphData[3], resGraphData[4], resGraphData[5],
            resGraphData[6], resGraphData[7], resGraphData[8],
            resGraphData[9], resGraphData[10], resGraphData[11],
            resGraphData[12], resGraphData[13], resGraphData[14],
            resGraphData[15], resGraphData[16], resGraphData[17],
            resGraphData[18], resGraphData[19])
        lowGraphData = lowGraphDatas[count]
        plotter(lowNames[count])
        plotter.plot(
            lowGraphData[0], lowGraphData[1], lowGraphData[2],
            lowGraphData[3], lowGraphData[4], lowGraphData[5],
            lowGraphData[6], lowGraphData[7], lowGraphData[8],
            lowGraphData[9], lowGraphData[10])
        highGraphData = highGraphDatas[count]
        plotter(highNames[count])
        plotter.plot(
            highGraphData[0], highGraphData[1], highGraphData[2],
            highGraphData[3], highGraphData[4], highGraphData[5],
            highGraphData[6], highGraphData[7], highGraphData[8])
      #limits on x dim
      outNames2 = (
          'set output "' + outputName + ".residues16" + '.png"',
          'set output "' + outputName + ".residues" + '.cumulative16.png"')
      for count, resGraphData in enumerate(resGraphDatas):
        plotter(outNames2[count])
        plotter('set key right bottom')
        plotter('set xrange [1:6]')
        plotter(ranges[count])
        plotter.plot(
            resGraphData[0], resGraphData[1], resGraphData[2],
            resGraphData[3], resGraphData[4], resGraphData[5],
            resGraphData[6], resGraphData[7], resGraphData[8],
            resGraphData[9], resGraphData[10], resGraphData[11],
            resGraphData[12], resGraphData[13], resGraphData[14],
            resGraphData[15], resGraphData[16], resGraphData[17],
            resGraphData[18], resGraphData[19])

def analyzePdbB(filenameList=False):
  residues = {}  # a dict of dicts where the sub-dicts are keyed on atom name
  eachRes = {}  # dict keyed on RESNAME-RESNUM
  if filenameList:
    for filename in filenameList:
      pdbD = pdb.pdbData(filename)
      for index, resName in enumerate(pdbD.resNames):
        if pdbD.radii[index] > 0.:
          if resName not in residues:
            residues[resName] = {}  # init sub-dict
          if string.strip(pdbD.atoms[index]) not in residues[resName]:
            residues[resName][string.strip(pdbD.atoms[index])] = []
          residues[resName][string.strip(pdbD.atoms[index])].append(
              pdbD.factors[index][1])
          resNum = pdbD.resNums[index]
          longName = str(resName) + str(resNum)
          if longName not in eachRes:
            eachRes[longName] = {}  # init sub-dict
          if string.strip(pdbD.atoms[index]) not in eachRes[longName]:
            eachRes[longName][string.strip(pdbD.atoms[index])] = []
          eachRes[longName][string.strip(pdbD.atoms[index])].append(
              pdbD.factors[index][1])
    #residues now contains all the b-factor (travelin) data
  makeResidueReport(residues)
  makeResidueReport(
      eachRes, outputFilename="individual.res.bfactor", runGraphs=False)
  makeAtomReport(residues)
  makeAtomReport(
      eachRes, outputFilename="individual.atom.bfactor", runGraphs=False)
  makeHistogramReport(residues)

def makeMeanPerProteinReport(pdbs, outName):
  '''pdbs is a list of dict of dicts, outname is filename'''
  outFile = open(outName, 'w')
  for pdb in pdbs:
    totalList = []
    for resList in pdbs[pdb].itervalues():
      for atomList in resList.itervalues():
        totalList.extend(atomList)
    avg = statistics.computeMean(totalList)
    outFile.write(pdb + "\t")
    outFile.write(str(avg) + "\n")
  outFile.close()

#this is main
if -1 != string.find(sys.argv[0], "analyzePdbB.py"):
  if len(sys.argv) > 1:
    analyzePdbB(sys.argv[1:])
