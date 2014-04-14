#!/usr/bin/env python

#ryan g. coleman ryan.g.coleman@gmail.com ryangc@mail.med.upenn.edu
#kim sharp lab http://crystal.med.upenn.edu

import string
import sys
import pdb

def calcRMSDprintTable(pdbFileNames):
  pdbDatas = []
  for pdbFileName in pdbFileNames:
    pdbDatas.append(pdb.pdbData(pdbFileName))
  matrix = []
  for index1 in xrange(len(pdbDatas)):
    tempRow = []
    for index2 in xrange(len(pdbDatas)):
      if index1 > index2:
        rmsd = pdbDatas[index1].calcRMSD(pdbDatas[index2], alphas=True)
        tempRow.append(rmsd)
    matrix.append(tempRow)
  print "matrix",
  for index1 in xrange(len(pdbDatas)):
    print pdbFileNames[index1],
  print " "  # title row done
  for index1 in xrange(len(pdbDatas)):
    print pdbFileNames[index1],
    for index2 in xrange(len(pdbDatas)):
      if index1 == index2:  # always 0
        print 0.,
      elif index1 > index2:
        print matrix[index1][index2],
      else:
        print matrix[index2][index1],
    print " "  # this row done

if -1 != string.find(sys.argv[0], "rmsdTable.py"):
  try:
    calcRMSDprintTable(sys.argv[1:])
  except IndexError:
    print "rmsdTable.py list_of_pdb_names"
    print "outputs to standard out"
    sys.exit(1)
