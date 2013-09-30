#!/usr/bin/env python2.5
#ryan g. coleman ryangc@mail.med.upenn.edu 
#usage: fixNMRpdb.py file

#file.pdb is copied to file.one.pdb and all but first model are removed

#grab system, string,  regular expression, and operating system modules
import sys, string, re, os, math
import pdb #for chain sorting ease

if -1 != string.find(sys.argv[0], "fixNMRpdb.py"): 
  fileName = sys.argv[1]
  pdbD = pdb.pdbData(fileName)
  modelNums = pdbD.getModelNumbers()
  newPdb = pdbD.getOneModel(modelNums[0])
  outputFileName = sys.argv[1][:-4] + ".one.pdb"
  newPdb.write(outputFileName)
  
