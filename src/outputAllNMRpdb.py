#!/usr/bin/env python
#ryan g. coleman ryangc@mail.med.upenn.edu

#outputs all NMR models as single pdb files

#grab system, string,  regular expression, and operating system modules
import sys
import string
import re
import os
import math
import pdb  # for chain sorting ease

if -1 != string.find(sys.argv[0], "outputAllNMRpdb.py"):
  fileName = sys.argv[1]
  pdbD = pdb.pdbData(fileName)
  modelNums = pdbD.getModelNumbers()
  for modelNum in modelNums:
    newPdb = pdbD.getOneModel(modelNum)
    outputFileName = sys.argv[1][:-4] + "." + str(modelNum) + ".pdb"
    newPdb.write(outputFileName)
