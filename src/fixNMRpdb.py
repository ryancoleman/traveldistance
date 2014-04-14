#!/usr/bin/env python
#ryan g. coleman ryangc@mail.med.upenn.edu
#usage: fixNMRpdb.py file

#file.pdb is copied to file.one.pdb and all but first model are removed

import sys
import string
import re
import os
import math
import pdb  # for chain sorting ease

if -1 != string.find(sys.argv[0], "fixNMRpdb.py"):
  fileName = sys.argv[1]
  pdbD = pdb.pdbData(fileName)
  modelNums = pdbD.getModelNumbers()
  newPdb = pdbD.getOneModel(modelNums[0])
  outputFileName = sys.argv[1][:-4] + ".one.pdb"
  newPdb.write(outputFileName)
