#!/usr/bin/env python

#ryan g. coleman ryan.g.coleman@gmail.com ryangc@mail.med.upenn.edu
#kim sharp lab http://crystal.med.upenn.edu

import string
import sys
import geometry
import pdb

if -1 != string.find(sys.argv[0], "pdbDistances.py"):
  try:
    for pdbName in sys.argv[1:]:
      pdbD = pdb.pdbData(pdbName)
      outputName = pdbName.replace("pdb", "").replace(".", "")
      longestDist, meanDist = geometry.longestAndMeanDist(
          pdbD.getHeavyAtomXYZ())
      print outputName, "\t", longestDist, "\t", meanDist
  except IndexError:
    print "pdbDistances.py pdbName [list of more pdbs]"
    print "outputs to standard out"
    sys.exit(1)
