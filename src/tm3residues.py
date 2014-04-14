#!/usr/bin/env python

#reports pockets with highest overlap with the set of residues

import tm3
import string
import sys
import statistics
import dot

if -1 != string.find(sys.argv[0], "tm3residues.py"):
  prefix = "temp"
  tmDataList = []
  nameOrNum = sys.argv[1]
  resList = sys.argv[2]
  for filename in sys.argv[3:]:
    tmDataList.append(tm3.tmTreeFromFile(filename))
  colName = "Residue Number List"
  if nameOrNum == "name":
    colName = "Residue Name List"
  dotData = dot.dotResidues(
      tmDataList, resList, colNameIn=colName, printHelpful=False)
  dotData.printHelpful(prefix)
  dotData.write(prefix + ".dot")
  dotData.writeGdl(prefix + ".gdl", force=False)
  print "dot -Tpng "+prefix+".dot > "+prefix+".png"
  print "output files written to "+prefix+".tab.txt and "+prefix+".animation.py"
