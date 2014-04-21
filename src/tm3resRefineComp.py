#!/usr/bin/env python

#reports pockets with highest overlap with the set of residues

import tm3
import string
import sys
import dot

if -1 != string.find(sys.argv[0], "tm3resRefineComp.py"):
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
  #print dotData.treeToBestNode
  #dotData.write("temp.dot")
  #dotData.writeGdl("temp.gdl", force=False)
  #print "dot -Tpng temp.dot > temp.png"
  matchList = tm3.findSimilarNodes(
      tmDataList, [
          "Surface Area", "Volume", "height",
          #"mean absolute Curvature",
          "mean Curvature", "mouths",
          "longest dimension", "middle dimension", "short dimension",
          "Area of Biggest Mouth", "Diameter of Biggest Mouth", "mean height"],
      10000000, 1000, 50,
      resCols=["Atom Name List"], sizeColName="Volume",
      sizeMin=-1., sizeMax=10000000000,
      doSelfScore=False, justNodes=dotData.treeToBestNode,
      mst=True, lineMst=True, lineMstEnds=False,
      refinePockets=True, possNodes=dotData.treeToPossNodes,
      dotResidues=dotData)
  dotData.printHelpful()
  #linemstends doesn't work that well so disabled for now
  names, matrix = [], {}
  for match in matchList:
    name1 = match[0].inputFileName + "-" + match[2].getId()
    name2 = match[1].inputFileName + "-" + match[3].getId()
    for name, otherName in [(name1, name2), (name2, name1)]:
      if name not in names:
        names.append(name)
      if name not in matrix:
        matrix[name] = {}
      matrix[name][otherName] = match[4]
      matrix[name][name] = 0.  # fill in the diagonal
  names.sort()
  print "matrix",
  for name in names:
    print name,
  print " "
  for name in names:
    print name,
    for otherName in names:
      print round(matrix[name][otherName], 3),
    print " "
