#!/usr/bin/env python2.5

#reads in tm3 files, ignores tree structure, just looks for nodes/pockets
#that are similar

import tm3, string, sys

def animateScript(tmData, tstData, sortCol, outName):
  '''sort by threshold, make a script to open and display each pocket in turn'''
  nodes = tmData.idNode.values()
  nodes.sort(lambda x,y: cmp(x.attributes[sortCol], y.attributes[sortCol]))
  nodes.reverse() #sorted high to low now
  print "from pymol import cmd"
  print 'cmd.do("delete all")'
  print 'cmd.do("tstOpen ' + tstData + '")'
  print 'cmd.do("tstDraw surf=pdb")'
  print 'cmd.do("color grey")'
  for index,node in enumerate(nodes):
    print 'cmd.do("turn x,1")'
    print 'cmd.do("tstDraw group=', str(node.getId()), '")'
    if index % 5 == 0:
      print 'cmd.do("ray")'
      print 'cmd.do("png ', str(outName + string.zfill(index, 6)+'.png') + '")'
  print 'cmd.do("ray")'
  print 'cmd.do("png ', str(outName + string.zfill(index+1, 6)+'.png') + '")'
  

if -1 != string.find(sys.argv[0], "tm3animatepockets.py"):
  tmData = tm3.tmTreeFromFile(sys.argv[1])
  tstData = sys.argv[1].replace(".tree.tm3", "")
  animateScript(tmData, tstData, 3, tstData+".movie.turn.") 
