#!/usr/bin/env python

#reads in what are assumed to be ligand pdb files
#for each one, cluster based on some threshold, break into distinct clusters,
#write each cluster

import string
import sys
import pdb

if -1 != string.find(sys.argv[0], "clusterLigands.py"):
  filenames = sys.argv[1:]
  for filename in filenames:
    pdbD = pdb.pdbData(filename)
    clusters = pdbD.clusterAtoms(distanceCutoff=5.0)
    #print filename, len(clusters) #debug
    width = len(str(len(clusters)))
    for clusterIndex, cluster in enumerate(clusters):
      outputNum = string.zfill(clusterIndex, width)
      outputName = outputNum + "_" + filename
      cluster.write(outputName)
