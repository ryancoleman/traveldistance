#!/usr/bin/env python2.5

#reads in a list of codes, finds pdb + ligand files, outputs 
#residues nearby each ligand in the pdb

import string, sys, pdb, glob

if -1 != string.find(sys.argv[0], "getResiduesNearbyLigands.py"):
  prefixes = sys.argv[1:]
  for prefix in prefixes:
    files = glob.glob("*" + prefix + "*pdb")
    mainPdb, ligandPdbs  = False, []
    for filename in files:
      if -1 == filename.find("ligand"): #is main
        mainPdb = filename
      elif -1 == filename.find("nearby"): #is not output from previous run
        ligandPdbs.append(filename)
    for ligandName in ligandPdbs:
      ligandPdbD = pdb.pdbData(ligandName)
      mainPdbD = pdb.pdbData(mainPdb)
      nearbyPdb = mainPdbD.getNearbyResidues(ligandPdbD.coords, 5.0)
      nearbyPdb.write("nearby_" + ligandName)
      nearbyPdb = pdb.pdbData("nearby_" + ligandName)
      justResString = pdb.turnListIntoString(nearbyPdb.getResidueNamesChains())
      outFile = open("nearby_" + ligandName + ".res", 'w')
      outFile.write(justResString)
      outFile.close()
