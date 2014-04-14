#!/usr/bin/env python

#copyright 2005-2009 Ryan G. Coleman, Kim A. Sharp, University of Pennsylvania
#cite the following papers if you use this code:
#
#Coleman, RG, Sharp, KA. Travel Depth, A New Shape Descriptior for
#Macromolecules: Application to Ligand Binding. Journal of Molecular
#Biology 362(3), pp. 441-458, 22 September 2006.
#http://dx.doi.org/10.1016/j.jmb.2006.07.022
#
#Coleman, RG, Sharp, KA. Shape and evolution of thermostable protein
#structure. Proteins: Structure, Function, and Bioinformatics. 2009.
#http://dx.doi.org/10.1002/prot.22558

import sys,string
import tstTravelDist

if -1 != string.find(sys.argv[0], "tstBurialDepth.py"):
  if len(sys.argv) >= 3:
    tstFile, phiFile = sys.argv[1:3]
    tstTravelDist.tstTravelSurfInsideMesh(tstFile, phiFile)
