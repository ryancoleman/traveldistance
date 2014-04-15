#!/usr/bin/env python

#copyright 2005-2009 Ryan G. Coleman, Kim A. Sharp, University of Pennsylvania
#cite the following papers if you use this code:
#
#Coleman, RG, Sharp, KA. Travel Depth, A New Shape Descriptior for
#Macromolecules: Application to Ligand Binding. Journal of Molecular
#Biology 362(3), pp. 441-458, 22 September 2006.
#http://dx.doi.org/10.1016/j.jmb.2006.07.022
#
#Coleman, RG, Sharp, KA. Finding and characterizing tunnels in
#macromolecules with application to ion channels and pores. Biophysical
#Journal. January 2009.


import sys
import string
import tstTravelDist

if -1 != string.find(sys.argv[0], "tstChunnel.py"):
  if len(sys.argv) >= 3:
    tstFile, phiFile = sys.argv[1:3]
    tstTravelDist.tstTravelFindHoles(tstFile, phiFile, debugOut=False)
