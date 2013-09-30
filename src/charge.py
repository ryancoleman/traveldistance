#!/usr/bin/env python2.5
# Ryan G. Coleman, Kim A. Sharp, http://crystal.med.upenn.edu, 
# ryan.g.coleman@gmail.com ryangc@mail.med.upenn.edu

#file takes care of reading .crg files and putting them into dictionary struct
#charge files are column based and look like this:

#NH1    ARG      -0.350
#0123456789012345678901234
#          1         2

import string, sys, os

class charge(object):
  '''reads in a .crg file, makes a usable data structure'''
  hydrophobicThreshold = 0.45

  def __init__(self, \
        chargeFileName="$TDHOME/src/charge_parameters/parse.crg",
        altChg="/disks/node18/coleman/src/charge_parameters/parse.crg"):
    '''reads in the file, sets up the structure, etc'''
    chargeFileName = os.path.expandvars(chargeFileName)
    altChg = os.path.expandvars(altChg)
    self.residues = {}
    self.atoms = {}
    try:
      chargeFile = open(chargeFileName, 'r')
    except IOError:
      chargeFile = open(altChg, 'r')
    try:
      for line in chargeFile:
        if line[0] == '!' or line[:22] == 'atom__resnumbc_charge_':
          pass #means a comment
        else:
          try:
            atom = string.strip(line[0:4]).upper()
            res = string.strip(line[5:14]).upper()
            ch = float(line[15:22])
            if len(res) == 0: #no residue, default for atoms if not found
              self.atoms[atom] = ch
            else:
              if res not in self.residues:
                self.residues[res] = {}
              self.residues[res][atom] = ch
          except TypeError:
            print "warning: error reading line: " + line
    except StopIteration:
      pass #EOF

  def getCharge(self, atomName, resName):
    '''given a residue and atom, find the charge'''
    atomUp = string.strip(atomName.upper())
    resUp = string.strip(resName.upper()) #don't want to mess with case sensitivity
    if resUp in self.residues and atomUp in self.residues[resUp]: #normal
      return self.residues[resUp][atomUp]
    elif atomUp in self.atoms: #trying to use default
      return self.atoms[atomUp]
    elif atomUp[0] in self.atoms: #trying to use just atom type
      return self.atoms[atomUp[0]]
    else: #can't do it
      return None 

  def getTrinaryCharge(self, atomName, resName):
    '''returns -1, 0 or +1 depending on charge and hydrophobic threshold'''
    tempCharge = self.getCharge(atomName, resName)
    if tempCharge is None:
      return None
    elif abs(tempCharge) < charge.hydrophobicThreshold:
      return 0
    elif tempCharge > 0:
      return +1
    else: #must be < 0
      return -1

#main, only runs if testing reading of files
if -1 != string.find(sys.argv[0], "charge.py"):
  if len(sys.argv) >= 2:
    chargeD = charge(sys.argv[1])
    print chargeD.residues, chargeD.atoms
    print chargeD.getCharge('CE1', 'HSE')
