#!/usr/bin/env python2.5
#module for importing/exporting pdb files
#rather primitive for now
#ryan g. coleman ryangc@mail.med.upenn.edu

import string,math,sys
import geometry  #distance function
import unionfind2

#declare this here
radii = {'C':1.9,'O':1.6,'N':1.65,'P':1.9,'S':1.9,'H':0.,'F':0.,'I':0.,'U':0.,'A':0.,'B':0.,'L':0.,'*':0.,'Z':0.} #F should really be about 1.5 but not in fortran so not here
#note that changing these breaks compatibility with trisrf/meshsrf surface 
#generation programs and does not actually affect the radii used in those
#processes. in other words don't change them. DON'T DO IT. it won't change 
#the radii used AT ALL, it will just break things.
aminoAcid3Codes = ["ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",\
        "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"]

aminoAcidCodes = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P',\
        'S','T','W','Y','V']

#following are for later features, automatic ANISOU removal and automatic
#nonstandard residue HETATM->ATOM replacement

nonStandardResidues = ['PTR','TPO','SEP','ACE','KCX','MSE','CME','CSD']

removeLinesStartingWith = ['ANISOU']


def turnListIntoString(resChainList, separator='+'):
  '''turns a list into a string without spaces, stupid helper'''
  if len(resChainList) > 0:
    resChainStr = ""
    for resChain in resChainList:
      resChainStr += separator + resChain
    returnStr = string.join(string.split(resChainStr), "")
    return returnStr #do this in case chain is " " 
  else:
    return separator #no residues nearby... this is weird

class pdbData(object):
  '''stores data for a pdb file consisting of atoms'''

  def __init__(self,filename=False):
    '''default constructor takes a pdb file name as input, other ways later'''
    self.__nonZeroRadiiCount = 0
    self.__nonZeroRadiiXYZ = None
    self.__nonZeroRadiiXYZr = None
    self.__factorsByRC = {}
    self.rawData = []
    self.coords = []
    self.radii = []
    self.charges = []
    self.hydroCharges = []
    self.factors = []
    self.atoms = []
    self.resNums = []
    self.resNames = []
    self.chains = []
    self.modelNums = [] #keeps track of NMR models if present
    self.atomToRaw = {}
    self.rawToAtom = {}
    if filename:
      pdbFile = open(filename, 'r')
      try:
        modelNum = 0
        for line in pdbFile:
          if (string.find(line, "MODEL", 0, 5) <> -1):
            modelNum = int(line.split()[1]) #[0] is MODEL, [1] is the number
          self.processLine(line, modelNum)
      except StopIteration:
        pass #read the end of the file
      pdbFile.close()

  def processLines(self, lines):
    '''alternate way to initialize if just given a list of lines'''
    modelNum = 0
    for line in lines:
      if (string.find(line, "MODEL", 0, 5) <> -1):
        modelNum = int(line.split()[1]) #[0] is MODEL, [1] is the number
      self.processLine(line, modelNum)

  def processLine(self,line,modelNumber=0):
    if (string.find(line, "ATOM", 0, 4) <> -1) or \
                     (string.find(line, "HETATM", 0, 6) <> -1):
      name = line[12:16]
      if name[0] == ' ': #because apparently hetatm entries can start one column
        name = name[1:]  #before atom entries (for the atom name)
      else:
        try:
          numberNotLetter = int(name[0])
          name = name[1:] #otherwise would have triggered exception
        except ValueError:
          pass #first character is not a number
      resName = line[17:20]
      if name != 'HOH' and resName != 'HOH': #ignore waters all the time for now
        self.rawData.append(line)
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        chain = line[21:22]
        self.coords.append((x,y,z))
        self.atoms.append(name)
        try:
          radius = radii[name[0]]
        except KeyError:
          radius = 0. #default is to ignore
        #sometimes hydrogens are formatted badly and don't have H in 1st column
        if radius > 0. and name.count('H') > 0: #any H's are bad
          radius = 0.
        self.radii.append(radius)
        if radius > 0.:
          self.__nonZeroRadiiCount += 1
        self.modelNums.append(modelNumber)
        factorStrings = (line[55:60], line[61:66])
        try:
          factors = (float(factorStrings[0]), float(factorStrings[1]))
        except ValueError: #in case no factors or occupancies present
          factors = (0.,0.) #set to zero for now
        self.resNums.append(int(line[22:26]))
        self.chains.append(chain)
        self.resNames.append(resName)
        self.factors.append(factors)
        self.atomToRaw[len(self.atoms)-1] = len(self.rawData)-1
        self.rawToAtom[len(self.rawData)-1] = len(self.atoms)-1

  def write(self, filename):
    '''opens the file, writes to it, closes it'''
    fileOut = open(filename, 'w')
    self.outputLines(fileOut)
    fileOut.close()

  def outputLines(self, fileOut):
    '''writes the lines to an already open file'''
    for rawDataLine in self.rawData:
      if rawDataLine:
        fileOut.write(rawDataLine)
        if not rawDataLine.endswith('\n'):
          fileOut.write('\n')

  def copy(self):
    newPdb = pdbData(False)
    newPdb.__nonZeroRadiiCount = 0
    newPdb.__nonZeroRadiiXYZ = False
    newPdb.__nonZeroRadiiXYZr = False
    newPdb.__factorsByRC = {}
    newPdb.rawData = self.rawData[:] #copy everything
    newPdb.coords = self.coords[:]
    newPdb.atoms = self.atoms[:]
    newPdb.radii = self.radii[:]
    newPdb.charges = self.charges[:]
    newPdb.hydroCharges = self.hydroCharges[:]
    newPdb.factors = self.factors[:]
    newPdb.resNums = self.resNums[:]
    newPdb.chains = self.chains[:]
    newPdb.resNames = self.resNames[:]
    newPdb.modelNums = self.modelNums[:]
    newPdb.atomToRaw = self.atomToRaw.copy()
    newPdb.rawToAtom = self.rawToAtom.copy()
    return newPdb

  def getOrderedRawIndices(self):
    '''gets a list of the ordered raw data indices'''
    indices = self.atomToRaw.values()
    indices.sort()
    return indices

  def updateNewXYZ(self, rawDataIndex, newX, newY, newZ):
    self.coords[self.rawToAtom[rawDataIndex]] = (newX,newY,newZ)
    newData = self.rawData[rawDataIndex][:31] 
    for data in (newX, newY, newZ):
      tempData = "%+6.3f" % data
      tempInt = float(tempData)
      if math.floor(abs(tempInt)) < 10:
        newData += " "
      if math.floor(abs(tempInt)) < 100:
        newData += "%+6.3f" % data + " "
      elif math.floor(abs(tempInt)) < 1000: # one less decimal place
        newData += "%+6.2f" % data + " "
      else: #two less decimal places... if this happens the pdb is huge
        newData += "%+6.1f" % data + " "
    newData += " " + self.rawData[rawDataIndex][56:]
    newDataNoPlus = string.replace(newData, "+", " ")
    self.rawData[rawDataIndex] = newDataNoPlus
    #print newDataNoPlus,

  def updateNewCoord(self, rawDataIndex, newCoord):
    self.updateNewXYZ(rawDataIndex, newCoord[0], newCoord[1], newCoord[2])

  def updateFactors(self, rawDataIndex, newFactors):
    self.factors[rawDataIndex] = newFactors
    #print self.rawData[rawDataIndex] #start
    newData = "" 
    for data in newFactors:
      tempData = "%+2.2f" % data
      newData += "%+2.2f" % data + " "
    newData =  self.rawData[rawDataIndex][:55] + newData + "\n"
    newDataNoPlus = string.replace(newData, "+", " ")
    self.rawData[rawDataIndex] = newDataNoPlus
    #print newDataNoPlus, #changed to

  def getModelNumbers(self):
    '''gets a list of unique model numbers (from NMR data)'''
    uniqueModNums = []
    for modelNum in self.modelNums:
      if modelNum not in uniqueModNums:
        uniqueModNums.append(modelNum)
    uniqueModNums.sort()
    return uniqueModNums

  def getOneModel(self, modelNumber):
    '''copies the pdb, removes all but one NMR model, returns new pdb'''
    newPdb = self.copy()
    markedForRemoval = []
    for index, thisModNum in enumerate(newPdb.modelNums):
      if thisModNum != modelNumber:
        markedForRemoval.append(index)
    for index in markedForRemoval:
      newPdb.removeLine(newPdb.atomToRaw[index])
    return newPdb

  def getAllResidues(self):
    '''returns 2 lists, one of numbers and one of names, for all unique ress'''
    nums, names = [],[] #cumulators
    curNum, curName = None, None
    for index, thisResNum in enumerate(self.resNums):
      thisResName = self.resNames[index]
      if curNum is None or curNum != thisResNum or curName != thisResName:
        curNum = thisResNum
        curName = thisResName
        nums.append(curNum)
        names.append(curName)
    return nums, names
 
  def getAllResidueNumbers(self):
    '''returns a list of all residue numbers'''
    return self.getAllResidues()[0]

  def getAllResidueNames(self):
    '''returns a list of all residue names'''
    return self.getAllResidues()[1]

  def getListResidues(self, residueNumbers):
    '''copies the pdb, removes all residues not in the list, returns new pdb'''
    newPdb = self.copy()
    markedForRemoval = []
    for index, thisResNum in enumerate(newPdb.resNums):
      if thisResNum not in residueNumbers:
        markedForRemoval.append(index)
    for index in markedForRemoval:
      newPdb.removeLine(newPdb.atomToRaw[index])
    return newPdb

  def getListResiduesChains(self, residueChainNumbers):
    '''copies the pdb, removes all residues not in the list, returns new pdb'''
    newPdb = self.copy()
    markedForRemoval = []
    for index, thisResNum in enumerate(newPdb.resNums):
      chain = newPdb.chains[index]
      thisCheck = str(thisResNum) + str(chain)
      if thisCheck not in residueChainNumbers:
        markedForRemoval.append(index)
    for index in markedForRemoval:
      newPdb.removeLine(newPdb.atomToRaw[index])
    return newPdb

  def getFactorsByResidueChain(self, resNum, chain):
    '''gets all the factors (burial depths) for a given resnum + chain combo'''
    hashString = str(resNum) + str(chain)
    if hashString in self.__factorsByRC:
      if len(self.__factorsByRC[hashString]) > 0:
        return self.__factorsByRC[hashString]
    #otherwise compute it
    factorsRet = []
    for index, thisChain in enumerate(self.chains):
      if thisChain == chain:
        if resNum == self.resNums[index]:
          if self.radii[index] > 0.:
            factorsRet.append(self.factors[index][1]) 
    self.__factorsByRC[hashString] = factorsRet
    return factorsRet

  def getOneChain(self, chainId):
    '''copies the pdb, removes all but one chain, returns it'''
    newPdb = self.copy()
    markedForRemoval = []
    for index, thisChain in enumerate(newPdb.chains):
      if thisChain != chainId:
        markedForRemoval.append(index)
    for index in markedForRemoval:
      newPdb.removeLine(newPdb.atomToRaw[index])
    return newPdb

  def removeLine(self, rawDataIndex):
    self.rawData[rawDataIndex] = False

  def removeAtomMatching(self, index, character):
    '''for removing all zeta-atoms for instance (call with 1, 'Z')'''
    markedForRemoval = []
    for atomIndex, atom in enumerate(self.atoms):
      if atom[index] == character:
        markedForRemoval.append(self.atomToRaw[atomIndex])
    for rawIndex in markedForRemoval:
      self.removeLine(rawIndex)
    #no return necessary, as self has been modified

  def updateOneResidue(self, rawDataIndex, newNumber):
    '''updates the data and raw line for a residue number'''
    self.resNums[self.rawToAtom[rawDataIndex]] = newNumber
    newData =  self.rawData[rawDataIndex][:22]
    tempData = "%4d" % newNumber
    newData += tempData + self.rawData[rawDataIndex][26:] 
    self.rawData[rawDataIndex] = newData

  def renumberResidues(self):
    '''renumbers residues to be unique as the pdb biological unit doesn't 
    do this properly amazingly enough'''
    chainsHere ={}
    for chainId in self.chains:
      chainsHere[chainId] = True
    for chainId in chainsHere.keys():
      self.renumberResiduesOneChain(chainId) #does all work

  def renumberResiduesOneChain(self, chainId):
    '''renumbers residues on one chain'''
    seenResNums,lastRes,highestResNum = {},-1000, -1000
    for atomIndex,residueNumber in enumerate(self.resNums):
      if self.chains[atomIndex] == chainId: #otherwise skip
        if residueNumber > highestResNum:
          highestResNum = residueNumber
        if lastRes != residueNumber:
          if residueNumber in seenResNums: #renumber
            newNum = seenResNums[residueNumber]
            if newNum:
              newNum = highestResNum + 1
              highestResNum += 1
            seenResNums[residueNumber] = newNum
            rawDataIndex = self.atomToRaw[atomIndex]
            self.updateOneResidue(rawDataIndex, newNum)
          else:        
            seenResNums[residueNumber] = True #indicate we've seen it before
          lastRes = residueNumber #always do this
        elif lastRes == residueNumber:
          if seenResNums[residueNumber]: #nothing to do actually
            pass
          else:
            newNum = seenResNums[residueNumber]
            rawDataIndex = self.atomToRaw[atomIndex]
            self.updateOneResidue(rawDataIndex, newNum)

  def getAverageCoords(self, listPdbFileNames):
    '''takes current pdb and adds the data from the list, averages the coords
    and makes a new pdb with the output. absolutely no error checking is done
    to ensure the atoms are in the same order in the files'''
    newPdb = self.copy()
    otherData = []
    for pdbFileName in listPdbFileNames:
      otherData.append(pdbData(pdbFileName))
    for index, coord in enumerate(self.coords):    
      badData = 0
      sums = list(coord) #instead of tuple
      for data in otherData:
        try:
          for sumIndex in xrange(len(sums)):
            sums[sumIndex] += data.coords[index][sumIndex]
        except IndexError: #somehow the model doesn't have the same #atoms
          badData += 1
      for sumIndex in xrange(len(sums)):
        sums[sumIndex] /= (len(otherData) + 1 - badData) #badData keeps avgs ok 
      newPdb.updateNewXYZ(index, sums[0], sums[1], sums[2])
    return newPdb

  def calcRMSDfile(self, otherFileName, alphas=False):
    '''wrapper for calcRMSD'''
    otherData = pdbData(otherFileName)
    return self.calcRMSD(otherData, alphas=alphas)

  def calcRMSD(self, other, alphas=False):
    '''calculates rmsd of all atoms between self and other, no checking or 
    alignment is done'''
    squaredSum = 0.0
    badData = 0 #keeps track of missing atoms (some NMR models wrong)
    for index, coord in enumerate(self.coords):
      if not alphas or self.atoms[index][0:2] == "CA":
        try: 
          otherCoord = other.coords[index]
          squaredSum += geometry.distL2Squared(coord, otherCoord)    
        except IndexError:
          badData += 1
      else:
        badData += 1 #also keep track of non-carbon alphas when looking for CA
    #print squaredSum, (len(self.coords) - badData)
    squaredSum /= (len(self.coords) - badData)
    return squaredSum**0.5

  def getHeavyAtomXYZ(self):
    '''returns the xyz of the nonzero radius atoms'''
    if self.__nonZeroRadiiXYZ is None: #recalculate
      self.__nonZeroRadiiXYZ = []
      for index,radius in enumerate(self.radii):
        if radius > 0.:
          self.__nonZeroRadiiXYZ.append(self.coords[index])
    return self.__nonZeroRadiiXYZ

  def getHeavyAtomXYZRadius(self):
    '''returns the x,y,z,radius of the nonzero radius atoms'''
    if self.__nonZeroRadiiXYZr is None: #recalculate
      self.__nonZeroRadiiXYZr = []
      for index,radius in enumerate(self.radii):
        if radius > 0.:
          newList = list(self.coords[index])
          newList.append(radius)
          self.__nonZeroRadiiXYZr.append(newList)
    return self.__nonZeroRadiiXYZr

  def getHeavyAtomCount(self):
    '''counts and returns number of heavy atoms'''
    if self.__nonZeroRadiiCount == 0: #recalculate
      for radius in self.radii:
        if radius > 0.:
          self.__nonZeroRadiiCount += 1
    return self.__nonZeroRadiiCount      

  def getIndexByResidueAtom(self, resNum, resCode, atomName):
    '''gets the index matching the input data, returns false if no match'''
    atomNameStr = string.strip(atomName)
    for index in xrange(len(self.atoms)):
      if resNum == self.resNums[index]:
        if resCode == self.resNames[index]:
          if atomNameStr == string.strip(self.atoms[index]):
            return index
    return False #not found

  def getListResidueNumberChain(self):
    '''outputs a unique list of residue number and chain combos'''
    listR = []
    for index, resNum in enumerate(self.resNums):
      chain = self.chains[index]
      if 0 == len(listR) or (resNum, chain) != listR[-1]:
        listR.append((resNum, chain))
    return listR

  def assignCharges(self, charge):
    '''assigns charges from the charge object, assumes caller really wants em'''
    self.charges = [] #delete old charges if present
    self.hydroCharges = []
    for index, resName in enumerate(self.resNames):
      atomName = self.atoms[index]
      self.charges.append(charge.getCharge(atomName, resName))
      self.hydroCharges.append(charge.getTrinaryCharge(atomName, resName))

  def clusterAtoms(self, distanceCutoff=2.0):
    '''breaks into distinct unions of atoms based on distance cutoff'''
    ligandClusters = unionfind2.unionFind()
    cutoffSquared = distanceCutoff ** 2. #faster comparisons
    for index,coord in enumerate(self.coords):
      for index2,coord2 in enumerate(self.coords):
        if index2 > index: #only do comparisons once each
          distBetweenSquared = geometry.distL2Squared(coord, coord2)
          if distBetweenSquared <= cutoffSquared:
            ligandClusters.union(index, index2)
    clusteredLists = ligandClusters.toLists()
    newPdbs = [] #list of pdbData objects to return
    for oneCluster in clusteredLists:
      newPdb = self.copy()
      markedForRemoval = []
      for index in xrange(len(self.coords)):
        if index not in oneCluster:
          markedForRemoval.append(index)
      for index in markedForRemoval:
        newPdb.removeLine(newPdb.atomToRaw[index])
      newPdbs.append(newPdb)
    return newPdbs

  def getNearbyAtoms(self, pointList, nearbyDistance=0.):
    '''returns the list of line numbers of atoms '''  
    lines = {}
    if nearbyDistance > 0.: #actually do distance cutoff
      nearbyDistanceSquared = nearbyDistance ** 2.
      for pt in pointList:
        tempSet = set()
        for index,coord in enumerate(self.coords):
          distanceBetween = geometry.distL2Squared(pt, coord)
          if distanceBetween < nearbyDistanceSquared:   
            tempSet.update([self.atomToRaw[index]])
        tempList = list(tempSet)
        tempList.sort()
        lines[tuple(pt)] = tempList
    else: #just find closese atom for each pt
      for pt in pointList:
        bestDist, bestPt = 10000000., False
        for index,coord in enumerate(self.coords):
          distanceBetween = geometry.distL2Squared(pt, coord)
          if distanceBetween < bestDist:   
            bestPt = self.atomToRaw[index]
            bestDist = distanceBetween
        lines[tuple(pt)] = [bestPt] #still needs to be a list
    return lines

  def getResidueChainsFromNums(self, atomList):
    '''returns a list of residuenumber+chain strings from the list of atoms'''
    resChainList = []
    for index in atomList:
      residueNumber = self.resNums[index]
      chain = self.chains[index]
      resChain = str(residueNumber) + str(chain)
      if resChain not in resChainList: #guarantee uniqueness
        resChainList.append(resChain)
    resChainList.sort()
    return resChainList

  def getResidueNamesChainsFromNums(self, atomList):
    '''returns a list residuename+number+chain strings from the list of atoms'''
    resChainList = []
    for index in atomList:
      residueName = self.resNames[index]
      residueNumber = self.resNums[index]
      chain = self.chains[index]
      resChain = str(residueName) + str(residueNumber) + str(chain)
      if resChain not in resChainList: #guarantee uniqueness
        resChainList.append(resChain)
    resChainList.sort() #is this useful anymore? 
    return resChainList

  def getResidueNamesChains(self):
    '''returns a list residuename+number+chain strings from all atoms'''
    resChainList = []
    for index in xrange(len(self.resNames)):
      residueName = self.resNames[index]
      residueNumber = self.resNums[index]
      chain = self.chains[index]
      resChain = str(residueName) + str(residueNumber) + str(chain)
      if resChain not in resChainList: #guarantee uniqueness
        resChainList.append(resChain)
    resChainList.sort() #is this useful anymore? 
    return resChainList

  def getAtomsFromNums(self, atomList):
    '''returns a list of atom number+name+chain string from the list of atoms'''
    atomChainList = []
    for index in atomList:
      atomName = self.atoms[index]
      chain = self.chains[index]
      atomChain = str(atomName) + str(index) + str(chain)
      if atomChain not in atomChainList: #guarantee uniqueness
        atomChainList.append(atomChain)
    return atomChainList

  def getNearbyResidues(self, pointList, nearbyDistance=0.):
    '''returns a new pdb of residues in this pdb near the points given'''
    nearbyDistanceSquared = nearbyDistance ** 2.
    residuesNearPts = []
    for pt in pointList:
      for index,coord in enumerate(self.coords):
        distanceBetweenSquared = geometry.distL2Squared(pt, coord)
        if distanceBetweenSquared < nearbyDistanceSquared:   
          residueNumber = self.resNums[index]
          chain = self.chains[index]
          resChain = str(residueNumber) + str(chain)
          if resChain not in residuesNearPts: #guarantee uniqueness
            residuesNearPts.append(resChain)
    residuesNearPts.sort()
    return self.getListResiduesChains(residuesNearPts)

  def countChains(self):
    '''counts the chains'''
    chains = []
    for chainId in self.chains:
      if chainId not in chains:
        chains.append(chainId)
    return len(chains)

  def getFasta(self):
    '''returns a string in FASTA format, one letter per residue'''
    fasta = ""
    lastResChain = "0X"
    for index in xrange(len(self.resNums)):
      residueName = self.resNames[index]
      residueNumber = self.resNums[index]
      chain = self.chains[index]
      resChain = str(residueNumber) + str(chain)
      if resChain != lastResChain: #new residue
        try:
          position = aminoAcid3Codes.index(residueName) 
          fasta += aminoAcidCodes[position]
        except ValueError:
          fasta += "X" #unknown amino acid
      lastResChain = resChain
    return fasta 

if -1 != string.find(sys.argv[0], "pdb.py"):
  for filename in sys.argv[1:]:
    pdbD = pdbData(filename)
    print filename, pdbD.getHeavyAtomCount()
