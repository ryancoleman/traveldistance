#script to read .pdb and make pdb-pymol objects
#ryan g. coleman ryangc@mail.med.upenn.edu

from pymol.cgo import *
from pymol import cmd
import math,  string

class displayObj(object):
  '''contains all methods and local variables used to display things'''
  pdbDatas = [] #global arrays basically
  pdbNames = []
  drawCount = 0
  cylRad = 0.3

  class pdbData(object):
    '''stores data for a pdb file consisting of atoms'''

    def __init__(self,filename=False):
      '''really really dumb that pymol can't import things sanely'''
      self.__nonZeroRadiiCount = 0
      self.__nonZeroRadiiXYZ = None
      self.__nonZeroRadiiXYZr = None
      self.__factorsByRC = {}
      self.rawData = []
      self.coords = []
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

  def pdbOpen(self,filename):
    '''this function is the opening to pymol to open a pdb file'''
    if len(self.pdbDatas) > 0:
      self.pdbDatas[0] = self.pdbData(filename)
      self.pdbNames[0] = filename[:4]
    else:
      self.pdbDatas.append(self.pdbData(filename))
      self.pdbNames.append(filename[:4])

  #follows is all color stuff
  def decideColor(self, value, min, med, max, gradient,cycle=False):
    #print value,min,med,max,gradient #useful debugging routine
    value = value[0]
    if cycle: #new case, make more like topo map...
      #ignore mediumValueColor
      #use the colors in the gradient as whole colors, switch based on mod value
      whichCol = int((value-min)/4) % (len(gradient)/3)
      red = gradient[int(whichCol*3+0)]
      gre = gradient[int(whichCol*3+1)]
      blu = gradient[int(whichCol*3+2)]
      return [red,gre,blu]
    elif len(gradient) == 9:
      if gradient[0]==gradient[3] and gradient[3]==gradient[6] and \
            gradient[1]==gradient[4] and gradient[4]==gradient[7] and \
            gradient[2]==gradient[5] and gradient[5]==gradient[8]: #all the same
        return [gradient[0],gradient[1],gradient[2]] #just output first
      if min > value:
        return [gradient[0],gradient[1],gradient[2]]
      elif min <= value and value < med:
        alpha = (value-min)/(med-min)
        red = (1.-alpha)*gradient[0] + (alpha)*gradient[3]
        gre = (1.-alpha)*gradient[1] + (alpha)*gradient[4]
        blu = (1.-alpha)*gradient[2] + (alpha)*gradient[5]
        return [red,gre,blu]
      elif med <= value and value <= max:
        alpha = (value-med)/(max-med)
        red = (1.-alpha)*gradient[3] + (alpha)*gradient[6]
        gre = (1.-alpha)*gradient[4] + (alpha)*gradient[7]
        blu = (1.-alpha)*gradient[5] + (alpha)*gradient[8]
        return [red,gre,blu]
      elif value > max:
        return [gradient[6],gradient[7],gradient[8]]
      else:
        print "something funny with colors"
        return [0.0,0.0,0.0]
    elif len(gradient) == 24:
      valueOct = int(value)-1
      #print valueOct
      red = gradient[valueOct*3]
      gre = gradient[valueOct*3+1]
      blu = gradient[valueOct*3+2]
      return [red,gre,blu]

  colorGrads = {'flat':[1.0,1.0,1.0, 1.0,1.0,1.0, 1.0,1.0,1.0],
                'black':[0.0,0.0,0.0, 0.0,0.0,0.0, 0.0,0.0,0.0],
                'mono':[1.0,1.0,1.0, 0.5,0.5,0.5, 0.0,0.0,0.0],
                'monorev':[0.0,0.0,0.0, 0.5,0.5,0.5, 1.0,1.0,1.0],
                'gwg':[0.5,0.5,0.5, 1.0,1.0,1.0, 0.0,1.0,0.0],
                'rwb':[1.0,0.0,0.0, 1.0,1.0,1.0, 0.0,0.0,1.0],
                'bwr':[0.0,0.0,1.0, 1.0,1.0,1.0, 1.0,0.0,0.0],
                'rgb':[1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0],
                'bgr':[0.0,0.0,1.0, 0.0,1.0,0.0, 1.0,0.0,0.0],
                'red':[1.0,0.0,0.0, 1.0,0.0,0.0, 1.0,0.0,0.0],
                'green':[0.0,1.0,0.0, 0.0,1.0,0.0, 0.0,1.0,0.0],
                'blue':[0.0,0.0,1.0, 0.0,0.0,1.0, 0.0,0.0,1.0],
                'yellow':[1.0,1.0,0.0, 1.0,1.0,0.0, 1.0,1.0,0.0],
                'magenta':[1.0,0.0,1.0, 1.0,0.0,1.0, 1.0,0.0,1.0],
                'cyan':[0.0,1.0,1.0, 0.0,1.0,1.0, 0.0,1.0,1.0],
                'octant':[1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0, 0.5,0.5,0.0,
                          0.5,0.0,0.5, 0.0,0.5,0.5, 1.0,1.0,1.0, 0.1,0.1,0.1]}

  def pdbMinMaxB(self):
    pdbData = self.pdbDatas[-1]
    pdbName = self.pdbNames[-1]
    minBvalue = None
    maxBvalue = None
    for factors in pdbData.factors:
      bFactor = factors[1]
      if minBvalue is None or bFactor < minBvalue:
        minBvalue = bFactor
      if maxBvalue is None or bFactor > maxBvalue:
        maxBvalue = bFactor
    return minBvalue, maxBvalue

  def pdbDisplayCyls(self, minB=None, maxB=None, masklow=None, maskhigh=None, gradient='bgr'):
    '''draws small spheres, cylinders requires connectivity, maybe later'''
    pdbData = self.pdbDatas[-1]
    pdbName = self.pdbNames[-1]
    minBvalue = minB
    maxBvalue = maxB
    if minB is None or maxB is None:
      for factors in pdbData.factors:
        bFactor = factors[1]
        if minB is None:
          if minBvalue is None or bFactor < minBvalue:
            minBvalue = bFactor
        if maxB is None:
          if maxBvalue is None or bFactor > maxBvalue:
            maxBvalue = bFactor
    medBvalue = (minBvalue + maxBvalue)/2.
    if masklow is None:
      masklow = minBvalue
    else:
      masklow = float(masklow)
    if maskhigh is None:
      maskhigh = maxBvalue
    else:
      maskhigh = float(maskhigh)
    pdbObj = []
    for index, coords in enumerate(pdbData.coords):
      bFactor = float(pdbData.factors[index][1])
      if bFactor >= masklow and bFactor <= maskhigh:
        color = self.decideColor([bFactor], minBvalue, medBvalue, maxBvalue, self.colorGrads[gradient])
        pdbObj.append(COLOR)
        pdbObj.extend(color)
        pdbObj.append(SPHERE)
        pdbObj.extend(coords)
        pdbObj.append(self.cylRad)
      else:
        print bFactor, masklow, maskhigh
    cmd.load_cgo(pdbObj, pdbName +"."+ str(self.drawCount))
    self.drawCount += 1

  def pdbDisplay(self):
    '''creates an actual pdb object of the protein or just the part near the
    surface/group of interest'''
    pdbText = "\n".join(self.pdbDatas[-1].rawData)
    cmd.read_pdbstr(pdbText, self.pdbNames[-1]+"."+str(self.drawCount)+".pdb")
    self.drawCount += 1

#this is really really dumb but i can't get it to work otherwise
displayInstance = displayObj()

def pdbOpen(filename):
  displayInstance.pdbOpen(filename)

def pdbDisplay():
  displayInstance.pdbDisplay()

def pdbDisplayCyls(minB=None, maxB=None, masklow=None, maskhigh=None):
  displayInstance.pdbDisplayCyls(minB=minB,maxB=maxB, \
          masklow=masklow,maskhigh=maskhigh)

def pdbOpenBurial(filename, outName):
  displayInstance.pdbOpen(filename)
  minb,maxb = displayInstance.pdbMinMaxB()
  curHi = minb+0.01
  cmd.do("delete all")
  count = 0
  while curHi <= maxb:
    displayInstance.pdbDisplayCyls(masklow=minb, maskhigh=curHi)
    curHi += 0.25
    prefix = str(count)
    if count < 100:
      prefix = "0" + str(count)
    if count < 10:
      prefix = "00" + str(count)
    cmd.do("ray")
    cmd.do("png " + prefix + "." + outName + ".down.png")
    cmd.do("delete all")
    count += 1
  curLo = maxb-0.01
  cmd.do("delete all")
  count = 0
  while curLo >= minb:
    displayInstance.pdbDisplayCyls(masklow=curLo, maskhigh=maxb)
    curLo -= 0.25
    prefix = str(count)
    if count < 100:
      prefix = "0" + str(count)
    if count < 10:
      prefix = "00" + str(count)
    cmd.do("ray")
    cmd.do("png " + prefix + "." + outName + ".up.png")
    cmd.do("delete all")
    count += 1


cmd.extend("pdbOpen",pdbOpen)
cmd.extend("pdbOpenBurial",pdbOpenBurial)
cmd.extend("pdbDisplay",pdbDisplay)
cmd.extend("pdbDisplayCyls",pdbDisplayCyls)
