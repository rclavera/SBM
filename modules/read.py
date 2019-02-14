# -*- coding: utf-8 -*-

import modules.var as v
import numpy as np


def meshNodes():
    try:
        fc = open('meshNodes.txt', 'r')
        fc.readline()
        fc.readline()
        line = fc.readline()
        words = str.split(line)
        v.nCol = int(words[0])
        v.nRow = int(words[1])
        v.nNode = v.nCol*v.nRow
        v.nElem = (v.nCol-1)*(v.nRow-1)*2
        v.readBasem = np.zeros(v.nNode, dtype=np.float64)
        v.xCo = np.zeros(v.nNode, dtype=np.float64)
        v.yCo = np.zeros(v.nNode, dtype=np.float64)
        v.headBC = np.zeros(v.nNode, dtype=np.float64)
        v.readQ = np.zeros(v.nNode, dtype=np.float64)
        v.readtypBC = np.zeros(v.nNode, int)
        fc.readline()
        for nod in range(v.nNode):
            line = fc.readline()
            words = str.split(line)
            if len(words) > 3:
                v.xCo[nod] = (float(words[0]))
                v.yCo[nod] = (float(words[1]))
                v.readBasem[nod] = float(words[2])
                v.readtypBC[nod] = int(words[3])
                if v.readtypBC[nod] == 0: v.headBC[nod] = float(words[4])
                if v.readtypBC[nod] == 3: v.readQ[nod] = float(words[4])*31557600.0
            else:
                print('There are problems reading nodes file, or nodes file format is not correct')
                fc.close()
                quit()
        print('mesh data read')
        fc.close()
    except IOError:
        print('There are problems opening nodes.txt file, or nodes.txt file does not exist')
        quit()


def globalParameters():
    try:
        fco = open('globalParameters.txt', 'r')
        fco.readline()
        fco.readline()
        line = fco.readline()
        words = str.split(line)
        if len(words) > 1:
            if words[1] == 'y':
                v.modelingTime = float(words[0])
            elif words[1] == 'ky':
                v.modelingTime = float(words[0])*1000.0
            elif words[1] == 'My':
                v.modelingTime = float(words[0])*1000000.0
            else:
                print('There are problems reading total modelling time', words[0], words[1])
                fco.close()
                quit()
        fco.readline()
        line = fco.readline()
        words = str.split(line)
        if len(words) > 2:
            if words[1] == 't':
                v.stepTime = float(words[0])
                print(v.modelingTime, v.stepTime)
                if v.modelingTime%v.stepTime == 0.0:
                    v.nti = int(v.modelingTime/v.stepTime)
                else:
                    v.nti = int(v.modelingTime/v.stepTime)+1
            elif words[1] == 'n':
                v.nti = int(words[0])
                v.stepTime = v.modelingTime/float(v.nti)
            else:
                print('There are problems reading steps to save data', words[0], words[1])
                fco.close()
                quit()
        fco.readline()
        line = fco.readline()
        words = str.split(line)
        if len(words) > 1:
            if words[1] == 'y': v.comTransport = True
        line = fco.readline()
        words = str.split(line)
        if len(words) > 1:
            if words[1] == 'y': v.comMassBalance = True
        line=fco.readline()
        words = str.split(line)
        if len(words) > 1:
            if words[1] == 'y': v.comClasticSedi = True
        line = fco.readline()
        words = str.split(line)
        if len(words) > 1:
            if words[1] == 'y': v.comCarboSedi = True
        line = fco.readline()
        words = str.split(line)
        if len(words) > 1:
            if words[1] == 'y': v.comErosion = True
        line = fco.readline()
        words = str.split(line)
        if len(words) > 1:
            if words[1] == 'y': v.comIsostasy = True
        line=fco.readline()
        words=str.split(line)
        if len(words)>1:
            if words[1] == 'y': v.comCompaction = True
        line=fco.readline()
        words=str.split(line)
        if len(words)>1:
            if words[1] == 'y': v.comSubsidence = True
        fco.readline()
        line = fco.readline()
        words = str.split(line)
        if len(words) > 1:
            v.nSiliMat = int(words[0])
            v.nCboMat = int(words[1])
        if v.nCboMat > 0:
            v.nClsMat = v.nSiliMat+v.nCboMat
            v.nTotalMat = v.nClsMat+v.nCboMat
        else:
            v.nClsMat = v.nSiliMat
            v.nTotalMat = v.nClsMat
        v.density = np.zeros((v.nTotalMat+2), dtype=np.float64)
        fco.readline()
        fco.readline()
        line = fco.readline()
        words = str.split(line)
        v.cn = float(words[0])
        v.tInitStep = float(words[1])
        fco.readline()
        line = fco.readline()
        words = str.split(line)
        v.wMin = float(words[0])
        fco.readline()
        line = fco.readline()
        words = str.split(line)
        v.density[v.nTotalMat] = float(words[0])
        line = fco.readline()
        words = str.split(line)
        v.density[v.nTotalMat+1] = float(words[0])
        fco.readline()
        fco.readline()
        line = fco.readline()
        words = str.split(line)
        v.waySL = int(words[0])
        v.sl = np.zeros((2, v.waySL), dtype=np.float64)
        fco.readline()
        for i in range(v.waySL):
            line = fco.readline()
            words = str.split(line)
            v.sl[0][i] = float(words[0])
            v.sl[1][i] = float(words[1])
        print('Global parameters read')
        fco.close()
    except IOError:
        print('There are problems opening globalParameters.txt file, or globalParameters.txt file do not exist')
        quit()


def sediParam():
  for i in range(v.nClsMat):
    sediClstFile='clastic'+str(i+1)+'.txt'
    fsp=open(sediClstFile,'r')
    fsp.readline()
    line=fsp.readline()
    words=str.split(line)
    v.sediName[i]=(words[0])
    fsp.readline()
    line=fsp.readline()
    words=str.split(line)
    v.setClast[i]=float(words[0])
    v.vCritic[i]=float(words[1])
    fsp.readline()
    line=fsp.readline()
    words=str.split(line)
    v.vCriticE[i]=float(words[1])
    v.setClastE[i]=float(words[0])
    fsp.readline()      
    line=fsp.readline()
    words=str.split(line)
    v.density[i]=float(words[0])    
    fsp.readline()
    line=fsp.readline()
    words=str.split(line)
    v.ss[i]=float(words[0]) 
    v.poroIni[i]=float(words[1]) 
    v.poroMin[i]=float(words[2])
    for k in range(v.nNode):
      for j in range(v.nti):
        v.porosity[k][j][i]=v.poroIni[i]
    fsp.readline()  
    line=fsp.readline()
    words=str.split(line)
    v.dispL[i]=float(words[0])   
    v.dispT[i]=float(words[1]) 
    v.diff[i]=float(words[2])*31557600.0 #rescale dimension of diff in m**2/y (diff comes in m**2/sec)
    v.wDisp[i]=float(words[3])
    fsp.readline()
    line=fsp.readline()
    words=str.split(line)
    op=int(words[0])
    fsp.readline()
    for ops in range(op):
      line=fsp.readline()
      words=str.split(line)
      node=int(words[0])
      v.readClastQ[node][i]=float(words[1])*31557600.0
    print ('Siliciclastic sediment',v.sediName[i],'read')
    fsp.close()
  print ('Siliciclastic sediments read')


  for i in range (v.nCboMat):
    sediCboFile='carbonate'+str(i+1)+'.txt'
    fsp=open(sediCboFile,'r')
    fsp.readline()
    line=fsp.readline()
    words=str.split(line)
    v.sediName[v.nClsMat+i]=(words[0])    
    fsp.readline()
    line=fsp.readline()
    words=str.split(line)
    v.tol=float(words[0])
    v.chicken=float(words[1])
    v.hIni=float(words[2])
    fsp.readline()
    line=fsp.readline()
    words=str.split(line)
    v.specPopIni[i]=float(words[0])
    for j in range (v.nNode):
      v.specPop[j][i]=v.specPopIni[i]
    v.specPopMin[i]=float(words[1]) 
    v.birth[i]=float(words[2])
    v.cboProd[i]=float(words[3])
    v.nutriCom[i]=float(words[4])
    fsp.readline()
    fsp.readline() 
    line=fsp.readline()
    words=str.split(line)
    for j in range(v.nCboMat): 
      v.interMatrix[i][j]=float(words[j])
    fsp.readline()
    line=fsp.readline()
    words=str.split(line)
    v.vCriticE[v.nClsMat+i]=float(words[1])
    v.setClastE[v.nClsMat+i]=float(words[0])
    fsp.readline()      
    line=fsp.readline()
    words=str.split(line)
    v.density[v.nClsMat+i]=float(words[0])    
    fsp.readline()
    line=fsp.readline()
    words=str.split(line)
    v.ss[v.nClsMat+i-1]=float(words[0]) 
    v.poroIni[v.nClsMat+i]=float(words[1]) 
    v.poroMin[v.nClsMat+i]=float(words[2])
    for k in range(v.nNode):
      for j in range(v.nti):
        v.porosity[k][j][v.nClsMat+i]=v.poroIni[v.nClsMat+i]
    fsp.readline()
    fsp.readline()
    line=fsp.readline()
    words=str.split(line)
    v.factorType=int(words[0])
    if (v.factorType>0):
      fsp.readline()
      fsp.readline()
      line=fsp.readline()
      words=str.split(line)
      num=int(words[0])
      v.deptFactorX[i]=np.zeros((num),dtype=np.float64)
      v.deptFactorY[i]=np.zeros((num),dtype=np.float64)
      fsp.readline()
      for j in range(num):
        line=fsp.readline()
        words=str.split(line)
        v.deptFactorX[i][j]=float(words[0])
        v.deptFactorY[i][j]=float(words[1])
      fsp.readline()
      fsp.readline()
      line=fsp.readline()
      words=str.split(line)
      num=int(words[0])
      v.slopeFactorX[i]=np.zeros((num),dtype=np.float64)
      v.slopeFactorY[i]=np.zeros((num),dtype=np.float64)
      fsp.readline()
      for j in range(num):
        line=fsp.readline()
        words=str.split(line)
        v.slopeFactorX[i][j]=float(words[0])
        v.slopeFactorY[i][j]=float(words[1])
      fsp.readline()
      fsp.readline()
      line=fsp.readline()
      words=str.split(line)
      num=int(words[0])
      v.flowFactorX[i]=np.zeros((num),dtype=np.float64)
      v.flowFactorY[i]=np.zeros((num),dtype=np.float64)
      fsp.readline()
      for j in range(num):
        line=fsp.readline()
        words=str.split(line)
        v.flowFactorX[i][j]=float(words[0])
        v.flowFactorY[i][j]=float(words[1])
      for k in range(v.nClsMat):
        fsp.readline()
        fsp.readline()
        line=fsp.readline()
        words=str.split(line)
        num=int(words[0])
        v.clsFactorX[i][k]=np.zeros((num),dtype=np.float64)
        v.clsFactorY[i][k]=np.zeros((num),dtype=np.float64)
        fsp.readline()
        for j in range(num):
          line=fsp.readline()
          words=str.split(line)
          v.clsFactorX[i][k][j]=float(words[0])
          v.clsFactorY[i][k][j]=float(words[1])
      fsp.readline()
      fsp.readline()
      line=fsp.readline()
      words=str.split(line)
      num=int(words[0])
      v.nutrFactorX[i]=np.zeros((num),dtype=np.float64)
      v.nutrFactorY[i]=np.zeros((num),dtype=np.float64)
      fsp.readline()
      for j in range(num):
        line=fsp.readline()
        words=str.split(line)
        v.nutrFactorX[i][j]=float(words[0])
        v.nutrFactorY[i][j]=float(words[1])
    print ('Carbonate sediment ',v.sediName[v.nClsMat+i],'read')
    fsp.close()
  print ('Carbonate sediments read')
      
  print ("sedimentation parameters read")



def subsidence():
  sFile='subsidence.txt'
  fsp=open(sFile,'r')

#  read subsidence data
  fsp.readline()
  line=fsp.readline()
  words=str.split(line)
  subsiType=int(words[0])
  fsp.readline()

  if subsiType==0: #different subsidence for each node
    for node in range(v.nNode):
      line=fsp.readline()
      words=str.split(line)
      cn.subsidence[int(words[0])]=float(words[1])
  elif subsiType==1:  # same subsidence for all nodes
    subsidence=float(words[0])
  else:    #case default
    cn.subsidence=0.0

  fsp.close()
  print ("subsidence parameters read")

