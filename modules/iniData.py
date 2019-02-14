
import math
from numpy import *
import modules.var as v
import modules.read,modules.compute


def initialData():
  #Read input data
  #Read nodes coordinates and depth from "meshNodes.txt" file
  modules.read.meshNodes()

  #Read initial parameters from "globalParameters.txt" file
  modules.read.globalParameters()

  v.head=zeros((v.nNode),dtype=float64)
  v.tt=zeros((v.nElem),dtype=float64)
  v.q=zeros((v.nNode),dtype=float64)
  v.wDepthE=zeros((v.nElem),dtype=float64)
  v.typBC=zeros((v.nNode),int)
  v.wDepth=zeros((v.nNode),dtype=float64)
  v.nodElems=ones((v.nNode,7),int)
  v.nodElems*=(-1)
  v.readConClastBC=zeros((v.nNode,v.nClsMat),dtype=float64)
  v.conClastBC=zeros((v.nNode,v.nClsMat),dtype=float64)
  v.veloNode=zeros((v.nNode),dtype=float64)
  v.readClastQ=zeros((v.nNode,v.nClsMat),dtype=float64)
  v.clastQ=zeros((v.nNode,v.nClsMat),dtype=float64)
  v.conClast=zeros((v.nNode,v.nClsMat),dtype=float64)
  v.grid=zeros((v.nNode,v.nti,v.nTotalMat),dtype=float64)
  v.ve3=zeros((v.nElem),dtype=float64)
  v.velDi3=zeros((v.nElem),dtype=float64)
  v.area=zeros((v.nElem),dtype=float64)
  v.inn=zeros((v.nElem,3),int)
  v.poroIni=zeros((v.nTotalMat),dtype=float64)
  v.vCritic=zeros((v.nClsMat),dtype=float64)
  v.setClast=zeros((v.nClsMat),dtype=float64)
  v.diff=zeros((v.nClsMat),dtype=float64)
  v.dispL=zeros((v.nClsMat),dtype=float64)
  v.dispT=zeros((v.nClsMat),dtype=float64)
  v.wDisp=zeros((v.nClsMat),dtype=float64)
  v.vCriticE=zeros((v.nTotalMat),dtype=float64)
  v.setClastE=zeros((v.nTotalMat),dtype=float64)
  v.subsidence=zeros((v.nNode),dtype=float64)
  v.slope=zeros((v.nti,v.nElem),dtype=float64)
  v.slopeNode=zeros((v.nti,v.nNode),dtype=float64)
  v.interMatrix=zeros((v.nCboMat,v.nCboMat),dtype=float64)
  v.birth=zeros((v.nCboMat),dtype=float64)
  v.specNutConsum=zeros((v.nCboMat),dtype=float64)
  v.cboProd=zeros((v.nCboMat),dtype=float64)
  v.specPopIni=zeros((v.nCboMat),dtype=float64)
  v.specPop=zeros((v.nNode,v.nti,v.nCboMat),dtype=float64)
  v.specPopMin=zeros((v.nCboMat),dtype=float64)
  v.nutriCom=zeros((v.nCboMat),dtype=float64)
  v.ss=zeros((v.nTotalMat),dtype=float64)
  v.porosity=zeros((v.nNode,v.nti,v.nTotalMat),dtype=float64)
  v.poroMin=zeros((v.nTotalMat),dtype=float64)
  v.totalVolSedSuspJtiMat=zeros((v.nTotalMat,v.nti))
  v.dist=zeros((v.nNode,v.nNode), dtype=float64)
  v.connectNodes=zeros((v.nNode,6),int)
  v.depo=zeros((v.nNode,v.nTotalMat),dtype=float64)
  #v.numEq=(v.nCboMat*3)+4
  v.times=ones((v.nti),dtype=float64)
  v.surface=zeros((v.nNode),dtype=float64)
  v.surf=zeros((v.nti,v.nNode),dtype=float64)
  v.basem=zeros((v.nNode),dtype=float64)
  v.prevTotTime=0.0
  v.deptFactorX=[None for i in range(v.nCboMat)]
  v.deptFactorY=[None for i in range(v.nCboMat)]
  v.slopeFactorX=[None for i in range(v.nCboMat)]
  v.slopeFactorY=[None for i in range(v.nCboMat)]
  v.flowFactorX=[None for i in range(v.nCboMat)]
  v.flowFactorY=[None for i in range(v.nCboMat)]
  v.clsFactorX=[[None for i in range(v.nClsMat)]for j in range(v.nCboMat)]
  v.clsFactorY=[[None for i in range(v.nClsMat)]for j in range(v.nCboMat)]
  v.nutrFactorX=[None for i in range(v.nCboMat)]
  v.nutrFactorY=[None for i in range(v.nCboMat)]
  v.factorWD=ones((v.nNode,v.nti,v.nCboMat),dtype=float64)
  v.factorVL=ones((v.nNode,v.nti,v.nCboMat),dtype=float64)
  v.factorSL=ones((v.nNode,v.nti,v.nCboMat),dtype=float64)
  v.factorNU=ones((v.nNode,v.nti,v.nCboMat),dtype=float64)
  v.factorCL=ones((v.nNode,v.nti,v.nCboMat,v.nClsMat),dtype=float64)
  v.envFactor=zeros((v.nNode,v.nti,v.nCboMat),dtype=float64)
  
  v.sediName=[None for i in range(v.nTotalMat)]
  
  for i in range(v.nNode):
    v.surface[i]=float(v.readBasem[i])
    v.basem[i]=float(v.readBasem[i])
    
  for i in range(v.nti):
    v.times[i]=round(v.modelingTime/float(v.nti),1)
  
  
  for j in range(v.nRow-1):
    for i in range(v.nCol-1):
      nElLo=(((v.nCol-1)*2)*(j)+(i+1))-1
      nElUp=(2*(v.nCol-1)*(j)+(v.nCol-1)+(i+1))-1
      v.inn[nElLo][0]=((i+1)+(j)*v.nCol)-1
      v.inn[nElLo][1]=((i+1)+1+(j)*v.nCol)-1
      v.inn[nElLo][2]=((i+1)+(j+1)*v.nCol)-1
      v.inn[nElUp][0]=((i+1)+(j+1)*v.nCol)-1
      v.inn[nElUp][2]=((i+1)+1+(j+1)*v.nCol)-1
      v.inn[nElUp][1]=(((i+1)+1)+(j)*v.nCol)-1
  for i in range(v.nNode):
    ncount=0
    for j in range(v.nElem):
      if(v.inn[j][0]==i):
        v.nodElems[i][ncount]=j
        ncount=ncount+1
      if(v.inn[j][1]==i):
        v.nodElems[i][ncount]=j
        ncount=ncount+1
      if(v.inn[j][2]==i):
        v.nodElems[i][ncount]=j
        ncount=ncount+1      
    v.nodElems[i][6]=ncount
  v.nl=0
  for l in range(v.nElem):
    n1=v.inn[l][2]
    for iv in range(3):
      n2=v.inn[l][iv]
      id=n1-n2
      if(id>v.nl): v.nl=id
      n1=n2
  v.bw=v.nl*2+1
  if(v.bw>(2*v.nCol+1)): print ("bandwidth too large: ",v.bw); quit()
  
  
# Calculate and store the distance among nodes (used in several subroutines)
  for i in range(v.nNode):
    for j in range(v.nNode):
      v.dist[i][j]=math.sqrt((v.xCo[i]-v.xCo[j])**2+(v.yCo[i]-v.yCo[j])**2)

  for i in range(v.nNode):
    cnTemp=zeros((18),int)
    n=0
    for j in range(6):
      for m in range(3):
        if(v.nodElems[i][j]>0): cnTemp[n]=v.inn[v.nodElems[i][j]][m]
        n=n+1
    for j in range(18):
      if (cnTemp[j]==i): cnTemp[j]=0
      for  k in range(j-1):
        if (cnTemp[j]==cnTemp[k]): cnTemp[j]=0
    n=0  
    for j in range(18):
      if (cnTemp[j]>0): 
        v.connectNodes[i][n]=cnTemp[j]
        n=n+1
  modules.compute.areaElem()
  modules.read.sediParam()
  if v.comSubsidence: read.subsidence()


# initialise head values to the maximum fixed head value
  headMax=-10000000.0

# find maximum value of all fixed head nodes
  for node in range(v.nNode):
    if (v.readtypBC[node]==1 and v.headBC[node]>headMax): headMax=v.headBC[node]

# set all nodes, which are not fixed head nodes to headmax
  for node in range(v.nNode):
    if (v.readtypBC[node]!=1): v.headBC[node]=headMax

  print ('Initial data dimensioned and initialized')
