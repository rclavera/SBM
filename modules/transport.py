# -*- coding: utf-8 -*-

import modules.var as v
import numpy as np
import math
#import matsolflow
import modules.compute as com


def transport():
#       subroutine for transport simulation in combination
#       with consolidation driven fluid flow
#       general program description see Kinzelbach 1986: 282 ff.
#       translation from BASIC with major changes in node and
#       element creation. version 6.1 from 2.11.2000

#       calculate dispersion tensor
  for material in range(v.nClsMat):
    xd=np.zeros((v.nElem),dtype=np.float64)
    yd=np.zeros((v.nElem),dtype=np.float64)
    xy=np.zeros((v.nElem),dtype=np.float64)
    cc=np.zeros((v.nNode),dtype=np.float64)
    qm=np.zeros((v.nNode),dtype=np.float64)
    itn=np.zeros((v.nNode),int)

    disptensor(xd,yd,xy,material)

#       define corresponding boundary conditions to transport
    setboundary(qm,material)

    #for node in range (v.nNode):
      #cc[node]=v.conClast[node][material]
    cc=np.copy(v.conClast[:,material])

    if any(qm)>0.0 or any(cc)>0.0:
      matsol(qm,cc,xd,yd,xy)

    for node in range(v.nNode):
      if cc[node]<0.0 or  v.wDepth[node]<=0.0: cc[node]=0.0

#     Calculate mass balance
    if (v.comMassBalance): com.massBalance(material,cc)

    for node in range (v.nNode):
      if v.wDepth[node]<=0.0: cc[node]=0.0
      v.conClast[node][material]=cc[node]

def disptensor(xd,yd,xy,material):

# It calculates the dispersion tensor of each element
  alphl=v.dispL[material]
  alpht=v.dispT[material]
  dumdif=v.diff[material]


  for n1 in range(v.nElem):
#       calculate dispersion tensor over all elements
#       vx ynd vy are velocity components in x and y in m/year

#       correct dispersivity and diffusivity for depth in order to include
#       increased mixing from wave (diffusion) and tidal flow in shallow areas
#       (dispersion) if wdepth of element is less than wdisp (matnum)
    if v.wDepthE[n1]<v.wDisp[material]:
      dumdif=v.diff[material]*(v.wDisp[material]-v.wDepthE[n1])
    um=v.ve3[n1]
    if(um==0.0): um=0.00000000001
#       calculate dispersivity from velocity and dip.coefficient
    diri=v.velDi3[n1]*0.017453292519943295
    vx=math.sin(diri)*um
    vy=math.cos(diri)*um
    dispx=(alphl*vx*vx+alpht*vy*vy)/um
    dispy=(alphl*vy*vy+alpht*vx*vx)/um
    xd[n1]=dumdif+dispx
    yd[n1]=dumdif+dispy
    xy[n1]=(alphl-alpht)*vx*vy/um

#       change xd and yd tensor if element has no flow node
  for n1 in range (v.nElem):
    if v.wDepthE[n1]<=0.0:
      xd[n1]=0.0
      yd[n1]=0.0
      xy[n1]=0.0
    
    
def setboundary(qm,material):

# It sets boundary conditions to transport calculation

#       assign boundary conditions to transport calculation of material
#       set all boundaries and qm to 0
# set bc identifyer, fixed concentration nodes and source rate nodes for material


  for i in range(v.nNode):
    if v.conClast[i][material]<0.0: v.conClast[i][material]=0.0
    if v.clastQ[i][material]>0.0 and v.wDepth[i]>0: qm[i]=v.clastQ[i][material]/v.wDepth[i]


def matsol(qm,cc,xd,yd,xy):
#       solution of system of equations from fe mesh
#       adapted from Kinzelbach 1985 p. 290

  bb=np.zeros((v.nNode),dtype=np.float64)
  a=np.zeros((v.nNode,v.bw+1),dtype=np.float64)
  b=np.zeros((3),dtype=np.float64)
  c=np.zeros((3),dtype=np.float64)
  re=np.ones((3,3),dtype=np.float64)
  pe=np.zeros((3,3),dtype=np.float64)
  ue=np.zeros((3,3),dtype=np.float64)

  rd=1.0/v.dt
  bb=np.copy(qm)

  #for i in range(v.nNode):
    #if itn[i]==1:
      #a[i][v.nl]=1.0
      #bb[i]=cc[i]

#       matrix assembly loop over elements
  for l in range(v.nElem):
    ux=v.ve3[l]*math.sin(v.velDi3[l]*0.017453292519943295)
    uy=v.ve3[l]*math.cos(v.velDi3[l]*0.017453292519943295)
    b[0]=v.yCo[v.inn[l][1]]-v.yCo[v.inn[l][2]]
    b[1]=v.yCo[v.inn[l][2]]-v.yCo[v.inn[l][0]]
    b[2]=v.yCo[v.inn[l][0]]-v.yCo[v.inn[l][1]]
    c[0]=v.xCo[v.inn[l][2]]-v.xCo[v.inn[l][1]]
    c[1]=v.xCo[v.inn[l][0]]-v.xCo[v.inn[l][2]]
    c[2]=v.xCo[v.inn[l][1]]-v.xCo[v.inn[l][0]]
    de=(b[0]*c[1]-b[1]*c[0])*0.5
    pe[0][0]=(b[0]*b[0]*xd[l]+c[0]*c[0]*yd[l]+(b[0]*c[0]+c[0]*b[0])*xy[l])/de*0.25
    pe[0][1]=(b[0]*b[1]*xd[l]+c[0]*c[1]*yd[l]+(b[0]*c[1]+c[0]*b[1])*xy[l])/de*0.25
    pe[0][2]=(b[0]*b[2]*xd[l]+c[0]*c[2]*yd[l]+(b[0]*c[2]+c[0]*b[2])*xy[l])/de*0.25
    pe[1][0]=(b[1]*b[0]*xd[l]+c[1]*c[0]*yd[l]+(b[1]*c[0]+c[1]*b[0])*xy[l])/de*0.25
    pe[1][1]=(b[1]*b[1]*xd[l]+c[1]*c[1]*yd[l]+(b[1]*c[1]+c[1]*b[1])*xy[l])/de*0.25
    pe[1][2]=(b[1]*b[2]*xd[l]+c[1]*c[2]*yd[l]+(b[1]*c[2]+c[1]*b[2])*xy[l])/de*0.25
    pe[2][0]=(b[2]*b[0]*xd[l]+c[2]*c[0]*yd[l]+(b[2]*c[0]+c[2]*b[0])*xy[l])/de*0.25
    pe[2][1]=(b[2]*b[1]*xd[l]+c[2]*c[1]*yd[l]+(b[2]*c[1]+c[2]*b[1])*xy[l])/de*0.25
    pe[2][2]=(b[2]*b[2]*xd[l]+c[2]*c[2]*yd[l]+(b[2]*c[2]+c[2]*b[2])*xy[l])/de*0.25
    ue[0][0]=(b[0]*ux+c[0]*uy)/6.0
    ue[0][1]=(b[1]*ux+c[1]*uy)/6.0
    ue[0][2]=(b[2]*ux+c[2]*uy)/6.0
    ue[1][0]=(b[0]*ux+c[0]*uy)/6.0
    ue[1][1]=(b[1]*ux+c[1]*uy)/6.0
    ue[1][2]=(b[2]*ux+c[2]*uy)/6.0
    ue[2][0]=(b[0]*ux+c[0]*uy)/6.0
    ue[2][1]=(b[1]*ux+c[1]*uy)/6.0
    ue[2][2]=(b[2]*ux+c[2]*uy)/6.0

    for i in range (3):
      for j in range (3):
        re[i][j]=de/12.0
    re[0][0]=de/6.0
    re[1][1]=de/6.0
    re[2][2]=de/6.0
# global matrix assembly
    ii=v.inn[l][0]
    bb[ii]=bb[ii]-(pe[0][0]+ue[0][0])*cc[v.inn[l][0]]
    jj=v.inn[l][0]-ii+v.nl
    if (jj>-1):
      a[ii][jj]=a[ii][jj]+(pe[0][0]+ue[0][0])*v.cn+ re[0][0]*rd
    bb[ii]=bb[ii]-(pe[0][1]+ue[0][1])*cc[v.inn[l][1]]
    jj=v.inn[l][1]-ii+v.nl
    if (jj>-1):
      a[ii][jj]=a[ii][jj]+(pe[0][1]+ue[0][1])*v.cn+ re[0][1]*rd
    bb[ii]=bb[ii]-(pe[0][2]+ue[0][2])*cc[v.inn[l][2]]
    jj=v.inn[l][2]-ii+v.nl
    if (jj>-1): 
      a[ii][jj]=a[ii][jj]+(pe[0][2]+ue[0][2])*v.cn+ re[0][2]*rd

    ii=v.inn[l][1]
    bb[ii]=bb[ii]-(pe[1][0]+ue[1][0])*cc[v.inn[l][0]]
    jj=v.inn[l][0]-ii+v.nl
    if (jj>-1):
      a[ii][jj]=a[ii][jj]+(pe[1][0]+ue[1][0])*v.cn+ re[1][0]*rd
    bb[ii]=bb[ii]-(pe[1][1]+ue[1][1])*cc[v.inn[l][1]]
    jj=v.inn[l][1]-ii+v.nl
    if (jj>-1):
      a[ii][jj]=a[ii][jj]+(pe[1][1]+ue[1][1])*v.cn+ re[1][1]*rd
    bb[ii]=bb[ii]-(pe[1][2]+ue[1][2])*cc[v.inn[l][2]]
    jj=v.inn[l][2]-ii+v.nl
    if (jj>-1): 
      a[ii][jj]=a[ii][jj]+(pe[1][2]+ue[1][2])*v.cn+ re[1][2]*rd

    ii=v.inn[l][2]
    bb[ii]=bb[ii]-(pe[2][0]+ue[2][0])*cc[v.inn[l][0]]
    jj=v.inn[l][0]-ii+v.nl
    if (jj>-1):
      a[ii][jj]=a[ii][jj]+(pe[2][0]+ue[2][0])*v.cn+ re[0][0]*rd
    bb[ii]=bb[ii]-(pe[2][1]+ue[2][1])*cc[v.inn[l][1]]
    jj=v.inn[l][1]-ii+v.nl
    if (jj>-1):
      a[ii][jj]=a[ii][jj]+(pe[2][1]+ue[2][1])*v.cn+ re[0][0]*rd
    bb[ii]=bb[ii]-(pe[2][2]+ue[2][2])*cc[v.inn[l][2]]
    jj=v.inn[l][2]-ii+v.nl
    if (jj>-1): 
      a[ii][jj]=a[ii][jj]+(pe[2][2]+ue[2][2])*v.cn+ re[0][0]*rd
  
# Solution of equation system
  kd=v.nl+1
# Factor matrix
  for ni in range(1,v.nNode):
    pi=-a[ni-1][kd-1]
    if (pi!=0.0):
      nj=ni+1
      ib=kd
      nk=ni+v.nl
      if (nk>v.nNode): nk=v.nNode
      for ll in range(nj,nk+1):
        ib=ib-1
        h=a[ll-1][ib-1]
        if (h!=0.0):
          h=h/pi
          a[ll-1][ib-1]=h
          jb=ib+1
          kb=ib+nk-ni
          lb=kd-ib
          for mb in range(jb,kb+1):
            a[ll-1][mb-1]=a[ll-1][mb-1]+(h*a[ni-1][lb+mb-1])

# Solution of equation system
# solution using factored matrix
  for ni in range(2,v.nNode+1):
    ko=ni-kd
    ib=1
    if (ko<=0): ib=1-ko
    nj=ib+ko
    su=0.0
    for jb in range(ib,v.nl+1):
      su=su+a[ni-1][jb-1]*bb[nj-1]
      nj=nj+1
    bb[ni-1]=bb[ni-1]+su


# back substitution

  bb[v.nNode-1]=bb[v.nNode-1]/a[v.nNode-1][kd-1]
  ll=kd+1
  ni=v.nNode
  for ib in range(2,v.nNode+1):
    ni=ni-1
    nj=ni
    mb=v.bw
    if (ib<=v.nl): mb=v.nl+ib
    su=0.0
    for jb in range(ll,mb+1):
      nj=nj+1
      su=su+a[ni-1][jb-1]*bb[nj-1]
    bb[ni-1]=(bb[ni-1]-su)/a[ni-1][kd-1]



# Equation solver finished
# Update concentrations
  for i in range(v.nNode):
    cc[i]=cc[i]+bb[i]
