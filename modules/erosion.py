#!/usr/bin/env python
# -*- coding: utf-8 -*-


def erosion():

  for mat in range(v.nTotalMat):
    for node in range(v.nNode):
      clastMass=0.0
      for i in range(v.jti-1):
        clastMass=clastMass+v.grid[node][i][mat]

    if ((v.veloNode[node]/365.25)>v.vCriticE[mat] and v.wDepth[node]>0.0 and clastMass>0.0):
## calculate erosion velocity (comes in m/day)
      factor=((v.veloNode[node]/365.25)-v.vCriticE[mat]/(v.velonode[node]/365.25)
      if (factor>1): factor=1.0
      if (factor<0): factor=0.0
      sediPortion=setClastE[mat]*365.25*factor*v.dt/v.wDepth[node]
      if (sediPortion>1.0) sediPortion=1.0
      v.eros[node][mat]=clastMass*sediPortion

  for mat in range(v.nClsMat):
    for node in range(v.nNode):   
      v.conClast[node][mat]=v.conClast[node][mat]+(v.eros[node][mat]*v.density[mat]/v.wDepth[node])
      resta=v.eros[node][mat]
      jt=v.jti
      er=v.eros[node][mat]
      while (resta>0.0 or jt>0):
        if (er>v.grid[node][jt][mat]):
          er=er-v.grid[node][jt][mat]
          v.grid[node][jt][mat]=0.0
        else:
          v.grid[node][jt][mat]=v.grid[node][jt][mat]-er
          er=0.0
        jt-=1

  for mat in range(v.nCboMat):
    for node in range(v.nNode):   
      v.conClast[node][v.nSiliMat+mat]=v.conClast[node][v.nSiliMat+mat]+
                                     (v.eros[node][v.nClsMat+mat]*v.density[v.nClsMat+mat]/v.wDepth[node])
      resta=v.eros[node][v.nClsMat+mat]
      jt=v.jti
      er=v.eros[node][v.nClsMat+mat]
      while (resta>0.0 or jt>0):
        if (er>v.grid[node][jt][v.nClsMat+mat]):
          er=er-v.grid[node][jt][v.nClsMat+mat]
          v.grid[node][jt][v.nClsMat+mat]=0.0
        else:
          v.grid[node][jt][v.nClsMat+mat]=v.grid[node][jt][v.nClsMat+mat]-er
          er=0.0
        jt-=1  