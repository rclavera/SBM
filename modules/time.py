# -*- coding: utf-8 -*-

import numpy as np
import modules.var as v


def dtcalc():

# It calculates time step for a dt time step

  timeCourant=np.copy(v.times[v.jti-1])
  dtClast=np.copy(v.times[v.jti-1])
  timeTemp=np.ones((v.nElem),dtype=np.float64)
  timeTemp*=1000000000.0
#     calculate maximum time steps to satisfy all stability criteria

#     maximum time step is duration of complete sedimentation period

#       check max timestep size from clastic sedimentation rate
#       if clastic sedims are considered
  if(v.comClasticSedi==True):
    dtClast=tstepsize(dtClast)

#     constrain timestep size from transport velocity (courant time)
#     only if transport is considered
  if(v.comTransport==True):
    for i in range(v.nElem):
      if v.wDepthE[i]>0.0:
        aix=[v.dist[v.inn[i][0]][v.inn[i][1]],v.dist[v.inn[i][0]][v.inn[i][2]],v.dist[v.inn[i][1]][v.inn[i][2]]]
        timeTemp[i]=min(aix)/v.ve3[i]

    timeCourant=np.amin(timeTemp)

#     if actual time + new time step size less than complete time step
#     then new time step size = complete time step - actual time
  aix=[dtClast,timeCourant,(v.times[v.jti-1]-v.cTime)]

  v.dt=min(aix)
  if v.cTime+v.dt>v.stepTime:
    v.dt=v.times[v.jti-1]-v.cTime
  if (v.dt<=0.0):
    print ("There was problems computing time. Program has been stopped",v.times, v.cTime)
    quit()    


def tstepsize(timest):

#       calculate maximum timestep from settling velocity
#       initialise timestep
#       loop for all nodes
  timestep=timest
  settletime=np.ones((v.nNode),dtype=np.float64)
  settletime*=10000000.0
  factor=np.ones((v.nClsMat),dtype=np.float64)
  for node in range(v.nNode):
#       calculate timestepsize only for cells above sealevel
    if ((v.wDepth[node]-v.head[node])>=v.wMin):
#       loop for all materials
      for i in range(v.nClsMat):
        if all(v.conClast[:][i])>0.0:
          factor[i]=((v.vCritic[i]-(v.veloNode[node]/365.25))/v.vCritic[i])
          if factor[i]>1.0: factor[i]=1.0
          if factor[i]<=0.0: factor[i]=0.00000001
          factor[i]=v.setClast[i]*365.25*factor[i]

      settletime[node]=v.wDepth[node]/np.amin(factor)
      if settletime[node]<=0.0: settletime[node]=1000000.0

  timestep=min(settletime)
  return timestep
