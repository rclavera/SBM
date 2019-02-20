# -*- coding: utf-8 -*-


from numpy import *
import numpy as np
import math
import modules.var as v


def area_elem():

    # calculate area of elements and relate nodes to element numbers
    # calculate area of each element
    # loop for elements
    x = np.zeros(3, dtype=np.float64)
    y = np.zeros(3, dtype=np.float64)
    v.totalArea = 0.0
  
    for l in range(v.nElem):
        # calculate area of element
        x[0] = v.xCo[v.inn[l][0]]
        y[0] = v.yCo[v.inn[l][0]]
        x[1] = v.xCo[v.inn[l][1]]
        y[1] = v.yCo[v.inn[l][1]]
        x[2] = v.xCo[v.inn[l][2]]
        y[2] = v.yCo[v.inn[l][2]]
        v.area[l] = ((y[0]+y[1])*(x[0]-x[1])+(y[1]+y[2])*(x[1]-x[2])+(y[2]+y[0])*(x[2]-x[0]))*0.5
        v.area[l] = abs(v.area[l])
        v.totalArea = v.totalArea+v.area[l]


def isostasy():
    # calculates isostatic compensation after sediment load
    # AIRY model see example given by Eisbacher 1991, p. 219
  
    subiso = np.zeros(v.nNode, dtype=np.float64)
    for node in range(v.nNode):
        dummy = 0.0
        dumdens = 0.0
        dumpor = 0.0
    # dummy is the total amount deposited in jti time step
    # calculate sediment density
    # average density of all sedim. types
        for m in range(v.nTotalMat):
            for t in range(v.jti):
                dummy = dummy+v.grid[node][t][m]
                dumdens = dumdens+v.grid[node][t][m]*v.density[m]
                dumpor = dumpor+v.grid[node][t][m]*v.porosity[node][t][m]

        if dummy > 0.00005:
            dumdens = dumdens/dummy
            dumpor = dumpor/dummy
        # average with water
            densitySed = dumpor*v.density[v.nTotalMat+1]+(1-dumpor)*dumdens
            para = (v.density[v.nTotalMat+1]-v.density[v.nTotalMat])/(v.density[v.nTotalMat+1]-densitySed)
            subiso[node] = dummy-(dummy/para)
            newBasem = v.readBasem[node]-subiso[node]-v.subsidence[node]*v.totalTime/1000
            diffe = v.basem[node]-newBasem
            v.basem[node] = newBasem
            v.surface[node] = v.readBasem[node]+dummy
            for t in range(v.jti-1):
                v.surf[t][node] = v.surf[t][node]-diffe


def subsidence():

    # calculate subsidence for a defined timestep dt
    # subsidence comes in m/1000 y
    # downward movement is positive
    # upward movement is negative

    temps = v.dt/1000.0
    for node in range(v.nNode):
        for j in range(v.jti-1):
            v.surf[j][node] = v.surf[j][node]-(v.subsidence[node]*temps)

    v.surface = v.surface-v.subsidence*temps
    v.basem = v.basem-v.subsidence*temps


def slope(te):
    punt = np.zeros((3, 3), dtype=np.float64)
    u = np.zeros(3, dtype=np.float64)
    b = np.zeros(3, dtype=np.float64)
    n = np.zeros((v.nElem, 3), dtype=np.float64)
    medNormal = np.zeros(3, dtype=np.float64)
    # this calculates the slope of each element and node
    for i in range(v.nElem):
        # coordenates x,y,z of punt 1
        punt[0][0] = v.xCo[v.inn[i][0]]
        punt[0][1] = v.yCo[v.inn[i][0]]
        punt[0][2] = v.surface[v.inn[i][0]]
        # coordenates x,y,z of punt 2
        punt[1][0] = v.xCo[v.inn[i][1]]
        punt[1][1] = v.yCo[v.inn[i][1]]
        punt[1][2] = v.surface[v.inn[i][1]]
        # coordenates x,y,z of punt 3
        punt[2][0] = v.xCo[v.inn[i][2]]
        punt[2][1] = v.yCo[v.inn[i][2]]
        punt[2][2] = v.surface[v.inn[i][2]]
        # calcula el vector del 'punt 1' al 'punt 2'
        # (x2-x1,y2-y1,z2-z1)
        u[0] = punt[1][0]-punt[0][0]
        u[1] = punt[1][1]-punt[0][1]
        u[2] = punt[1][2]-punt[0][2]
        # calcula el vector del 'punt 1' al 'punt 3'
        # (x3-x1,y3-y1,z3-z1)
        b[0] = punt[2][0]-punt[0][0]
        b[1] = punt[2][1]-punt[0][1]
        b[2] = punt[2][2]-punt[0][2]
        # calcula el vector normal al pla format pels vectors 'u' i 'b'
        # n=u x b (producte escalar de dos vectors)
        n[i][0] = (u[1]*b[2])-(u[2]*b[1])
        n[i][1] = -((u[0]*b[2])-(u[2]*b[0]))
        n[i][2] = (u[0]*b[1])-(u[1]*b[0])
        if n[i][2] < 0.0:
            n[i][0] = -n[i][0]
            n[i][1] = -n[i][1]
            n[i][2] = -n[i][2]
        v.slope[te][i] = (math.acos(n[i][2]/(math.sqrt((n[i][0]**2.0)+(n[i][1]**2.0)+(n[i][2]**2.0))))) * \
                                                                                             57.29577951308232  # 180/PI

    for i in range(v.nNode):
        if v.nodElems[i][6] == 1:
            medNormal[0] = n[v.nodElems[i][0]][0]
            medNormal[1] = n[v.nodElems[i][0]][1]
            medNormal[2] = n[v.nodElems[i][0]][2]
        elif v.nodElems[i][6] == 2:
            medNormal[0] = (n[v.nodElems[i][0]][0]+n[v.nodElems[i][1]][0])/2
            medNormal[1] = (n[v.nodElems[i][0]][1]+n[v.nodElems[i][1]][1])/2
            medNormal[2] = (n[v.nodElems[i][0]][2]+n[v.nodElems[i][1]][2])/2
        elif v.nodElems[i][6] == 3:
            medNormal[0] = (n[v.nodElems[i][0]][0] +
                           n[v.nodElems[i][1]][0] +
                           n[v.nodElems[i][2]][0])/3
            medNormal[1] = (n[v.nodElems[i][0]][1] +
                           n[v.nodElems[i][1]][1] +
                           n[v.nodElems[i][2]][1])/3
            medNormal[2] = (n[v.nodElems[i][0]][2] +
                           n[v.nodElems[i][1]][2] +
                           n[v.nodElems[i][2]][2])/3
        elif v.nodElems[i][6] == 6:
            medNormal[0] = (n[v.nodElems[i][0]][0] +
                          n[v.nodElems[i][1]][0] +
                          n[v.nodElems[i][2]][0] +
                          n[v.nodElems[i][3]][0] +
                          n[v.nodElems[i][4]][0] +
                          n[v.nodElems[i][5]][0])/6
            medNormal[1] = (n[v.nodElems[i][0]][1] +
                          n[v.nodElems[i][1]][1] +
                          n[v.nodElems[i][2]][1] +
                          n[v.nodElems[i][3]][1] +
                          n[v.nodElems[i][4]][1] +
                          n[v.nodElems[i][5]][1])/6
            medNormal[2] = (n[v.nodElems[i][0]][2] +
                          n[v.nodElems[i][1]][2] +
                          n[v.nodElems[i][2]][2] +
                          n[v.nodElems[i][3]][2] +
                          n[v.nodElems[i][4]][2] +
                          n[v.nodElems[i][5]][2])/6

        v.slopeNode[te][i] = (math.acos(medNormal[2]/(math.sqrt((medNormal[0]**2.0)+(medNormal[1]**2.0) +
                                                                (medNormal[2]**2.0)))))*57.29577951308232  # 180/PI


def sea_level():
    v.seaLevel = np.interp(v.totalTime, v.sl[0][:], v.sl[1][:])


def waderdepth():
    v.wDepth[:] = v.head[:]+v.seaLevel-v.surface[:]
    mask = v.wDepth < 0.0
    v.wDepth[mask] = v.wDepth[mask]*(-0.0)
    v.wDepthE = np.zeros(v.nElem, dtype=float64)
    for i in range(v.nElem):
        w = np.zeros(3, float)
        w[0] = v.wDepth[v.inn[i][0]]
        w[1] = v.wDepth[v.inn[i][1]]
        w[2] = v.wDepth[v.inn[i][2]]
        if min(w) > 0.0:
            v.wDepthE[i] = (w[0]+w[1]+w[2])/3.0


def boundary_conditions():
    import random
    # It checks the boun#  print 'cTime:',v.cTime,'dt:',v.dt, 'totalTime:',v.totalTime, 'stepTime',v.stepTimedary
    # conditions to move them when sealevel changes

    nCLE = np.zeros(v.nElem, int)
    proxNod = np.zeros(2, int)
    nLAN = np.zeros(v.nElem, int)
    dista = np.zeros(v.nNode, dtype=np.float64)
    distaT = np.zeros(v.nNode, dtype=np.float64)
    headd = np.zeros(v.nNode, dtype=np.float64)

    # Add sealevel change to initial head boundary condition
    headd = np.copy(v.head)
    # for i in range(v.nNode):
    # if v.typBC[i]==0:
    # v.head[i]=v.headBC[i]
    mask = v.typBC == 0
    v.head[mask] = v.headBC[mask]
    # Copy boundary conditions to local for move boundary

    v.typBC = np.copy(v.readtypBC)
    v.q = np.copy(v.readQ)
    v.clastQ = np.copy(v.readClastQ)
    v.conClastBC = np.copy(v.readConClastBC)

    someElem = 0
    dista = 10000000000.0
    proxNod = 0
    someLand = 0

    # look for coast line and land elements

    for el in range(v.nElem):
        wa = np.zeros(3, float)
        wa[0] = v.wDepth[v.inn[el][0]]-v.head[v.inn[el][0]]
        wa[1] = v.wDepth[v.inn[el][1]]-v.head[v.inn[el][1]]
        wa[2] = v.wDepth[v.inn[el][2]]-v.head[v.inn[el][2]]
        if max(wa) > 0.0 and min(wa) <= 0.0:
            nCLE[el] = el
            someElem += 1

        if (wa[0] <= 0.0) and (wa[1] <= 0.0) and (wa[2] <= 0.0):
            nLAN[el] = el
            someLand += 1

    if someElem > 0:
        for node in range(v.nNode):
            # if a node is over sea level and have BC look for the nearest nodes in the coastline
            if v.wDepth[node] <= 0.0:
                tempc = 0
                for j in range(v.nClsMat):
                    if v.clastQ[node][j] > 0.0:
                        tempc = tempc+1
                if v.q[node] > 0.0 or v.typBC[node] < 1 or tempc > 0:
                    dista = np.zeros(v.nNode, float)
                    # search nodes near input node over sea level
                    for el in range(v.nElem):
                        if nCLE[el] > 0:
                            for i in range(3):
                                if (v.wDepth[v.inn[el][i]]-v.head[v.inn[el][i]] > 0.0) and (v.dist[node][v.inn[el][i]]) > 0.0:
                                    dista[v.inn[el][i]] = v.dist[node][v.inn[el][i]]
                    for i in range(v.nNode):
                        if dista[i] > 0.0:
                            canvi = 0
                            for ii in range(6):
                                if v.connectNodes[i][ii] == 0.0:
                                    canvi = canvi+1
                                else:
                                    if v.wDepth[v.connectNodes[i][ii]]-v.head[v.connectNodes[i][ii]] <= 0.0:
                                        canvi = canvi+1
                            if canvi >= 3:
                                dista[i] = 0.0
                    for i in range(v.nNode):
                        if dista[i] <= 0.0:
                            dista[i] = 1000000000.0
                    aixxx1 = np.argmin(dista)
                    distaT = dista
                    distaT[aixxx1] = 1000000000.0
                    aixxx2 = np.argmin(distaT)

                    if dista[aixxx1] == distaT[aixxx2]:
                        rannum = random.random()
                        if rannum<0.5:
                            aixxx1 = aixxx2

                    v.q[aixxx1] = v.q[node]+v.q[aixxx1]
                    v.q[node] = 0.0
                    v.typBC[aixxx1] = v.typBC[node]
                    for j in range(v.nClsMat):
                        v.clastQ[aixxx1][j] = v.clastQ[node][j]+v.clastQ[aixxx1][j]
                        v.conClastBC[aixxx1][j] = v.conClastBC[node][j]+v.conClastBC[aixxx1][j]
                        v.clastQ[node][j] = 0.0
                        v.conClastBC[node][j] = 0.0

    if someLand > 0:
        for el in range(v.nElem):
            if nLAN[el] > 0:
                v.head[v.inn[el][0]] = 0.0
                v.head[v.inn[el][1]] = 0.0
                v.head[v.inn[el][2]] = 0.0
                for mat in range(v.nClsMat):
                    v.conClast[v.inn[el][0]][mat] = 0.0
                    v.conClast[v.inn[el][1]][mat] = 0.0
                    v.conClast[v.inn[el][2]][mat] = 0.0


def mass_balance(material, cc):

    # function to correct mass balance due to hyperbolic dominance in
    # transport equation from advective transport term and eulerian approach
    # in the transport algorithm (non moving grid)

  j = np.ones(v.nElem, dtype=np.float64)
  j *= -1.0
  elemTag = np.zeros(v.nElem, int)

  multiplicationFactor = 1
  additionFactor = 0
  j = -1.0
  newMassTotalCc = 0.0
  oldMassTotalConc = 0.0
  sedInFlow = 0.0
  sediOutFlow = 0.0
  sediOutFlow_2 = 0.0
  volTotal = 0.0

#  inflow mass
  for node in range(v.nNode):
    sedInFlow = sedInFlow+(v.clastQ[node][material]*v.dt)

  for elem in range(v.nElem):
    # mass in suspension before transport
    oldMassTotalConc = oldMassTotalConc+((v.conClast[v.inn[elem][0]][material]+
                     v.conClast[v.inn[elem][1]][material]+v.conClast[v.inn[elem][2]][material])/
                     3.0)*(v.area[elem]*v.wDepthE[elem])
  
#      new mass in suspension after trasnport
    newMassTotalCc = newMassTotalCc+((cc[v.inn[elem][0]]+cc[v.inn[elem][1]]+cc[v.inn[elem][2]])/
                                    3.0)*(v.area[elem]*v.wDepthE[elem])

#       Calculation of outflow mass in a dt	 

  for node in range(v.nNode):
    if v.typBC[node] == 1:
      for numElems in range(v.nodElems[node][6]):
        if elemTag[v.nodElems[node][numElems]] != 1:
          for k in range(3):
            alter = v.inn[v.nodElems[node][numElems]][k]
            if alter != node:
              if v.typBC[alter] == 1:
                elemTag[v.nodElems[node][numElems]] = 1
                x1 = v.xCo[node]
                x2 = v.xCo[alter]
                y1 = v.yCo[node]
                y2 = v.yCo[alter]
                z1 = v.wDepth[node]
                z2 = v.wDepth[alter]
                vx = x2-x1
                vy = y2-y1
                beta = math.atan2(vy,vx)*180.0/3.1415927
                dist12 = math.sqrt(vx*vx+vy*vy)
                gamma = v.velDi3[v.nodElems[node][numElems]]-beta
                gamma = gamma*3.1415927/180.0
                r1 = math.sin(gamma)
                dista = r1*dist12
                if dista < 0.0:
                    dista = -dista
                crossection = dista*(z1+z2)/2.0
                qq = crossection*v.ve3[v.nodElems[node][numElems]]
                sediOutFlow = sediOutFlow+qq*v.dt*((cc[node]+cc[alter])/2.0)

  if sediOutFlow < 0.0:
        sediOutFlow = -sediOutFlow
  
  mBalance = sedInFlow-sediOutFlow-newMassTotalCc+oldMassTotalConc

  if newMassTotalCc == 0.0 and sediOutFlow == 0.0:
      factor = 0.0
  else:
      factor = (sedInFlow+oldMassTotalConc)/(newMassTotalCc+sediOutFlow)

  cc *= factor

# Tornem a calcular la massa total del sistema per comprovar que el balan\8D de massa ha funcionat per un dt. Aquest calcul es nom\8Es de control, es pot el.liminar
  newMassTotalCc = 0.0;
  for elem in range(v.nElem):
#       calculated mass in suspension
    newMassTotalCc = newMassTotalCc+((cc[v.inn[elem][0]]+cc[v.inn[elem][1]]+cc[v.inn[elem][2]])/
                                     3.0)*(v.area[elem]*v.wDepthE[elem])
 
# We have to re-calculated the new mass that flows out of the model 
# Calculation of outflow mass	 
  elemTag = np.zeros(v.nElem, int)
  sediOutFlow = 0.0
  for node in range(v.nNode):
    if v.typBC[node] == 1:
      for numElems in range(v.nodElems[node][6]):
        if elemTag[v.nodElems[node][numElems]] != 1:
          for k in range(3):
            alter = v.inn[v.nodElems[node][numElems]][k]
            if alter != node:
              if v.typBC[alter] == 1:
                elemTag[v.nodElems[node][numElems]] = 1
                x1 = v.xCo[node]
                x2 = v.xCo[alter]
                y1 = v.yCo[node]
                y2 = v.yCo[alter]
                z1 = v.wDepth[node]
                z2 = v.wDepth[alter]
                vx = x2-x1
                vy = y2-y1
                beta = math.atan2(vy,vx)*180.0/3.1415927
                dist12 = math.sqrt(vx*vx+vy*vy)
                gamma = v.velDi3[v.nodElems[node][numElems]]-beta
                gamma = gamma*3.1415927/180.0
                r1 = math.sin(gamma)
                dista = r1*dist12
                if dista < 0.0:
                    dista = -dista
                crossection = dista*(z1+z2)/2.0
                qq = crossection*v.ve3[v.nodElems[node][numElems]]
                sediOutFlow = sediOutFlow+qq*v.dt*((cc[node]+cc[alter])/2.0)

# New loop to calculate in a diferent way the material that flows out of the model
# The material that can flow out of the model is the material located in the element that has nodes with ouflow of sediment

  sediOutFlow_2 = 0.0;
  for elem in range(v.nElem):
    if (((v.typBC[v.inn[elem][0]] == 1)and(v.typBC[v.inn[elem][1]] == 1))or\
        ((v.typBC[v.inn[elem][0]] == 1)and(v.typBC[v.inn[elem][2]] == 1))or\
        ((v.typBC[v.inn[elem][1]] == 1)and(v.typBC[v.inn[elem][2]] == 1))):

      sediOutFlow_2 = sediOutFlow_2+((cc[v.inn[elem][0]]+cc[v.inn[elem][1]]+cc[v.inn[elem][2]])/3.0)*(v.area[elem]*v.wDepthE[elem])

  if sediOutFlow < 0.0:
        sediOutFlow = -sediOutFlow
  v.totalVolSedOut = v.totalVolSedOut+sediOutFlow/v.density[material]
  v.totalVolSedOut_2 = v.totalVolSedOut_2+sediOutFlow_2/v.density[material]
