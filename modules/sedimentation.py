# -*- coding: utf-8 -*-

import modules.var as v
import numpy as np
from scipy.integrate import ode
import matplotlib as plt

def update():

# calculate water depth as difference between hydraulic head and surface
  depoN=np.zeros((v.nNode,v.nClsMat),dtype=np.float64)

  for node in range (v.nNode):
    totalDepo=0.0
    aix=v.wDepth[node]-v.wMin
    for mat in range (v.nTotalMat):
      totalDepo=totalDepo+v.depo[node][mat]
    if (totalDepo<=0.0 or aix<=0.0):
      continue
    elif (totalDepo>aix):
      x=aix/totalDepo
      for mat in range (v.nTotalMat):
        if (mat<v.nClsMat):
          depoN[node][mat]=v.depo[node][mat]*x
        #s'ha d'actualitzar conClast
          retorn=v.depo[node][mat]-depoN[node][mat]
          v.conClast[node][mat]=v.conClast[node][mat]+retorn
          v.depo[node][mat]=depoN[node][mat]
        else:
          v.depo[node][mat]=v.depo[node][mat]*x

    for mat in range(v.nTotalMat):
# reduce surface
      v.surface[node]=v.surface[node]+v.depo[node][mat]
# update thickness of clastic sed deposited (grid array)
      v.grid[node][v.jti-1][mat]=v.grid[node][v.jti-1][mat]+v.depo[node][mat]

  v.depo=np.zeros((v.nNode,v.nTotalMat),dtype=np.float64)


def dn_dt(t,y,eF,iM,cboP,imp):
  ans=0.0
  eq=[]
  for i in range(v.nCboMat):
    for j in range(v.nCboMat):
      ans+=iM[i][j]*y[i]*y[j]
    eq.append(eF[i]*y[i]+ans)
  for i in range(v.nCboMat):
    if imp[i]!=0.0:
      eq.append(cboP[i]*y[i]/imp[i])
    else:
      eq.append(cboP[i]*y[i]*imp[i])
  return eq


def carbonate():
  ti =0.0
  tf =v.dt
  iM= v.interMatrix  
  #backend = 'vode'
  backend = 'dopri5'
  #backend = 'dop853'
  
# calculate growthrate of species as function of environmental parameters
  for node in range(v.nNode):
    if v.wDepth[node]>v.wMin:
      eF=np.zeros((v.nCboMat),dtype=np.float64)
      invmp=np.zeros((v.nCboMat),dtype=np.float64)
      for i in range (v.nCboMat):
        v.factorWD[node][v.jti-1][i]=np.interp(v.wDepth[node],v.deptFactorX[i][:],v.deptFactorY[i][:])
        v.factorVL[node][v.jti-1][i]=np.interp(v.veloNode[node]/365.25,v.flowFactorX[i][:],v.flowFactorY[i][:])
        v.factorSL[node][v.jti-1][i]=np.interp(v.slopeNode[v.jti-1][node],v.slopeFactorX[i][:],v.slopeFactorY[i][:])
        v.factorNU[node][v.jti-1][i]=1.0

        #factorNU=np.interp(v.nutriConc[node],v.nutrFactorX[i][:],v.nutrFactorY[i][:])
        for j in range(v.nClsMat):
          v.factorCL[node][v.jti-1][i][j]=np.interp(v.conClast[node][j],v.clsFactorX[i][j][:],v.clsFactorY[i][j][:])
        if v.factorType==2:
          fcl=np.amin(v.factorCL)
          v.envFactor[node][v.jti-1][i]=min([v.factorWD[node][v.jti-1][i],v.factorVL[node][v.jti-1][i],v.factorSL[node][v.jti-1][i],v.factorNU[node][v.jti-1][i],fcl])
        elif v.factorType==1:
          fcl=1.0
          for j in range(v.nClsMat):
            fcl=fcl*v.factorCL[node][i][j]
          v.envFactor[node][v.jti-1][i]=v.factorWD[node][v.jti-1][i]*v.factorVL[node][v.jti-1][i]*v.factorSL[node][v.jti-1][i]*v.factorNU[node][v.jti-1][i]*fcl
        else:
          v.envFactor[node][v.jti-1][i]=1.0

        eF[i]=v.birth[i]*v.envFactor[node][v.jti-1][i]

        if v.specPop[node][v.jti-1][i]<v.specPopMin[i]:
          v.specPop[node][v.jti-1][i]=v.specPopMin[i]

        ixx=0.0
        for j in range(v.nCboMat):
          if i!=j:
            ixx+=v.interMatrix[i][j]
        if ixx!=0.0 and eF[i]>0.0:
          ainv=np.linalg.inv(v.interMatrix)
          invmp=np.dot(ainv,-eF)
        elif ixx==0.0 and eF[i]>0.0:
          invmp[i]=-eF[i]/v.interMatrix[i][i]
        else: 
          invmp[i]=0.0
      if any(invmp)>0.0:
        y0=np.zeros((2*v.nCboMat),dtype=np.float64)
        for ix in range(v.nCboMat):
          y0[ix]=v.specPop[node][v.jti-1][ix]
          y0[v.nCboMat+ix]=0.0
        solver = ode(dn_dt)
        solver.set_integrator(backend)
        solver.set_initial_value(y0, ti)
        solver.set_f_params(eF,iM,v.cboProd,invmp)
        solt = []
        soly = []
        while solver.t < tf:
          solver.integrate(tf, step=True)
          solt.append(solver.t)
          soly.append(solver.y)
        #print soly
        solyp=np.asarray(soly)
        #if (solyp.shape[0])>1:
          #print solyp
          #print solt        
          #plt.plot(solt, solyp[:,0])
          #plt.plot(solt, solyp[:,1])
          #plt.show()
          #quit()
        for i in range(v.nCboMat):
          v.specPop[node][v.jti-1][i]=solyp[solyp.shape[0]-1][i]
          if solyp[solyp.shape[0]-1][v.nCboMat+i]>0.0:
            v.depo[node][v.nClsMat+i]=solyp[solyp.shape[0]-1][v.nCboMat+i]
      else:
        for i in range(v.nCboMat):
          v.specPop[node][v.jti-1][i]=0.0
    else:
      for i in range (v.nCboMat):
        v.specPop[node][v.jti-1][i]=0.0


def clastic():

# It calculates the clastic sedimentation
  for mat in range(v.nClsMat):
    for node in range(v.nNode):
#       check if velocity is below critical value for sedim type
#       velo comes in m/y,  vCritic im m/day, set velo to m/day

      if((v.veloNode[node]/365.25)<=v.vCritic[mat] and (v.wDepth[node]-v.head[node])>v.wMin and v.conClast[node][mat]>0.0):

#       clastmass is suspended mass (in m thickness of sediment)
        clastMass=(v.conClast[node][mat]*v.wDepth[node])/v.density[mat]
#       calculate settle velocity (comes in m/day)
        factor=(v.vCritic[mat]-(v.veloNode[node]/365.25))/v.vCritic[mat]
        if (factor>1): factor=1.0
        if (factor<0): factor=0.0
#       calculate settle distance for dt
        sediPortion=v.setClast[mat]*365.25*factor*v.dt/v.wDepth[node]
        if(sediPortion>1.0): sediPortion=1.0
#       calculate mass deposited in dt
        v.depo[node][mat]=clastMass*sediPortion
#       calculate new concentration after deposition
        clastMass=clastMass-v.depo[node][mat]
        v.conClast[node][mat]=(v.density[mat]*clastMass)/(v.wDepth[node]-v.depo[node][mat])
#        print 'depo[',node,'][',mat,']=',v.depo[node][mat]



def carbonate2():
  a=np.zeros((6),dtype=np.float64)
  b=np.zeros((6,5),dtype=np.float64)
  ch=np.zeros((6),dtype=np.float64)
  ct=np.zeros((6),dtype=np.float64)
  a[1]=0.222222222222222      #2.0/9.0
  a[2]=0.333333333333333      #1.0/3.0
  a[3]=0.75                   #3.0/4.0
  a[4]=1.0
  a[5]=0.833333333333334      #5.0/6.0
  b[1][0]=0.222222222222222    #2.0/9.0
  b[2][0]=0.083333333333333    #1.0/12.0
  b[2][1]=0.25                 #1.0/4.0
  b[3][0]=0.5390625            #69.0/128.0
  b[3][1]=-1.8984375           #-243.0/128.0
  b[3][2]=2.109375             #135.0/64.0
  b[4][0]=-1.4166666666666667  #-17.0/12.0
  b[4][1]=6.75                 #27.0/4.0
  b[4][2]=-5.4                 #-27.0/5.0
  b[4][3]=1.0666666666666667   #16.0/15.0
  b[5][0]=0.15046296296296297  #65.0/432.0
  b[5][1]=-0.3125              #-5.0/16.0
  b[5][2]=0.8125               #13.0/16.0
  b[5][3]=0.14814814814814814  #4.0/27.0
  b[5][4]=0.034722222222222224 #5.0/144.0
  ch[0]=0.10444444444444445   #47.0/450.0
  ch[1]=0.0
  ch[2]=0.48                  #12.0/25.0
  ch[3]=0.14222222222222222   #32.0/225.0
  ch[4]=0.03333333333333333   #1.0/30.0
  ch[5]=0.24                  #6.0/25.0
  ct[0]=-0.006666666666666667 #-1.0/150.0
  ct[1]=0.0
  ct[2]=0.03                  #3.0/100.0
  ct[3]=-0.21333333333333335  #-16.0/75.0
  ct[4]=-0.05                 #-1.0/20.0
  ct[5]=0.24                  #6.0/25.0

  for node in range(v.nNode):
    if v.wDepth[node]>v.wMin:
      neq=v.nCboMat*2
      y=np.zeros((neq),dtype=np.float64)
      ytemp=np.zeros((neq),dtype=np.float64)
      z=np.zeros((neq),dtype=np.float64)
      f=np.zeros((6,neq),dtype=np.float64)
      te=np.zeros((neq),dtype=np.float64)
      for i in range (v.nCboMat):
        if v.specPop[node][v.jti-1][i]>=v.specPopMin[i]:
          y[i]=v.specPop[node][v.jti-1][i]
        else:
          y[i]=v.specPopMin[i]
      teMax=10000.0
      hold=v.hIni
      h=v.hIni
      x=0.0 #x is time
      xMax=x+v.dt #xmax is end carbonate model time(dt)
      while(x<xMax):
        xTemp=x
        fun(y,z,node)
        if (v.wDepth[node]<v.wMin):
          z[:]=0.0
        #for i in range (neq):
        f[0][:]=z[:]
        ytemp=np.copy(y)
        teMax=10000.0
        while (teMax>v.tol):
          for k in range (1,6):
            x=xTemp+a[k]*h
            for n in range (neq):
              y[n]=np.copy(ytemp[n])
              k1=k-1
              for l in range (k1):
                y[n]=y[n]+h*b[k][l]*f[l][n]
            fun(y,z,node)
            for n in range(neq):
              f[k][n]=z[n]
          for n in range(neq):
            te[n]=0.0
            y[n]=ytemp[n]
            for k in range(6):
              y[n]=y[n]+h*ch[k]*f[k][n]
              te[n]=te[n]+h*ct[k]*f[k][n]
            te[n]=np.absolute(te[n])
  #       check if sealevel below surface
          if (v.wDepth[node]<v.wMin):
            z[:]=0.0
            y[:]=0.0

  #       keep population above lower limit
          for i in range(v.nCboMat):
            if y[i]<=0.0:
              y[i]=v.specPopMin[i]
          teMax=np.amax(te)
          hold=h
          h=0.9*h*(v.tol/teMax)**0.2
  #       set an upper limit t in case y"=0 t gets infinity
          if(h>xMax*0.5): h=v.hIni
          if(h>10000.0): h=10000.0
          if(h>(xMax-xTemp)): h=xMax-xTemp
          if(h<v.hIni):
            h=v.hIni
            teMax=v.tol*0.9
            hold=h
  #     step is done
        x=xTemp+hold
      v.hIni=hold
      for i in range(v.nCboMat):
        v.specPop[node][i]=y[i]
  #     Update depo due to Carbonate sedimentation
        v.depo[node][v.nClsMat+i]=y[v.nCboMat+i]
    else:
      for i in range(v.nCboMat):
        v.specPop[node][i]=0.0
        v.depo[node][v.nClsMat+i]=0.0
def fun(y,z,node):

# in case of emersion species die out and always keep population above lower limit
  for i in range(v.nCboMat):
    if v.wDepth[node]<=v.wMin or y[i]<v.specPopMin[i]:
      y[i]=v.specPopMin[i]

# change parameters through time
# calculate growthrate of species as function of environmental parameters

  eF=np.zeros((v.nCboMat),dtype=np.float64)
  invmp=np.zeros((v.nCboMat),dtype=np.float64)
  for i in range (v.nCboMat):
    v.factorWD[node][v.jti-1][i]=np.interp(v.wDepth[node],v.deptFactorX[i][:],v.deptFactorY[i][:])
    v.factorVL[node][v.jti-1][i]=np.interp(v.veloNode[node]/365.25,v.flowFactorX[i][:],v.flowFactorY[i][:])
    v.factorSL[node][v.jti-1][i]=np.interp(v.slopeNode[v.jti-1][node],v.slopeFactorX[i][:],v.slopeFactorY[i][:])
    v.factorNU[node][v.jti-1][i]=1.0

    #factorNU=np.interp(v.nutriConc[node],v.nutrFactorX[i][:],v.nutrFactorY[i][:])
    for j in range(v.nClsMat):
      v.factorCL[node][v.jti-1][i][j]=np.interp(v.conClast[node][j],v.clsFactorX[i][j][:],v.clsFactorY[i][j][:])
    if v.factorType==2:
      fcl=np.amin(v.factorCL)
      v.envFactor[node][v.jti-1][i]=np.amin([v.factorWD[node][v.jti-1][i],v.factorVL[node][v.jti-1][i],v.factorSL[node][v.jti-1][i],v.factorNU[node][v.jti-1][i],fcl])
    elif v.factorType==1:
      fcl=1.0
      for j in range(v.nClsMat):
        fcl=fcl*v.factorCL[node][i][j]
      v.envFactor[node][v.jti-1][i]=v.factorWD[node][v.jti-1][i]*v.factorVL[node][v.jti-1][i]*v.factorSL[node][v.jti-1][i]*v.factorNU[node][v.jti-1][i]*fcl
    else:
      v.envFactor[node][v.jti-1][i]=1.0

    eF[i]=v.birth[i]*v.envFactor[node][v.jti-1][i]
    ixx=0.0
    for j in range(v.nCboMat):
      if i!=j:
       ixx+=v.interMatrix[i][j]
    if ixx!=0.0 and eF[i]>0.0:
      ainv=np.linalg.inv(v.interMatrix)
      invmp=np.dot(ainv,-eF)
    elif ixx==0.0 and eF[i]>0.0:
      invmp[i]=-eF[i]/v.interMatrix[i][i]
    else: 
      invmp[i]=0.0
# differential equations (z = population species)
  ans=0.0
  for i in range (v.nCboMat):
    for j in range(v.nCboMat):
      ans+=v.interMatrix[i][j]*y[i]*y[j]
    z[i]=eF[i]*y[i]+ans
  for i in range (v.nCboMat):
    if invmp[i]!=0.0:
      z[v.nCboMat+i]=v.cboProd[i]*y[i]/invmp[i]
    else:
      z[v.nCboMat+i]=0.0

#       check if critical waterdepth is reached
  if (v.wDepth[node]<=v.wMin):
    z[:]=0.0
