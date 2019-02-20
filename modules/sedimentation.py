# -*- coding: utf-8 -*-

import modules.var as v
import numpy as np
from scipy.integrate import ode


def update():

    # calculate water depth as difference between hydraulic head and surface
    depoN = np.zeros((v.nNode, v.nClsMat), dtype=np.float64)

    for node in range(v.nNode):
        totalDepo = 0.0
        aix = v.wDepth[node]-v.wMin
        for mat in range(v.nTotalMat):
            totalDepo = totalDepo+v.depo[node][mat]
        if totalDepo <= 0.0 or aix <= 0.0:
            continue
        elif totalDepo > aix:
            x = aix/totalDepo
            for mat in range(v.nTotalMat):
                if mat < v.nClsMat:
                    depoN[node][mat] = v.depo[node][mat]*x
                    # s'ha d'actualitzar conClast
                    retorn = v.depo[node][mat]-depoN[node][mat]
                    v.conClast[node][mat] = v.conClast[node][mat]+retorn
                    v.depo[node][mat] = depoN[node][mat]
                else:
                    v.depo[node][mat] = v.depo[node][mat]*x

        for mat in range(v.nTotalMat):
            # reduce surface
            v.surface[node] = v.surface[node]+v.depo[node][mat]
            # update thickness of clastic sed deposited (grid array)
            v.grid[node][v.jti-1][mat] = v.grid[node][v.jti-1][mat]+v.depo[node][mat]

    v.depo = np.zeros((v.nNode, v.nTotalMat), dtype=np.float64)


def dn_dt(t, y, eF, iM, cboP, imp):
    ans = 0.0
    eq = []
    for i in range(v.nCboMat):
        for j in range(v.nCboMat):
            ans += iM[i][j]*y[i]*y[j]
        eq.append(eF[i]*y[i]+ans)
    for i in range(v.nCboMat):
        if imp[i] != 0.0:
            eq.append(cboP[i]*y[i]/imp[i])
        else:
            eq.append(cboP[i]*y[i]*imp[i])
    return eq


def carbonate():
    ti = 0.0
    tf = v.dt
    iM = v.interMatrix
    # backend = 'vode'
    backend = 'dopri5'
    # backend = 'dop853'

    # calculate growthrate of species as function of environmental parameters
    for node in range(v.nNode):
        if v.wDepth[node] > v.wMin:
            eF = np.zeros(v.nCboMat, dtype=np.float64)
            invmp = np.zeros(v.nCboMat, dtype=np.float64)
            for i in range(v.nCboMat):
                v.factorWD[node][v.jti-1][i] = np.interp(v.wDepth[node], v.deptFactorX[i][:], v.deptFactorY[i][:])
                v.factorVL[node][v.jti-1][i] = np.interp(v.veloNode[node]/365.25, v.flowFactorX[i][:], v.flowFactorY[i][:])
                v.factorSL[node][v.jti-1][i] = np.interp(v.slopeNode[v.jti-1][node], v.slopeFactorX[i][:],
                                                         v.slopeFactorY[i][:])
                v.factorNU[node][v.jti-1][i] = 1.0

                # factorNU=np.interp(v.nutriConc[node],v.nutrFactorX[i][:],v.nutrFactorY[i][:])
                for j in range(v.nClsMat):
                    v.factorCL[node][v.jti-1][i][j] = np.interp(v.conClast[node][j], v.clsFactorX[i][j][:],
                                                                v.clsFactorY[i][j][:])
                if v.factorType == 2:
                    fcl = np.amin(v.factorCL)
                    v.envFactor[node][v.jti-1][i] = min([v.factorWD[node][v.jti-1][i], v.factorVL[node][v.jti-1][i],
                                                         v.factorSL[node][v.jti-1][i], v.factorNU[node][v.jti-1][i],
                                                         fcl])
                elif v.factorType == 1:
                    fcl = 1.0
                    for j in range(v.nClsMat):
                        fcl = fcl*v.factorCL[node][i][j]
                    v.envFactor[node][v.jti-1][i] = v.factorWD[node][v.jti-1][i]*v.factorVL[node][v.jti-1][i]\
                                                    * v.factorSL[node][v.jti-1][i]*v.factorNU[node][v.jti-1][i]*fcl
                else:
                    v.envFactor[node][v.jti-1][i] = 1.0

                eF[i] = v.birth[i]*v.envFactor[node][v.jti-1][i]

                if v.specPop[node][v.jti-1][i] < v.specPopMin[i]:
                    v.specPop[node][v.jti-1][i] = v.specPopMin[i]

                ixx = 0.0
                for j in range(v.nCboMat):
                    if i != j:
                        ixx += v.interMatrix[i][j]
                if ixx != 0.0 and eF[i] > 0.0:
                    ainv = np.linalg.inv(v.interMatrix)
                    invmp = np.dot(ainv, -eF)
                elif ixx == 0.0 and eF[i] > 0.0:
                    invmp[i] = -eF[i]/v.interMatrix[i][i]
                else:
                    invmp[i] = 0.0
            if any(invmp) > 0.0:
                y0 = np.zeros((2*v.nCboMat), dtype=np.float64)
                for ix in range(v.nCboMat):
                    y0[ix] = v.specPop[node][v.jti-1][ix]
                    y0[v.nCboMat+ix] = 0.0
                solver = ode(dn_dt)
                solver.set_integrator(backend)
                solver.set_initial_value(y0, ti)
                solver.set_f_params(eF, iM, v.cboProd, invmp)
                solt = []
                soly = []
                while solver.t < tf:
                    solver.integrate(tf, step=True)
                    solt.append(solver.t)
                    soly.append(solver.y)
                # print soly
                solyp = np.asarray(soly)
                # if (solyp.shape[0])>1:
                #     print solyp
                #     print solt
                #     plt.plot(solt, solyp[:,0])
                #     plt.plot(solt, solyp[:,1])
                #     plt.show()
                #     quit()
                for i in range(v.nCboMat):
                    v.specPop[node][v.jti-1][i] = solyp[solyp.shape[0]-1][i]
                    if solyp[solyp.shape[0]-1][v.nCboMat+i] > 0.0:
                        v.depo[node][v.nClsMat+i] = solyp[solyp.shape[0]-1][v.nCboMat+i]
            else:
                for i in range(v.nCboMat):
                    v.specPop[node][v.jti-1][i] = 0.0
        else:
            for i in range(v.nCboMat):
                v.specPop[node][v.jti-1][i] = 0.0


def clastic():
    # It calculates the clastic sedimentation
    for mat in range(v.nClsMat):
        for node in range(v.nNode):
            # check if velocity is below critical value for sedim type
            # velo comes in m/y,  vCritic im m/day, set velo to m/day

            if v.veloNode[node]/365.25 <= v.vCritic[mat] and (v.wDepth[node]-v.head[node]) > v.wMin and v.conClast[node][mat] > 0.0:

                # clastmass is suspended mass (in m thickness of sediment)
                clastMass = (v.conClast[node][mat]*v.wDepth[node])/v.density[mat]
                # calculate settle velocity (comes in m/day)
                factor = (v.vCritic[mat]-(v.veloNode[node]/365.25))/v.vCritic[mat]
                if factor > 1:
                    factor = 1.0
                if factor < 0:
                    factor = 0.0
                # calculate settle distance for dt
                sediPortion = v.setClast[mat]*365.25*factor*v.dt/v.wDepth[node]
                if sediPortion > 1.0:
                    sediPortion = 1.0
                # calculate mass deposited in dt
                v.depo[node][mat] = clastMass*sediPortion
                # calculate new concentration after deposition
                clastMass = clastMass-v.depo[node][mat]
                v.conClast[node][mat] = (v.density[mat]*clastMass)/(v.wDepth[node]-v.depo[node][mat])
