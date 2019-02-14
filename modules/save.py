#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import numpy as np
from vtk import *
from vtk.util.numpy_support import *
import modules.var as v


def jtiRes():
  
#saving sealevelData  
  saveSeaData()
  saveBasinData()
  
def saveSeaData():
  sedi_vtk=[]
  for i in range(v.nClsMat):
    sedi=np.zeros((v.nNode),dtype=np.float64)
    sedi[:]= v.conClast[:,i]
    sedi_vtk.append(numpy_to_vtk(sedi,deep=True,array_type=vtk.VTK_DOUBLE))
    sedi_vtk[i].SetName((v.sediName[i]+'_conc'))
  head_vtk=numpy_to_vtk(v.head)
  head_vtk.SetName('head')
  wd_vtk=numpy_to_vtk(v.wDepth)
  wd_vtk.SetName('Water_depth')
  veloNode_vtk = numpy_to_vtk(v.veloNode/365.25)
  veloNode_vtk.SetName('velocity_node')
  io_vtk = numpy_to_vtk(v.typBC) 
  io_vtk.SetName('i/o')
  veloCell_vtk = numpy_to_vtk(v.ve3/365.25)
  veloCell_vtk.SetName('velocity_cell')
  Points=vtkPoints()
  Points.SetNumberOfPoints(v.nNode)
  for i in range(v.nNode):
    Points.SetPoint(i,(v.xCo[i],v.yCo[i],v.seaLevel))
  Cells = vtkCellArray()
  Cells.Allocate(v.nElem)
  for ne in range (v.nElem):
    Cells.InsertNextCell(3)
    Cells.InsertCellPoint(v.inn[ne][0])
    Cells.InsertCellPoint(v.inn[ne][1])
    Cells.InsertCellPoint(v.inn[ne][2])

  G = vtkUnstructuredGrid() 
  G.GetPointData().AddArray(veloNode_vtk)
  G.GetPointData().AddArray(head_vtk)
  G.GetPointData().AddArray(wd_vtk)
  for i in range (v.nClsMat):
    G.GetPointData().AddArray(sedi_vtk[i])
  G.GetPointData().AddArray(io_vtk)
  G.GetCellData().AddArray(veloCell_vtk)
  G.SetPoints(Points)
  G.SetCells(vtk.VTK_TRIANGLE,Cells)

  writer = vtkXMLUnstructuredGridWriter()
  outputname='results/seaLevel'+str(v.jti)+'.vtu'
  writer.SetFileName(outputname)
  writer.SetInputData(G)
  writer.SetDataModeToBinary()
  writer.SetCompressorTypeToZLib()
  writer.Write()

def saveBasinData():
  if v.jti==v.nti:
    ff=open('results/population1.txt','w')
    for n in range(v.nNode):
      ff.write(str(n)+' ')
      for j in range(v.nti):
        ff.write(str(v.specPop[n][j][0])+' ')
      ff.write('\n')
    ff.close
    ff=open('results/factorWD_pop1.txt','w')
    for n in range(v.nNode):
      ff.write(str(n)+' ')
      for j in range(v.nti):
        ff.write(str(v.factorWD[n][j][0])+' ')
      ff.write('\n')
    ff.close
    ff=open('results/factorSL_pop1.txt','w')
    for n in range(v.nNode):
      ff.write(str(n)+' ')
      for j in range(v.nti):
        ff.write(str(v.factorSL[n][j][0])+' ')
      ff.write('\n')
    ff.close
    ff=open('results/factorVL_pop1.txt','w')
    for n in range(v.nNode):
      ff.write(str(n)+' ')
      for j in range(v.nti):
        ff.write(str(v.factorVL[n][j][0])+' ')
      ff.write('\n')
    ff.close
    ff=open('results/factorCL1_pop1.txt','w')
    for n in range(v.nNode):
      ff.write(str(n)+' ')
      for j in range(v.nti):
        ff.write(str(v.factorCL[n][j][0][0])+' ')
      ff.write('\n')
    ff.close
    ff=open('results/factorCL2_pop1.txt','w')
    for n in range(v.nNode):
      ff.write(str(n)+' ')
      for j in range(v.nti):
        ff.write(str(v.factorCL[n][j][0][1])+' ')
      ff.write('\n')
    ff.close

#saving new basinData
  Points=vtkPoints()
  Points.SetNumberOfPoints(2*v.nNode*v.jti)
  NewNumNodes=2*v.nNode*v.jti
  pp=0
  for i in range(v.nNode):
    x=float(v.xCo[i])
    y=float(v.yCo[i])
    z=float(v.basem[i])
    Points.SetPoint(pp,(x,y,z))
    pp+=1
  for j in range (v.jti-1):
    for cf in range (2):
      for i in range(v.nNode):
        x=float(v.xCo[i])
        y=float(v.yCo[i])
        z=float(v.surf[j][i])
        Points.SetPoint(pp,(x,y,z))
        pp+=1
  for i in range(v.nNode):
    x=float(v.xCo[i])
    y=float(v.yCo[i])
    z=float(v.surface[i])
    Points.SetPoint(pp,(x,y,z))
    pp+=1
    
  Cells = vtkCellArray()
  Cells.Allocate(v.nElem*v.jti)
  for i in range(v.jti):
    ij=2*i
    for ne in range (v.nElem):
      Cells.InsertNextCell(6)
      Cells.InsertCellPoint((v.nNode*ij)+v.inn[ne][0])
      Cells.InsertCellPoint((v.nNode*ij)+v.inn[ne][1])
      Cells.InsertCellPoint((v.nNode*ij)+v.inn[ne][2])
      Cells.InsertCellPoint((v.nNode*(ij+1))+v.inn[ne][0])
      Cells.InsertCellPoint((v.nNode*(ij+1))+v.inn[ne][1])
      Cells.InsertCellPoint((v.nNode*(ij+1))+v.inn[ne][2])

  G = vtkUnstructuredGrid() 
  slN=np.zeros((NewNumNodes),float)
  cont=0
  for jj in range (v.jti):
    jjj=jj+1
    for i in range(v.nNode):
      slN[cont]=v.slopeNode[jj][i]
      slN[v.nNode+cont]=v.slopeNode[jj][i]
      cont+=1
    cont=(2*jjj*v.nNode)

  slopeNode_vtk=numpy_to_vtk(slN)
  slopeNode_vtk.SetName('Node_slope')
  G.GetPointData().AddArray(slopeNode_vtk)  
  slC=np.zeros((v.nElem*v.jti),float)
  pp=0
  for i in range(v.jti):
    for node in range(v.nElem):
      slC[pp]=v.slope[i][node]
      pp+=1
  slopeElem_vtk=numpy_to_vtk(slC)
  slopeElem_vtk.SetName('cell_slope')
  G.GetCellData().AddArray(slopeElem_vtk)

  sedi_vtk=[]
  for mat in range(v.nTotalMat):
    sedi=np.zeros((NewNumNodes),dtype=np.float64)
    cont=0
    for jj in range (v.jti):
      jjj=jj+1
      for i in range(v.nNode):
        sedi[cont]=v.grid[i][jj][mat]
        sedi[v.nNode+cont]=v.grid[i][jj][mat]
        cont+=1
      cont=(2*jjj*v.nNode)
    sedi_vtk.append(numpy_to_vtk(sedi,deep=True,array_type=vtk.VTK_DOUBLE))
    sedi_vtk[mat].SetName(v.sediName[mat])
    G.GetPointData().AddArray(sedi_vtk[mat])

  
  sediPC=np.zeros((v.nNode,v.jti,v.nTotalMat),dtype=np.float64)
  for i in range(v.nNode):
    for jj in range(v.jti):
      total=0.0
      for mat in range(v.nTotalMat):
        total+=v.grid[i][jj][mat]
      if total>0.0:
        for mat in range(v.nTotalMat):
          sediPC[i][jj][mat]=v.grid[i][jj][mat]/total*100.0

  
  sedi_vtk=[]
  for mat in range(v.nTotalMat):
    sedi=np.zeros((NewNumNodes),dtype=np.float64)
    cont=0
    for jj in range (v.jti):
      jjj=jj+1
      for i in range(v.nNode):
        sedi[cont]=sediPC[i][jj][mat]
        sedi[v.nNode+cont]=sediPC[i][jj][mat]
        cont+=1
      cont=(2*jjj*v.nNode)
    sedi_vtk.append(numpy_to_vtk(sedi,deep=True,array_type=vtk.VTK_DOUBLE))
    sedi_vtk[mat].SetName('%'+v.sediName[mat])
    G.GetPointData().AddArray(sedi_vtk[mat])

  sedi_vtk=[]
  for mat in range(v.nCboMat):
    sedi=np.zeros((NewNumNodes),dtype=np.float64)
    cont=0
    for jj in range (v.jti):
      jjj=jj+1
      for i in range(v.nNode):
        sedi[cont]=v.specPop[i][jj][mat]
        sedi[v.nNode+cont]=v.specPop[i][jj][mat]
        cont+=1
      cont=(2*jjj*v.nNode)
    sedi_vtk.append(numpy_to_vtk(sedi,deep=True,array_type=vtk.VTK_DOUBLE))
    sedi_vtk[mat].SetName('pop_'+v.sediName[v.nClsMat+mat])
    G.GetPointData().AddArray(sedi_vtk[mat])



  enf_vtk=[]
  for mat in range(v.nCboMat):
    enf=np.zeros((NewNumNodes),dtype=np.float64)
    cont=0
    for jj in range (v.jti):
      jjj=jj+1
      for i in range(v.nNode):
        enf[cont]=v.envFactor[i][jj][mat]
        enf[v.nNode+cont]=v.envFactor[i][jj][mat]
        cont+=1
      cont=(2*jjj*v.nNode)
    enf_vtk.append(numpy_to_vtk(enf,deep=True,array_type=vtk.VTK_DOUBLE))
    enf_vtk[mat].SetName(v.sediName[v.nClsMat+mat]+'_Factor_ENV')  
    G.GetPointData().AddArray(enf_vtk[mat])  

  enf_vtk=[]
  for mat in range(v.nCboMat):
    enf=np.zeros((NewNumNodes),dtype=np.float64)
    cont=0
    for jj in range (v.jti):
      jjj=jj+1
      for i in range(v.nNode):
        enf[cont]=v.factorWD[i][jj][mat]
        enf[v.nNode+cont]=v.factorWD[i][jj][mat]
        cont+=1
      cont=(2*jjj*v.nNode)
    enf_vtk.append(numpy_to_vtk(enf,deep=True,array_type=vtk.VTK_DOUBLE))
    enf_vtk[mat].SetName(v.sediName[v.nClsMat+mat]+'_Factor_WD')  
    G.GetPointData().AddArray(enf_vtk[mat])  

  enf_vtk=[]
  for mat in range(v.nCboMat):
    enf=np.zeros((NewNumNodes),dtype=np.float64)
    cont=0
    for jj in range (v.jti):
      jjj=jj+1
      for i in range(v.nNode):
        enf[cont]=v.factorVL[i][jj][mat]
        enf[v.nNode+cont]=v.factorVL[i][jj][mat]
        cont+=1
      cont=(2*jjj*v.nNode)
    enf_vtk.append(numpy_to_vtk(enf,deep=True,array_type=vtk.VTK_DOUBLE))
    enf_vtk[mat].SetName(v.sediName[v.nClsMat+mat]+'_Factor_VL')  
    G.GetPointData().AddArray(enf_vtk[mat])  


  enf_vtk=[]
  for mat in range(v.nCboMat):
    enf=np.zeros((NewNumNodes),dtype=np.float64)
    cont=0
    for jj in range (v.jti):
      jjj=jj+1
      for i in range(v.nNode):
        enf[cont]=v.factorSL[i][jj][mat]
        enf[v.nNode+cont]=v.factorSL[i][jj][mat]
        cont+=1
      cont=(2*jjj*v.nNode)
    enf_vtk.append(numpy_to_vtk(enf,deep=True,array_type=vtk.VTK_DOUBLE))
    enf_vtk[mat].SetName(v.sediName[v.nClsMat+mat]+'_Factor_SL')  
    G.GetPointData().AddArray(enf_vtk[mat])

  
  for mat2 in range(v.nClsMat):
    enf_vtk=[]
    for mat in range(v.nCboMat):
      enf=np.zeros((NewNumNodes),dtype=np.float64)
      cont=0
      for jj in range (v.jti):
        jjj=jj+1
        for i in range(v.nNode):
          enf[cont]=v.factorCL[i][jj][mat][mat2]
          enf[v.nNode+cont]=v.factorCL[i][jj][mat][mat2]
          cont+=1
        cont=(2*jjj*v.nNode)
      enf_vtk.append(numpy_to_vtk(enf,deep=True,array_type=vtk.VTK_DOUBLE))
      enf_vtk[mat].SetName(v.sediName[v.nClsMat+mat]+'_Factor_CL_'+v.sediName[mat2])
      G.GetPointData().AddArray(enf_vtk[mat])
      

  G.SetPoints(Points)
  G.SetCells(vtk.VTK_WEDGE,Cells)

  writer = vtkXMLUnstructuredGridWriter()
  outputname='results/surface'+str(v.jti)+'.vtu'
  writer.SetFileName(outputname)
  writer.SetInputData(G)
  writer.SetDataModeToBinary()
  writer.SetCompressorTypeToZLib()
  writer.Write()

def initialData():

#  It saves initial basement and surface data in numeric and vtk files

  dominant=np.zeros((v.nNode),float)
  sedi=np.zeros((v.nTotalMat,v.nNode),float)
  total=np.zeros((v.nNode),float)
  percent=np.zeros((v.nNode),float)
  
  dominant_vtk = numpy_to_vtk(dominant)
  dominant_vtk.SetName('dominant')
  total_vtk = numpy_to_vtk(total) 
  total_vtk.SetName('sed_total')
  percent_vtk = numpy_to_vtk(percent) 
  percent_vtk.SetName('percentage')
  slope_vtk = numpy_to_vtk(v.slopeNode) 
  slope_vtk.SetName('slope')

  Points=vtkPoints()
  Points.SetNumberOfPoints(v.nNode)
  for i in range(v.nNode):
    Points.SetPoint(i,(v.xCo[i],v.yCo[i],v.basem[i]))
  Cells = vtkCellArray()
  Cells.Allocate(v.nElem)
  for ne in range (v.nElem):
    Cells.InsertNextCell(3)
    Cells.InsertCellPoint(v.inn[ne][0])
    Cells.InsertCellPoint(v.inn[ne][1])
    Cells.InsertCellPoint(v.inn[ne][2])

  G = vtkUnstructuredGrid() 
  G.GetPointData().AddArray(dominant_vtk)
  for i in range(v.nTotalMat):
    sedi_vtk=numpy_to_vtk(sedi[i][:])
    sedi_vtk.SetName(v.sediName[i])
    G.GetPointData().AddArray(sedi_vtk)
  G.GetPointData().AddArray(total_vtk)
  G.GetPointData().AddArray(percent_vtk)
  G.GetPointData().AddArray(slope_vtk)
  G.SetPoints(Points)
  G.SetCells(vtk.VTK_TRIANGLE,Cells)

  writer = vtkXMLUnstructuredGridWriter()
  outputname='results/surf0.vtu'
  writer.SetFileName(outputname)
  writer.SetInputData(G)
  writer.SetDataModeToBinary()
  writer.SetCompressorTypeToZLib()
  writer.Write()


  veloNode=np.zeros((v.nNode),float)
  veloCell=np.zeros((v.nElem),float)
  sediConc=np.zeros((v.nTotalMat,v.nNode),float)
  
  veloNode_vtk = numpy_to_vtk(veloNode)
  veloNode_vtk.SetName('velocity_node')
  io_vtk = numpy_to_vtk(v.typBC) 
  io_vtk.SetName('i/o')
  veloCell_vtk = numpy_to_vtk(veloCell)
  veloCell_vtk.SetName('velocity_cell')

  Points=vtkPoints()
  Points.SetNumberOfPoints(v.nNode)
  for i in range(v.nNode):
    Points.SetPoint(i,(v.xCo[i],v.yCo[i],v.seaLevel))
  Cells = vtkCellArray()
  Cells.Allocate(v.nElem)
  for ne in range (v.nElem):
    Cells.InsertNextCell(3)
    Cells.InsertCellPoint(v.inn[ne][0])
    Cells.InsertCellPoint(v.inn[ne][1])
    Cells.InsertCellPoint(v.inn[ne][2])

  G = vtkUnstructuredGrid() 
  G.GetPointData().AddArray(veloNode_vtk)
  for i in range(v.nTotalMat):
    sedi_vtk=numpy_to_vtk(sediConc[i][:])
    sedi_vtk.SetName((v.sediName[i]+'_conc'))
    G.GetPointData().AddArray(io_vtk)
  G.GetCellData().AddArray(veloCell_vtk)
  G.SetPoints(Points)
  G.SetCells(vtk.VTK_TRIANGLE,Cells)

  writer = vtkXMLUnstructuredGridWriter()
  outputname='results/seaLevel0.vtu'
  writer.SetFileName(outputname)
  writer.SetInputData(G)
  writer.SetDataModeToBinary()
  writer.SetCompressorTypeToZLib()
  writer.Write()
  import matplotlib.pyplot as plt
  for i in range (v.nCboMat):
    
    xWD = np.linspace(0.0, max(v.deptFactorX[i][:])+10, num=1001)
    xVL = np.linspace(0.0, max(v.flowFactorX[i][:])+10, num=1001)
    xSL = np.linspace(0.0, max(v.slopeFactorX[i][:])+10, num=1001)
    
    fWD=np.interp(xWD,v.deptFactorX[i][:],v.deptFactorY[i][:])
    fVL=np.interp(xVL,v.flowFactorX[i][:],v.flowFactorY[i][:])
    fSL=np.interp(xSL,v.slopeFactorX[i][:],v.slopeFactorY[i][:])
    f1= plt.figure()
    plt.ylim(-0.1,1.1)
    plt.xlim(0.0,max(v.deptFactorX[i][:])+10)
    plt.grid()
    plt.xlabel('water depth (m)')
    plt.ylabel('f_wd')
    plt.plot(v.deptFactorX[i][:], v.deptFactorY[i][:], 'o', mew=2)
    plt.plot(xWD, fWD)
    leg = plt.legend(['points', 'interpolation'])
    leg.get_frame().set_facecolor('#fafafa')
    f1.savefig('results/fig_'+v.sediName[v.nClsMat+i]+'_WD.eps')
    
    f2= plt.figure()
    plt.ylim(-0.1,1.1)
    plt.xlim(0.0,max(v.flowFactorX[i][:]))
    plt.grid()
    plt.xlabel('flow velocity (m/d)')
    plt.ylabel('f_vl')
    plt.plot(v.flowFactorX[i][:], v.flowFactorY[i][:], 'o', mew=2)
    plt.plot(xVL, fVL)
    leg = plt.legend(['points', 'interpolation'])
    leg.get_frame().set_facecolor('#fafafa')
    f2.savefig('results/fig_'+v.sediName[v.nClsMat+i]+'_VL.eps')

    f3= plt.figure()
    plt.ylim(-0.1,1.1)
    plt.xlim(0.0,90.0)
    plt.grid()
    plt.xlabel('slope')
    plt.ylabel('f_sl')
    plt.plot(v.slopeFactorX[i][:], v.slopeFactorY[i][:], 'o', mew=2)
    plt.plot(xSL, fSL)
    leg = plt.legend(['points', 'interpolation'])
    leg.get_frame().set_facecolor('#fafafa')
    f3.savefig('results/fig_'+v.sediName[v.nClsMat+i]+'_SL.eps') 
  
  print ('initial surface and sea level saved')

