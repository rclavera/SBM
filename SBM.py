# -*- coding: utf-8 -*-

#  Sedimentary Basin Modelling v.0.01
#
#***********************************************************************
# SBM is a program designed and developed by:
#
#           Dr. Roger Clavera-Gispert
#           Dr. Oscar Gratacos-Torra
#
#***********************************************************************
#   v.0.01-beta1
#***********************************************************************
# SBM is a 2D1/2 stratigraphic forward numerical model
# to simulate a basin infilling.
# It is originally based on SIMSAFADIM-CLASTIC program.
#****************************************************************************


#import sys
#import numpy as np
#from numpy import *
import modules.var as v
import modules.save as s
import modules.compute as com
import modules.sedimentation as sedi
import modules.time as ti
import modules.flow as f
import modules.transport as tr
import modules.iniData as i


i.initialData()
s.initialData()
v.dt = v.tInitStep
v.jti = 1

print ("Starting time step",v.jti,"of",v.nti)

# -----------------------------------------  Start loop over dt time
while v.totalTime<v.modelingTime:

# Computing sealevel variation through time
    com.seaLevel()

# Computing basin subsidence
    if v.comSubsidence: com.subsidence()

# Computing isostacy
    if v.comIsostasy: com.isostasy()
  
# Computing slope of all elements
    com.slope(v.jti-1)

# Computing water depth for all elements and nodes
    com.waderdepth()

# Computing Boundary Conditions movement
    com.boundaryConditions()

# Computing flow field
    f.flow()

# Computing transport of clastic sediment
    if v.comTransport: tr.transport()

# Computing sedimentation of clastic sediments
    if v.comClasticSedi: sedi.clastic()

# Computing sedimentation of carbonate sediments
    if v.comCarboSedi: sedi.carbonate()

# updating sediment thickness due to sedimentation
    sedi.update()

# Computing timestepping
