# -*- coding: utf-8 -*-

import modules.var as v
import numpy as np
from scipy.linalg import lu, inv
#import matsolflow as fort
#from dolfin import *


def flow():

#     Subroutine for instationary flow simulation
#     General program description see Kinzelbach 1986: 282 ff.
#     Version 1.1 from 21.10.1997

  w=np.zeros((3),dtype=np.float64)

  v.tt=np.copy(v.wDepthE)*31557600.0 #meters/year

  for n in range(v.nElem):
    w[0] = v.wDepth[v.inn[n][0]]-v.head[v.inn[n][0]]
    w[1] = v.wDepth[v.inn[n][1]]-v.head[v.inn[n][1]]
    w[2] = v.wDepth[v.inn[n][2]]-v.head[v.inn[n][2]]
    if (max(w)>0.0 and min(w)<=0.0): v.tt[n]=0.0
  
  
  #print fort.matsolflow.__doc__
  #for i in range (len(v.head)):
    #headF[i]=v.head[i]
    #headd[i]=v.head[i]

  #fort.matsolflow(v.dt,v.q,v.bw,v.typBC,v.head,v.inn,v.tt,v.cn,v.nl,v.xCo,v.yCo,v.nNode,v.nElem)

  #headF=v.head
  #v.head=headd
  #fms=open('matsol.txt','w')
  
  matsolKinzelbach()
 
  #matsolDolfin()
  
  for node in range (v.nNode):
    if v.surface[node]-v.seaLevel>0.0:
      v.head[node]=0.0

  vector()
  veloTrans()


def veloTrans():

# change velocity from element property to nodal property

  for i in range(v.nNode):
    vlN=0.0
    for j in range (v.nodElems[i][6]):
      vlN=vlN+v.ve3[v.nodElems[i][j]]
    v.veloNode[i]=vlN/float(v.nodElems[i][6])

def vector():
  
  w=np.zeros((3),dtype=np.float64)
# This subroutine calculates flow vectors from a given potential field
  xTemp=np.zeros((3),dtype=np.float64)
  yTemp=np.zeros((3),dtype=np.float64)
  b=np.zeros((3),dtype=np.float64)
  c=np.zeros((3),dtype=np.float64)

# Calculate flow vectors vx and vy, flow direction and velocity

  for n in range(v.nElem):
    hTemp=np.zeros((3),dtype=np.float64)
    xTemp[0]=v.xCo[v.inn[n][0]]
    yTemp[0]=v.yCo[v.inn[n][0]]
    if (v.wDepth[v.inn[n][0]]-v.head[v.inn[n][0]]>0.0): hTemp[0]=v.head[v.inn[n][0]]

    xTemp[1]=v.xCo[v.inn[n][1]]
    yTemp[1]=v.yCo[v.inn[n][1]]
    if (v.wDepth[v.inn[n][1]]-v.head[v.inn[n][1]]>0.0): hTemp[1]=v.head[v.inn[n][1]]

    xTemp[2]=v.xCo[v.inn[n][2]]
    yTemp[2]=v.yCo[v.inn[n][2]]
    if (v.wDepth[v.inn[n][2]]-v.head[v.inn[n][2]]>0.0): hTemp[2]=v.head[v.inn[n][2]]

    b[0]=yTemp[1]-yTemp[2]
    b[1]=yTemp[2]-yTemp[0]
    b[2]=yTemp[0]-yTemp[1]
    c[0]=xTemp[2]-xTemp[1]
    c[1]=xTemp[0]-xTemp[2]
    c[2]=xTemp[1]-xTemp[0]

    de=abs(xTemp[0]*yTemp[1]-xTemp[0]*yTemp[2]+xTemp[1]*yTemp[2]-xTemp[1]*yTemp[0]+xTemp[2]*yTemp[0]-xTemp[2]*yTemp[1])

# Calculate flow for each triangular element velocity per length of time step
    vx1=(-1)*v.tt[n]*(hTemp[0]*b[0]+hTemp[1]*b[1]+hTemp[2]*b[2])/de
    vy1=(-1)*v.tt[n]*(hTemp[0]*c[0]+hTemp[1]*c[1]+hTemp[2]*c[2])/de

    v.velDi3[n]=np.arctan2(vx1,vy1)*57.29577951308232 #180/PI
    if(v.wDepthE[n]>0.0):
      v.ve3[n]=np.sqrt(vx1*vx1+vy1*vy1)/v.wDepthE[n]
      if(v.ve3[n]<=0.0): v.ve3[n]=0.0000000001
    else:
      v.ve3[n]=0.0000000001

    w[0] = v.wDepth[v.inn[n][0]]-v.head[v.inn[n][0]]
    w[1] = v.wDepth[v.inn[n][1]]-v.head[v.inn[n][1]]
    w[2] = v.wDepth[v.inn[n][2]]-v.head[v.inn[n][2]]
    if (max(w)>0.0 and min(w)<=0.0): v.ve3[n]=0.0000000001

def matsolKinzelbach():
#       solution of system of equations from fe mesh
#       adapted from Kinzelbach 1985 p. 290

  bb=np.zeros((v.nNode),dtype=np.float64)
  a=np.zeros((v.nNode,v.bw+1),dtype=np.float64)
  b=np.zeros((3),dtype=np.float64)
  c=np.zeros((3),dtype=np.float64)
  re=np.ones((3,3),dtype=np.float64)
  pe=np.zeros((3,3),dtype=np.float64)

  rd=1.0/v.dt
  for i in range(v.nNode):
    bb=np.copy(v.q)

  for i in range(v.nNode):
    if v.typBC[i]==0:
      a[i][v.nl]=1.0
      bb[i]=v.head[i]
#       matrix assembly loop over elements
  for l in range(v.nElem):
    b[0]=v.yCo[v.inn[l][1]]-v.yCo[v.inn[l][2]]
    b[1]=v.yCo[v.inn[l][2]]-v.yCo[v.inn[l][0]]
    b[2]=v.yCo[v.inn[l][0]]-v.yCo[v.inn[l][1]]
    c[0]=v.xCo[v.inn[l][2]]-v.xCo[v.inn[l][1]]
    c[1]=v.xCo[v.inn[l][0]]-v.xCo[v.inn[l][2]]
    c[2]=v.xCo[v.inn[l][1]]-v.xCo[v.inn[l][0]]
    de=(b[0]*c[1]-b[1]*c[0])*0.5
    pe[0][0]=(b[0]*b[0]*v.tt[l]+c[0]*c[0]*v.tt[l]+b[0]*c[0])/de*0.25
    pe[0][1]=(b[0]*b[1]*v.tt[l]+c[0]*c[1]*v.tt[l]+b[0]*c[1])/de*0.25
    pe[0][2]=(b[0]*b[2]*v.tt[l]+c[0]*c[2]*v.tt[l]+b[0]*c[2])/de*0.25
    pe[1][0]=(b[1]*b[0]*v.tt[l]+c[1]*c[0]*v.tt[l]+b[1]*c[0])/de*0.25
    pe[1][1]=(b[1]*b[1]*v.tt[l]+c[1]*c[1]*v.tt[l]+b[1]*c[1])/de*0.25
    pe[1][2]=(b[1]*b[2]*v.tt[l]+c[1]*c[2]*v.tt[l]+b[1]*c[2])/de*0.25
    pe[2][0]=(b[2]*b[0]*v.tt[l]+c[2]*c[0]*v.tt[l]+b[2]*c[0])/de*0.25
    pe[2][1]=(b[2]*b[1]*v.tt[l]+c[2]*c[1]*v.tt[l]+b[2]*c[1])/de*0.25
    pe[2][2]=(b[2]*b[2]*v.tt[l]+c[2]*c[2]*v.tt[l]+b[2]*c[2])/de*0.25
    for i in range (3):
      for j in range (3):
        re[i][j]=de/12.0
    re[0][0]=de/6.0
    re[1][1]=de/6.0
    re[2][2]=de/6.0
# global matrix assembly
    ii=v.inn[l][0]
    if (v.typBC[ii]!=0):
      bb[ii]=bb[ii]-(pe[0][0]*v.head[v.inn[l][0]])
      if (v.typBC[v.inn[l][0]])!=0:
        jj=v.inn[l][0]-ii+v.nl
        if (jj>-1):
          a[ii][jj]=a[ii][jj]+pe[0][0]*v.cn+ re[0][0]*rd
      bb[ii]=bb[ii]-(pe[0][1]*v.head[v.inn[l][1]])
      if (v.typBC[v.inn[l][1]]!=0):
        jj=v.inn[l][1]-ii+v.nl
        if (jj>-1):
          a[ii][jj]=a[ii][jj]+(pe[0][1])*v.cn+ re[0][1]*rd
      bb[ii]=bb[ii]-(pe[0][2])*v.head[v.inn[l][2]]
      if (v.typBC[v.inn[l][2]]!=0):
        jj=v.inn[l][2]-ii+v.nl
        if (jj>-1): 
          a[ii][jj]=a[ii][jj]+(pe[0][2])*v.cn+ re[0][2]*rd

    ii=v.inn[l][1]
    if (v.typBC[ii]!=0):
      bb[ii]=bb[ii]-(pe[1][0]*v.head[v.inn[l][0]])
      if (v.typBC[v.inn[l][0]])!=0:
        jj=v.inn[l][0]-ii+v.nl
        if (jj>-1):
          a[ii][jj]=a[ii][jj]+pe[1][0]*v.cn+ re[1][0]*rd
      bb[ii]=bb[ii]-(pe[1][1]*v.head[v.inn[l][1]])
      if (v.typBC[v.inn[l][1]]!=0):
        jj=v.inn[l][1]-ii+v.nl
        if (jj>-1):
          a[ii][jj]=a[ii][jj]+(pe[1][1])*v.cn+ re[1][1]*rd
      bb[ii]=bb[ii]-(pe[1][2])*v.head[v.inn[l][2]]
      if (v.typBC[v.inn[l][2]]!=0):
        jj=v.inn[l][2]-ii+v.nl
        if (jj>-1): 
          a[ii][jj]=a[ii][jj]+(pe[1][2])*v.cn+ re[1][2]*rd

    ii=v.inn[l][2]
    if (v.typBC[ii]!=0):
      bb[ii]=bb[ii]-(pe[2][0]*v.head[v.inn[l][0]])
      if (v.typBC[v.inn[l][0]])!=0:
        jj=v.inn[l][0]-ii+v.nl
        if (jj>-1):
          a[ii][jj]=a[ii][jj]+pe[2][0]*v.cn+ re[2][0]*rd
      bb[ii]=bb[ii]-(pe[2][1]*v.head[v.inn[l][1]])
      if (v.typBC[v.inn[l][1]]!=0):
        jj=v.inn[l][1]-ii+v.nl
        if (jj>-1):
          a[ii][jj]=a[ii][jj]+(pe[2][1])*v.cn+ re[2][1]*rd
      bb[ii]=bb[ii]-(pe[2][2])*v.head[v.inn[l][2]]
      if (v.typBC[v.inn[l][2]]!=0):
        jj=v.inn[l][2]-ii+v.nl
        if (jj>-1): 
          a[ii][jj]=a[ii][jj]+(pe[2][2])*v.cn+ re[2][2]*rd



  #prova= gausselim(a,bb)
  #print prova
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

  #print bb-prova

#  quit()
# Equation solver finished
# Update head
  for i in range(v.nNode):
    if v.typBC[i]!=0:
      v.head[i]=v.head[i]+bb[i]



def gausselim(A,B):
    """
    Solve Ax = B using Gaussian elimination and LU decomposition.
    A = LU   decompose A into lower and upper triangular matrices
    LUx = B  substitute into original equation for A
    Let y = Ux and solve:
    Ly = B --> y = (L^-1)B  solve for y using "forward" substitution
    Ux = y --> x = (U^-1)y  solve for x using "backward" substitution
    :param A: coefficients in Ax = B
    :type A: numpy.ndarray of size (m, n)
    :param B: dependent variable in Ax = B
    :type B: numpy.ndarray of size (m, 1)
    """
    
    c= np.reshape(B,(B.size,1))
    # LU decomposition with pivot
    pl, u = lu(A, permute_l=True)
    # forward substitution to solve for Ly = B
    y = np.zeros(c.size)
    print (pl.shape)
    for m, b in enumerate(c.flatten()):
        y[m] = b
        # skip for loop if m == 0
        if m:
            for n in xrange(m):
                y[m] -= y[n] * pl[m,n]
        y[m] /= pl[m, m]

    # backward substitution to solve for y = Ux
    x = np.zeros(c.size)
    lastidx = c.size - 1  # last index
    for midx in xrange(c.size):
        m = c.size - 1 - midx  # backwards index
        x[m] = y[m]
        if midx:
            for nidx in xrange(midx):
                n = c.size - 1  - nidx
                x[m] -= x[n] * u[m,n]
        x[m] /= u[m, m]
    return x






def matsolDolfin():

  fb=open('mesh_mesh.xml','w')
  fb.write('<dolfin>\n<mesh celltype=\"triangle\" dim=\"2\">\n<vertices size=\"'+str(v.nNode)+'\">\n')
  for i in range(v.nNode):
    fb.write('<vertex index="'+str(i)+'\" x=\"'+str(v.xCo[i])+'\" y=\"'+str(v.yCo[i])+'\" z=\"'+str(v.wDepth[i])+'\"/>\n')
  fb.write('</vertices>\n<cells size="'+str(v.nElem)+'\">\n')
  for i in range(v.nElem):
    fb.write('<triangle index="'+str(i)+'\" v0=\"'+str(v.inn[i][0])+'\" v1=\"'+str(v.inn[i][1])+'\" v2=\"'+str(v.inn[i][2])+'\"/>\n')
  fb.write('</cells>\n</mesh>\n</dolfin>')
    
    
  fb=open('mesh_subdomain.xml','w')
  fb.write('<dolfin>\n<mesh_function>\n<mesh_value_collection type="uint" dim="1" size=\"'+str(v.nElem*3)+'\">\n')
  for i in range(v.nElem):
    for j in range(3):
      fb.write('<value cell_index="'+str(i)+'\" local_entity=\"'+str(j)+'\" value=\"'+str(int(v.typBC[v.inn[i][j]]))+'\"/>\n')
  fb.write('</mesh_value_collection>\n</mesh_function>\n</dolfin>')
  fb.close()



  # Load mesh and subdomains
  mesh = Mesh("mesh_mesh.xml")
  sub_domains = MeshFunction("size_t", mesh, "mesh_subdomain.xml")

  
# Define function spaces
  scalar = FunctionSpace(mesh, "CG", 1)
  vector = VectorFunctionSpace(mesh, "CG", 1)
  system = vector * scalar

# Create functions for boundary conditions
  noslip = Constant((0, 0))
  bc0 = DirichletBC(system.sub(0), noslip, sub_domains, 2)

  inflow = Constant((100, 0))#Expression(("(x[2])", "2"))
  bc1 = DirichletBC(system.sub(0), inflow, sub_domains, 3)

  aix =  Expression(("(x[2])", "0"))
  bc2 = DirichletBC(system.sub(0), aix, sub_domains, 1)

  ups = Constant((0, 0))
  bc3 = DirichletBC(system.sub(0), ups, sub_domains, 0)

# Collect boundary conditions
  bcs = [bc0,bc1,bc3]

# Define variational problem
  (vv, q) = TestFunctions(system)
  (u, p) = TrialFunctions(system)
  f = Expression(('x[2]', '0'))
  h = CellSize(mesh)
  beta  = 0.2
  delta = beta*h*h
  a = (inner(grad(vv), grad(u)) - div(vv)*p + q*div(u) + delta*inner(grad(q), grad(p)))*dx
  L = inner(vv + delta*grad(q), f)*dx

# Compute solution
  w = Function(system)
  solve(a == L, w, bcs)
  u, p = w.split()

# Save solution in VTK format
  ufile_pvd = File("velocity.pvd")
  ufile_pvd << u
  ufile=File("velo.xml")
  ufile << w
# Plot solution
  plot(u)
  interactive()

  quit()
  
def matsolDolfin2():
  fb=open('mesh_mesh.xml','w')
  fb.write('<dolfin>\n<mesh celltype=\"triangle\" dim=\"2\">\n<vertices size=\"'+str(v.nNode)+'\">\n')
  for i in range(v.nNode):
    fb.write('<vertex index="'+str(i)+'\" x=\"'+str(v.xCo[i])+'\" y=\"'+str(v.yCo[i])+'\" z=\"'+str(v.wDepth[i])+'\"/>\n')
  fb.write('</vertices>\n<cells size="'+str(v.nElem)+'\">\n')
  for i in range(v.nElem):
    fb.write('<triangle index="'+str(i)+'\" v0=\"'+str(v.inn[i][0])+'\" v1=\"'+str(v.inn[i][1])+'\" v2=\"'+str(v.inn[i][2])+'\"/>\n')
  fb.write('</cells>\n</mesh>\n</dolfin>')
    
    
  fb=open('mesh_subdomain.xml','w')
  fb.write('<dolfin>\n<mesh_function>\n<mesh_value_collection type="uint" dim="1" size=\"'+str(v.nElem*3)+'\">\n')
  for i in range(v.nElem):
    for j in range(3):
      fb.write('<value cell_index="'+str(i)+'\" local_entity=\"'+str(j)+'\" value=\"'+str(int(v.typBC[v.inn[i][j]]))+'\"/>\n')
  fb.write('</mesh_value_collection>\n</mesh_function>\n</dolfin>')
  fb.close()



  # Load mesh and subdomains
  mesh = Mesh("mesh_mesh.xml")
  sub_domains = MeshFunction("size_t", mesh, "mesh_subdomain.xml")

  
# Define function spaces
  V = FunctionSpace(mesh, "CG", 1)
# Create functions for boundary conditions
  noslip = Constant(0.0)
  aix =  Expression("(x[2]*31557600.0)")
  inflow = Constant(1.0)
  zero   = Constant(0.0)
  ups = Constant(0.0)

# No-slip boundary condition for velocity
  bc0 = DirichletBC(V, noslip, sub_domains, 2)
# Inflow boundary condition for velocity
  bc1 = DirichletBC(V, inflow, sub_domains, 3)
  bc2 = DirichletBC(V, aix, sub_domains, 1)
  bc3 = DirichletBC(V, ups, sub_domains, 0)

# Collect boundary conditions
  bcs = [bc1,bc3]#[bc0, bc1, bc2]

# Define variational problem
  u = TrialFunction(V)
  vv = TestFunction(V)
  f = Expression("x[2]*31557600.0")
  a = dot(grad(u), grad(vv))*dx
  L = f*vv*dx
  if has_petsc():
    PETScOptions.set("mat_mumps_icntl_14", 40.0)
# Compute solution
  u = Function(V)
  solve(a == L, u, bcs)

# Plot solution and solution gradient
  plot(u, title="Solution")
  plot(grad(u), title="Solution gradient")
  interactive()

# Save solution in VTK format
  file = File("nonlinear_poisson.xml")
  file << u


  quit()
