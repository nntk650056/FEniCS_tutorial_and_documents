
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 21:01:33 2020

@author: root
"""

from __future__ import print_function
import matplotlib.pyplot as plt
from fenics import *
import time
import numpy as np
import copy as cp

# Optimization options for the form compiler
parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["representation"] = "uflacs"
parameters["form_compiler"]["quadrature_degree"] = 3

# Create mesh and define function space
nx = 100
ny = 52
#mesh = RectangleMesh.create(MPI.comm_world, [Point(0, 0), Point(20, 20)], [nx, ny], CellType.Type.quadrilateral)
mesh = RectangleMesh(Point(0.0,0.0),Point(400.0,52.0),nx,ny)

#define global parameters and initial parameters
#Initialize C0 Nv mu0
D = 1 #diffusion coefficience
tol = 1E-14
chi = 0.2
C0_0 = 10
C0_1 = 0.02
C0 = Expression('x[1]<=50.0 + tol ? C0_0 : C0_1', degree=0, tol = tol, C0_0 = C0_0, C0_1 = C0_1)
Nv_0 = 0.001
Nv_1 = 1 
Nv = Expression('x[1]<=50.0 + tol ? Nv_0 : Nv_1', degree=0, tol = tol, Nv_0 = Nv_0, Nv_1 = Nv_1)   
         
lambda0 = (1+C0)**(1/3)

lambda01 = (1+C0_0)**(1/3)
lambda02 = (1+C0_1)**(1/3)  
mu_0 = ln(C0_0/(1+C0_0))+1/(1+C0_0)+chi/(1+C0_0)**2+Nv_0*(1/lambda01-1/(1+C0_0))
mu_1 = ln(C0_1/(1+C0_1))+1/(1+C0_1)+chi/(1+C0_1)**2+Nv_1*(1/lambda02-1/(1+C0_1))
mu0 = Expression('x[1]<=50.0 + tol ? mu_0 : mu_1', degree=0, tol = tol, mu_0 = mu_0, mu_1 = mu_1)
print(mu_0)
print(mu_1)
print(lambda01)
print(lambda02)


#define different C0 Nv mu0 in defferent domains

# Define chemical potential boundary conditions
#mu_D = Expression("mu_1 + (mu_0-mu_1)*t/90.0", degree = 1, mu_1 = mu_1, mu_0 = mu_0, t = 0)

# Create mesh and define function space
V = VectorElement('P',triangle,1)
Q = FiniteElement('P',triangle,1)
TH = V*Q
V1 = FunctionSpace(mesh, V)
Q1 = FunctionSpace(mesh, Q)
W = FunctionSpace(mesh,TH)

# Mark boundary subdomians
left =  CompiledSubDomain("near(x[0], side) && on_boundary", side = 0.0)
right = CompiledSubDomain("near(x[0], side) && on_boundary", side = 400.0)
top = CompiledSubDomain("near(x[1], side) && on_boundary", side = 52.0)
bottom = CompiledSubDomain("near(x[1], side) && on_boundary", side = 0.0)

# Define chemical potential boundary conditions


#fixed left and bottom sides and Define chemical potential boundary conditions
bcl = DirichletBC(W.sub(0).sub(0), Constant(0), left)
bcr = DirichletBC(W.sub(0).sub(0), Constant(0), right)
bcb = DirichletBC(W.sub(0).sub(1), Constant(0), bottom)
#bcu_mu = DirichletBC(W.sub(1), mu_D, top)

bcs = [bcl,bcr,bcb]

# Define test functions
w = Function(W)
u, mu = split(w)
wt = TestFunction(W)
v, mu_ = split(wt)
dw = TrialFunction(W)
w_n  = Function(W)
u_n, mu_n = split(w_n)


# Define deformation gradient tensor
F = grad(u) + Identity(len(u))
I_ = Identity(len(u))
Fv = grad(v)
F_mu = grad(mu)
F_muv = grad(mu_)
H = inv(F)
C = F.T*F                   # Right Cauchy-Green tensor
C_2 = C*C
# Invariants of deformation tensors
Ic = tr(C)
Ic2 = 0.5*(tr(C)*tr(C)-tr(C_2))
J  = det(F)
J_v = det(F)*H*Fv

C1 = lambda0**3*det(F)-1

# Save solution in VTK format
vtkfile_u = File("2D_gel_D1_wrinkle/displacement.pvd");
vtkfile_mu = File("2D_gel_D1_wrinkle/chemical_potential.pvd");
vtkfile_phi = File("2D_gel_D1_wrinkle/polumer_volume_fraction.pvd");
vtkfile_stress = File("2D_gel_D1_wrinkle/cauchy_stress.pvd");

#initial chemical potential value in the domian
e_u0 = Expression(('0', '0'), degree=1)
#e_mu0 = Expression('mu0', mu0=mu0, degree=1)
assign(w, [interpolate(e_u0, V1), interpolate(mu0, Q1)])


# p1 = plot(mu)
# plt.colorbar(p1)
# plt.show()
# Newton solver parameters
parameters={"newton_solver":{"relative_tolerance":1e-4, "absolute_tolerance":1e-6, "maximum_iterations":50}}

dt = Expression('dt1', degree=0, dt1=0.1)
gradv = grad((u-u_n)/dt)
dv2dt = gradv.T*gradv
W_FC = 0.5*Nv*(lambda0**2*Ic-3-2*ln(lambda0**3*det(F))) - C1*(ln(1+1/C1))-chi/(1+C1)
vis = 0.01
W_vis = vis*det(F)*(0.5*inner(gradv*H, (gradv*H+H.T*gradv.T)) - 1/3* tr(gradv*H)*tr(gradv*H))
#W_vis = vis*det(F)*(inner(H.T*gradv, H.T*grad((v)/dt)) + 0.5*inner(H.T*grad((v)/dt),H*gradv) + 0.5*inner(H.T*gradv,H*grad((v)/dt)) - 2/3* dot(H.T*div(grad((v)/dt)),H.T*div(gradv)))

Pi0 = -W_FC - W_vis*dt
a0 = (derivative(Pi0, w, wt) + mu*derivative(C1, w, wt))*dx + (-C1*D*dot(H.T*grad(mu),H.T*grad(mu_)) - lambda0**3*det(F)*inner(gradv,H.T)*mu_)*dx

J0 = derivative(a0, w, dw)



t = 0.0

for n in range(50+95*2+80):
    # Update previous solution
    w_n.assign(w)
    
    # Update current time range(0,0.1,1) range(2,1,50) range(60,10,90) range(100,50,500)
    if t < 0.1:
        dt.dt1 = 0.01
    elif t < 1:
        dt.dt1 = 0.1
    elif t < 100:
        dt.dt1 = 1
    elif t < 499.9:
        dt.dt1 = 5
    t += dt.dt1
        
    print("dt")
    print(dt)
    
    print(t)

    solve(a0==0, w, bcs,  J = J0, solver_parameters=parameters)

    p1 = plot(w.sub(1),title=str(t))
    plt.colorbar(p1)
    plt.show()
    p1 = plot(w.sub(0), mode = "displacement",title = 'disp t='+str(t))
    #plt.colorbar(plot(p1))
    plt.show()
    
    vtkfile_u << (w.sub(0), t)
    vtkfile_mu << (w.sub(1), t)
    
