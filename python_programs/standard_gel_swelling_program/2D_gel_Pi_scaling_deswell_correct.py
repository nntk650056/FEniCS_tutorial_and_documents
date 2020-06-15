
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
nx = 40
ny = 40
#mesh = RectangleMesh.create(MPI.comm_world, [Point(0, 0), Point(20, 20)], [nx, ny], CellType.Type.quadrilateral)
mesh = RectangleMesh(Point(0.0,0.0),Point(20.0,20.0),nx,ny)

#define global parameters and initial parameters
C0 = 0.5
lambda0 = (1+C0)**(1/3)
Nv = 0.001
chi = 0.2
mu0 = ln(C0/(1+C0))+1/(1+C0)+chi/(1+C0)**2+Nv*(1/lambda0-1/(1+C0)) # mu0 = -0.81943 
vis = 0.01
v_const = 0.588
vtotal = -(3*v_const**2-2.5*v_const+0.5)
phi0v = (1/(C0+1))**vtotal
print(phi0v)

# Create mesh and define function space
V = VectorElement('P',triangle,2)
Q = FiniteElement('P',triangle,1)
TH = V*Q
V1 = FunctionSpace(mesh, V)
Q1 = FunctionSpace(mesh, Q)
W = FunctionSpace(mesh,TH)

# Mark boundary subdomians
left =  CompiledSubDomain("near(x[0], side) && on_boundary", side = 0.0)
right = CompiledSubDomain("near(x[0], side) && on_boundary", side = 20.0)
top = CompiledSubDomain("near(x[1], side) && on_boundary", side = 20.0)
bottom = CompiledSubDomain("near(x[1], side) && on_boundary", side = 0.0)

# Define chemical potential boundary conditions
mu_D = Expression("mu00*(1+t/100.0)", degree = 1, mu00 = mu0, t = 0)

#fixed left and bottom sides and Define chemical potential boundary conditions
bcl = DirichletBC(W.sub(0).sub(0), Constant(0), left)
bcb = DirichletBC(W.sub(0).sub(1), Constant(0), bottom)
bct_mu = DirichletBC(W.sub(1), mu_D,top)
bcr_mu = DirichletBC(W.sub(1), mu_D, right)
bcs = [bcl,bcb,bct_mu,bcr_mu]

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
phi = 1/(1+C1)
print(phi)
lambda_ref = (phi**vtotal)/phi0v





# Save solution in VTK format
vtkfile_u = File("2D_gel_D_deswelling/displacement.pvd");
vtkfile_mu = File("2D_gel_D_deswelling/chemical_potential.pvd");
vtkfile_phi = File("2D_gel_D_deswelling/polumer_volume_fraction.pvd");
vtkfile_stress = File("2D_gel_D_deswelling/cauchy_stress.pvd");

#initial chemical potential value in the domian
e_u0 = Expression(('0', '0'), degree=1)
e_mu0 = Expression('mu0', mu0=mu0, degree=1)
assign(w, [interpolate(e_u0, V1), interpolate(e_mu0, Q1)])


p1 = plot(mu)
plt.colorbar(p1)
plt.show()
# Newton solver parameters
parameters={"newton_solver":{"relative_tolerance":1e-4, "absolute_tolerance":1e-6, "maximum_iterations":10}}

dt = Expression('dt1', degree=0, dt1=0.1)

gradv = grad((u-u_n)/dt)
W_vis = vis*det(F)*(0.5*inner(gradv*H, (gradv*H+H.T*gradv.T)) - 1/3* tr(gradv*H)*tr(gradv*H))
W_FC = 0.5*Nv*(phi*1.2)**1.07*(Ic/(lambda_ref**2)-3-2*ln(det(F)/( lambda_ref**3))) - C1*(ln(1+1/C1))-chi/(1+C1)

#W_vis = vis*det(F)*(inner(H.T*gradv, H.T*grad((v)/dt)) + 0.5*inner(H.T*grad((v)/dt),H*gradv) + 0.5*inner(H.T*gradv,H*grad((v)/dt)) - 2/3* dot(H.T*div(grad((v)/dt)),H.T*div(gradv)))

Pi0 = -W_FC - W_vis*dt
a0 = (derivative(Pi0, w, wt) + mu*derivative(C1, w, wt))*dx + (-C1*dot(H.T*grad(mu),H.T*grad(mu_)) - lambda0**3*det(F)*inner(gradv,H.T)*mu_)*dx
J1 = derivative(a0, w, dw)

# Newton solver parameters
parameters={"newton_solver":{"relative_tolerance":1e-4, "absolute_tolerance":1e-6, "maximum_iterations":10}}

#define the intial and Maximum time increment and allowed Minimum time increment 
ini_time_int = 0.1
max_time_int = 1
min_time_int = 0.0001
total_time = 200
max_steps = 1000
#initial time
t= ini_time_int

for n in range(max_steps):
    # Update previous solution
    w_n.assign(w)
    
    # Update current time range(0,0.1,1) range(2,1,50) range(60,10,90) range(100,50,500)
    
    
    
    
    str1 = "time = " + str(t)
    print(str1)
    
    try:
        if t >89.9:
            mu_D.t = 90.0
        else:
            mu_D.t = t
        solve(a0==0, w, bcs,  J = J1, solver_parameters=parameters)
        if dt.dt1 <=1-1E-6:
            dt.dt1 = dt.dt1*2
        else:
            dt.dt1 = max_time_int
        str1 = "time increment = " + str(dt.dt1)
        print(str1)
        t = t + dt.dt1
        
    except RuntimeError:
        print("time increment is too large")
        t = t - dt.dt1
        dt.dt1 = dt.dt1/5
        t = t + dt.dt1
        
        if t >89.9:
            mu_D.t = 90.0
        else:
            mu_D.t = t
        
        solve(a0==0, w, bcs,  J = J1, solver_parameters=parameters)
        if dt.dt1 < 1:
            dt.dt1 = dt.dt1*2
        
        if dt.dt1 < min_time_int:
            break
        
        str1 = "time increment = " + str(dt.dt1)
        print(str1)
        t = t + dt.dt1
    
    except RuntimeError:
        print("time increment is too large")
        t = t - dt.dt1
        dt.dt1 = dt.dt1/5
        t = t + dt.dt1
        
        if t >89.9:
            mu_D.t = 90.0
        else:
            mu_D.t = t
        
        solve(a0==0, w, bcs,  J = J1, solver_parameters=parameters)
        if dt.dt1 < 1:
            dt.dt1 = dt.dt1*2
        
        if dt.dt1 < min_time_int:
            break
        
        str1 = "time increment = " + str(dt.dt1)
        print(str1)
        t = t + dt.dt1
        
    except RuntimeError:
        print("time increment is too large")
        t = t - dt.dt1
        dt.dt1 = dt.dt1/5
        t = t + dt.dt1
        
        if t >89.9:
            mu_D.t = 90.0
        else:
            mu_D.t = t
        
        solve(a0==0, w, bcs,  J = J1, solver_parameters=parameters)
        if dt.dt1 < 1:
            dt.dt1 = dt.dt1*2
        
        if dt.dt1 < min_time_int:
            break
        
        str1 = "time increment = " + str(dt.dt1)
        print(str1)
        t = t + dt.dt1
    
    

    

    p1 = plot(w.sub(1),mode = "displacement", title='swell mu t='+str(t))
    plt.colorbar(plot(w.sub(1)))
    plt.show()
    p1 = plot(w.sub(0), mode = "displacement",title = 'swell displacement t='+str(t))
    plt.colorbar(p1)
    plt.xlim(-1.5,27.5)
    plt.ylim(-1.5,27.5)
    plt.show()
    
    vtkfile_u << (w.sub(0), t)
    vtkfile_mu << (w.sub(1), t)
    
    if t>=total_time:
        break
    
