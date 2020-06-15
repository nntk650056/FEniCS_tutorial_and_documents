#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 21:01:33 2020

@author: root
"""


import matplotlib.pyplot as plt
from fenics import *
#from mshr import *
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
mesh = RectangleMesh(Point(0.0,0.0),Point(100.0,52.0),nx,ny)
#define global parameters and initial parameters



vis = 0.01
D = 5   #diffusion constant

class K(Expression):
    def set_k_values(self, k_0, k_1):
        self.k_0, self.k_1 = k_0, k_1
    def eval(self , value , x):
        "set value[0] to value at point x"
        tol = 1E-14
        if x[1] <= 20 + tol:
            value[0] = self.k_0
        else:
            value[0] = self.k_1
            
#Initialize C0 Nv mu0
tol = 1E-14
chi = 0.2
C0_0 = 0.2
C0_1 = 0.2
Nv_0 = 0.001
Nv_1 = 1.0    
lambda01 = (1+C0_0)**(1/3)
lambda02 = (1+C0_1)**(1/3)  
mu_0 = ln(C0_0/(1+C0_0))+1/(1+C0_0)+chi/(1+C0_0)**2+Nv_0*(1/lambda01-1/(1+C0_0)) # mu10 = -0.81943 
mu_1 = ln(C0_1/(1+C0_1))+1/(1+C0_1)+chi/(1+C0_1)**2+Nv_1*(1/lambda02-1/(1+C0_1)) # mu10 = -0.81943       

  
#define different C0 Nv mu0 in defferent domains
C0 = Expression('x[1]<=50.0 + tol ? C0_0 : C0_1', degree=0, tol = tol, C0_0 = C0_0, C0_1 = C0_1)
Nv = Expression('x[1]<=50.0 + tol ? Nv_0 : Nv_1', degree=0, tol = tol, Nv_0 = Nv_0, Nv_1 = Nv_1)
mu0 = Expression('x[1]<=50.0 + tol ? mu_0 : mu_1', degree=0, tol = tol, mu_0 = mu_0, mu_1 = mu_1)
lambda0 = Expression('x[1]<=50.0 + tol ? lambda01: lambda02', degree=0, tol = tol, lambda01 = lambda01, lambda02 = lambda02)

#lambda0 = (1+C0)**(1/3)

#chi = 0.2
#mu10 = ln(C0_1/(1+C0_1))+1/(1+C0_1)+chi/(1+C0_1)**2+Nv*(1/lambda0-1/(1+C0_1)) # mu10 = -0.81943 


# Create mesh and define function space
V = VectorElement('P',triangle,2)
Q = FiniteElement('P',triangle,2)
plot(mesh,title='mesh')
plt.show()

V1 = VectorFunctionSpace(mesh, 'P', 2)
Q1 = FunctionSpace(mesh, 'P', 2)
TH = V*Q
W = FunctionSpace(mesh,TH)

# Mark boundary subdomians
left =  CompiledSubDomain("near(x[0], side) && on_boundary", side = 0.0)
right = CompiledSubDomain("near(x[0], side) && on_boundary", side = 100.0)
top = CompiledSubDomain("near(x[1], side) && on_boundary", side = 52.0)
bottom = CompiledSubDomain("near(x[1], side) && on_boundary", side = 0.0)

# Define initial value
#mu0 = Expression('mu00+a*x[0] + a*x[1]', degree=1, mu00 = mu10,a=0.0)
u0 = Expression('a*x[0] + a*x[1]', degree=1, a=0.0)
# Define chemical potential boundary conditions
#mu_D1 = Expression("mu00*(1-t/200.0)", degree = 1, mu00 = mu10, t = 0)
mu_D = Expression("mu00*(1-t/100.0)", degree = 1, mu00 = mu_1, t = 0)
#mu_Dbc = Expression("mu00+(-0.221-mu00)*((x[0]+x[1])*(x[0]+x[1])+1)*pow((t),0.5)", degree = 1, mu00 = mu10, t = 0)

#fixed left and bottom sides and Define chemical potential boundary conditions
#bcl = DirichletBC(W.sub(0).sub(0), Constant(0),left)
bcl = DirichletBC(W.sub(0).sub(0), Constant(0), left)
bcr = DirichletBC(W.sub(0).sub(0), Constant(0), right)
bcb = DirichletBC(W.sub(0).sub(1), Constant(0), bottom)
#bcl_mu = DirichletBC(W.sub(1), mu_D1,left)
#bcb_mu = DirichletBC(W.sub(1), mu_D1, bottom)
bct_mu = DirichletBC(W.sub(1), mu_D,top)
#bcr_mu = DirichletBC(W.sub(1), mu_D, right)
bcs = [bcl,bcr,bcb,bct_mu]

# Define test functions
w = Function(W)
#print(w)
w = Function(W)
u, mu = split(w)
wt = TestFunction(W)
v, mu_ = split(wt)
dw = TrialFunction(W)
w_n  = Function(W)
u_n, mu_n = split(w_n)

#initial chemical potential value in the domian
#initial chemical potential value in the domian
e_u0 = Expression(('0', '0'), degree=1)
#e_mu0 = Expression('mu0', mu0=mu0, degree=1)
assign(w, [interpolate(e_u0, V1), interpolate(Constant(mu_0), Q1)])




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
c1 = C1/det(F)/lambda0**3

#define the time increment
dt = Expression('dt1', degree=0, dt1=0.1)


#define governing equation
gradv = grad((u-u_n)/dt)

W_vis = vis*det(F)*(0.5*inner(gradv*H, (gradv*H+H.T*gradv.T)) -\
        1/3* tr(gradv*H)*tr(gradv*H))

W_FC = 0.5*Nv*(lambda0**2*Ic-3-2*ln(lambda0**3*det(F))) - \
        C1*(ln(1+1/C1))-chi/(1+C1)
Pi0 = -W_FC - W_vis*dt
a0 = (derivative(Pi0, w, wt) + mu*derivative(C1, w, wt))*dx +\
    (-C1*dot(H.T*grad(mu),H.T*grad(mu_)) -\
     lambda0**3*det(F)*inner(gradv,H.T)*mu_)*dx
J1 = derivative(a0, w, dw)




# Save solution in VTK format
vtkfile_u = File("wrinkle_test1/displacement.pvd");
vtkfile_mu = File("wrinkle_test1/chemical_potential.pvd");
vtkfile_w = File("wrinkle_test1/w.pvd");
vtkfile_phi = File("wrinkle_test1/polumer_volume_fraction.pvd");
vtkfile_stress = File("wrinkle_test1/cauchy_stress.pvd");
vtkfile_modulus = File("wrinkle_test1/elastic_modulus.pvd");



# Newton solver parameters
parameters={"newton_solver":{"relative_tolerance":1e-4, "absolute_tolerance":1e-6, "maximum_iterations":10}}

#define the intial and Maximum time increment and allowed Minimum time increment 
ini_time_int = 0.1
max_time_int = 1
min_time_int = 0.0001
total_time = 20
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
    
    
    
        
    #postprocessing
    p1 = plot(w.sub(1), title='chemical potential t='+str(t))
    plt.colorbar(p1)
    plt.show()
    p1 = plot(w.sub(0), mode = "displacement",title = 'displacement t='+str(t))
    plt.colorbar(p1)
    plt.show()
    
    vtkfile_u << (w.sub(0), t)
    vtkfile_mu << (w.sub(1), t)
    
    if t>=total_time:
        break
    
    
