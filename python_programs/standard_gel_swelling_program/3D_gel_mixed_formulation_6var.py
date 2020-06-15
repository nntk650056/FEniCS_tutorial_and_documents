
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
#mesh = RectangleMesh(Point(0.0,0.0),Point(20.0,20.0),nx,ny)
mesh = Mesh("3D_gel_size0_5_old.xml")

#define global parameters and initial parameters

C10 = 0.2
C20 = 0.001
C1_final = 0.001
C2_final = 1
J0 = 1+C10+C20
lambda0 = (J0)**(1/3)
J_final = 1+C1_final+C2_final
lambda_final = J_final**(1/3)
Nv = 0.01
chi = 0.2
chi1p = 0
chi2p = 0
chi12 = 0
mu10 = Nv*(1/lambda0-1/J0)+ln(C10/J0)+(1+chi1p)/J0-(chi1p*C10+chi2p*C20)/J0**2+chi12*C20/J0-chi12*C10*C20/J0**2
mu20 = Nv*(1/lambda0-1/J0)+ln(C20/J0)+(1+chi2p)/J0-(chi1p*C10+chi2p*C20)/J0**2+chi12*C10/J0-chi12*C10*C20/J0**2
mu1_final = Nv*(1/lambda_final-1/J_final)+ln(C1_final/J_final)+(1+chi1p)/J_final-(chi1p*C1_final+chi2p*C2_final)/J_final**2+chi12*C2_final/J_final-chi12*C1_final*C2_final/J_final**2
mu2_final = Nv*(1/lambda_final-1/J_final)+ln(C2_final/J_final)+(1+chi1p)/J_final-(chi1p*C1_final+chi2p*C2_final)/J_final**2+chi12*C2_final/J_final-chi12*C1_final*C2_final/J_final**2
v_const = 0.588
vtotal = -(3*v_const**2-2.5*v_const+0.5)
phi0v = (1/(C10+C20+1))**vtotal
print(phi0v)

vis = 0.01
D1 = 0.005   #diffusion constant
D2 = 0.1
# Create mesh and define function space

# plot(mesh,title='mesh')
# plt.show()

# Create mesh and define function space
V = VectorElement('P',tetrahedron,1)
Q = FiniteElement('P',tetrahedron,1)
TH = MixedElement([V, Q, Q, Q, Q, Q])
V1 = FunctionSpace(mesh, V)
Q1 = FunctionSpace(mesh, Q)
W = FunctionSpace(mesh,TH)
# unknown fields: u, Pi, Cs1, mu1, Cs2, mu2

# Mark boundary subdomians
# left =  CompiledSubDomain("near(x[2], side) && on_boundary", side = 0.0)
# right = CompiledSubDomain("near(x[2], side) && on_boundary", side = 10.0)
# #up1 = CompiledSubDomain(" (x[0]*x[0]+x[1]*x[1])>=0.2499 && (x[0]*x[0]+x[1]*x[1])<=0.2501 ", side = 0.25)
# up1 = CompiledSubDomain(" x[2]>1E-10 && x[0]>1E-10 && x[1]>0.25+1E-10 && x[2]< 10-(1E-10) && on_boundary", side = 0.0)

# left_side = CompiledSubDomain("x[2]<=side1 && near(x[1], side2) && on_boundary", side1 = 0.3, side2 = 0.25)
# bottom = CompiledSubDomain("near(x[1], side) && on_boundary", side = 0.25)
# front = CompiledSubDomain("near(x[0], side) && on_boundary", side = 0.0)


left =  CompiledSubDomain("near(x[0], side) && on_boundary", side = 0.0)
right = CompiledSubDomain("near(x[0], side) && on_boundary", side = 50.0)
up1 = CompiledSubDomain(" x[0]>1E-10 && x[1]>1E-10 && x[2]>1E-10 && x[0]< 48-(1E-10) && on_boundary", side = 0.0)
up2 = CompiledSubDomain("near(x[1], side1) && near(x[2], side2) && on_boundary", side1 = 0.0, side2 = 4.0)
bottom = CompiledSubDomain("near(x[2], side) && on_boundary", side = 0.0)
front = CompiledSubDomain("near(x[1], side) && on_boundary", side = 0.0) 
left_side = CompiledSubDomain("x[0]<=side1 && near(x[2], side2) && on_boundary", side1 = 1, side2 = 0.0)

# Define chemical potential boundary conditions
C1_bc_RU = Expression('C10+(C1_final-C10)*t/100', degree=1, C10 = C10,C1_final = C1_final,t=0)
C2_bc_RU = Expression('C20+(C2_final-C20)*t/100', degree=1, C20 = C20,C2_final = C2_final,t=0)


J_bc = Expression('1+C1_bc_RU+C2_bc_RU', degree=1, C1_bc_RU = C1_bc_RU, C2_bc_RU = C2_bc_RU)

#Nv*(1/J_bc^(1/3)-1/J_bc)+log(C1_bc/J_bc)+(1+chi1p)/J_bc-(chi1p*C1_bc+chi2p*C2_bc)/J_bc^2+chi12*C2_bc/J_bc-chi12*C1_bc*C2_bc/J_bc^2
#Nv*(1/J_bc^(1/3)-1/J_bc)+log(C2_bc/J_bc)+(1+chi2p)/J_bc-(chi1p*C1_bc+chi2p*C2_bc)/J_bc^2+chi12*C1_bc/J_bc-chi12*C1_bc*C2_bc/J_bc^2

mu1_bc_RU = Expression('Nv*(1/pow(J_bc,1/3)-1/J_bc) + log10(C1_bc/J_bc)/log10(2.71828) + 1/J_bc', degree=1, Nv = Nv, J_bc = J_bc, C1_bc = C1_bc_RU, C2_bc = C2_bc_RU, chi1p = chi1p)
mu2_bc_RU = Expression('Nv*(1/pow(J_bc,1/3)-1/J_bc) + log10(C2_bc/J_bc)/log10(2.71828) + 1/J_bc', degree=1, Nv = Nv, J_bc = J_bc, C1_bc = C1_bc_RU, C2_bc = C2_bc_RU, chi2p = chi2p)

#mu1_D = Expression("mu00+(-mu_final-mu00)*(t/1000.0)", degree = 1, mu00 = mu10, mu_final=mu1_final)
#mu2_D = Expression("mu00+(-mu_final-mu00)*(t/100.0)", degree = 1, mu00 = mu20, mu_final=mu2_final)

#fixed left and bottom sides and Define chemical potential boundary conditions

bcl1 = DirichletBC(W.sub(0).sub(2), Constant(0), left_side)
bcl2 = DirichletBC(W.sub(0).sub(0), Constant(0), left)
bcf = DirichletBC(W.sub(0).sub(1), Constant(0), front)
bcb_mu1 = DirichletBC(W.sub(3), mu1_bc_RU, bottom)
bcb_mu2 = DirichletBC(W.sub(5), mu2_bc_RU, bottom)
bcu_mu1 = DirichletBC(W.sub(3), mu1_bc_RU, up1)
bcu_mu2 = DirichletBC(W.sub(5), mu2_bc_RU, up1)
bcs = [bcl1,bcl2,bcf,bcb_mu1,bcb_mu2,bcu_mu1,bcu_mu2]
#bcs = [bcl1,bcl2,bcf,bcb_mu1,bcb_mu2]


# Define test functions
w = Function(W)
u, Pi, Cs1, mu1, Cs2, mu2 = split(w)
wt = TestFunction(W)
v, Pi_, Cs1_, mu1_, Cs2_, mu2_ = split(wt)
dw = TrialFunction(W)
w_n  = Function(W)
u_n, Pi_n, Cs1_n, mu1_n, Cs2_n, mu2_n  = split(w_n)


# Define deformation gradient tensor
F = grad(u) + Identity(len(u))
I_ = Identity(len(u))
Fv = grad(v)

H = inv(F)
C = F.T*F                   # Right Cauchy-Green tensor
C_2 = C*C
# Invariants of deformation tensors
Ic = tr(C)
Ic2 = 0.5*(tr(C)*tr(C)-tr(C_2))
J  = det(F)
J_v = det(F)*H*Fv

# Save solution in VTK format
vtkfile_u = File("3D_two_test5/displacement.pvd");
vtkfile_C1 = File("3D_two_test5/C1.pvd");
vtkfile_C2 = File("3D_two_test5/C2.pvd");
vtkfile_mu1 = File("3D_two_test5/chemical_potential1.pvd");
vtkfile_mu2 = File("3D_two_test5/chemical_potential2.pvd");
vtkfile_pi = File("3D_two_test5/osmotic_pressure.pvd");
vtkfile_stress = File("3D_two_test5/cauchy_stress.pvd");

#initial chemical potential value in the domian
e_u0 = Expression(('0', '0','0'), degree=1)
e_Pi0 = Expression('0', degree=1)
e_C1 = Expression('C0', C0 = C10, degree = 1)
e_mu1 = Expression('mu0', mu0=mu10, degree=1)
e_C2 = Expression('C0', C0 = C20, degree = 1)
e_mu2 = Expression('mu0', mu0=mu20, degree=1)
assign(w, [interpolate(e_u0, V1), interpolate(e_Pi0, Q1), interpolate(e_C1, Q1), interpolate(e_mu1, Q1), interpolate(e_C2, Q1), interpolate(e_mu2, Q1)])



cs1 = Cs1/(det(F)*lambda0**3)
cs2 = Cs2/(det(F)*lambda0**3)

dt = Expression('dt1', degree=0, dt1=0.01)
gradv = grad((u-u_n)/dt)
dv2dt = gradv.T*gradv
C1t = (Cs1-Cs1_n)/dt
C2t = (Cs2-Cs2_n)/dt
W_FC = 0.5*Nv*(lambda0**2*Ic-3-2*ln(lambda0**3*det(F))) + Cs1*ln(cs1)+Cs2*ln(cs2)+chi1p*cs1+chi2p*cs2+chi12*Cs1*cs2
#W_vis = 10*vis*det(F)*(0.5*inner(gradv*H, (gradv*H+H.T*gradv.T)) - 1/3*tr(gradv*H)*tr(gradv*H))
Pi0 = -W_FC #- W_vis*dt
IMCOMPRE = Pi*(det(F)-(1+Cs1+Cs2)/J0)
a0 = (derivative(Pi0, w, wt) + mu1*Cs1_ + mu2*Cs2_)/J0*dx + derivative(IMCOMPRE, w, wt)*dx+ (-Cs1/J0*D1*dot(H.T*grad(mu1),H.T*grad(mu1_)) - C1t/J0*mu1_)*dx + (-Cs2/J0*D2*dot(H.T*grad(mu2),H.T*grad(mu2_)) - C2t/J0*mu2_)*dx
J1 = derivative(a0, w, dw)

# Newton solver parameters
parameters={"newton_solver":{"relative_tolerance":1e-4, "absolute_tolerance":1e-6, "maximum_iterations":10}}

#define the intial and Maximum time increment and allowed Minimum time increment 
ini_time_int = 0.001
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
        if t >=96.0:
            C1_bc_RU.t = 96
            C2_bc_RU.t = 96
        else:
            C1_bc_RU.t  = t
            C2_bc_RU.t  = t
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
        
        if t >=96.0:
            C1_bc_RU.t = 96
            C2_bc_RU.t = 96
        else:
            C1_bc_RU.t  = t
            C2_bc_RU.t  = t
        
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
        
        if t >=96.0:
            C1_bc_RU.t = 96
            C2_bc_RU.t = 96
        else:
            C1_bc_RU.t  = t
            C2_bc_RU.t  = t
        
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
        
        if t >=96.0:
            C1_bc_RU.t = 96
            C2_bc_RU.t = 96
        else:
            C1_bc_RU.t  = t
            C2_bc_RU.t  = t
        
        solve(a0==0, w, bcs,  J = J1, solver_parameters=parameters)
        if dt.dt1 < 1:
            dt.dt1 = dt.dt1*2
        
        if dt.dt1 < min_time_int:
            break
        
        str1 = "time increment = " + str(dt.dt1)
        print(str1)
        t = t + dt.dt1
    
    
    







    p1 = plot(w.sub(3),title=str(t))
    plt.colorbar(p1)
    plt.show()
    p1 = plot(w.sub(0), mode = "displacement",title = 'disp t='+str(t))
    plt.colorbar(plot((w.sub(0).sub(0)**2+w.sub(0).sub(1)**2+w.sub(0).sub(2)**2)**0.5))
    plt.show()
    
    vtkfile_u << (w.sub(0), t)
    vtkfile_C1 << (w.sub(2), t)
    vtkfile_mu1 << (w.sub(3), t)
    vtkfile_C2 << (w.sub(4), t)
    vtkfile_mu2 << (w.sub(5), t)
    
    if t>=total_time:
        break
