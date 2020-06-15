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
mesh = RectangleMesh(Point(0.0,0.0),Point(20.0,20.0),nx,ny)  
#define global parameters and initial parameters  
C0 = 0.2  
J0 = 1+C0
lambda0 = (1+C0)**(1/3)  
Nv = 0.001  
chi = 0.2  
mu0 = ln(C0/(1+C0))+1/(1+C0)+chi/(1+C0)**2+Nv*(1/lambda0-1/(1+C0)) # mu0 = -0.81943   
vis = 0.01  
# Create mesh and define function space  
V = VectorElement('P',triangle,2)  
Q = FiniteElement('P',triangle,1)  
TH = MixedElement([V, Q, Q, Q]) 
V1 = FunctionSpace(mesh, V)  
Q1 = FunctionSpace(mesh, Q)  
W = FunctionSpace(mesh,TH)  
# Mark boundary subdomians  
left =  CompiledSubDomain("near(x[0], side) && on_boundary", side = 0.0)  
right = CompiledSubDomain("near(x[0], side) && on_boundary", side = 20.0)  
top = CompiledSubDomain("near(x[1], side) && on_boundary", side = 20.0)  
bottom = CompiledSubDomain("near(x[1], side) && on_boundary", side = 0.0)  
# Define chemical potential boundary conditions  
mu_D = Expression("mu00*(1-t/100.0)", degree = 1, mu00 = mu0, t = 0)  
#fixed left and bottom sides and Define chemical potential boundary conditions  
bcl = DirichletBC(W.sub(0).sub(0), Constant(0), left)  
bcb = DirichletBC(W.sub(0).sub(1), Constant(0), bottom)  
bct_mu = DirichletBC(W.sub(3), mu_D,top)  
bcr_mu = DirichletBC(W.sub(3), mu_D, right)  
bcs = [bcl,bcb,bct_mu,bcr_mu]  
# Define test functions  
w = Function(W)  
u, Pi, Cs, mu = split(w)  
wt = TestFunction(W)  
v, Pi_, Cs_, mu_ = split(wt)  
dw = TrialFunction(W)  
w_n  = Function(W)  
u_n, Pi_n, Cs_n, mu_n = split(w_n)  
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
vtkfile_u = File("2D_gel_D1_t1_n5/displacement.pvd");  
vtkfile_mu = File("2D_gel_D1_t1_n5/chemical_potential.pvd");  
vtkfile_phi = File("2D_gel_D1_t1_n5/polumer_volume_fraction.pvd");  
vtkfile_stress = File("2D_gel_D1_t1_n5/cauchy_stress.pvd");  
#initial chemical potential value in the domian  
e_u0 = Expression(('0', '0'), degree=1)  
e_Pi0 = Expression('0', degree=1)
e_C0 = Expression('C0', C0 = C0, degree = 1)
e_mu0 = Expression('mu0', mu0=mu0, degree=1)  
assign(w, [interpolate(e_u0, V1), interpolate(e_Pi0, Q1),interpolate(e_C0, Q1), interpolate(e_mu0, Q1)]) 
#define govern equation 
dt = Expression('dt1', degree=0, dt1=0.1)  
Ct = (Cs-Cs_n)/dt 

W_FC= 0.5*Nv*(lambda0**2*Ic-3-2*ln(lambda0**3*det(F))) - Cs*(ln(1+1/Cs))-chi/(1+Cs) + Pi*(det(F)-(1+Cs)/J0) 
Pi0 = -W_FC
a0 =(derivative(Pi0, w, wt) + mu*Cs_)*dx +(-Cs*dot(H.T*grad(mu),H.T*grad(mu_)) - Ct*mu_)*dx  
J1 = derivative(a0, w, dw)  
# Newton solver parameters  
parameters={"newton_solver":{"relative_tolerance":1e-4, "absolute_tolerance":1e-6, "maximum_iterations":10}}  
t = 0.0  
for n in range(10+89+82):  
	    # Update previous solution  
    w_n.assign(w)  
    # Update current time range(0,0.1,1) range(2,1,90) range(90,5,500)  
    if t < 0.99:  
        dt.dt1 = 0.1  
    elif t < 89.99:  
        dt.dt1 = 1
    elif t<499.99:
        dt.dt1 = 5  
    t += dt.dt1   
    if t >89.9:  
        mu_D.t = 90.0  
    else:  
        mu_D.t = t  
    print(t)  
  
    solve(a0==0, w, bcs,  J = J1, solver_parameters=parameters)  
	  
    p1 = plot(w.sub(1),title=str(t))  
    plt.colorbar(p1)  
    plt.show()  
    p1 = plot(w.sub(0), mode = "displacement",title = 'displacement t='+str(t))  
    plt.colorbar(p1)  
    plt.show()  
    vtkfile_u << (w.sub(0), t)  
    vtkfile_mu << (w.sub(1), t) 
