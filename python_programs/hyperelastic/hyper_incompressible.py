import matplotlib.pyplot as plt
from fenics import *
#from mshr import *
import time
import numpy as np
import copy as cp
T = 1.0
num_steps = 10
dt = T/num_steps
# Optimization options for the form compiler
parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["representation"] = "uflacs"


# Create mesh and define function space
nx = nz = 5
ny = 100
#domain = Box(Point(0.0,0.0,0.0),Point(5.0,100.0,5.0))
#mesh = generate_mesh(domain, 16)

mesh = BoxMesh(Point(0.0,0.0,0.0),Point(5.0,100.0,5.0),nx,ny,nz)
V = VectorElement('P', tetrahedron, 2)
Q = FiniteElement('P', tetrahedron, 1)
#V = VectorElement('P', tetrahedron, 2)
#Q = FiniteElement('P', tetrahedron, 1)
TH = V*Q
W = FunctionSpace(mesh,TH)
#print(W)

# Mark boundary subdomians
left =  CompiledSubDomain("near(x[1], side) && on_boundary", side = 0.0)
right = CompiledSubDomain("near(x[1], side) && on_boundary", side = 100.0)
top = CompiledSubDomain("near(x[2], side) && on_boundary", side = 5.0)
bottom = CompiledSubDomain("near(x[2], side) && on_boundary", side = 0.0)
front = CompiledSubDomain("near(x[0], side) && on_boundary", side = 5.0)
behind = CompiledSubDomain("near(x[0], side) && on_boundary", side = 0.0)

# Define boundary conditions
u1 = Expression("300*t", degree = 1,t = 0)


bcl = DirichletBC(W.sub(0).sub(1), Constant(0),left)
bcr = DirichletBC(W.sub(0).sub(1), u1,right)
bcb = DirichletBC(W.sub(0).sub(2), Constant(0), bottom)
bcd = DirichletBC(W.sub(0).sub(0), Constant(0), behind)
bcs = [bcl, bcr,bcb,bcd]

# Define test functions
w = Function(W)
u, p = split(w)
(v, q) = TestFunctions(W)
wt = TestFunction(W)
dw = TrialFunction(W)

# Define deformation gradient tensor
F = grad(u) + Identity(len(u))
Fv = grad(v)


B  = Constant((0.0, 0.0, 0.0))  # Body force per unit volume
T  = Constant((0.0, 0.0, 0.0))  # Traction force on the boundary



"""
# Define functions
du = TrialFunction(W)            # Incremental displacement
v  = TestFunction(W)             # Test function
u  = Function(W)            # Displacement from previous iteration

u_k = interpolate(u,W)#previous (known) u
B  = Constant((0.0, 0.0, 0.0))  # Body force per unit volume
T  = Constant((0.0, 0.0, 0.0))  # Traction force on the boundary


# Kinematics
d = len(u)
I = Identity(d)             # Identity tensor
F = I + grad(u)             # Deformation gradient
"""
C = F.T*F                   # Right Cauchy-Green tensor
C_2 = C*C
# Invariants of deformation tensors
Ic = tr(C)
Ic2 = 0.5*(tr(C)*tr(C)-tr(C_2))
J  = det(F)

# Elasticity parameters
E, nu = 10.0, 0.3
mu, lmbda = Constant(E/(2*(1 + nu))), Constant(E*nu/((1 + nu)*(1 - 2*nu)))
#Neo-Hookean model
Gn1 = 2.290
#Mooney-Rivlin model
Cm1 = 1.030
Cm2 = 0.114
#Yeoh model
Cy1 = 1.202
Cy2 = -0.057
Cy2 = 0.004
#Gent model
Gg = 2.290
Jm = 30
#Arruda-Boyce model
Ca1 = 2.078
lambda_m = 2.8



# Define the weak form (incompressible neo-Hookean)
#a = mu*inner(F, Fv)*dx - p*det(F)*inner(inv(F).T, Fv)*dx - (det(F)-1)*q*dx

# Stored strain energy density (incompressible neo-Hookean model)
psi = (mu/2)*(tr(F.T*F) - 3) - p*(det(F)-1)
#psi = (mu/2)*(Ic - 3) - p*(det(F)-1)
# Total potential energy
Pi = psi*dx - dot(B, u)*dx - dot(T, u)*ds

# Stored strain energy density (incompressible Mooney-Rivlin model)
#a = 2*Cm1*inner(F, Fv)*dx -2*Cm2*inner(F, Fv)*dx - p*det(F)*inner(inv(F).T, Fv)*dx - (det(F)-1)*q*dx  #  a = derivative(Pi, w, wt)
#psi = Cm1*(Ic - 3) + Cm2*(Ic2 - 3) - p*(det(F)-1)
#Pi = psi*dx - dot(B, u)*dx - dot(T, u)*ds


# Stored strain energy density (incompressible Gent model)
#psi = -0.5*Gg*Jm*ln(1-(Ic-3)/Jm) - p*(det(F)-1)
#Pi = psi*dx - dot(B, u)*dx - dot(T, u)*ds
# Stored strain energy density (incompressible Yeoh model)
#psi = Cy1*(Ic - 3)+Cy2*(Ic - 3)*(Ic - 3)+Cy3*(Ic - 3)*(Ic - 3)*(Ic - 3) - p*(det(F)-1)
#Pi = psi*dx - dot(B, u)*dx - dot(T, u)*ds
# Stored strain energy density (incompressible Arruda-Boyce model)
#psi = Ca1*(0.5*(Ic-3) + 1/(20*lambda_m**2)*(Ic*Ic-9) + 11/(1050*lambda_m**4)*(Ic**3-27) + 19/(7000*lambda_m**6)*(Ic**4-81) + 519/(673750*lambda_m**8)*(Ic**5-243)) - p*(det(F)-1)
#Pi = psi*dx - dot(B, u)*dx - dot(T, u)*ds


a = derivative(Pi, w, wt)
J0 = derivative(a, w, dw)


# Save solution in VTK format
vtkfile = File("hyper_incompressible/displacement0107.pvd");
#file << u;

# Plot solution
#plot(u)
#plot(mesh)
#plt.show()

eps = 1.0
tol = 1.0E-5
iter = 0
maxiter = 10
t = 0.0
print(t)
for n in range(num_steps):

    # Update current time
    t += dt
    u1.t = t
    print(t)
    solve(a == 0, w, bcs,J=J0)
    # Compute solution
    #while eps>tol and iter< maxiter:
    #    iter +=1
    #    solve(F == 0, u, bcs, J=J)
    #    diff = u.vector().array() - u_k.vector().array()
    #    eps = np.linalg.norm(diff, ord=np.Inf)
    #    print("iter=%d: norm=%g"%(iter,eps))
    #    u_k.assign(u)
    # Save to file and plot solution
    #print(u)
    #print("w")
    
    vtkfile << (w.sub(0), t)
    
    plot(u)
    #plt.colorbar(fig)
    plt.show()
    

    # Update previous solution
    #u_n.assign(u)

# Hold plot
#interactive()
