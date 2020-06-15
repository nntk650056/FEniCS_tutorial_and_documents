import matplotlib.pyplot as plt
from fenics import *
#from mshr import *
import time
import numpy as np
import copy as cp
T = 1.0
num_steps = 50
dt = T/num_steps
# Optimization options for the form compiler
#parameters["form_compiler"]["cpp_optimize"] = True
#parameters["form_compiler"]["representation"] = "uflacs"


# Create mesh and define function space
nx = 5
ny = 5
#domain = Box(Point(0.0,0.0,0.0),Point(5.0,100.0,5.0))
#mesh = generate_mesh(domain, 16)
#V = VectorElement('P', tetrahedron, 3)
#Q = FiniteElement('P', tetrahedron, 2)
#TH = V*Q
#W = FunctionSpace(mesh, TH)
#V = VectorFunctionSpace(mesh, "Lagrange", 1)
#W = FunctionSpace(mesh, V)
mesh = RectangleMesh(Point(0.0,0.0),Point(10.0,10.0),nx,ny)
V = VectorElement('P',triangle,1)
Q = FiniteElement('P',triangle,1)
#V = VectorElement('P', tetrahedron, 2)
#Q = FiniteElement('P', tetrahedron, 1)
TH = V*Q
W = FunctionSpace(mesh,TH)
print(W)

# Mark boundary subdomians
left =  CompiledSubDomain("near(x[0], side) && on_boundary", side = 0.0)
right = CompiledSubDomain("near(x[0], side) && on_boundary", side = 10.0)
top = CompiledSubDomain("near(x[1], side) && on_boundary", side = 10.0)
bottom = CompiledSubDomain("near(x[1], side) && on_boundary", side = 0.0)


# Define boundary conditions
u1 = Expression("10*t", degree = 1,t = 0)


bcl = DirichletBC(W.sub(0).sub(0), Constant(0),left)
bcr = DirichletBC(W.sub(0).sub(0), u1,right)
bcb = DirichletBC(W.sub(0).sub(1), Constant(0), bottom)
bcs = [bcl, bcr,bcb]

# Define test functions
w = Function(W)
u, p = split(w)
(v, q) = TestFunctions(W)
wt = TestFunction(W)
dw = TrialFunction(W)

# Define deformation gradient tensor
F = grad(u) + Identity(len(u))
Fv = grad(v)


B  = Constant((0.0, 0.0))  # Body force per unit volume
T  = Constant((0.0, 0.0))  # Traction force on the boundary



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
E, nu = 10.0, 0.49
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
print('a')
print(a)

J = derivative(a, w, dw)
print('J')
print(J)
# Stored strain energy density (compressible neo-Hookean model)
#psi = (mu/2)*(Ic - 3) - mu*ln(J) + (lmbda/2)*(ln(J))**2

# Total potential energy
#Pi = psi*dx - dot(B, u)*dx - dot(T, u)*ds

# Stored strain energy density (compressible neo-Hookean model)
#psi = (mu/2)*(Ic - 3) - mu*ln(J) + (lmbda/2)*(ln(J))**2

# Total potential energy
#Pi = psi*dx - dot(B, u)*dx - dot(T, u)*ds

# Stored strain energy density (incompressible neo-Hookean model)
#psi = Gn1*(Ic - 3)
#Pi = psi*dx - dot(B, u)*dx - dot(T, u)*ds

# Stored strain energy density (incompressible Mooney-Rivlin model)
#psi = Cm1*(Ic - 3) + Cm2*(Ic2 - 3)
#Pi = psi*dx - dot(B, u)*dx - dot(T, u)*ds
# Stored strain energy density (incompressible Gent model)
#psi = -0.5*Gg*Jm*ln(1-(Ic-3)/Jm)
#Pi = psi*dx - dot(B, u)*dx - dot(T, u)*ds
# Stored strain energy density (incompressible Yeoh model)
#psi = Cy1*(Ic - 3)+Cy2*(Ic - 3)*(Ic - 3)+Cy3*(Ic - 3)*(Ic - 3)*(Ic - 3)
#Pi = psi*dx - dot(B, u)*dx - dot(T, u)*ds
# Stored strain energy density (incompressible Arruda-Boyce model)
#psi = Ca1*(0.5*(Ic-3) + 1/(20*lambda_m**2)*(Ic*Ic-9) + 11/(1050*lambda_m**4)*(Ic**3-27) + 19/(7000*lambda_m**6)*(Ic**4-81) + 519/(673750*lambda_m**8)*(Ic**5-243))
#Pi = psi*dx - dot(B, u)*dx - dot(T, u)*ds
# Stored strain energy density (incompressible Ogden model)
#psi = (mu/2)*(Ic - 3)


# Compute first variation of Pi (directional derivative about u in the direction of v)
#F = derivative(Pi, u, v)

# Compute Jacobian of F
#J = derivative(F, u_k, du)

# The complete variational problem can now be solved by a single call to
# :py:func:`solve <dolfin.fem.solving.solve>`::

# Solve variational problem
#solve(F == 0, u, bcs, J=J)

# Finally, the solution ``u`` is saved to a file named
# ``displacement.pvd`` in VTK format, and the deformed mesh is plotted
# to the screen::

# Save solution in VTK format
vtkfile = File("hyperiter2D_NH/displacementNH.pvd");
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
    solve(a == 0, w, bcs, J=J)
    #plot(w.sub(0),title='disp'+str(n+1))
    #plt.show()
    #plot(w.sub(1),title='imp'+str(n+1))
    #plt.show()
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

    # Update previous solution
    #u_n.assign(u)

# Hold plot
#interactive()
