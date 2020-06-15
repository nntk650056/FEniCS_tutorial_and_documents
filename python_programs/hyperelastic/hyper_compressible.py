#import matplotlib.pyplot as plt
from dolfin import *
import time
import numpy as np
T = 1.0
num_steps = 50
dt = T/num_steps
# Optimization options for the form compiler
parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["representation"] = "uflacs"


# Create mesh and define function space
nx = nz = 5
ny = 100
mesh = BoxMesh(Point(0.0,0.0,0.0),Point(5.0,100.0,5.0),nx,ny,nz)


V = VectorFunctionSpace(mesh, "Lagrange", 1)
#print(V.sub(0))
#print(V.sub(1))

# Mark boundary subdomians
left =  CompiledSubDomain("near(x[1], side) && on_boundary", side = 0.0)
right = CompiledSubDomain("near(x[1], side) && on_boundary", side = 100.0)
top = CompiledSubDomain("near(x[2], side) && on_boundary", side = 5.0)
bottom = CompiledSubDomain("near(x[2], side) && on_boundary", side = 0.0)
front = CompiledSubDomain("near(x[0], side) && on_boundary", side = 5.0)
behind = CompiledSubDomain("near(x[0], side) && on_boundary", side = 0.0)

# Define Dirichlet boundary (x = 0 or x = 1)
# Define boundary conditions
u1 = Expression("300*t", degree = 1,t = 0)


bcl = DirichletBC(V.sub(1), Constant(0),left)
bcr = DirichletBC(V.sub(1), u1,right)
bcb = DirichletBC(V.sub(2), Constant(0), bottom)
bcd = DirichletBC(V.sub(0), Constant(0), behind)
bcs = [bcl, bcr,bcb,bcd]

# Define test functions



# Define functions
du = TrialFunction(V)            # Incremental displacement
v  = TestFunction(V)             # Test function
u  = Function(V)                 # Displacement from previous iteration
B  = Constant((0.0, 0.0, 0.0))  # Body force per unit volume
T  = Constant((0.0, 0.0, 0.0))  # Traction force on the boundary


# Kinematics
d = len(u)
I = Identity(d)             # Identity tensor
F = I + grad(u)             # Deformation gradient
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



# Stored strain energy density (compressible neo-Hookean model)
#psi = (mu/2)*(Ic - 3) - mu*ln(J) + (lmbda/2)*(ln(J))**2

# Total potential energy
#Pi = psi*dx - dot(B, u)*dx - dot(T, u)*ds

# Stored strain energy density (compressible neo-Hookean model)
psi = (mu/2)*(Ic - 3) - mu*ln(J) + (lmbda/2)*(ln(J))**2

# Total potential energy
Pi = psi*dx - dot(B, u)*dx - dot(T, u)*ds

# Stored strain energy density (incompressible neo-Hookean model)
#psi = Gn1*(Ic - 3)
#Pi = psi*dx - dot(B, u)*dx - dot(T, u)*ds

# Stored strain energy density (incompressible Mooney-Rivlin model)
#psi = Cm1*(Ic - 3) + Cm2*(Ic2 - 3)
#Pi = psi*dx - dot(B, u)*dx - dot(T, u)*ds
# Stored strain energy density (incompressible Gent model)
#psi = (mu/2)*(Ic - 3)

# Stored strain energy density (incompressible Ogden model)
#psi = (mu/2)*(Ic - 3)

# Stored strain energy density (incompressible Yeoh model)
#psi = (mu/2)*(Ic - 3)

# Stored strain energy density (incompressible Arruda-Boyce model)
#psi = Ca1*(0.5*(Ic-3) + 1/(20*lambda_m**2)*(Ic*Ic-9) + 11/(1050*lambda_m**4)*(Ic**3-27) + 19/(7000*lambda_m**6)*(Ic**4-81) + 519/(673750*lambda_m**8)*(Ic**5-243))
#Pi = psi*dx - dot(B, u)*dx - dot(T, u)*ds


# Compute first variation of Pi (directional derivative about u in the direction of v)
F = derivative(Pi, u, v)

# Compute Jacobian of F
J = derivative(F, u, du)

# The complete variational problem can now be solved by a single call to
# :py:func:`solve <dolfin.fem.solving.solve>`::

# Solve variational problem
#solve(F == 0, u, bcs, J=J)

# Finally, the solution ``u`` is saved to a file named
# ``displacement.pvd`` in VTK format, and the deformed mesh is plotted
# to the screen::

# Save solution in VTK format
vtkfile = File("hyper_compressible/displacement.pvd");
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
for n in range(50):

    # Update current time
    t += dt
    u1.t = t
    print(t)
    solve(F == 0, u, bcs, J=J)
    # Compute solution
    #"""while eps>tol and iter< maxiter:
    #    iter +=1
    #    solve(F == 0, u, bcs, J=J)
    #    diff = u.ve"""
    print(u)

    # Save to file and plot solution
    vtkfile << (u, t)
    #plot(u)

    # Update previous solution
    #u_n.assign(u)

# Hold plot
#interactive()
