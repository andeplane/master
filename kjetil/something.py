from dolfin import *
parameters["linear_algebra_backend"] = "uBLAS";

import numpy as np
import time

# Set parameter values
rho = 1.0   # Assumed to 1.0, not included in the calculations
dt = 0.01   # Timestep
nu = 0.01    # Kinematic viscosity

x0 = -2.0
x1 = 14.0

PA = 1.0    # Pressure at x=0
PB = 0.0    # Pressure at x=Lx

start_time = time.time()

mesh = Mesh("meshtofile.xml")
boundaries = MeshFunction('size_t', mesh, "meshtofile.xml")

elapsed_time = time.time() - start_time

print "File read in ",elapsed_time, " seconds"

# Define function spaces (P2-P1)
V = VectorFunctionSpace(mesh, "CG", 2)
Q = FunctionSpace(mesh, "CG", 1)

# Define trial and test functions
u = TrialFunction(V) # Trial function for u
p = TrialFunction(Q) # Trial function for p

v = TestFunction(V) # Test function for u
q = TestFunction(Q) # Test function for p

# Normal vector
n = FacetNormal(mesh)

p_in = Constant(PA) # Pressure at x=0
p_out = Constant(PB) # Pressure at x=1

def u0_noslip(x,on_boundary):
    return on_boundary

# No slip at the walls
noslip = DirichletBC(V, (0,0), boundaries, 3)

# Pressure boundary conditions
inflow = DirichletBC(Q, p_in, boundaries, 1)
outflow = DirichletBC(Q, p_out, boundaries, 2)

bcu = [noslip]
bcp = [inflow, outflow]

# Create functions to save solutions
u0 = Function(V)
u1 = Function(V)
p0 = Function(Q)
p1 = Function(Q)
F  = Function(Q)

# Define coefficients
k = Constant(dt)
nu_const = Constant(nu)

f = Constant((0, 0))

U = 0.5*(u0 + u)

# Tentative velocity step
F1 = (1/k)*inner(u - u0, v)*dx \
    + nu_const*inner(nabla_grad(U), nabla_grad(v))*dx \
    - p0*div(v)*dx\
    - inner(f, v)*dx\
    + inner(dot(u0,nabla_grad(u0)), v)*dx \
    + inner(p0*n,v)*ds\

a1 = lhs(F1)
L1 = rhs(F1)

# Calculation of phi
a2 = inner(grad(p), grad(q))*dx
L2 = inner(grad(p0),grad(q)) *dx\
    -(1.0/k)*div(u1)*q*dx

# Velocity update
a3 = inner(u, v)*dx
L3 = inner(u1, v)*dx - k*inner(grad(p1 - p0), v)*dx

# Assemble matrices
A1 = assemble(a1)
A2 = assemble(a2)
A3 = assemble(a3)

# Time-stepping
t = 0.0

def step():
    global t

    # Compute tentative velocity step
    b1 = assemble(L1)
    [bc.apply(A1, b1) for bc in bcu]
    solve(A1, u1.vector(), b1, "lu", "none")
    
    # Pressure correction
    b2 = assemble(L2)
    [bc.apply(A2, b2) for bc in bcp]
    solve(A2, p1.vector(), b2, "lu", "none")
    
    # Velocity correction
    b3 = assemble(L3)
    [bc.apply(A3, b3) for bc in bcu]
    solve(A3, u1.vector(), b3, "lu", "none")
    
    ds_2 = Measure("ds")[boundaries]

    m1 = p1*n[0]*ds_2(4)
    m2 = p1*n[1]*ds_2(4)
    
    v1 = assemble(m1)
    v2 = assemble(m2)
    #print "Force = ", v1.array().max()
    print "F_x = ", v1
    print "F_y = ", v2

    
    # Move to next time step
    u0.assign(u1)
    p0.assign(p1)
    t += dt
    
    return t,u0,p1