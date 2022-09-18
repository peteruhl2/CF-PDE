"""
Going to be an example from a paper about turing patterns
8/28/22

FEniCS tutorial demo program: Convection-diffusion-reaction for a system
describing the concentration of three species A, B, C undergoing a simple
first-order reaction A + B --> C with first-order decay of C. The velocity
is given by the flow field w from the demo navier_stokes_cylinder.py.
  u_1' + w . nabla(u_1) - div(eps*grad(u_1)) = f_1 - K*u_1*u_2
  u_2' + w . nabla(u_2) - div(eps*grad(u_2)) = f_2 - K*u_1*u_2
  u_3' + w . nabla(u_3) - div(eps*grad(u_3)) = f_3 + K*u_1*u_2 - K*u_3
"""

from __future__ import print_function
from fenics import *

T = 5.0            # final time
num_steps = 500    # number of time steps
dt = T / num_steps # time step size
eps = 0.1         # diffusion coefficient
K = 2.0           # reaction rate

Da = 0.1
Db = 0.03
rho = 0.0001
rhoa = 0.01
rhob = 0.02
ka = 0.02
kb = 0.03
mua = 0.1

a0 = 0.1
b0 = 0.05

# Read mesh from file
# mesh = Mesh('navier_stokes_cylinder/cylinder.xml.gz')
nx = ny = 50
mesh = RectangleMesh(Point(0, 0), Point(1, 1), nx, ny)

# Define function space for velocity
W = VectorFunctionSpace(mesh, 'P', 2)

# Define function space for system of concentrations
P1 = FiniteElement('P', triangle, 1)
element = MixedElement([P1, P1])
V = FunctionSpace(mesh, element)

# Define test functions
v_1, v_2 = TestFunctions(V)

# Define functions for velocity and concentrations
w = Function(W)
u = Function(V)
# u_n = Function(V)

u_0 = Expression(('1.1','1.05'), degree = 1)
u_n = interpolate(u_0, V)

# Split system functions to access components
u_1, u_2 = split(u)
u_n1, u_n2 = split(u_n)

# Define source terms
f_1 = Expression('pow(x[0]-0.1,2)+pow(x[1]-0.1,2)<0.05*0.05 ? 0.5 : 0',
                 degree=1)
f_2 = Expression('pow(x[0]-0.1,2)+pow(x[1]-0.3,2)<0.05*0.05 ? 0.5 : 0',
                 degree=1)
f_3 = Constant(0)

# Define expressions used in variational forms
k = Constant(dt)
K = Constant(K)
eps = Constant(eps)
Da = Constant(Da)
Db = Constant(Db)
rho = Constant(rho)
rhoa = Constant(rhoa)
rhob - Constant(rhob)
ka = Constant(ka)
kb = Constant(kb)
mua = Constant(mua)

# # Define variational problem
# F = ((u_1 - u_n1) / k)*v_1*dx  \
#   + Da*dot(grad(u_1), grad(v_1))*dx + rho*u_1*u_1/( u_2*( 1+mua*u_1*u_1 ) )*v_1*dx - ka*u_1*v_1*dx + rhoa*v_1*dx  \
#   + ((u_2 - u_n2) / k)*v_2*dx  \
#   + Db*dot(grad(u_2), grad(v_2))*dx + rho*u_1*u_1*v_2*dx - kb*u_2*v_2*dx + rhob*v_2*dx
  
# Define variational problem
F = ((u_1 - u_n1) / k)*v_1*dx  \
  + Da*dot(grad(u_1), grad(v_1))*dx - rhoa*v_1*dx + ka*u_1*v_1*dx - ( (u_1+u_2)/( (u_2 + u_1) ) )*v_1*dx \
  + ((u_2 - u_n2) / k)*v_2*dx  \
  + Db*dot(grad(u_2), grad(v_2))*dx - rhob*v_2*dx + kb*u_2*v_2*dx - rho*u_1*u_1*v_2*dx


# Create VTK files for visualization output
vtkfile_u_1 = File('turing_ex/u_1.pvd')
vtkfile_u_2 = File('turing_ex/u_2.pvd')

# Create progress bar
# progress = Progress('Time-stepping')
# set_log_level(PROGRESS)
progress = Progress('Time-stepping', num_steps)

# Time-stepping
t = 0
for n in range(num_steps):

    # Update current time
    t += dt

    # Read velocity from file
    # timeseries_w.retrieve(w.vector(), t)

    # Solve variational problem for time step
    solve(F == 0, u)

    # Save solution to file (VTK)
    _u_1, _u_2 = u.split()
    vtkfile_u_1 << (_u_1, t)
    vtkfile_u_2 << (_u_2, t)

    # Update previous solution
    u_n.assign(u)

    # Update progress bar
    set_log_level(LogLevel.PROGRESS)
    progress += 1

# Hold plot
#interactive()