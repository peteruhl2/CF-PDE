"""
CF-PDE system on a ciruclar domain
Model is:

ct = beta*w/(b + w)*c*(1-c-f) - dc*c + Dc Lap(c)
ft = beta(1 - w/(b + w))*f*(1-c-f) - df*f + Df Lap(f)
wt = lambda - mu*w - eta*c*w + Dw Lap(w)

Want to have Dirichlet BC for oxygen and no-flow
for c and f

"""

from __future__ import print_function
from fenics import *
from mshr import *
from dolfin import *
import numpy as np
import sys

T = 20.0            # final time
num_steps = 1000    # number of time steps
dt = T / num_steps # time step size
k = Constant(dt)

# CF params
beta = 5.0
b = 1.
n = 1.0
dc = 0.1
df = 0.1
q = 0.5
eta = .5
Dc = 4e-6
Df = 4e-4
Dw = 0.5e-1
w_r = 1.0

beta = Constant(beta)
b = Constant(b)
n = Constant(n)
dc = Constant(dc)
df = Constant(df)
q = Constant(q)
eta = Constant(eta)
Dc = Constant(Dc)
Df = Constant(Df)
Dw = Constant(Dw)
w_r = Constant(w_r)

# Read mesh from file
# mesh = Mesh('navier_stokes_cylinder/cylinder.xml.gz')
# L = .2
L = float(sys.argv[1])
nx = ny = 40
# mesh = RectangleMesh(Point(-L, -L), Point(L, L), nx, ny)
domain = Circle(Point(0, 0), L)
mesh = generate_mesh(domain, nx)

# Define function space for system of concentrations
P1 = FiniteElement('P', triangle, 1)
element = MixedElement([P1, P1, P1])
V = FunctionSpace(mesh, element)

# Define boundary condition for oxygen
def boundary(x, on_boundary):
    return on_boundary

# v.sub(2) puts a dirichlet BC only on the
# third component (oxygen)
bc = DirichletBC(V.sub(2), w_r, boundary)

# Define test functions
v_1, v_2, v_3 = TestFunctions(V)

# Define functions for velocity and concentrations
u = Function(V)
u_n = Function(V)

# Guasian ICs for C and F
# u_0 = Expression(('0.4*exp(-1*pow(x[0], 2) - 1*pow(x[1], 2))','0.2*exp(-1*pow(x[0], 2) - 1*pow(x[1], 2))','1.0'), degree = 2, L=L)
# pow(x[0],2) + pow(x[1],2) <= L ? 0 : w_r
u_0 = Expression(('0.4*exp(-1*pow(x[0], 2) - 1*pow(x[1], 2))','0.2*exp(-1*pow(x[0], 2) - 1*pow(x[1], 2))',\
      'pow(x[0],2) + pow(x[1],2) <= pow(L,2) ? 0 : w_r'), degree = 2, L=L, w_r = w_r)

u_n = interpolate(u_0, V)

# # Constant initial conditions for checking against ODE
# u_0 = Expression(('0.4','0.3','0.1'), degree = 2)
# u_n = interpolate(u_0, V)

# Split system functions to access components
u_1, u_2, u_3 = split(u)
u_n1, u_n2, u_n3 = split(u_n)


# Define variational problem
F = ((u_1 - u_n1) / k)*v_1*dx  \
  + Dc*dot(grad(u_1), grad(v_1))*dx - (beta/( b + u_3 ))*u_3*u_1*(1-u_1-u_2)*v_1*dx + dc*u_1*v_1*dx  \
  + ((u_2 - u_n2) / k)*v_2*dx  \
  + Df*dot(grad(u_2), grad(v_2))*dx - beta*(1 - u_3/( b + u_3 ))*u_2*(1-u_1-u_2)*v_2*dx + df*u_2*v_2*dx + q*u_2*u_3*v_2*dx  \
  + ((u_3 - u_n3) / k)*v_3*dx  \
  + Dw*dot(grad(u_3), grad(v_3))*dx + eta*u_1*u_3*v_3*dx

# # Define variational problem
# F = ((u_1 - u_n1) / k)*v_1*dx  \
#   + Dc*dot(grad(u_1), grad(v_1))*dx - (pow(beta,n)/( pow(b,n) + pow(u_3,n) ))*pow(u_3,n)*u_1*(1-u_1-u_2)*v_1*dx + dc*u_1*v_1*dx  \
#   + ((u_2 - u_n2) / k)*v_2*dx  \
#   + Df*dot(grad(u_2), grad(v_2))*dx + df*u_2*v_2*dx + q*u_2*u_3*v_2*dx  - pow(beta,n)*(1 - pow(u_3,n)/( pow(b,n) + pow(u_3,n) ))*u_2*(1-u_1-u_2)*v_2*dx  \
#   + ((u_3 - u_n3) / k)*v_3*dx  \
#   + Dw*dot(grad(u_3), grad(v_3))*dx - _lambda*v_3*dx + mu*u_3*v_3*dx + eta*u_1*u_3*v_3*dx



# Create VTK files for visualization output
vtkfile_u_1 = File('cf_sys/u_1.pvd')
vtkfile_u_2 = File('cf_sys/u_2.pvd')
vtkfile_u_3 = File('cf_sys/u_3.pvd')

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
    solve(F == 0, u, bc)

    # Save solution to file (VTK)
    _u_1, _u_2, _u_3 = u.split()
    vtkfile_u_1 << (_u_1, t)
    vtkfile_u_2 << (_u_2, t)
    vtkfile_u_3 << (_u_3, t)

    # Update previous solution
    u_n.assign(u)

    # Update progress bar
    set_log_level(LogLevel.PROGRESS)
    progress += 1

# Hold plot
#interactive()

### write out function values for anaerobic population
u_1, u_2, u_3 = split(u)

# outfile = open('cf_sys/u out.txt','w')
      
# outfile.close()


print("Argv[1] is ", sys.argv[1])
print("Type of argv[1] is ", type(sys.argv[1]))
# print(len(sys.argv))