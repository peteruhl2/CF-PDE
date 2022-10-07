"""
Attempt at CF pde system
8/28/22

TODO: Write desciption

"""

from __future__ import print_function
from fenics import *
import numpy as np
import sys

T = 10.0            # final time
num_steps = 500    # number of time steps
dt = T / num_steps # time step size
k = Constant(dt)

# CF params
beta = 2.2
b = 1.
n = 1.0
dc = 9e-5
df = 9e-1
q = 4.48e-0
mu = 0.01312
_lambda = 0.01312
eta = 10.54
Dc = 1.32e-8
Df = 1.32e-2
Dw = 1.32e-3

beta = Constant(beta)
b = Constant(b)
n = Constant(n)
dc = Constant(dc)
df = Constant(df)
q = Constant(q)
mu = Constant(mu)
_lambda = Constant(_lambda)
eta = Constant(eta)
Dc = Constant(Dc)
Df = Constant(Df)
Dw = Constant(Dw)

# Read mesh from file
# mesh = Mesh('navier_stokes_cylinder/cylinder.xml.gz')
L = 0.5
nx = ny = 20
mesh = RectangleMesh(Point(-L, -L), Point(L, L), nx, ny)

# Define function space for system of concentrations
P1 = FiniteElement('P', triangle, 1)
element = MixedElement([P1, P1, P1])
V = FunctionSpace(mesh, element)

# Define test functions
v_1, v_2, v_3 = TestFunctions(V)

# Define functions for velocity and concentrations
u = Function(V)
u_n = Function(V)

# Guasian ICs for C and F
u_0 = Expression(('0.1*exp(-100*pow(x[0]-0.4, 2) - 100*pow(x[1]-0.5, 2))','0.4*exp(-100*pow(x[0]- L/2, 2) - 100*pow(x[1]-0.6, 2))','0.1'), degree = 2, L=L)
u_n = interpolate(u_0, V)

# # Constant initial conditions for checking against ODE
# u_0 = Expression(('0.4','0.3','0.1'), degree = 2)
# u_n = interpolate(u_0, V)

# Split system functions to access components
u_1, u_2, u_3 = split(u)
u_n1, u_n2, u_n3 = split(u_n)

# Define source terms
f_1 = Expression('pow(x[0]-0.4,2)+pow(x[1]-0.4,2)<0.05*0.05 ? 0.5 : 0',
                 degree=1)
f_2 = Expression('pow(x[0]-0.1,2)+pow(x[1]-0.3,2)<0.05*0.05 ? 0.5 : 0',
                 degree=1)
f_3 = Constant(0)


# Define variational problem
F = ((u_1 - u_n1) / k)*v_1*dx  \
  + Dc*dot(grad(u_1), grad(v_1))*dx - (beta/( b + u_3 ))*u_3*u_1*(1-u_1-u_2)*v_1*dx + dc*u_1*v_1*dx  \
  + ((u_2 - u_n2) / k)*v_2*dx  \
  + Df*dot(grad(u_2), grad(v_2))*dx - beta*(1 - u_3/( b + u_3 ))*u_2*(1-u_1-u_2)*v_2*dx + df*u_2*v_2*dx + q*u_2*u_3*v_2*dx  \
  + ((u_3 - u_n3) / k)*v_3*dx  \
  + Dw*dot(grad(u_3), grad(v_3))*dx - _lambda*v_3*dx + eta*u_1*u_3*v_3*dx + mu*u_3*v_3*dx

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
    solve(F == 0, u)

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

outfile = open('cf_sys/u out.txt','w')

x = np.linspace(-L,L,100)
y = np.linspace(-L,L,100)
X,Y = np.meshgrid(x,y)
results = np.zeros((len(x),len(y)))

### note that the solution funciton u can be called
print("Writing anaerobe function to file...")
for i in range(len(x)):
      for j in range(len(y)):
            results[i,j] = u_2( [x[i], y[j]])
            
results = np.reshape(results,((len(x)*len(y))))

for p in results:
      print(p, file = outfile)
      
outfile.close()


print("Argv[1] is ", sys.argv[1])
print("Type of argv[1] is ", type(sys.argv[1]))
# print(len(sys.argv))