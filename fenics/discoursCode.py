from fenics import *
import matplotlib.pyplot as plt
#final time
T = 85
#number of time steps
num_steps = 70

#time step size
dt = T / num_steps

#Generate Mesh
x_lim = 100
y_lim = 100
mesh = RectangleMesh (Point(0,0),Point(x_lim,y_lim), x_lim, y_lim)
plot(mesh)
plt.show()

#Define function space for given system of pde
P1 = FiniteElement('P', triangle, 1)
element = MixedElement([P1, P1])
V = FunctionSpace(mesh, element)

#Define test functions
v_1, v_2 = TestFunctions(V)

#Define functions for system
u = Function(V)
u_n = Function(V)

#Define Initial Condition
u_0 = Expression (('6.69-(1.25pow(10,-18))(x[0]-1600)(x[1]-800)-(0.1875pow(10,-10))(x[0]-2)(x[1]-300)', '1.19-(0.3pow(10,-8))(x[0]-1500)-(1.5pow(10,-8))(x[1]-50)'), degree = 1)
u_n = project(u_0,V)

#Split system functions to access components
u_1, u_2 = split(u)
u_n1, u_n2= split(u_n)

#Define source terms
f_1 = Constant(0)
f_2 = Constant(0)

#Parameter values
d = 1.2
alpha = 1.30039
beta = 1.3
rho = 0.24
kappa = 42
xi = 0.110013
mu = 0.4614
k = Constant(dt)

F = ((u_1 - u_n1) / k)v_1dx -rhou_1(1-(u_1/kappa))v_1dx+((betau_1u_2)/(1+u_1))v_1dx+dot( grad(u_1),grad(v_1))*dx\
((u_2 - u_n2) / k)v_2dx - mu*(((alphau_1u_2)/(1+u_1))-u_2-xiu_2u_2)v_2dx+d*dot(grad(u_2),grad(v_2))*dx \
f_1v_1dx - f_2v_2dx

#Create VTK files for visualization output
vtkfile_u_1 = File('Pattern formation_system/u_1.pvd')
vtkfile_u_2 = File('Pattern_formation/u_2.pvd')

#Create progress bar
progress = Progress('Time-stepping')
set_log_level(LogLevel.PROGRESS)

#Time-stepping
t = 0

for n in range(num_steps):
#Update current time
t += dt

#Solve variational problem for time step
solve(F == 0, u)

print(t)
plot(u_1)
plt.show()
plot(u_2)
plt.show()

#Save solution to file (VTK)
_u_1, _u_2= u.split()
vtkfile_u_1 << (_u_1, t)
vtkfile_u_2 << (_u_2, t)

#Update previous solution
u_n.assign(u)