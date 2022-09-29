### post processing for CF pde fenics data
### started 9/26/22

using Plots, CSV, DataFrames
using DelimitedFiles, PyPlot

mesh = readdlm("cf_sys/meshout.csv")
udat = readdlm("cf_sys/u out.txt")

### rearrange udat
# udat = reshape(udat, (3,length(mesh[:,1])))
# udat = reshape(udat, (3,:))
# cdat = udat[1,:]
# fdat = udat[2,:]
# wdat = udat[3,:]

L = 2.0
x = collect(range(-L,L,length=100))
y = collect(range(-L,L,length=100))

udat = reshape(udat,length(x),length(y))


# p = surface(udat)
p = Plots.heatmap(udat)

display(p)

### plot with pyplot
pygui(true)
fig = plt.figure()
ax = plt.contourf(udat)
# ax = fig.add_subplot(projection="3d")
# ax.plot_surface(x, y, udat)
