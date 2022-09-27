### post processing for CF pde fenics data
### started 9/26/22

using Plots, CSV, DataFrames
using DelimitedFiles

mesh = readdlm("cf_sys/meshout.csv")
udat = readdlm("cf_sys/u out.csv")