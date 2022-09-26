#!/bin/bash
# Execute this file to recompile locally
c++ -Wall -shared -fPIC -std=c++11 -O3 -fno-math-errno -fno-trapping-math -ffinite-math-only -I/Users/peteruhl/opt/anaconda3/envs/fenicsproject/include -I/Users/peteruhl/opt/anaconda3/envs/fenicsproject/include/eigen3 -I/Users/peteruhl/opt/anaconda3/envs/fenicsproject/.cache/dijitso/include dolfin_expression_9ba4f92c645d99df7459f1972a666a1d.cpp -L/Users/peteruhl/opt/anaconda3/envs/fenicsproject/lib -L/Users/peteruhl/opt/anaconda3/envs/fenicsproject/Users/peteruhl/opt/anaconda3/envs/fenicsproject/lib -L/Users/peteruhl/opt/anaconda3/envs/fenicsproject/.cache/dijitso/lib -Wl,-rpath,/Users/peteruhl/opt/anaconda3/envs/fenicsproject/.cache/dijitso/lib -lpmpi -lmpi -lmpicxx -lpetsc -lslepc -lhdf5 -lboost_timer -ldolfin -Wl,-install_name,/Users/peteruhl/opt/anaconda3/envs/fenicsproject/.cache/dijitso/lib/libdijitso-dolfin_expression_9ba4f92c645d99df7459f1972a666a1d.so -olibdijitso-dolfin_expression_9ba4f92c645d99df7459f1972a666a1d.so