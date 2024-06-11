'''
MIT License

Copyright (c) 2024 Jorge Pérez Zerpa

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
'''

import sys
from os.path import dirname
sys.path.append(dirname('../src/'))

import numpy as np
from math import pi

from cablenets import solve, plot

# scalar parameters
L      = 2 
youngs = [ 20 ]
areas  = [ pi*.01**2 ]
nelems = 100

nodes  = np.zeros((nelems+1,3))
for i in range(nelems+1):
    nodes[i,:] = [i*L/nelems, 0.0, 0.0]

connec = np.zeros((nelems,2),dtype=int)
for i in range(nelems):
    connec[i,:] = [i, i+1 ]

ks_vec = np.ones((nelems))* k

disp_mat = np.array([ [          0, 0   , 0, 0   ],
                      [   nelems/2, L*.5, 0, L*.05],
                      [     nelems, L   , 0, 0   ]])

fext_mat  = np.zeros((nelems-1,4))
for i in range(nelems-1):
    fext_mat[i,:] = [ i+1, 0, 0.0, -1.0e-3] # node fx fy fz

nodes_def, normal_forces = solve( nodes, connec, youngs, areas, disp_mat, fext_mat )

plot( nodes, connec, nodes_def, normal_forces )

