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

# example 
import sys
sys.path.append('../src')

import numpy as np
import cablenets as cn

Lx = 2
Ly = 2
k = 2

nodes  = np.array([ [  0,  0, 0 ],
                    [ Lx,  0, 0 ],
                    [  0, Ly, 0 ] ])

connec = np.array([ [0,1],
                    [1,2],
                    [2,0] ])

print("nodes", nodes)
print("conec", connec)

ks_vec = np.ones(np.shape(connec)[0])* k

disp_mat = np.array([ [0, 0,      0 , 0],
                      [1, Lx*1.1, 0  , 0] ])

fext_mat  = np.array([ [2, 1,1,0 ]]) # node fx fy fz

print("f", fext_mat)

nodes_def, normal_forces = cn.solve( nodes, connec, ks_vec, disp_mat, fext_mat )
print("normal forces", normal_forces)

cn.plot( nodes, connec, nodes_def, normal_forces )
