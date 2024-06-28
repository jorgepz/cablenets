# cablenets: A Python Module for Structural Cable Nets Analysis

[![codecov](https://codecov.io/gh/jorgepz/cablenets/graph/badge.svg?token=4HE7F6GB1Y)](https://codecov.io/gh/jorgepz/cablenets)

_cablenets_ is a tiny python module for solving structural cable nets analysis problems using [Convex Optimization](https://en.wikipedia.org/wiki/Convex_optimization) formulations, in particular [CVXOPT](https://cvxopt.org/) Python library.

The solver implements an optimization formulation presented in [this paper](https://doi.org/10.1016/S0020-7683(03)00215-4), published in the International Journal of Solids and Structures.

![image saddle net](https://github.com/jorgepz/cablenets/blob/main/docs/assets/docs/assets/saddle_net.png?raw=true)

## Installing :crossed_fingers:

You can install `cablenets` and its dependencies using `pip` with this command:
```
pip install git+https://github.com/jorgepz/cablenets.git#egg=cablenets
```

## Using :muscle:

Include the `import cablenets` in your code. By now the interface consists in two functions: `solve` and `plot`. You are encouraged to explore the examples folder to see how to use `cablenets`.

The function `solve` receives numpy arrays with: the reference nodes coordinates, a connectivity matrix (see examples), the young moduli, the cross-section areas, a matrix with the imposed-deformed positions of certain nodes and a matrix with the applied forces.
```python
nodes_def, normal_forces = solve( nodes, connec, youngs, areas, disp_mat, fext_mat )
```
The output consists in: the matrix of the deformed nodes and a vector with the normal forces of the elements.

The function `plot` receives the nodes, the connectivity, the deformed positions of the nodes and the normal forces.
```python
plot( nodes, connec, nodes_def, normal_forces )
```

## Connecting :call_me_hand:

Email: jorgepz _AT_ fing edu uy

## Contributing :handshake:

Contributions are welcome. Open an issue and let's discuss there. The [JOSS authorship guidelines](https://joss.readthedocs.io/en/latest/submitting.html#authorship) are considered as a starting point.