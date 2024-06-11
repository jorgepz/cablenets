# cablenets: A Python Module for Structural Cable Nets Analysis

[![codecov](https://codecov.io/gh/jorgepz/cablenets/graph/badge.svg?token=4HE7F6GB1Y)](https://codecov.io/gh/jorgepz/cablenets)

_cablenets_ is a tiny python module for solving structural cable nets analysis problems using [Convex Optimization](https://en.wikipedia.org/wiki/Convex_optimization) formulations, in particular [CVXOPT](https://cvxopt.org/) Python library.


## Installing

```
pip install git+https://github.com/jorgepz/cablenets.git#egg=cablenets
```

## Using

By now the interface consists in two simple functions: `solve` and `plot`.

```python
nodes_def, normal_forces = solve( nodes, connec, youngs, areas, disp_mat, fext_mat )
```

