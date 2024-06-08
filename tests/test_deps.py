import sys
sys.path.append('./src')

import numpy as np
import matplotlib as mpl

import cablenets

def test_answer():
    b = cablenets.solve(2)
    assert b == 3
