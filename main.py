from FEM import *
from Optimization import *
from Plots import *
import sys
import numpy as np
from numpy.linalg import inv, det
from math import floor
import matplotlib.pyplot as plt
from scipy.sparse import csc_matrix, linalg
import gif

defaultInputs = 3

if __name__ == '__main__':
	# user input
	if len(sys.argv) == defaultInputs+1:
		nx = sys.argv[1]
		ny = sys.argv[2]
		v = sys.argv[3]
	else:
		nx = 1
		ny = 1
		v = 1

	it = 0
	change = 0
	tol = 1e-4
	while change > tol:
		pass