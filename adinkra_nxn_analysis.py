# ******************************************************************************
# Name:    Adinkra NxN Pattern Analysis
# Author:  Vadim Korotkikh
# Email:   va.korotki@gmail.com
# Date:	   February 2017
#
# Description: Script for analyzing generated NxN Adinkras of k colours for
# repeating pattern

from numpy.linalg import inv
import numpy as np
import itertools
import time
import sys

def main():
	pass


# ********************************
def inverse_tally(matrix_list):
	""" Function for finding all the inverse pairs in a i sized list of
	sign permutation matrices """

	temp	= []

	for i, mat in enumerate(matrix_list):
		imat_inv	= inv(mat)
		for j, mat in enumerate(matrix_list):
			if np.array_equal(imat_inv, mat):
				temp.append((i,j)
