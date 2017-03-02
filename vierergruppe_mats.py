# ******************************************************************************
# Name:  Generate the Vierergruppe
# Author:  Vadim Korotkikh
# Email:   va.korotki@gmail.com
# Date:    February 2017
# Version:
#
# Description:
#
# ******************************************************************************
# Begin Imports

import math
import sys
import time
import numpy as np
import numpy.matlib
import itertools
from numpy import array


# ********************************
def main():

	assemble_tetrads()

# ********************************
# Defining the six Vierergruppe representations
def vierergruppe_sets():

	vp1 	= np.matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
	vp2 	= np.matrix([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]])
	vp3 	= np.matrix([[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]])
	vp4 	= np.matrix([[0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0]])
	vprime 	= [vp1, vp2, vp3, vp4]

	"""Elements for flopping bosonic fields"""
	b12 	= np.matrix([[0,1,0,0], [1,0,0,0], [0,0,1,0], [0,0,0,1]])
	b13 	= np.matrix([[0,0,1,0], [0,1,0,0], [1,0,0,0], [0,0,0,1]])
	b23 	= np.matrix([[1,0,0,0], [0,0,1,0], [0,1,0,0], [0,0,0,1]])
	b123	= np.matrix([[0,0,1,0], [1,0,0,0], [0,1,0,0], [0,0,0,1]])
	b132	= np.matrix([[0,1,0,0], [0,0,1,0], [1,0,0,0], [0,0,0,1]])

	vgruppe	= 	{'()': vprime,
				'(12)': [np.dot(b12, i) for i in vprime],
				'(13)': [np.dot(b13, i) for i in vprime],
				'(23)': [np.dot(b23, i) for i in vprime],
				'(123)':[np.dot(b123, i) for i in vprime],
				'(132)':[np.dot(b132, i) for i in vprime]
				}

	return vgruppe

# ********************************
# Defining the binary multiplication matrices
def binaries(bin_code):

	binaries_lt	= [(0, [1,1,1,1]), (2, [1,-1,1,1]), (4, [1,1,-1,1]),
					(6, [1,-1,-1,1]), (8, [1,1,1,-1]), (10, [1,-1,1,-1]),
					(12, [1,1,-1,-1]), (14, [1,-1,-1,-1])]

	for btuple in binaries_lt:
		if bin_code == btuple[0]:
			tarray = np.array(btuple[1])
			return np.diag(tarray)


# ********************************
# Defining the elle binary representations for the Vierergruppe
def vgrp_ellebin():

	vgrp_elle			= {}
	vgrp_elle['()']		= [[14,8,2,4], [2,4,14,8], [4,2,8,14],[8,14,4,2],
							[6,0,10,12], [10,12,6,0], [12,10,0,6], [0,6,12,10]]

	vgrp_elle['(12)'] 	= [[14,4,2,8], [2,8,14,4], [4,14,8,2], [8,2,4,14],
							[6,12,10,0], [10,0,6,12], [12,6,0,10], [0,10,12,6]]

	vgrp_elle['(13)'] 	= [[14,2,8,4], [2,14,4,8], [4,8,2,14], [8,4,14,2],
							[6,10,0,12], [10,6,12,0], [12,0,10,6], [0,12,6,10]]

	vgrp_elle['(23)'] 	= [[2,4,8,14], [14,8,4,2], [8,14,2,4], [4,2,14,8],
							[10,12,0,6], [6,0,12,10], [0,6,10,12], [12,10,6,0]]

	vgrp_elle['(123)']	= [[14,4,8,2], [2,8,4,14], [4,14,2,8], [8,2,14,4],
							[6,12,0,10], [10,0,12,6], [12,6,10,0], [0,10,6,12]]

	vgrp_elle['(132)']	= [[14,2,4,8], [2,14,8,4], [4,8,14,2], [8,4,2,14],
							[6,10,12,0], [10,6,0,12], [12,0,6,10], [0,12,10,6]]

	return vgrp_elle


# ********************************
# Defining the tilde-elle binary representations for the Vierergruppe
def vgrp_tildebin():

	vgrp_tilde			= {}

	vgrp_tilde['()']	= [[14,8,2,4], [2,4,14,8], [4,2,8,14], [8,14,4,2],
							[6,0,10,12], [10,12,6,0], [12,10,0,6], [0,6,12,10]]

	vgrp_tilde['(12)']	= [[14,4,2,8], [2,8,14,4], [4,14,8,2], [8,2,4,14],
							[6,12,10,0], [10,0,6,12], [12,6,0,10], [0,10,12,6]]

	vgrp_tilde['(13)']	= [[14,2,8,4], [2,14,4,8], [4,8,2,14], [8,4,14,2],
							[6,10,0,12], [10,6,12,0], [12,0,10,6], [0,12,6,10]]

	vgrp_tilde['(23)']	= [[2,4,8,14], [14,8,4,2], [8,14,2,4], [4,2,14,8],
							[10,12,0,6], [6,0,12,10], [0,6,10,12], [12,10,6,0]]

	vgrp_tilde['(123)']	= [[14,4,8,2], [2,8,4,14], [4,14,2,8], [8,2,14,4],
							[6,12,0,10], [10,0,12,6], [12,6,10,0], [0,10,6,12]]

	vgrp_tilde['(132)'] = [[14,2,4,8], [2,14,8,4], [4,8,14,2], [8,4,2,14],
							[6,10,12,0], [10,6,0,12], [12,0,6,10], [0,12,10,6]]

	return vgrp_tilde


# ********************************
# Compiling the tetrads from predfined Adinkras
def assemble_tetrads():

	"""Vierergrupe dictionaries with binary quadsets for ells and tilde-ells
	"""

	main_tetrad			= []

	vgruppe_sets		= vierergruppe_sets()
	vierergruppe_elle	= vgrp_ellebin()
	vierergruppe_tilde	= vgrp_tildebin()

	for vgrp, binaries_list in vierergruppe_elle.items():
		vbasis	= vgruppe_sets[vgrp]
		temp 	= lmat_flipping(vbasis, binaries_list)
		main_tetrad.extend(temp)

	for vgrp, binaries_list in vierergruppe_tilde.items():
		vbasis	= vgruppe_sets[vgrp]
		temp 	= lmat_flipping(vbasis, binaries_list)
		# for i, tet in enumerate(temp):
		# 	print("Length of tet:", len(tet), "Type", type(tet))
		# print("Length lmat_flipping", len(temp), vgrp, binaries_list)
		main_tetrad.extend(temp)

	# print(main_tetrad)
	for i, tetrad in enumerate(main_tetrad):
		print("Tetrad #", i, "type:", type(tetrad), "Ls", tetrad)


# ********************************
# Use the binary representation info to perform flips on L mats in each tetrad
def lmat_flipping(vbasis, binaries_list):

	lmat_list		= []

	for xbin in binaries_list:
		binmats = [binaries(b) for b in xbin]
		temp	= [np.dot(vbasis[i], binmats[i]) for i in range(0, len(binmats))]
		lmat_list.append(temp)

	return lmat_list

# ********************************
# Alpha and Beta matrices hardcoded
def pauli_byproducts():

	""" These are the alpha and beta matrices multiplied by 2i
		alpha and beta are results of outerproducts inbetween
		Pauli matrices and Identity matrix They are defined in equtions (4.5)
		in Isaac Chappell II, S. James Gates, Jr - 2012
	"""

	alpha1i = np.matrix([[0, 0, 0, 2], [0, 0, 2, 0], [0, -2, 0, 0], [-2, 0, 0, 0]])
	alpha2i = np.matrix([[0, 2, 0, 0], [-2, 0, 0, 0], [0, 0, 0, 2], [0, 0, -2, 0]])
	alpha3i = np.matrix([[0, 0, 2, 0], [0, 0, 0, -2], [-2, 0, 0, 0], [0, 2, 0, 0]])

	# Betas
	beta1i = np.matrix([[0, 0, 0, 2], [0, 0, -2, 0], [0, 2, 0, 0], [-2, 0, 0, 0]])
	beta2i = np.matrix([[0, 0, 2, 0], [0, 0, 0, 2], [-2, 0, 0, 0], [0, -2, 0, 0]])
	beta3i = np.matrix([[0, 2, 0, 0], [-2, 0, 0, 0], [0, 0, 0, -2], [0, 0, 2, 0]])


# ********************************
# Run main()
start_time = time.time()

main()
print("-- Execution time --")
print("---- %s seconds ----" % (time.time() - start_time))
