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
# Library Imports

import math
import sys
import time
import numpy as np
import numpy.matlib
import itertools
from numpy import array

import vij_holoraumy_prime
import vij_holoraumy_4x4

# ********************************
def main():

	elle_class_calc()


# ********************************
# Performing the coefficient calculation
def elle_class_calc():

	pieslices	= pieslicing()
	vflops		= vierergruppe_flops()

	elle_bin	= flip_ellebin()

	# for i in vflops:
	for vset, flop_ops in vflops.items():
		print(vset,"V")
		vset_flips	= elle_bin[vset]
		print("Bin length", len(vset_flips))
		print("Flips:", vset_flips)
		print("")
		for flip_ops in vset_flips:
			for x in range(0, len(flop_ops)):
				print("Flop x:", flop_ops[x], "Flip x:", flip_ops[x])
				for pie in pieslices:
					temp_flop	= colorspace_flop(pie, flop_ops[x])
					temp_flip	= colorspace_flip(temp_flop, binaries(flip_ops[x]))
					vijset 		= vij_holoraumy_4x4.calculate_vijmatset(temp_flip)
				# temp_adink	= colorspace_flip(colorspace_flop(flop[x], pieslices), binaries[x])


# ********************************
# Defining the Vierergruppe representations for flop operations
def vierergruppe_flops():

	vprime 	= [ "()", "(12)(34)", "(13)(24)", "(14)(23)" ]
	# vprime	= [ "1234", "2143", "3412", "4321"]
	vprime	= [ [1,2,3,4], [2,1,4,3], [3,4,1,2], [4,3,2,1] ]


	vgrpv 		= [ ("()", [1,2,3,4]), ("(12)(34)", [2,1,4,3]),
					("13)(24)", [3,4,1,2]), ("(14)(23)", [4,3,2,1])]

	vgrp12v		= [ ("(12)", [2,1,3,4]), ("(34)", [1,2,4,3]),
					("(1324)", [3,4,2,1]), ("(1423)", [4,3,1,2]) ]

	vgrp13v		= [ ("(13)", [3,2,1,4]), ("(1234)", [2,3,4,1]),
					("(24)", [1,4,3,2]), ("(1432)", [4,1,2,3])	]

	vgruppe	= 	{
				'()': vgrpv,
				'(12)': vgrp12v,
				# '(13)': vgrp13v
				}
	# return [ vgrpv, vgrp12v, vgrp13v ]
	return vgruppe


# ********************************
# Defining the Pizza slices
def pieslicing():

	p1	= [np.matrix([[1, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1], [0, 1, 0, 0]]),
			np.matrix([[0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0], [1, 0, 0, 0]]),
			np.matrix([[0, 0, 1, 0], [1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1]]),
			np.matrix([[0, 0, 0, 1], [0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 1, 0]])
			]

	p2	= [np.matrix([[1, 0, 0, 0], [0, 0, 0, 1], [0, 1, 0, 0], [0, 0, 1, 0]]),
			np.matrix([[0, 1, 0, 0], [0, 0, 1, 0], [1, 0, 0, 0], [0, 0, 0, 1]]),
			np.matrix([[0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1], [1, 0, 0, 0]]),
			np.matrix([[0, 0, 0, 1], [1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0]])
			]

	""" Just one pie piece for now """
	# return [ p1 ]
	return [p1, p2]


# ********************************
# Perform flop operation over Adinkra color space
def colorspace_flop(adinkra, flop_op):
	print("Flop Operation", flop_op[1])
	print("Adinkra")
	print(adinkra)
	print("")

	new_adinkra	= [ adinkra[(ind - 1)] for ind in flop_op[1] ]

	# test_adinkra	= []
	# for ind in flop_op[1]:
	# 	test_adinkra.append(adinkra[(ind - 1)])

	return new_adinkra


# ********************************
# Perform flip operation over Adinkra color space
def colorspace_flip(adinkra, flip_op):

	print("Flip Operation", flip_op)
	new_adinkra	= []
	new_adinkra	= [(adinkra[i] * flip_op[i]) for i in range(0, len(adinkra))]

	print("Flipped Adinkra")
	print(new_adinkra)
	print("")
	return new_adinkra


# ********************************
# Defining the binary multiplication matrices
def binaries(bin_code):

	binaries_lt	= [(0, [1,1,1,1]), (2, [1,-1,1,1]), (4, [1,1,-1,1]),
					(6, [1,-1,-1,1]), (8, [1,1,1,-1]), (10, [1,-1,1,-1]),
					(12, [1,1,-1,-1]), (14, [1,-1,-1,-1])]

	for btuple in binaries_lt:
		if bin_code == btuple[0]:
			# tarray = np.array(btuple[1])
			# temp   = np.diag(tarray)
			return btuple[1]

# ********************************
# Defining the elle binary representations for the Vierergruppe
def flip_ellebin():

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
def flip_tildebin():

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
	vierergruppe_elle	= flip_ellebin()
	vierergruppe_tilde	= flip_tildebin()

	for vgrp, binaries_list in vierergruppe_elle.items():
		vbasis	= vgruppe_sets[vgrp]
		print("")
		print("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx ")
		print("Calculating Vij elle coefficients")
		print("							")
		print("Vierergruppe flop: ",vgrp)
		temp 	= lmat_flipping(vbasis, binaries_list)
		# print("Flip sets:", binaries_list)
		# vij_holoraumy_prime.calculate_vij_matrices(temp)
		calculate_vgruppe_sets(temp, binaries_list)

		print("<<<>>>")
		main_tetrad.extend(temp)

	# for vgrp, binaries_list in vierergruppe_tilde.items():
	# 	vbasis	= vgruppe_sets[vgrp]
	# 	temp 	= lmat_flipping(vbasis, binaries_list)
	# 	print("")
	# 	print("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx ")
	# 	print("Calculating Vij tilde elle coefficients")
	# 	print("							")
	# 	print("Vierergruppe flop: ",vgrp)
	# 	# print("Flip sets:", binaries_list)
	# 	# vij_holoraumy_prime.calculate_vij_matrices(temp)
	# 	calculate_vgruppe_sets(temp, binaries_list)
	#
	# 	print("<<<>>>")
	# 	main_tetrad.extend(temp)

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
# Function for calling calculate_vij_matrices
def calculate_vgruppe_sets(gruppe_adinkras, gruppe_binaries):
	"""
	Function for printing out the details of each Adinkra - Vij matrix
	sixset calculation, including the binary representation and the
	corresponding resulting Vij matrices and their elles/tilde elles
	Coefficients.
	"""

	for i in range(0, len(gruppe_adinkras)):

		print("")
		print("Calculating for binary flip:", gruppe_binaries[i])

		vijset = vij_holoraumy_4x4.calculate_vijmatset(gruppe_adinkras[i])
		for i in vijset:
			print(i)

	# for vgrp, binaries_list in vierergruppe_tilde.items():
	# 	vbasis	= vgruppe_sets[vgrp]
	# 	temp 	= lmat_flipping(vbasis, binaries_list)
	# 	# for i, tet in enumerate(temp):
	# 	# 	print("Length of tet:", len(tet), "Type", type(tet))
	# 	# print("Length lmat_flipping", len(temp), vgrp, binaries_list)
	# 	main_tetrad.extend(temp)

	# print(main_tetrad)
	# for i, tetrad in enumerate(main_tetrad):
		# print("Tetrad #", i, "type:", type(tetrad), "Ls", tetrad)

	# vij_holoraumy_prime.calculate_vij_matrices(main_tetrad)


# ********************************
# Use the binary representation info to perform flips on L mats in each tetrad
def lmat_flipping(vbasis, binaries_list):

	lmat_list		= []

	for xbin in binaries_list:
		print("Flip:", xbin)
		binmats = [binaries(b) for b in xbin]
		temp	= [np.dot(binmats[i], vbasis[i]) for i in range(0, len(binmats))]
		lmat_list.append(temp)

	return lmat_list




# ********************************
# Alpha and Beta matrices hardcoded
def alphas_betas():

	""" These are the alpha and beta matrices multiplied by 2i
		alpha and beta are tensor products of Pauli spin matrices
 		Identity matrix They are defined in equtions (4.5)
		in Isaac Chappell II, S. James Gates, Jr - 2012
	"""

	alpha1i = np.matrix([[0, 0, 0, 2], [0, 0, 2, 0], [0, -2, 0, 0], [-2, 0, 0, 0]])
	alpha2i = np.matrix([[0, 2, 0, 0], [-2, 0, 0, 0], [0, 0, 0, 2], [0, 0, -2, 0]])
	alpha3i = np.matrix([[0, 0, 2, 0], [0, 0, 0, -2], [-2, 0, 0, 0], [0, 2, 0, 0]])

	# Betas
	beta1i = np.matrix([[0, 0, 0, 2], [0, 0, -2, 0], [0, 2, 0, 0], [-2, 0, 0, 0]])
	beta2i = np.matrix([[0, 0, 2, 0], [0, 0, 0, 2], [-2, 0, 0, 0], [0, -2, 0, 0]])
	beta3i = np.matrix([[0, 2, 0, 0], [-2, 0, 0, 0], [0, 0, 0, -2], [0, 0, 2, 0]])

	return [alpha1i, alpha2i, alpha3i, beta1i, beta2i, beta3i]

# ********************************
# Run main()
start_time = time.time()

main()
print("-- Execution time --")
print("---- %s seconds ----" % (time.time() - start_time))
