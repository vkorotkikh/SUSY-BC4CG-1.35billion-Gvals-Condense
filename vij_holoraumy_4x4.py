# ******************************************************************************
# Name:    Calculate Vij matrices and the elle & tilde~elle Coefficients
# Author:  Vadim Korotkikh
# Email:   va.korotki@gmail.com
# Date:    November 2016
# Version: 1.3
#
# Description: Calculates the corresponding Vij Holoraumy matrices for Adinkras
# Using the Holoraumy matrices, finds the corresponding elle or tilde`elle
# coefficients for ALpha or Beta matrices
#
# ******************************************************************************


# ******************************************************************************
# Begin Imports

import sys
import math
import time
import itertools
import numpy as np
from numpy import array
from numpy.linalg import inv

# ******************************************************************************
# Main() function.
def main():
	# gen_signm(4)
	pass

# ******************************************************************************
# Do the final Vij calculation
def calculate_vij_matrices(main_tetrad_list):

	""" Remember that the main_tetrad_ark is a list of lists,
		with each list containing four tuples, with tuples being
		matrix number and the matrices itself. """


	vij_possibilities 	= alphas_betas()
	vij_sixset 			= []

	ij_indices			= list(itertools.combinations([0,1,2,3], 2))

	# print("							")
	# print("Calculating Vij matrices")
	# print("							")
	vij_alphas 		= []
	vij_betas  		= []
	calc_check		= []

	anomaly_switch  = 0
	debug			= 0

	for ti, teti in enumerate(main_tetrad_list):
		if debug:
			print("# ********************************")
			print("								     ")
			print("Tetrad i: ", ti)

		alpha_temp	= []
		beta_temp   = []
		vij_tempset = []
		""" Store 6 Vij matrices in temp_vijmat"""
		# temp_vijmat		= []

		for ijtup in ij_indices:
			limat 		= teti[ijtup[0]]
			ljmat 		= teti[ijtup[1]]
			ij_temp		= str(ijtup[0] + 1) + str(ijtup[1] + 1)
			ijstr		= ij_temp
			tr_limat	= np.transpose(limat)
			tr_ljmat	= np.transpose(ljmat)
			""" Vij eq from 1601.00 (3.2) """
			temp_mat	= np.dot(tr_limat, ljmat) - np.dot(tr_ljmat, limat)
			""" Compare against the 6 possible matrix solutions """
			tf_bool = 0

			for xi, ijx in enumerate(vij_possibilities):
				ijx_neg = np.multiply(ijx, -1)
				# print(xi)
				if np.array_equal(temp_mat, ijx):
					tf_bool = 1
					if debug:
						print("*************$$$$$$$$$$$$$$$$$$ ")
						print("l-solution found:")
						print(ijx)
					tmint = np.int(1)
					if xi < 3:
						tmp_str = "alpha" + str((xi + 1))
						# print(tmp_str)
						vij_tempset.append([tmp_str, ijstr, tmint])
						alpha_temp.append([tmp_str, ijstr, tmint])
					elif xi >= 3:
						tmp_str = "beta" + str((xi - 2))
						vij_tempset.append([tmp_str, ijstr, tmint])
						beta_temp.append([tmp_str, ijstr, tmint])
				elif np.array_equal(temp_mat, ijx_neg):
					tf_bool = 1
					if debug:
						print("*************$$$$$$$$$$$$$$$$$$ ")
						print("l-solution found:")
						print(ijx_neg)
					# xint = (xi + 1) * ( -1)
					tmint = np.int(-1)
					if xi < 3:
						tmp_str = "alpha" + str((xi + 1))
						# print(tmp_str)
						vij_tempset.append([tmp_str, ijstr, tmint])
						alpha_temp.append([tmp_str, ijstr, tmint])
					elif xi >= 3:
						tmp_str = "beta" + str((xi - 2))
						vij_tempset.append([tmp_str, ijstr, tmint])
						beta_temp.append([tmp_str, ijstr, tmint])
				else:
					if tf_bool == 0 and xi >= 5:
						if not(np.array_equal(temp_mat, ijx)) or not np.array_equal(temp_mat, ijx_neg):
							print("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx ")
							print("Anomaly found:",ijstr)
							print(temp_mat)
							anomaly_switch = 1
			tf_bool = 0

		calc_check.append(vij_tempset)


		if alpha_temp:
			vij_alphas.append(alpha_temp)
		elif beta_temp:
			vij_betas.append(beta_temp)
		beta_temp 	= []
		alpha_temp 	= []

	print("*************$$$$$$$$$$$$$$$$$$ ")
	print("Vij Matrix Coefficients Results:")
	print("")
	for mvals in calc_check:
		if any(x for x in mvals if x[0].startswith('alpha')) and any(x for x in mvals if x[0].startswith('beta')):
			print("MIXED ALPHA_BETA ERROR")
			print(mvals)
		else:
			print(mvals)

	print("Length Vij alphas adinkras:", len(vij_alphas))
	print("Length Vij beta adikras:", len(vij_betas))

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


		# 	""" This section does a double loop over the same tetrad to calculate
		# 	the set of 6 Vij matrices for the tetrad.
		# 	So for each matrix in the tetrad its checked against all the possible others,
		#  	bypassing the duplicate calculations
		# 	"""
		# 	for i, li in enumerate(teti):
		# 		# print(li[1])
		# 		bigli = li[1]
		# 		tr_bigli = np.transpose(bigli)
		#
		# 		for j, lj in enumerate(teti):
		# 			biglj = lj[1]
		# 			ij_temp = [i, j]
		# 			ij_temp.sort()
		# 			ir = i + 1
		# 			jr = j + 1
		# 			ijstr = str(ir) + str(jr)
		# 			if ij_temp not in temp_combos and i != j:
		# 				# print("Vij matrix i-j vals:", ij_temp)
		# 				# print("Vij matrix i-j vals:", ijstr)
		# 				temp_combos.append(ij_temp)
		# 				tr_biglj = np.transpose(biglj)
		# 				# temp_mat = np.dot(tr_bigli, biglj) - np.dot(tr_biglj, bigli)
		# 				""" Vij eq from 1601.00 (3.2) """
		# 				# temp_mat = np.matmul(tr_biglj, bigli) - np.matmul(tr_bigli, biglj)
		# 				temp_mat = np.dot(tr_bigli, biglj) - np.dot(tr_biglj, bigli)
		# 				""" Compare against the 6 possible matrix solutions """
		# 				tf_bool = 0
		# 				for xi, ijx in enumerate(vij_possibilities):
		# 					ijx_neg = np.multiply(ijx, -1)
		# 					# print(xi)
		# 					if np.array_equal(temp_mat, ijx):
		# 						tf_bool = 1
		# 						temp_vijmat.append(temp_mat)
		# 						if debug:
		# 							print("*************$$$$$$$$$$$$$$$$$$ ")
		# 							print("l-solution found:")
		# 							print(ijx)
		# 						tmint = np.int(1)
		# 						if xi < 3:
		# 							tmp_str = "alpha" + str((xi + 1))
		# 							# print(tmp_str)
		# 							vij_tempset.append([tmp_str, ijstr, tmint])
		# 							alpha_temp.append([tmp_str, ijstr, tmint])
		# 						elif xi >= 3:
		# 							tmp_str = "beta" + str((xi - 2))
		# 							vij_tempset.append([tmp_str, ijstr, tmint])
		# 							beta_temp.append([tmp_str, ijstr, tmint])
		# 					elif np.array_equal(temp_mat, ijx_neg):
		# 						tf_bool = 1
		# 						temp_vijmat.append(temp_mat)
		# 						if debug:
		# 							print("*************$$$$$$$$$$$$$$$$$$ ")
		# 							print("l-solution found:")
		# 							print(ijx_neg)
		# 						# xint = (xi + 1) * ( -1)
		# 						tmint = np.int(-1)
		# 						if xi < 3:
		# 							tmp_str = "alpha" + str((xi + 1))
		# 							# print(tmp_str)
		# 							vij_tempset.append([tmp_str, ijstr, tmint])
		# 							alpha_temp.append([tmp_str, ijstr, tmint])
		# 						elif xi >= 3:
		# 							tmp_str = "beta" + str((xi - 2))
		# 							vij_tempset.append([tmp_str, ijstr, tmint])
		# 							beta_temp.append([tmp_str, ijstr, tmint])
		# 					else:
		# 						if i != j and tf_bool == 0 and xi >= 5:
		# 							if not(np.array_equal(temp_mat, ijx)) or not np.array_equal(temp_mat, ijx_neg):
		# 								print("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx ")
		# 								print("Anomaly found:",i,j)
		# 								print(temp_mat)
		# 								anomaly_switch = 1
		#
		# 				tf_bool = 0
		#
		# 	vij_matrices.append(temp_vijmat)
		# 	calc_check.append(vij_tempset)
		#
		#
		# 	if alpha_temp:
		# 		vij_alphas.append(alpha_temp)
		# 	elif beta_temp:
		# 		vij_betas.append(beta_temp)
		# 	beta_temp 	= []
		# 	alpha_temp 	= []
		#
		# print("*************$$$$$$$$$$$$$$$$$$ ")
		# print("Vij Matrix Coefficients Results:")
		# print("")
		# for mvals in calc_check:
		# 	if any(x for x in mvals if x[0].startswith('alpha')) and any(x for x in mvals if x[0].startswith('beta')):
		# 		print("MIXED ALPHA_BETA ERROR")
		# 		print(mvals)
		# 	else:
		# 		print(mvals)
		#
		# print("Length Vij alphas tetrads:", len(vij_alphas))
		# print("length Vij beta tetrads:", len(vij_betas))
