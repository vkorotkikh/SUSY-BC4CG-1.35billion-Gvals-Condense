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

##************************************
def main():

	elle_class_calc()

##************************************
# Performing the coefficient calculation
def elle_class_calc():

	vflops		= vierergruppe_flops()
	# elle_bin	= flip_ellebin()
	for vset in vflops:
		# print(vset[0], "V set")
		print("Vierergruppe:", vset[0],"V")
		print("Flop Ops:", vset[1])
		print("")
		# print("Flips:",temp_flips)
		# print("")
		""" For every populated cis Pie slice - 8 Adinkras per slice"""
		for i in range(0, 6):
			temp_pie	= cis_seed_pies(i)
			# print("Pie i:", i)

			""" For i Adinkra out of selected Pie slice	"""
			pie_vijres	=	[]
			pie_newrep	=	[]
			for itet, adinkra in enumerate(temp_pie):
				""" temp_flips contains a Group of flips associated with a Vierergruppe
				so 8 sets of 4. """
				temp_flips	= flip_ellebin(vset[0])
				print("Flips:",temp_flips)
				print("")
				print("Pie #", i, "Adinkra", itet)
				for flip_ops in temp_flips:
					print("Current Flip:", flip_ops)
					print("")
					vijres_temp	= []
					newrep_temp = []
					for x in range(0, len(vset[1])):
						flop_tup 	=	vset[1][x]
						flop_op		=	vset[1][x][1]
						assoc_flip	=	flip_ops[x]
						# print("Current Flop:", flop_op, "	Flip:", flip_ops[x], binaries(flip_ops[x]))
						temp_flop		=	colorspace_flop(adinkra, flop_op)
						temp_flip		=	colorspace_flip(temp_flop, binaries(flip_ops[x]))
						vijset, newrep	=	vij_holoraumy_4x4.calculate_vijmatset_nicely(temp_flip)
						if vijset not in vijres_temp:
							vijres_temp.append(vijset)
							if len(vijres_temp) > 1:
								print("Current Flop:", flop_op, "	Flip:", flip_ops[x])
						else:
							pass

						if newrep not in newrep_temp:
							newrep_temp.append(newrep)
						else:
							pass
					print("vijres_temp")
					for zz in range(0, len(vijres_temp)):
						print(vijres_temp[zz])
					print("newrep_temp")
					for zz in range(0, len(newrep_temp)):
						print(newrep_temp[zz])
					print("		")
					pie_vijres	=	[vij for vij in vijres_temp]
					pie_newrep	=	[newrep for newrep in newrep_temp]


##************************************
# Performing the coefficient calculation
def elle_class_calc():

	vflops		= vierergruppe_flops()
	# elle_bin	= flip_ellebin()
	for vset in vflops:
		# print(vset[0], "V set")
		print("Vierergruppe:", vset[0],"V")
		print("Flop Ops:", vset[1])
		print("")
		# print("Flips:",temp_flips)
		# print("")
		""" For every populated cis Pie slice - 8 Adinkras per slice"""
		for i in range(0, 6):
			temp_pie	= cis_seed_pies(i)
			# print("Pie i:", i)

			""" For i Adinkra out of selected Pie slice	"""
			for itet, adinkra in enumerate(temp_pie):
				""" temp_flips contains a Group of flips associated with a Vierergruppe
				so 8 sets of 4. """
				temp_flips	= flip_ellebin(vset[0])
				print("Flips:",temp_flips)
				print("")
				print("Pie #", i, "Adinkra", itet)
				for flip_ops in temp_flips:
					print("Current Flip:", flip_ops)
					print("")
					vijres_temp	= []
					newrep_temp = []
					for x in range(0, len(vset[1])):
						flop_tup 	=	vset[1][x]
						flop_op		=	vset[1][x][1]
						assoc_flip	=	flip_ops[x]
						# print("Current Flop:", flop_op, "	Flip:", flip_ops[x], binaries(flip_ops[x]))
						temp_flop		=	colorspace_flop(adinkra, flop_op)
						temp_flip		=	colorspace_flip(temp_flop, binaries(flip_ops[x]))
						vijset, newrep	=	vij_holoraumy_4x4.calculate_vijmatset_nicely(temp_flip)
						if vijset not in vijres_temp:
							vijres_temp.append(vijset)
							if len(vijres_temp) > 1:
								print("Current Flop:", flop_op, "	Flip:", flip_ops[x])
						else:
							pass

						if newrep not in newrep_temp:
							newrep_temp.append(newrep)
						else:
							pass
					print("vijres_temp")
					for zz in range(0, len(vijres_temp)):
						print(vijres_temp[zz])
					print("newrep_temp")
					for zz in range(0, len(newrep_temp)):
						print(newrep_temp[zz])
					print("		")
				# for x in range(0, len(vset[1])):
				# 	flop_tup 	= vset[1][x]
				# 	flop_op		= vset[1][x][1]
				# 	# print("vset[1]:", flop_tup)
				# 	for flip_ops in temp_flips:
				# 		assoc_flip	=	flip_ops[x]
				# 		print("")
				# 		print("Pie i", i, "Adinkra", itet)
				# 		print("Current Flop:", flop_op, "	Flip:", flip_ops[x], binaries(flip_ops[x]))
				# 		temp_flop	=	colorspace_flop(adinkra, flop_op)
				# 		temp_flip	=	colorspace_flip(temp_flop, binaries(flip_ops[x]))
				# 		vijset		=	vij_holoraumy_4x4.calculate_vijmatset(temp_flip)
				# 		print(vijset)
	# for vset, flop_ops in vflops.items():
	# 	print(vset,"V")
	# 	vset_flips	= elle_bin[vset]
	# 	print("Bin length", len(vset_flips))
	# 	print("Flips:", vset_flips)
	# 	print("")
	# 	for flip_ops in vset_flips:
	# 		for x in range(0, len(flop_ops)):
	# 			print("Flop x:", flop_ops[x], "Flip x:", flip_ops[x])
	# 			for pie in pieslices:
	# 				temp_flop	= colorspace_flop(pie, flop_ops[x])
	# 				temp_flip	= colorspace_flip(temp_flop, binaries(flip_ops[x]))
	# 				vijset 		= vij_holoraumy_4x4.calculate_vijmatset(temp_flip)
	# 				print(vijset)


##************************************
# Defining the Vierergruppe representations for flop operations
def vierergruppe_flops():

	vprime 	= [ "()", "(12)(34)", "(13)(24)", "(14)(23)" ]
	# vprime	= [ "1234", "2143", "3412", "4321"]
	vprime	= [ [1,2,3,4], [2,1,4,3], [3,4,1,2], [4,3,2,1] ]


	vgrpv 		= [ ("()", [1,2,3,4]), ("(12)(34)", [2,1,4,3]),
					("13)(24)", [3,4,1,2]), ("(14)(23)", [4,3,2,1])]

	vgrp12v		= [ ("(12)", [2,1,3,4]), ("(34)", [1,2,4,3]),
					# ("(1324)", [3,4,2,1]), ("(1423)", [4,3,1,2]) ]
					("(1324)", [4,3,1,2]), ("(1423)", [3,4,2,1]) ]

	vgrp13v		= [ ("(13)", [3,2,1,4]), ("(1234)", [4,1,2,3]),
					("(24)", [1,4,3,2]), ("(1432)", [2,3,4,1])	]

	vgrp23v		= [ ("(23)", [1,3,2,4]), ("(1342)", [2,4,1,3]),
					("(1243)", [3,1,4,2]), ("(14)", [4,2,3,1])	]

	vgrp123v	= [ ("(132)", [3,1,2,4]), ("(132)", [4,2,1,3]),
					("(243)", [1,3,4,2]), ("(142)", [2,4,3,1])	]

	vgrp132v	= [ ("(132)", [2,3,1,4]), ("(234)", [1,4,2,3]),
					("(124)", [4,1,3,2]), ("(143)", [3,2,4,1])	]

	# vgrp13v		= [ ("(13)", [3,2,1,4]), ("(1234)", [2,3,4,1]),
	# 				("(24)", [1,4,3,2]), ("(1432)", [4,1,2,3])	]
	# vgruppe 	= [ ('()', vgrpv), ('(12)',vgrp12v), ('(13)', vgrp13v) ]
	# vgruppe 	= [ ('()', vgrpv), ('(12)',vgrp12v)]
	vgruppe	   = [ ('(23)', vgrp23v), ('(123)', vgrp123v), ('(132)', vgrp132v) ]

	return vgruppe

##************************************
# cis seed pie slices - elle coefficients
def cis_seed_pies(pie_index):

	""" Defining the words (aka collection of four boolean factors) that
		when applied to corresponding Pie slices "promote" the Pie slice
		to Adinkras
	"""
	p1plus	= 	[	[0,12,10,6], [2,14,8,4], [4,8,14,2], [6,10,12,0],
					[8,4,2,14], [10,6,0,12], [12,0,6,10], [14,2,4,8]
				]
	p2plus	=	[	[0,6,12,10], [2,4,14,8], [4,2,8,14], [6,0,10,12],
					[8,14,4,2], [10,12,6,0], [12,10,0,6], [14,8,2,4]
				]
	p3plus	=	[	[12,0,10,6], [14,2,8,4], [8,4,14,2], [10,6,12,0],
					[4,8,2,14],	[6,10,0,12], [0,12,6,10], [2,14,4,8]
				]
	p4plus	=	[	[0,10,12,6], [2,8,14,4], [4,14,8,2], [6,12,10,0],
					[8,2,4,14], [10,0,6,12], [12,6,0,10], [14,4,2,8]
				]
	p5plus	=	[	[0,6,10,12], [2,4,8,14], [4,2,14,8], [6,0,12,10],
					[8,14,2,4], [10,12,0,6], [12,10,6,0], [14,8,4,2]
				]
	p6plus	=	[	[0,10,6,12], [2,8,4,14], [4,14,2,8], [6,12,0,10],
					[8,2,14,4], [10,0,12,6], [12,6,10,0], [14,4,8,2]
				]

	cis_promotions	=	[p1plus, p2plus, p3plus, p4plus, p5plus, p6plus]

	""" Import pie slice definitions before applying the binaries	"""
	# pie_slices	=	pieslices()

	# for i in range(0, len(cis_promotions)):

	pslice			=	pieslices(pie_index)
	pslice_words	=	cis_promotions[pie_index]

	""" List to hold all the promoted Pie slice Adinkras using corresponding
		words """
	pie_adinkras	=	[]
	for word in pslice_words:
		temp_adinkra	=	[]
		bool_list		=	[]
		for bins in word:
			bin_list 	= binaries(bins)
			temp		= np.array(bin_list)
			bool_mat	= np.diag(temp)
			bool_list.append(bool_mat)
		temp_adinkra	= [(np.dot(bool_list[x],pslice[x])) for x in range(0,len(pslice))]
		# pie_adinkras.append(temp_adinkra)
		pie_adinkras.append(temp_adinkra)

	return pie_adinkras

##************************************
# trans seed pie slices - tilde~elle coefficients
def trans_seed_pies(pie_index):

	""" Defining the words (aka collection of four boolean factors) that
		when applied to corresponding Pie slices "promote" the Pie slice
		to Adinkras
	"""
	p1neg	=	[	[0,10,6,12], [2,8,4,14], [4,14,2,8], [6,12,0,10],
					[8,2,14,4], [10,0,12,6], [12,6,10,0], [14,4,8,2]
				]

	p2neg	=	[	[0,12,10,6], [2,14,8,4], [4,8,14,2], [6,10,12,0],
					[8,4,2,14], [10,6,0,12], [12,0,6,10], [14,2,4,8]
				]

	p3neg	=	[	[6,0,12,10], [4,2,14,8], [2,4,8,14], [0,6,10,12],
					[14,8,4,2], [12,10,6,0], [10,12,0,6], [8,14,2,4]
				]

	p4neg	=	[	[0,12,6,10], [2,4,14,8], [4,8,2,14], [6,10,0,12],
					[8,4,14,2],	[10,6,12,0], [12,0,10,6], [14,2,8,4]
				]

##************************************
# Defining the Pizza slices
def pieslices(pie_index):

	p1	= 	[np.matrix([[1, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1], [0, 1, 0, 0]]),
			np.matrix([[0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0], [1, 0, 0, 0]]),
			np.matrix([[0, 0, 1, 0], [1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1]]),
			np.matrix([[0, 0, 0, 1], [0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 1, 0]])
			]

	p2	= 	[np.matrix([[1, 0, 0, 0], [0, 0, 0, 1], [0, 1, 0, 0], [0, 0, 1, 0]]),
			np.matrix([[0, 1, 0, 0], [0, 0, 1, 0], [1, 0, 0, 0], [0, 0, 0, 1]]),
			np.matrix([[0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1], [1, 0, 0, 0]]),
			np.matrix([[0, 0, 0, 1], [1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0]])
			]

	p3	=	[np.matrix([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]]),
			np.matrix([[0, 1, 0, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 0, 1, 0]]),
			np.matrix([[0, 0, 1, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 1, 0, 0]]),
			np.matrix([[0, 0, 0, 1], [0, 1, 0, 0], [0, 0, 1, 0], [1, 0, 0, 0]])
			]

	p4	=	[np.matrix([[1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0]]),
			np.matrix([[0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0]]),
			np.matrix([[0, 0, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1]]),
			np.matrix([[0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0]])
			]

	p5	=	[np.matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]]),
			np.matrix([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]),
			np.matrix([[0, 0, 1, 0], [0, 0, 0, 1], [0, 1, 0, 0], [1, 0, 0, 0]]),
			np.matrix([[0, 0, 0, 1], [0, 0, 1, 0], [1, 0, 0, 0], [0, 1, 0, 0]])
			]

	p6	=	[np.matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]),
			np.matrix([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]]),
			np.matrix([[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]]),
			np.matrix([[0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0]])
			]

	pie_complete = [ p1, p2, p3, p4, p5, p6 ]


	""" Just one pie piece for now """
	# return [ p1 ]
	return pie_complete[pie_index]


##************************************
# Perform flop operation over Adinkra color space
def colorspace_flop(adinkra, flop_op):
	# print("	")
	# print("Executing colorspace_flip", flop_op)
	# print("Adinkra")
	# print(adinkra)
	# print("")
	new_adinkra	= [ adinkra[(ind - 1)] for ind in flop_op]

	return new_adinkra


##************************************
# Perform flip operation over Adinkra color space
def colorspace_flip(adinkra, flip_op):

	print("Executing colorspace_flip", flip_op)
	""" Weird bug here if you do the algorith this way """
	# new_adinkra	= []
	# for i in range(0, len(adinkra)):
	# 	tmat 		= adinkra[i]
	# 	tmat[1:]	= tmat[1:] * flip_op[i]
	# 	new_adinkra.append(tmat)
	""" Normal algorithm """
	new_adinkra	= []
	for i in range(0, len(flip_op)):
		# new_adinkra[i][1:]	= new_adinkra[i][1:] * flip_op[i]
		temp_mat	= adinkra[i] * flip_op[i]
		new_adinkra.append(temp_mat)
	# print("Flipped Adinkra")
	# print(new_adinkra)
	# print("")
	return new_adinkra


##************************************
# Defining the binary multiplication arrays
def binaries(bin_code):

	binaries_lt	= [(0, [1,1,1,1]), (2, [1,-1,1,1]), (4, [1,1,-1,1]),
					(6, [1,-1,-1,1]), (8, [1,1,1,-1]), (10, [1,-1,1,-1]),
					(12, [1,1,-1,-1]), (14, [1,-1,-1,-1])]

	for btuple in binaries_lt:
		if bin_code == btuple[0]:
			# tarray = np.array(btuple[1])
			# temp   = np.diag(tarray)
			return btuple[1]

##************************************
# Defining the elle binary representations for the Vierergruppe
def flip_ellebin(flip_set):

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

	return vgrp_elle[flip_set]


##************************************
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
