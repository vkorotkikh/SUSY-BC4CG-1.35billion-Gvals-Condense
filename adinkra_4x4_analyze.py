# ******************************************************************************
# Name:    Adinkra 4 x 4 analysis
# Author:  Vadim Korotkikh
# Email:   va.korotki@gmail.com
# Date:    June 2017
#
# Description: Added later
#
# ******************************************************************************
# Library Imports

from numpy.linalg import inv
import numpy as np
import itertools
import time
import sys

# Load scripts
import adinkra_nxn_constructor

def main():
	matlist		= adinkra_nxn_constructor.gen_product_matrices(4)
	print(len(matlist))
	# pass

def skew_sym():
	pass



main()
