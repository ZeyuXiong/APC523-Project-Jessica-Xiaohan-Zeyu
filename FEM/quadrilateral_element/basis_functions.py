import numpy as np

from utils import lagrange_basis, d_lagrange_basis

def hexahedral_base( poly_order, base_index, x ):
	'''
	Parameters
	----------
	poly_order: int 
		the order of the polynomial interpolant
	base_index: int
		the index of the basis function

	Returns
	-------
	The value of the base function corresponding to 
	base_index defined sover parameteric space 
	$[-1,1]\times[-1,1]$ evaluated at point x

	The bases are numbered consequentially. Eg. for
	quadratic elements 
	6--7--8
	|     |
	3  4  5
	|     |
	0--1--2
	'''

	# The nodes of the one dimensional polynomial 
	# over the parametric interval [-1,1]
	param_one_dim = np.linspace(-1., 1., poly_order + 1 )

	# Get the index in the xi_1 direction 
	base_index_xi_1 = base_index%(poly_order + 1 )

	# Get the index in the xi_2 direction 
	base_index_xi_2 = base_index//(poly_order + 1)

	l0 = lagrange_basis(param_one_dim,base_index_xi_1,poly_order,x[0])
	l1 = lagrange_basis(param_one_dim,base_index_xi_2,poly_order,x[1])

	return l0*l1
						

def grad_hexahedral_base( poly_order, base_index,x ):
	'''
	Parameters
	----------
	poly_order: int 
		the order of the polynomial interpolant
	base_index: int
		the index of the basis function
	
	Returns
	-------
	The value of the gradient of the base function corresponding to 
	base_index defined sover parameteric space 
	$[-1,1]\times[-1,1]$ evaluated at point x
	
	The bases are numbered consequentially. Eg. for
	quadratic elements 
	6--7--8
	|     |
	3  4  5
	|     |
	0--1--2
	'''

	# The nodes of the one dimensional polynomial 
	# over the parametric interval [-1,1]
	param_one_dim = np.linspace(-1., 1., poly_order + 1 )

	# Get the index in the xi_1 direction 
	base_index_xi_1 = base_index%(poly_order + 1 )

	# Get the index in the xi_2 direction 
	base_index_xi_2 = base_index//(poly_order + 1)

	l0 = lagrange_basis(param_one_dim,base_index_xi_1,poly_order,x[0])
	l1 = lagrange_basis(param_one_dim,base_index_xi_2,poly_order,x[1])
	dl0 = d_lagrange_basis(param_one_dim,base_index_xi_1,poly_order,x[0])
	dl1 = d_lagrange_basis(param_one_dim,base_index_xi_2,poly_order,x[1])

	# Return the gradient
	return np.array([ dl0*l1, l0*dl1 ])
