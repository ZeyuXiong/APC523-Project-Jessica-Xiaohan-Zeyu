import numpy as np
import matplotlib.pyplot as plt
# @param[in] nodes the coordinates of the nodes
# @param[in] index the index of the basis function
# @param[in] order the order of the polynomial basis
# @param[in] x the coordinate where to evaluate the basis
# @return the value of the index-th basis function at point x
def lagrange_basis(nodes, index, order, x): 
	
	# Check that we have the right number of nodes
	assert len(nodes) == (order + 1)

	# Create the initial value of the function
	ell = 1.

	# Loop over all nodes
	for i in range(0,order+1):
		
		# If this node is the same as the 
		# support node of the basis function, skip it
		if i == index:
			continue

		# Otherwise perform the multiplication
		ell *= ( x - nodes[i])/(nodes[index] - nodes[i])
	
	return ell	

# @param[in] nodes the coordinates of the nodes
# @param[in] index the index of the basis function
# @param[in] order the order of the polynomial basis
# @param[in] x the coordinate where to evaluate the basis
# @return the value of the index-th basis function at point x
def d_lagrange_basis(nodes, index, order, x): 
	
	# Check that we have the right number of nodes
	assert len(nodes) == (order + 1)

	# Create the initial value of the function
	d_ell = 0

	# Loop over all nodes
	for i in range(0,order+1):
		
		# If this node is the same as the 
		# support node of the basis function, skip it
		if i == index:
			continue
	
		coeff = 1./(nodes[index] - nodes[i]) 

		for j in range(0,order+1):

			if j == i or j == index:
				continue
			
			coeff *= (x - nodes[ j ] )/(nodes[index] - nodes[j]) 

		d_ell += coeff

	return d_ell	


if __name__ == "__main__":
	
	poly_order = 3
	nodes = np.linspace(-1,1,poly_order+1)
	
	xs = np.linspace(-1,1,100)
	
	for i in range( poly_order + 1 ):
		#y = lagrange_basis( nodes, i , poly_order, x )
		dy = np.array([ d_lagrange_basis( nodes, i , poly_order, x ) for x in xs ])
		plt.plot(xs,dy,'-r')

	plt.show()
		
