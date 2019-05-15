import os, sys
sys.path.append(os.getcwd()+'./../')

from element_operations import *
import numpy as np

class poisson(element_operation):

    def __init__(self, f = 0 ):

        self.f = f 

    ## @brief the implementation of the element array for the poisson problem
    def get_element_arrays(self,element):

        # Get the element local degress of freedom
        num_dofs = element.get_num_dofs( )

        # Create the local stiffness matrix
        ke = np.zeros( (num_dofs, num_dofs ) )

        # Create the local force vector
        fe = np.zeros( num_dofs )

        # Get the quadrature rule 
        gauss_points, gauss_weights = element.get_quadrature(quadrature_order=10)

        # Fill in the stiffness matrix
        for q in range(len(gauss_points)):

            # Get the jacobian at the quadrature point 
            jacobian = element.get_dmap( gauss_points[q], jacobian=True )

            # Get the image of the gauss point on the physical domain 
            x = element.get_map( gauss_points[q] ) 

            # Get the value of the gradient of all basis at the quadrature point
            grad_phi = [] 

            for k in range(num_dofs):
                # Get the base function gradients (already with respect to physical coordinates)
                grad_phi.append( element.get_base_function_grad( k, gauss_points[q] ) )

            # Loop over all degree of freedoms
            for i in range(num_dofs):

                for j in range(num_dofs):

                    # Add the contribution to the stiffness matrix
                    ke[i,j] +=  np.dot( grad_phi[i] , grad_phi[j] )*jacobian*gauss_weights[q]

                # If we specified a source term
                if self.f:
                    fe[i] += gauss_weights[q]*element.get_base_function_val(i,gauss_points[q])*\
                                    self.f( x )*jacobian

        return ke, fe

