import numpy as np
import time

# @brief assembles the local arrays into the global arrays
# @param[in] local_to_global_map the local to global map
# @param[in] element the element object
# @param[in] ke the element stiffness matrix
# @param[in] fe the element stiffness matrix
# @param[out] K the reference to the global stiffness matrix
def assemble_local(local_to_global_map, element, ke, fe, K, F ):

	# Loop over all degrees of freedom 
	for i in range(element.num_dofs):

		# Get the global index of the local i dof
		i_global = local_to_global_map.get_global_dof( element.element_index , i )

		# Add the contribution to the source vector
		F[i_global] += fe[i]

		for j in range(element.num_dofs):

			# Get the global index of the local j dof
			j_global = local_to_global_map.get_global_dof( element.element_index , j )

			# Add the contribution of the local stiffness to the global stiffness matrix		
			K[i_global, j_global] += ke[ i, j]

	return


# @brief assembles the local arrays into the global arrays
# @param[in] local_to_global_map the local to global map
# @param[in] element_operation the element operation 
# @param[in] elements the aray containing the lement objects 
# @param[in] ke the element stiffness matrix
# @param[in] fe the element stiffness matrix
# @param[out] K the reference to the global stiffness matrix
# @param[out] F the reference to the global source vector
def assemble_global(local_to_global_map, element_operation, elements, K, F ):

	time_elt_stiff = 0
	time_assemble = 0

	# Loop over all elements to assemble the global stiffness matrix
	for element in elements:

		# The element stiffness matrix and source vector
		t0 = time.time()
		ke, fe = element_operation.get_element_arrays(element)
		time_elt_stiff += time.time() - t0

		# Assemble them 
		t0 = time.time()
		assemble_local(local_to_global_map, element, ke, fe, K, F )
		time_assemble += time.time() - t0
	print('Time for computation of element stiffnesses: %.2e s'%time_elt_stiff ,\
			'Time for assembly: %.2e s'%time_assemble )
			

	return

