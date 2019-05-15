import os, sys
sys.path.append(os.getcwd()+'/../')

# @brief the parent class for element operations
class element_operation: 
    
	# @brief the method should compute the element arrays 
	# @param[in] the element object
	def get_element_arrays(self, element ):

		raise NotImplementedError()
