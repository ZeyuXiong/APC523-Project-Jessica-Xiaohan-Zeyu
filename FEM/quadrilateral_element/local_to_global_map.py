import numpy as np

# @brief The local to global map for abitrary order quadrilateral elements
class local_to_global_map:

	#	Parameters
	#	----------
	#	poly_order: int
	#		the order of the polynomial interpolant
	#	connectivity: array
	#		the numpy array of dimension num_elements \times vertex_per_element	
	def __init__(self, poly_order, connectivity ):

		# Save the connectivity 
		self.connectivity = connectivity

		# Save the polynomial order
		self.poly_order = poly_order 

		# The number of vertices
		self.vertices_element = 4
		
		# The number of elements in the mesh
		self.num_elements = len( connectivity )

		# Get the total number of vertex in the mesh 
		self.num_vertex = len(set(connectivity.flatten()))

		# The vertex dof per element
		self.num_vertex_dof = 4

		# The dofs per edge 
		self.num_edge_dof = poly_order - 1

		# The number of edges per element 
		self.num_edges_element = 4 

		# The interior dofs 
		self.num_interior_dof = pow( poly_order + 1, 2)  -\
				self.num_edge_dof*self.num_edges_element - self.num_vertex_dof 

		# Construct an edge list 
		self.edges = dict()

		# Construct an edge connectivity list
		self.edge_connectivity = dict()

		# Loop over all elements
		for e,conn_e in enumerate(connectivity):

			# Add the edge to the edge dictionary 
			for i in range( len(conn_e) ):

				edge = tuple(np.sort(( conn_e[ i ] , conn_e[ (i+1)%self.vertices_element ] )) )

				if not edge in self.edges:
					self.edge_connectivity[ edge ] = [ (e,i) ]

					self.edges[ edge ] = len( self.edges )
				else :

					self.edge_connectivity[edge].append( (e,i) )

		# The total number of unique edges in the mesh
		self.num_edges = len( self.edges )

		# The total number of degrees of freedom 
		self.total_dof = self.num_elements*self.num_interior_dof +\
				self.num_edges*self.num_edge_dof +\
				self.num_vertex

	def get_global_dof(self, element, local_dof ):
		'''
		Parameters
		----------
		element: int
			the element number
		local_dof: int 
			the local degree of freedom
		Returns
		-------
			The global degree of freedom
		'''

		nx = local_dof%(self.poly_order + 1) 
		ny = local_dof//(self.poly_order + 1)

		# If a vertex node
		if (not nx%self.poly_order) and (not ny%self.poly_order ) : 
			vertex_num = int(abs( ny//self.poly_order*3 - nx//self.poly_order ))
			global_dof = self.connectivity[element][ vertex_num ]
			return global_dof

		# If an edge node 
		elif (not nx%self.poly_order) or (not ny%self.poly_order ):

			# Local edge number for this dof
			local_edge_number =  ( 1 if  nx//self.poly_order else 0 ) +\
					( 2 if ny//self.poly_order else 0 ) + \
					( 3 if not nx else 0 ) 
					
			# Local dof number on edge 
			edge_dof = ny - 1 if ny%self.poly_order else nx - 1 
			
			if local_edge_number > 1:
				edge_dof = self.num_edge_dof - 1 - edge_dof
      
			# Get connectivity of this edge
			vertex_element = 4
			conn0 = self.connectivity[element][local_edge_number]
			conn1 = self.connectivity[element][(local_edge_number+1)%self.vertices_element]
		
			edge_num = self.edges[ tuple( np.sort([ conn0,conn1])) ]
			
			#print 'Conn / edge ',conn0,conn1, '\t', edge_num
			if conn0 > conn1 :
				edge_dof = self.num_edge_dof - 1 - edge_dof 

			#print 'local nx ny', local_dof, nx, ny
			return self.num_vertex + self.num_edge_dof*edge_num + edge_dof

		# If an interior node
		else:

			# The number of the interior dof 
			interior_dof = (ny - 1)*self.num_edge_dof + (nx - 1)

			return self.num_vertex + self.num_edges*self.num_edge_dof +\
					self.num_interior_dof*element + interior_dof

	## @brief assessor for the total degress of freedom
	def get_total_dof(self):
		return self.total_dof

	## @brief returns an array with the degrees of freedom on the boundary
	#  Only works for linear finite elements
	def get_boundary_dofs(self ):

		# Initialize the list of elements and local dof 
		boundary_dofs = { }

		# Loop over the edge connectivity array
		for edge, elements in self.edge_connectivity.items() :
			
			# If only one element shares this edge than this is on the boundary of the domain
			if len(elements) == 1 :
				
				# The edge number in the element
				edge_num = elements[0][1]

				# The element number
				element = elements[0][0] 

				# Append the element
				boundary_dofs[ element ] = [] 
				
				# Loop over the two vertices of the edge
				for i in range(2):
					
					# Get the local degrees of freedom
					# The 2 is added in there because according to the numbering scheme 
					# the edge attached to the 1st and second node is the last edge 
					local_dof = (edge_num + i )%self.num_vertex_dof

					# Get the global degree of freedom
					global_dof = self.get_global_dof( element, local_dof )

					# Append it
					boundary_dofs[element].append( local_dof )
		
		global_boundary_dofs = [] 

		for e, local_dofs in  boundary_dofs.items() :
			for local_dof in local_dofs:
				
				global_boundary_dofs.append( self.get_global_dof( e, local_dof ) ) 
			
		return np.array(global_boundary_dofs)
