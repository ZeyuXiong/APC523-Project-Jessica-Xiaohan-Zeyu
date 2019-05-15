import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.tri as tri
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

def plot_contour( coordinates, values, space_dim= 2 ):

	assert space_dim == 2, 'only works for two space dim'

	x = coordinates[:,0] 
	y = coordinates[:,1] 
	
	# Create the Triangulation just for plotting
	triang = tri.Triangulation(x, y)

	#	Plot contour
	plt.figure()
	plt.gca().set_aspect('equal')
	plt.tricontourf(triang, values)
	
	return


def plot_contour_3d( coordinates, values, space_dim= 2 ):

	assert space_dim == 2, 'only works for two space dim'

	x = coordinates[:,0] 
	y = coordinates[:,1] 
	z = values
	# Create the Triangulation just for plotting
	triang = tri.Triangulation(x, y)

	#	Plot contour
	fig = plt.figure(figsize=plt.figaspect(1.))
	ax = fig.add_subplot(1, 1, 1, projection='3d')
	ax.plot_trisurf(triang, values,cmap=plt.cm.Spectral)	
# 	ax.set_zlim(-0.25, 0.25)   
	ax.auto_scale_xyz([np.min(x), np.max(x)],[np.min(y), np.max(y)],[np.min(z), np.max(z)])
	ax.set_aspect('equal')
    


# 	plt.colorbar()

	return

        
def plot_quad_mesh( coordinates, connectivity ):
	
	for conn_e in connectivity: 
		
		crds = coordinates[ conn_e ] 

		plt.plot( crds[:,0], crds[:,1] , '-k',linewidth=1)

	return

