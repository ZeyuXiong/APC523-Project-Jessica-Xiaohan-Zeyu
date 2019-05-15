import numpy as np 
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D



def f(x,y):
	#return np.sin(np.pi*x)*np.sin(np.pi*y)
	t = 0.01
	return -1./2/np.pi/t/t*np.exp(-( (x-0.5)**2+(y-0.5)**2 )/2/t/t)
# Nx = Ny = N

N = 40
h = 1./N
print(h)
K = np.zeros(((N+1)*(N+1),(N+1)*(N+1)))
F = np.zeros((N+1)*(N+1))


#F[(N//2+1)*N+(N//2+1)] = 1./h/h
#K is coefficient matrix
for i in range(N+1):
	for j in range(N+1):
		F[i*(N+1)+j] = f(i*h,j*h)
		if i != 0:
			K[i*(N+1)+j][(i-1)*(N+1)+j] = 1./h/h
		if i != N:
			K[i*(N+1)+j][(i+1)*(N+1)+j] = 1./h/h
		if j != 0:
			K[i*(N+1)+j][i*(N+1)+j-1] = 1./h/h
		if j != N:
			K[i*(N+1)+j][i*(N+1)+j+1] = 1./h/h
		K[i*(N+1)+j][i*(N+1)+j] = -4./h/h

#print(K)

# apply bc
for i in range(N+1):
	K[i*(N+1)] = np.zeros((N+1)*(N+1))
	K[i*(N+1)+N] = np.zeros((N+1)*(N+1))
	K[i*(N+1)][i*(N+1)] = 1
	F[i*(N+1)] = 0
	K[i*(N+1)+N][i*(N+1)+N] = 1
	F[i*(N+1)+N] = 0

for j in range(N+1):
	K[j] = np.zeros((N+1)*(N+1))
	K[(N+1)*N+j] = np.zeros((N+1)*(N+1))
	K[j][j] = 1
	F[j] = 0
	K[N*(N+1)+j][N*(N+1)+j]=1
	F[N*(N+1)+j] = 0

#print(K)
#print("rank ", np.linalg.matrix_rank(K))

u = np.linalg.solve(K,F)
x = h*np.arange(N+1)

fig = plt.figure()
ax = fig.gca(projection='3d')

u = u.reshape(N+1,N+1)
xx,yy = np.meshgrid(x,x)
# ax.plot_surface(xx,yy,u)
# plt.savefig('N=180.png')
# print(np.max(u))
print(type(u))
print(type(K))
#plt.pcolormesh(xx,yy,u)
#plt.contour(xx,yy,u,levels=[0.01,0.02,0.03,0.04])
#plt.show()
#plt.figure(3)
#plt.plot(x,u[N//2])
# tot = 0
# m = np.asmatrix(u)
# for i in range(N+1):
# 	n = N//4
# 	tot = tot + m[n, i] * h
# #
# error = np.abs(tot-0.068184116)
# #
tot = 0
for i in range(N+1):
	n = N//4
	tot = tot + u[n, i] * h
error = np.abs(tot-0.068184116)
print(0.01, error)
print(np.max(u))
# num_modes = 40
# u_exact = 0*xx
# for m in range(1,num_modes+1):
# 	for n in range(1,num_modes+1):
# 		u_mn = 4*np.sin(m*np.pi/2)*np.sin(n*np.pi/2)/(m*m+n*n)/np.pi**2
# 		#print("(m,n,u_mn) =", m,n,u_mn)
# 		u_exact += u_mn*np.sin(m*np.pi*xx)*np.sin(n*np.pi*yy)
# plt.figure(2)
# plt.contour(xx,yy,u_exact,levels=[0.01,0.02,0.03,0.04])
# plt.show()
#
# Diff = u - u_exact
# Diff[0,0:N+1] = Diff[0,0:N+1]/2
# Diff[N,0:N+1] = Diff[N,0:N+1]/2
# Diff[0:N+1,0] = Diff[0:N+1,0]/2
# Diff[0:N+1,N] = Diff[0:N+1,N]/2

# L1_err = h*h*np.sum(np.sum(abs(Diff)))
# print("L1_err:", L1_err)