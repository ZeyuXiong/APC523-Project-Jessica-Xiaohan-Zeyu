import numpy as np 
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import axes3d, Axes3D
from cheb import chebcoeffs
from cheb import diff

# define the guass function
def f(x,y):
	t = 0.5
	return -1./2/np.pi/t/t*np.exp(-( (x-0.5)**2+(y-0.5)**2 )/2/t/t)

# M is the polynomial order (should be even otherwise will cause problems)
# N should be num of elements on one side
M = 10
N = 10

# node_X is Gauss Lobatto points
j = np.arange(M+1)
node_X = np.cos(j*np.pi/M)

# define connectivity and coordinates array
# label the elements as 
# -------------
# | 6 | 7 | 8 |
# -------------
# | 3 | 4 | 5 |
# -------------
# | 0 | 1 | 2 |
# -------------
# label the edges as
# - 9 - 10 - 11 -
# |   |    |    |
# - 6 - 7  - 8  -
# |   |    |    |
# - 3 - 4  - 5  -
# |   |    |    |
# - 0 - 1  - 2  -
connectivity = []
coordinates = []
for i in range(N*N):
	I = i//N
	J = i%N
	connectivity.append([I*N+J,I*N+J+N, J*N+I, J*N+I+N])
	coordinates.append(float(I)/N)
for i in range(N):
	coordinates.append(1.0)



# doing Chebyshev expansion and save the coefficients in f_hat
f_hat = np.array([])
for e in range(N*N):
	xa = coordinates[connectivity[e][2]]
	xb = coordinates[connectivity[e][3]]
	ya = coordinates[connectivity[e][0]]
	yb = coordinates[connectivity[e][1]]
	node_x = xb*(node_X+1)/2 + xa*(-node_X+1)/2
	node_y = yb*(node_X+1)/2 + ya*(-node_X+1)/2
	f_m = np.zeros((M+1,M+1))
	for i in range(M+1):
		f_m[i] = chebcoeffs(f(node_x,node_y[i]))

	f_m = np.transpose(f_m)
	fe = np.array([])
	for i in range(M+1):
		fe = np.append(fe, chebcoeffs(f_m[i]))
	f_hat = np.append(f_hat, fe)


# define differential matrix
D = diff(M)
In = np.eye(M+1) 
K = np.zeros(( (M+1)*(M+1)*N*N, (M+1)*(M+1)*N*N ))
KD2 = np.zeros(( (M+1)*(M+1)*N*N, (M+1)*(M+1)*N*N ))

left_bc1 = (-1)**j
right_bc1 = np.ones(M+1)
right_bc2 = j*j
left_bc2 = -left_bc1*right_bc2

# assenmble the matrix K with conectivity
for e in range(N*N):
	I = e//N
	J = e%N

	xa = coordinates[connectivity[e][2]]
	xb = coordinates[connectivity[e][3]]
	ya = coordinates[connectivity[e][0]]
	yb = coordinates[connectivity[e][1]]

	Dx = np.kron(D*2/(xb-xa),In)
	Dy = np.kron(In,D*2/(yb-ya))
	D2 = np.dot(Dx,Dx)+np.dot(Dy,Dy)
	K[e*(M+1)*(M+1):(e+1)*(M+1)*(M+1), e*(M+1)*(M+1):(e+1)*(M+1)*(M+1)] = D2
	KD2[e*(M+1)*(M+1):(e+1)*(M+1)*(M+1), e*(M+1)*(M+1):(e+1)*(M+1)*(M+1)] = D2

	if I != N-1:
		for i in range(1,(M+1)//2+1):
			K[e*(M+1)*(M+1)+i*(M+1)+M-1] = np.zeros( (M+1)*(M+1)*N*N )
			K[e*(M+1)*(M+1)+i*(M+1)+M] = np.zeros( (M+1)*(M+1)*N*N )
			K[e*(M+1)*(M+1)+i*(M+1)+M-1][e*(M+1)*(M+1)+(i-1)*(M+1):e*(M+1)*(M+1)+i*(M+1)] = -right_bc1
			K[e*(M+1)*(M+1)+i*(M+1)+M][e*(M+1)*(M+1)+(i-1)*(M+1):e*(M+1)*(M+1)+i*(M+1)] = -right_bc2
			K[e*(M+1)*(M+1)+i*(M+1)+M-1][(e+N)*(M+1)*(M+1)+(i-1)*(M+1):(e+N)*(M+1)*(M+1)+i*(M+1)] = left_bc1
			K[e*(M+1)*(M+1)+i*(M+1)+M][(e+N)*(M+1)*(M+1)+(i-1)*(M+1):(e+N)*(M+1)*(M+1)+i*(M+1)] = left_bc2
			f_hat[e*(M+1)*(M+1)+i*(M+1)+M-1] = 0
			f_hat[e*(M+1)*(M+1)+i*(M+1)+M] = 0
		i = (M+1)//2+1
		K[e*(M+1)*(M+1)+M] = np.zeros( (M+1)*(M+1)*N*N )
		f_hat[e*(M+1)*(M+1)+M] = 0
		K[e*(M+1)*(M+1)+M][e*(M+1)*(M+1)+(i-1)*(M+1):e*(M+1)*(M+1)+i*(M+1)] = right_bc1
		K[e*(M+1)*(M+1)+M][(e+N)*(M+1)*(M+1)+(i-1)*(M+1):(e+N)*(M+1)*(M+1)+i*(M+1)] = -left_bc1

	if J != N-1:
		for i in range((M+1)//2+1, M+1):
			K[e*(M+1)*(M+1)+i*(M+1)+M-1] = np.zeros( (M+1)*(M+1)*N*N )
			K[e*(M+1)*(M+1)+i*(M+1)+M] = np.zeros( (M+1)*(M+1)*N*N )
			f_hat[e*(M+1)*(M+1)+i*(M+1)+M-1] = 0
			f_hat[e*(M+1)*(M+1)+i*(M+1)+M] = 0
			L = i - (M+1)//2 - 1
			for k in range(M+1):
				K[e*(M+1)*(M+1)+i*(M+1)+M-1][e*(M+1)*(M+1)+L+k*(M+1)] = -1
				K[e*(M+1)*(M+1)+i*(M+1)+M-1][(e+1)*(M+1)*(M+1)+L+k*(M+1)] = (-1)**k
				K[e*(M+1)*(M+1)+i*(M+1)+M][e*(M+1)*(M+1)+L+k*(M+1)] = k*k
				K[e*(M+1)*(M+1)+i*(M+1)+M][(e+1)*(M+1)*(M+1)+L+k*(M+1)] = (-1)**k * k*k

# apply bc
for I in range(N):
	e = I*N + N-1
	for i in range( (M+1)//2+1, M+1 ):
		K[e*(M+1)*(M+1)+i*(M+1)+M-1] = np.zeros( (M+1)*(M+1)*N*N )
		K[e*(M+1)*(M+1)+i*(M+1)+M] = np.zeros( (M+1)*(M+1)*N*N )
		f_hat[e*(M+1)*(M+1)+i*(M+1)+M-1] = 0
		f_hat[e*(M+1)*(M+1)+i*(M+1)+M] = 0
		L = i - (M+1)//2 - 1
		for k in range(M+1):
			K[e*(M+1)*(M+1)+i*(M+1)+M-1][e*(M+1)*(M+1)+L+k*(M+1)] = 1
			K[e*(M+1)*(M+1)+i*(M+1)+M][I*N*(M+1)*(M+1)+L+k*(M+1)] = (-1)**k
for J in range(N):
	e = N*(N-1) + J
	for i in range(1, (M+1)//2+1):
		K[e*(M+1)*(M+1)+i*(M+1)+M-1] = np.zeros( (M+1)*(M+1)*N*N )
		K[e*(M+1)*(M+1)+i*(M+1)+M] = np.zeros( (M+1)*(M+1)*N*N )
		K[e*(M+1)*(M+1)+i*(M+1)+M-1][J*(M+1)*(M+1)+(i-1)*(M+1):J*(M+1)*(M+1)+i*(M+1)] = left_bc1
		K[e*(M+1)*(M+1)+i*(M+1)+M][e*(M+1)*(M+1)+(i-1)*(M+1):e*(M+1)*(M+1)+i*(M+1)] = right_bc1
		f_hat[e*(M+1)*(M+1)+i*(M+1)+M-1] = 0
		f_hat[e*(M+1)*(M+1)+i*(M+1)+M] = 0
	i = (M+1)//2+1
	K[e*(M+1)*(M+1)+M] = np.zeros( (M+1)*(M+1)*N*N )
	f_hat[e*(M+1)*(M+1)+M] = 0
	K[e*(M+1)*(M+1)+M][e*(M+1)*(M+1)+(i-1)*(M+1):e*(M+1)*(M+1)+i*(M+1)] = left_bc1



u_hat = np.linalg.solve(K,f_hat)
K = 0
Ns = 201
x = np.linspace(0,1,Ns)
y = np.linspace(0,1,Ns)
xx, yy = np.meshgrid(x,y)
u = np.zeros((Ns,Ns))
for i in range(Ns):
	for j in range(Ns):
		for e in range(N*N):
			xa = coordinates[connectivity[e][2]]
			xb = coordinates[connectivity[e][3]]
			ya = coordinates[connectivity[e][0]]
			yb = coordinates[connectivity[e][1]]

			if (x[i] <= xb and x[i] >= xa) and (y[j] <= yb and y[j] >= ya):
				break
		Xe = (x[i]-(xa+xb)/2)*2/(xb-xa)/(1.+1e-12)
		Ye = (y[j]-(ya+yb)/2)*2/(yb-ya)/(1.+1e-12)
		
		for m in range(M+1):
			for n in range(M+1):
				u[i][j] += u_hat[e*(M+1)*(M+1)+m*(M+1)+n]*np.cos(m*np.arccos(Xe))*np.cos(n*np.arccos(Ye))

plt.figure(1)
num_modes = 100
u_exact = 0*xx
for m in range(1,num_modes+1):
	for n in range(1,num_modes+1):
		u_mn = 4*np.sin(m*np.pi/2)*np.sin(n*np.pi/2)/(m*m+n*n)/np.pi**2
		u_exact += u_mn*np.sin(m*np.pi*xx)*np.sin(n*np.pi*yy)
plt.contour(xx,yy,u_exact,levels=[0.01,0.2,0.3,0.4])


plt.figure(2)
plt.contour(xx,yy,u,levels=[0.01,0.2,0.3,0.4])
#plt.figure(3)
#plt.plot(xx[Ns//2], u[Ns//2], xx[Ns//2], u_exact[Ns//2])

fig = plt.figure(3)
ax = Axes3D(fig)
ax.plot_surface(xx,yy,u)

plt.show()


#dx = 5
#Diff = u - u_exact
#Diff[0,0:Ns] = Diff[0,0:Ns]/2
#Diff[Ns-1,0:Ns] = Diff[Ns-1,0:Ns]/2
#Diff[0:Ns,0] = Diff[0:Ns,0]/2
#Diff[0:Ns,Ns-1] = Diff[0:Ns,Ns-1]/2
#Diff[Ns//2-dx:Ns//2+dx,Ns//2-dx:Ns//2+dx] = np.zeros((2*dx,2*dx))
#h = 1./(Ns-1)
#L1_err = h*h*np.sum(np.sum(abs(Diff)))
#print("L1_err:", L1_err)


h = 1./(Ns-1)
u025 = u[Ns//4]
u025[0] /=2
u025[Ns-1] /=2

u_int = np.sum(u025)*h
print("N =", N)
print("integral:",u_int)
#print(x[Ns//4])
print("max:", np.max(np.max(u)))
