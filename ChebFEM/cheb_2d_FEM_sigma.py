import numpy as np 
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import axes3d, Axes3D
from cheb import chebcoeffs
from cheb import diff


def f(x,y):
	#return np.sin(np.pi*x)*np.sin(np.pi*y)
	t = 0.05
	return -1./2/np.pi/t/t*np.exp(-( (x-0.5)**2+(y-0.5)**2 )/2/t/t)
M = 10
N = 10
j = np.arange(M+1)
node_X = np.cos(j*np.pi/M)


connectivity = []
coordinates = []
for i in range(N*N):
	I = i//N
	J = i%N
	connectivity.append([I*N+J,I*N+J+N, J*N+I, J*N+I+N])
	coordinates.append(float(I)/N)
for i in range(N):
	coordinates.append(1.0)

#print(connectivity,coordinates)


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
#print(f_hat)


D = diff(M)
In = np.eye(M+1) 
K = np.zeros(( (M+1)*(M+1)*N*N, (M+1)*(M+1)*N*N ))
KD2 = np.zeros(( (M+1)*(M+1)*N*N, (M+1)*(M+1)*N*N ))

left_bc1 = (-1)**j
right_bc1 = np.ones(M+1)
right_bc2 = j*j
left_bc2 = -left_bc1*right_bc2

for e in range(N*N):
	I = e//N
	J = e%N

	#if I == 0 or I == N-1 or J == 0 or J == N-1:
	#	continue
	xa = coordinates[connectivity[e][2]]
	xb = coordinates[connectivity[e][3]]
	ya = coordinates[connectivity[e][0]]
	yb = coordinates[connectivity[e][1]]

	#print("element", e, ", (xa,xb,ya,yb) = (", xa,xb,ya,yb, ")")

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

#print("rank ", np.linalg.matrix_rank(K))
#for e in range(N*N):
#	print(K[e*(M+1)*(M+1):(e+1)*(M+1)*(M+1), e*(M+1)*(M+1):(e+1)*(M+1)*(M+1)])

#for i in range(len(K)):
#	for j in range(len(K)):
#		print(K[i][j], end = "\t")
#	print()

u_hat = np.linalg.solve(K,f_hat)
K = 0
Ns = 10001
x = np.linspace(0,1,Ns)
y = np.linspace(0,1,Ns)
y0 = 0.25
xx, yy = np.meshgrid(x,y)
u = 0*x #np.zeros((Ns,Ns))
for i in range(Ns):
	for e in range(N*N):
		xa = coordinates[connectivity[e][2]]
		xb = coordinates[connectivity[e][3]]
		ya = coordinates[connectivity[e][0]]
		yb = coordinates[connectivity[e][1]]

		if (x[i] <= xb and x[i] >= xa) and (y0 <= yb and y0 >= ya):
			break
	Xe = (x[i]-(xa+xb)/2)*2/(xb-xa)/(1.+1e-12)
	Ye = (y0-(ya+yb)/2)*2/(yb-ya)/(1.+1e-12)
		
	for m in range(M+1):
		for n in range(M+1):
			u[i] += u_hat[e*(M+1)*(M+1)+m*(M+1)+n]*np.cos(m*np.arccos(Xe))*np.cos(n*np.arccos(Ye))


h = 1./(Ns-1)
#u025 = u[Ns//4]
u[0] /=2
u[Ns-1] /=2

u_int = np.sum(u)*h
print("N =", N)
print("integral:",u_int)
#print(x[Ns//4])
#print("max:", np.max(np.max(u)))
