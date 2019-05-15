import numpy as np 
import numpy.linalg as la 
from scipy.fftpack import dct 


def chebcoeffs(fval):
	N = len(fval)-1
	coeffs = dct(fval,type = 1)/N
	coeffs[0] /= 2
	coeffs[N] /= 2
	return coeffs

def diff(N):
	invD = np.eye(N)
	invD[0][0] = 2
	for i in range(N-2):
		invD[i][i+2] = -1
	upperright = np.dot(la.inv(invD),np.diag(np.arange(1,N+1)*2))
	Mat = np.block([ [np.zeros((N,1)), upperright], [ np.array([[0]]), np.zeros((1,N)) ] ] )
	return Mat

def get_val(f_hat,x):
	N = len(f_hat)-1
	f = 0*x
	for i in range(N+1):
		f += f_hat[i]*np.cos(i*np.arccos(x))
	return f
