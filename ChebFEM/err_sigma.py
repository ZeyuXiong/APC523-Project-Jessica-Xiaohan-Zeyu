import numpy as np 
import matplotlib.pyplot as plt

sigma = np.array([1, 0.75, 0.5,0.25,0.1,0.05, 0.025])
err = np.array([0.006016951068882793, 0.010272016345216019, 0.0206793402,0.0501246045,0.06798370862,0.0681841134,0.0681840042 ])
exact = 0.0681841165
plt.semilogy(sigma,abs(err-exact),'.b')
plt.semilogy(sigma,abs(err-exact),'--r')
plt.xlabel("sigma")
plt.ylabel("error with analytical soln")

plt.show()