import numpy as np 
import matplotlib.pyplot as plt

#N = np.array([2,4,5,6,7,8,9,10,12])
#L1_err = np.array([0.092023,0.0078232,0.0043006,0.0013221,0.00083341,0.00066108,0.00040118,0.00039229,0.00034661])


# sigma = 0.5
N = np.array([4,8,10,16])
L1_err = 10**np.array([-5.286727919,-6.25703933,-6.654038458,-7.590066877])
plt.subplot(321)
plt.loglog(N,L1_err,'.')
plt.loglog(N,1/N**6)
plt.xlabel("N")
plt.ylabel("L1_err @ y=0.25")
plt.legend(["sigma=0.5"])

# sigma = 0.25
N = np.array([4,8,10,16])
L1_err = 10**np.array([-5.249352777,-6.91649738,-7.2644011,-8.337242168])
plt.subplot(322)
plt.loglog(N,L1_err,'.')
plt.loglog(N,1/N**6)
plt.xlabel("N")
plt.ylabel("L1_err @ y=0.25")
plt.legend(["sigma=0.25","6th-order ref"])

# sigma = 0.1
N = np.array([13,14,15,16,17])
L1_err = 10**np.array([-7.214670165,-7.48148606,-7.991399828,-8.102372909,-7.764471553])
plt.subplot(323)
plt.loglog(N,L1_err,'.')
plt.loglog(N,1/N**6)
plt.xlabel("N")
plt.ylabel("L1_err @ y=0.25")
plt.legend(["sigma=0.1","6th-order ref"])


# sigma = 0.05
N = np.array([10,12,14,16,18])
L1_err = 10**np.array([-4.269730624,-5.176564898,-6.393618635,-8.167491088,-9.301029996])
plt.subplot(324)
plt.loglog(N,L1_err,'.')
plt.loglog(N,1/N**6)
plt.xlabel("N")
plt.ylabel("L1_err @ y=0.25")
plt.legend(["sigma=0.05","6th-order ref"])


# sigma = 0.01
N = np.array([18,20,22,24,26])
L1_err = 10**np.array([-2.479393297,-2.744811442,-3.132169789,-3.849233739,-4.013819219])
plt.subplot(325)
plt.loglog(N,L1_err,'.')
plt.loglog(N,1/N**6)
plt.xlabel("N")
plt.ylabel("L1_err @ y=0.25")
plt.legend(["sigma=0.01","6th-order ref"])

#plt.legend(["sigma=0.5","sigma=0.25","sigma=0.1","sigma=0.05","sigma=0.01"])
plt.show()