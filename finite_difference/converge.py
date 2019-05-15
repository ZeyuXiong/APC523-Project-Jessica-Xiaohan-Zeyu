import numpy as np
import matplotlib.pyplot as plt

def f(x,y):
    t = 0.01
    return -1./2/np.pi/t/t*np.exp(-( (x-0.5)**2+(y-0.5)**2 )/2/t/t)

lis =[10, 20, 40, 80, 160]
record=np.zeros(len(lis))
st=0
for N in lis:
    h = 1. / N
    print(h)
    K = np.zeros(((N + 1) * (N + 1), (N + 1) * (N + 1)))
    F = np.zeros((N + 1) * (N + 1))

    # F[(N//2+1)*N+(N//2+1)] = 1./h/h
    # K is coefficient matrix
    for i in range(N + 1):
        for j in range(N + 1):
            F[i * (N + 1) + j] = f(i * h, j * h)
            if i != 0:
                K[i * (N + 1) + j][(i - 1) * (N + 1) + j] = 1. / h / h
            if i != N:
                K[i * (N + 1) + j][(i + 1) * (N + 1) + j] = 1. / h / h
            if j != 0:
                K[i * (N + 1) + j][i * (N + 1) + j - 1] = 1. / h / h
            if j != N:
                K[i * (N + 1) + j][i * (N + 1) + j + 1] = 1. / h / h
            K[i * (N + 1) + j][i * (N + 1) + j] = -4. / h / h

    # print(K)

    # apply bc
    for i in range(N + 1):
        K[i * (N + 1)] = np.zeros((N + 1) * (N + 1))
        K[i * (N + 1) + N] = np.zeros((N + 1) * (N + 1))
        K[i * (N + 1)][i * (N + 1)] = 1
        F[i * (N + 1)] = 0
        K[i * (N + 1) + N][i * (N + 1) + N] = 1
        F[i * (N + 1) + N] = 0

    for j in range(N + 1):
        K[j] = np.zeros((N + 1) * (N + 1))
        K[(N + 1) * N + j] = np.zeros((N + 1) * (N + 1))
        K[j][j] = 1
        F[j] = 0
        K[N * (N + 1) + j][N * (N + 1) + j] = 1
        F[N * (N + 1) + j] = 0

    # print(K)
    # print("rank ", np.linalg.matrix_rank(K))

    u = np.linalg.solve(K, F)
    x = h * np.arange(N + 1)

    #fig = plt.figure()
    #ax = fig.gca(projection='3d')

    u = u.reshape(N + 1, N + 1)
    xx, yy = np.meshgrid(x, x)
    # ax.plot_surface(xx,yy,u)
    # plt.savefig('N=180.png')
    # print(np.max(u))
    error = abs(np.max(u)-0.6264281590815105)/0.6264281590815105
    record[st] = error
    st = st+1
    print(st)

print(record)
