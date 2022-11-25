#!/usr/bin/env python3
"""
Utilidades varias para curso Teoría y Algoritma de Optimización
autor: Ignacio Ramírez (nacho@fing.edu.uy)
"""

import numpy as np

def proj_l1(x):
    """
    proyecta los puntos de entrada sobre la bola unitaria de norma 1
    Este algoritmo es casi idéntico al que se utiliza para proyectar
    sobre un simplex de probabilidad en el examen de 2021, ya que la
    bola L1 se puede dividir en simplexes, uno por cuadrante del espacio.
    :param x: vector de entrada
    :return: vector proyectado
    """
    if np.sum(np.abs(x)) < 1:
        return np.copy(x)
    n = len(x)
    a = np.abs(x)
    s = np.sign(x)
    idx = np.argsort(-a) # indices de valores ordenados: orden decreciente de a = orden creciente de -a
    y = a[idx]           # vectores ordenados
    # cumsum: suma acumulativa de y: w[0] = y[0], w[1] = y[0] + y[1], etc.
    # w_j = (sum_{i=0}^{i=j} y_j - 1) / j
    w = (np.cumsum(y) - 1) / np.arange(1, n + 1)
    for j in range(len(w)):
        yj = np.maximum(y - w[j], 0)
        if np.abs(np.sum(yj) - 1) < 1e-8:
            xp = s * yj[idx]
            return xp
    print('no deberiamos haber llegado aca.')
    return None

#
# si se ejecuta como script se corre una pruebita
#
if __name__ == '__main__':
    import numpy.random as rand
    import matplotlib.pyplot as plt
    import mpl_toolkits.mplot3d as plt3

    n = 3
    N = 1000
    X = rand.uniform(low=0,high=1.2,size=(N,n))
    #print('x',x)
    Xp = np.empty(X.shape)
    for i in range(N):
        Xp[i,:] = proj_l1(X[i,:])

    fig = plt.figure(figsize=(15,15))
    ax = fig.add_subplot(projection='3d')
    ax.scatter(X[:,0],X[:,1],X[:,2],color='black',label='orig',alpha=0.2)
    ax.scatter(Xp[:,0],Xp[:,1],Xp[:,2],color='blue',label='proj',alpha=0.2)
    ax.legend()
    plt.show()
