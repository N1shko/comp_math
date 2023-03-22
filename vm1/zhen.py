import random
import matplotlib.pyplot as plt
import numpy as np
from base_1 import *
class function:
    def __init__(self):
        n = np.random.randint(7, 15)
        m = np.random.randint(7, 15)
        j = np.random.rand(m)
        k = np.random.rand(n)
        self.n = n
        self.m = m
        self.j = j
        self.k = k
    def read(self):
        return self.n, self.m, self.j, self.k


def generation():
    n = random.randint(7, 15)

    m = random.randint(7, 15)
    j = [random.randint(0, 1000) / 1000 for i in range(0, m + 1)]
    k = [random.randint(0, 1000) / 1000 for i in range(0, n + 1)]
    return n, m, j, k
def func_random(x,m, n, j, k):
    a = []
    b = []
    for i in range(len(j)):
        a.append(j[i]*x**i)
    for i in range(len(k)):
        b.append(k[i]*x**i)
    A = np.array(a)
    B = np.array(b)
    return np.sum(A)/(1+(np.sum(B)))
def tsk1_2(flag):
    functions = [function() for i in range(100)]
    N = 20
    # x_node = np.linspace(-1, 1, n)
    # y_node = f(x_node)

    fig, axes = plt.subplots(4, 2, figsize=(24, 12))
   # ax = []
    l = 0
    funk_res = []
    for ax in axes:
        i = 25*l
        l += 1
        n,m,j,k = functions[i].read()
        funk_res.append(functions[i])
        #ax.append(fig.add_subplot())
        # определение точек f(x)
        x_plt = np.linspace(-1, 1, 100)
        y_plt = np.array([func_random(x, n, m, j, k) for x in x_plt])
        x_node = np.linspace(-1,1,N)#np.array([x_plt[x] for x in range(1, len(x_plt)+1,100//N)])
        y_node = np.array([func_random(x, n, m, j, k) for x in x_node])
        x_chebishev = np.array([np.cos((2 * i - 1) / (2 * N) * np.pi) for i in range(1, N + 1)])
        y_chebishev = np.array([func_random(x, n, m, j, k) for x in x_chebishev])
        if(flag):
            #ax[0].set(title = (str(l)+'узлы равномерные'))
            ax[0].plot(x_plt, [L(x, np.array(x_node), np.array(y_node)) for x in np.array(x_plt)], 'blue', label = 'интерполянт')
            ax[0].plot(x_node, y_node, 'ro', markersize=5, label = 'узлы')
            ax[0].plot(x_plt, y_plt, 'grey', label = 'график')
            ax[0].grid()
            ax[0].legend()

            #ax[1].set(title=(str(l)+'узлы чебышева'))
            ax[1].plot(x_plt, [L(x, np.array(x_chebishev), np.array(y_chebishev)) for x in np.array(x_plt)], 'blue', label = 'интерполянт')
            ax[1].plot(x_chebishev, y_chebishev, 'ro', markersize=3, label = 'узлы')
            ax[1].plot(x_plt, y_plt, 'grey', label = 'график')
            ax[1].grid()
            ax[1].legend()




    #if(flag):
        #plt.show()
    #print(func_random(1, n, m, j, k))
    return funk_res
def tsk3(functions):
    N = 30
    #N_i = 0
    l = -1

    fig, axes = plt.subplots(4, 3, figsize=(100, 50))
    for ax in axes:
        l += 1
        n, m, j, k = functions[l].read()
        x_plt = np.linspace(-1, 1, 200)
        y_plt = np.array([func_random(x, n, m, j, k) for x in x_plt])
        x_node = np.linspace(-1,1,N)#np.array([x_plt[x] for x in range(0, len(x_plt), 100 // N)])
        y_node = np.array([func_random(x, n, m, j, k) for x in x_node])
        x_chebishev = np.array([np.cos((2 * i - 1) / (2 * N) * np.pi) for i in range(1, N + 1)])
        y_chebishev = np.array([func_random(x, n, m, j, k) for x in x_chebishev])

        ax[0].plot(x_plt, [L(x, np.array(x_node), np.array(y_node)) for x in np.array(x_plt)], 'blue',
                   label='интерполянт равномерный')
        ax[0].plot(x_plt, [L(x, np.array(x_chebishev), np.array(y_chebishev)) for x in np.array(x_plt)], 'green',
                   label='интерполянт чебышева')
        ax[0].plot(x_node, y_node, 'ro', markersize=5, label='узлы')
        ax[0].plot(x_plt, y_plt, 'grey', label='график')
        ax[0].grid()
        ax[0].legend()

        x_ro = np.array([N_i for N_i in range(1,31)])
        y_ro_1 = []
        y_ro_2 = []
        for N_i in x_ro:
            x_node = np.linspace(-1,1,N_i)#np.array([x_plt[x] for x in range(0, len(x_plt), 100 // N_i)])
            #print([x_plt[x] for x in range(0, len(x_plt), 100 // N_i)])
            y_node = np.array([func_random(x, n, m, j, k) for x in x_node])

            x_chebishev = np.array([np.cos((2 * i - 1) / (2 * N_i) * np.pi) for i in range(1, N_i + 1)])
            y_chebishev = np.array([func_random(x, n, m, j, k) for x in x_chebishev])

            y = max([abs(L(x, np.array(x_node), np.array(y_node))-func_random(x, n, m, j, k)) for x in np.array(x_plt)])
            y_ro_1.append(y)
            y_ch = max([abs(L(x, np.array(x_chebishev), np.array(y_chebishev))-func_random(x, n, m, j, k)) for x in np.array(x_plt)])
            y_ro_2.append(y_ch)
        ax[1].semilogy(x_ro, y_ro_1, 'blue', label='график')
        ax[2].semilogy(x_ro, y_ro_2, 'green', label='график')
    #plt.show()

    return 0