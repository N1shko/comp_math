import numpy as np
import matplotlib.pyplot as plt
import intvalpy as ip
from functools import lru_cache
# a = 3
# b = 0.002
# g = 0.5
# d = 0.0006
k = 1 
a = -0.9
b = 0.5
c = 0.7
d = -0.8
#k = ip.Interval(0.999, 1.001)

#1st ode system



def get_cmap(n, name='hsv'):
    '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct 
    RGB color; the keyword argument name must be a standard mpl colormap name.'''
    return plt.cm.get_cmap(name, n)


def func(x):
    #return np.array([-k * x[0], 2*k*x[0]])
    return np.array([a * x[0] + b * x[0] * x[1], c * x[1] + d * x[0] * x[1]])


def rk4(x_0, t_n, f, h):
    t = np.arange(0, t_n, h)
    n = t.size
    kolvo = len(x_0)
    x = np.zeros((kolvo, n))
    x[:, 0] = x_0
    print(x)
    for i in range(n - 1):
        k1 = h * f(t[i], x[:, i])
        k2 = h * f(t[i] + h / 2, x[:, i] + k1 / 2)
        k3 = h * f(t[i] + h / 2, x[:, i] + k2 / 2)
        k4 = h * f(t[i] + h, x[:, i] + k3)
        koef = (k1 + 2 * k2 + 2 * k3 + k4) / 6
        x[:, i + 1] = x[:, i] + koef
    print(x)
    return x


def rk4_iteration(x_0, t, f, h):
    k1 = h * f(t, x_0)
    k2 = h * f(t + h / 2, x_0 + k1 / 2)
    k3 = h * f(t + h / 2, x_0 + k2 / 2)
    k4 = h * f(t + h, x_0 + k3)
    koef = (k1 + 2 * k2 + 2 * k3 + k4) / 6
    return x_0 + koef


def lagrange_basis(i_offset, cur_values, value0_bonds, p=2, m=2):
    res=1
    for i in range(p+1):
        for j in range(m):
            if i_offset[j] != i:
                res *= (p*(cur_values[j]-value0_bonds[j][0])/(value0_bonds[j][1]-value0_bonds[j][0]) - i) / (i_offset[j] - i)
    return res


def lagrange_interpolant(y, cur_values, value0_bounds, p=2, m=2):
    res = 0
    for i in range(p+1):
        for j in range(p+1):
            res += y[i][j] * lagrange_basis([i, j], cur_values, value0_bounds, p, m)
    return res


if __name__ == '__main__':
    f = lambda t, x: func(x)
    t_n = 5
    h = 0.5
    time = np.arange(0, t_n, h)
    koef1 = np.arange(0.5, 2.5 + h, h)
    koef2 = np.arange(0.5, 2.5 + h, h)
    t = np.arange(0, t_n, h)
    n = t.size
    print(koef1)
    #x_0 = [koef1[0] + (koef1[1] - koef1[0]) * 0 / 2, koef1[0] + (koef2[1] - koef2[0]) * 1 / 2]
    #x_0 = [0.5, 0.5]
            # ax.legend()
    #x_0 = [[0.5, 2.5], [0.5, 2.5]]
    # print('begining here')
    # x_0 = np.array([200, 200])
    # print('res = ', rk4(x_0, t_n, f, h))
            # ax.plot(t, res[0, :], color='blue')
    # plt.show()
    #x_0 = ip.Interval([[0, 0], [150, 150], [0, 0]])
    # h_int = 0.1
    fig, ax = plt.subplots(1, 2, figsize=(20, 6))
    # for i in range(20):
    #     for j in range(20):
    #         x_0 = [0.5 + h_int*i, 0.5 + h_int*j]
    temp_t = []
    temp_x = []
    temp_y = []
    gor_grid_x = []
    gor_grid_y = []
    ax[1].set_xlim(0, 3)
    ax[1].set_ylim(0, 6)
    ax[1].set_box_aspect(1)
    for i in range(len(koef1)):
        temp = []
        gor_grid_x.append([])
        gor_grid_y.append([])
        for j in range(len(koef2)):

            res1 = rk4([koef1[i], koef2[j]], t_n, f, h)
                #ax.plot(res1[0, :], res1[1, :], '-', color='steelblue')
            temp_x.append(res1[0])
            temp_y.append(res1[1])
            ax[0].plot(time, res1[0, :], color='cornflowerblue')
            #ax[0].plot(time, res1[1, :], color='indigo')
            ax[1].scatter(res1[0][-1], res1[1][-1], color='red')
            #ax[1].plot()
            #print(res1[0][1])
            gor_grid_x[i].append(res1[0][-1])
            gor_grid_y[i].append(res1[1][-1])
        #ax[1].plot(vert_grid_x, vert_grid_y)
    #print('res[0, :]', res1[0, :])
    #print('res[1, :]', res1[1, :])
    #x_0 = [0.5, 0.5]
    #res1 = rk4(x_0, t_n, f, h)
    #ax.plot(res1[0, :], res1[1, :], '-', color='steelblue')
    ax[0].grid()
    ax[1].grid()
    print(gor_grid_x)
    for i in range(len(gor_grid_x)):
        ax[1].plot(gor_grid_x[i], gor_grid_y[i])
    for i in range(len(gor_grid_x)):
        x_res = [sublist[i] for sublist in gor_grid_x]
        y_res = [sublist[i] for sublist in gor_grid_y]
        ax[1].plot(x_res, y_res)
    #X, Y = np.meshgrid(gor_grid_x, gor_grid_y)
    #ax[1].plot(X, Y)
    plt.show()
    #print(temp_x)
    print(len(temp_x[0]))
    print(len(temp_x))
    res_x_min=np.array(temp_x[0])
    res_x_max=np.array(temp_x[0])
    for i in range(len(temp_x)):
        #res_x_min.append(temp_x[j][0])
        #res_x_max.append(temp_x[j][0])
        for j in range(len(temp_x[i])):
            if temp_x[i][j] < res_x_min[j]:
                res_x_min[j] = temp_x[i][j]
            if temp_x[i][j] > res_x_max[j]:
                res_x_max[j] = temp_x[i][j]
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))

    ax.plot(time, res_x_min, color ="red")
    ax.plot(time, res_x_max, color = "blue")
    ax.grid()
    plt.show()
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    print(gor_grid_x)
    ax.scatter(gor_grid_x, gor_grid_y, color ="red")
    ax.grid()

    p = int((koef1[-1] - koef1[0])//h)
    print(p)
    k1_for_plotting = np.linspace(0.5, 2.5, 101)
    k2_for_plotting = np.linspace(0.5, 2.5, 101)
    print(k1_for_plotting)
    for k1 in k1_for_plotting:
        x_for_plotting = []
        y_for_plotting = []
        for k2 in k2_for_plotting:
            if k1 in koef1:
                x_for_plotting.append(lagrange_interpolant(gor_grid_x, [k1, k2], [[0.5, 2.5], [0.5, 2.5]], p, 2))
                y_for_plotting.append(lagrange_interpolant(gor_grid_y, [k1, k2], [[0.5, 2.5], [0.5, 2.5]], p, 2))
                ax.plot(x_for_plotting, y_for_plotting)
    for k2 in k2_for_plotting:
        x_for_plotting = []
        y_for_plotting = []
        for k1 in k1_for_plotting:
            if k2 in koef2:
                x_for_plotting.append(lagrange_interpolant(gor_grid_x, [k1, k2], [[0.5, 2.5], [0.5, 2.5]], p, 2))
                y_for_plotting.append(lagrange_interpolant(gor_grid_y, [k1, k2], [[0.5, 2.5], [0.5, 2.5]], p, 2))
                ax.plot(x_for_plotting, y_for_plotting)
    plt.show()
#     ax.grid()
#     ax.set_xlabel('t')
#     ax.set_ylabel(r'$x,~y$')
#     # ax.plot(time, res1[0, :], color='cornflowerblue', label=r'$x(t)$')
#     # ax.plot(time, res1[1, :], color='indigo', label=r'$y(t)$')
#    # ax.plot(time, res1[2, :], color='red', label=r'$y(t)$')
#     ax.legend()
#     plt.savefig('xy_graf.png', dpi=300)

#     plt.show()