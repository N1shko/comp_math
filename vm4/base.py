import numpy as np
import matplotlib.pyplot as plt
import intvalpy as ip
from self_made_kd_tree import KDTree
from calc_func import *



k = 1 
a = -0.9
b = 0.5
c = 0.7
d = -0.8


def func(x):
    #return np.array([-k * x[0], 2*k*x[0]])
    return np.array([a * x[0] + b * x[0] * x[1], c * x[1] + d * x[0] * x[1]])


if __name__ == '__main__':
    p = 4
    m = 2
    tree = KDTree([[0.5, 2.5], [0.5, 2.5]], m, p)
    f = lambda t, x: func(x)
    t_n = 15
    h = 0.5
    t_h = 0.5
    t = np.arange(0, t_n, t_h)
    n = t.size
    fig, ax = plt.subplots(1, 1, figsize=(20, 20))
    gor_grid_x = []
    gor_grid_y = []
    max_val = []
    min_val = []
    is_set = False
    for dot in tree.nodes[0].dots:
        ttemp_x = [dot[0]]
        ttemp_y = [dot[1]]
        temp0_x = dot[0]
        temp0_y = dot[1]
        for k in range(len(t)):
            res = rk4_iteration([temp0_x, temp0_y], t[k], f, t_h)
            ttemp_x.append(res[0])
            ttemp_y.append(res[1])
            temp0_x = res[0]
            temp0_y = res[1]
            if not is_set:
                min_val.append([res[0], res[1]])
                max_val.append([res[0], res[1]])
            else:
                min_val[k][0] = min(min_val[k][0], res[0])
                min_val[k][1] = min(min_val[k][1], res[1])
                max_val[k][0] = max(max_val[k][0], res[0])
                max_val[k][1] = max(max_val[k][1], res[1])
        is_set = True
        ax.scatter(ttemp_x[n-1::n], ttemp_y[n-1::n])
        gor_grid_x.append(ttemp_x[n-1])
        gor_grid_y.append(ttemp_y[n-1])

    # ax.plot(t, [x[0] for x in min_val])
    # ax.plot(t, [x[0] for x in max_val])
    ax.grid()
    for i in range(p + 1):
        ax.plot(gor_grid_x[i*(p+1):(i+1)*(p+1)], gor_grid_y[i*(p+1):(i+1)*(p+1)])
        ax.plot(gor_grid_x[i::p+1], gor_grid_y[i::p+1])
    x_for_plotting = []
    y_for_plotting = []
    print(len(gor_grid_x))
    i = 0
    for dot in tree.nodes[0].plot_dots:        
        x_for_plotting.append(lagrange_interpolant(gor_grid_x, dot, tree.nodes[0].borders, p, m))
        y_for_plotting.append(lagrange_interpolant(gor_grid_y, dot, tree.nodes[0].borders, p, m))
    num = tree.nodes[0].plot_dots_num
    h_plot = num // p
    for i in range(0, num+1, h_plot):
        ax.plot(x_for_plotting[i*(num+1):(i+1)*(num+1)], y_for_plotting[i*(num+1):(i+1)*(num+1)])
        ax.plot(x_for_plotting[i::num+1], y_for_plotting[i::num+1])
    plt.show()