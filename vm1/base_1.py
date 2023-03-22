import random
import matplotlib.pyplot as plt
import numpy as np

def l_i(i, x, x_nodes):
    x_nodes_i = np.r_[x_nodes[:i],x_nodes[i+1:]]
    return np.prod((x-x_nodes_i)/(x_nodes[i]-x_nodes_i))

def L(x, x_nodes, y_nodes):
    return np.sum(y_nodes*np.array([l_i(i,x,x_nodes) for i in range(x_nodes.size)]))

def plot_data(x_node, f):
    y_node = f(x_node)
    fig, ax = plt.subplots(1, 1, figsize = (12, 6))
    #определение точек f(x)
    x_plt = np.linspace(-1,1,200)
    ax.plot(x_node, f(x_node), 'ro', markersize = 10)
    ax.plot(x_plt, f(x_plt), 'grey')

    #определение точек интерполянта ~f(x)
    ax.plot(x_plt, [L(i,x_node,y_node) for i in x_plt])

    ax.grid()
    plt.show()


def uniform():
    f = lambda x: 1. / (1. + 25 * x ** 2)
    for n in range(5, 24, 3):
        x = np.linspace(-1,1,n)
        plot_data(x,f)

def chebishev():
    f = lambda x: 1. / (1. + 25 * x ** 2)
    for n in range(5, 24, 3):
        x = np.array([np.cos((2*i-1)/(2*n)*np.pi) for i in range(1, n+1)])
        plot_data(x,f)
