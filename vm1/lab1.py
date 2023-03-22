import numpy as np
import matplotlib.pyplot as plt


def lagrange_basis(i, x, x_nodes):
    x_nodes_except_i = np.delete(x_nodes, i)
    return np.prod((x - x_nodes_except_i) / (x_nodes[i] - x_nodes_except_i))


def lagrange_interpolant(x, x_nodes, y_nodes):
    n = len(x_nodes)
    return np.sum(y_nodes* np.array([lagrange_basis(i, x , x_nodes) for i in range(n)]))


def plot_data_and_interpolant(x_nodes, f, i):
    y_nodes = f(x_nodes)
    fig, ax = plt.subplots(1, 1, figsize=(12, 6))
    x_for_plotting = np.linspace(-1, 1, 200)
    ax.plot(x_nodes, f(x_nodes), 'ro', markersize=10) #узлы
    ax.plot(x_for_plotting, f(x_for_plotting), '-', color='#aaa', linewidth=4) #о
    ax.plot(x_for_plotting, [lagrange_interpolant(x, x_nodes, y_nodes) for x in x_for_plotting], 'g-', linewidth=4)
    ax.grid()

if __name__ == '__main__':
    f = lambda x: 1./(1. + 25*x**2)
    n = [5, 8, 11, 14, 17, 20, 23]
    for i in range(7):
        x_equi_nodes = np.linspace(-1, 1, n[i])
        plot_data_and_interpolant(x_equi_nodes, f, i)
    plt.show()

