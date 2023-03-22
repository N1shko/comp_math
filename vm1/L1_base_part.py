import numpy as np
import matplotlib.pyplot as plt


def runge(x) -> float:
     return 1. / (1. + 25*x**2)


def chebish(n) -> list:
    return([np.cos((2 * i - 1) / (2 * n) * np.pi) for i in range(1, n+1)])


def l_i(i, x, x_nodes) -> float:
    iter = np.delete(x_nodes, i)
    return np.prod((x - iter) / (x_nodes[i] - iter))


def lagrange_interpolant(x, x_nodes, y_nodes) -> float:
    n = len(x_nodes)
    return np.sum(y_nodes * np.array([l_i(i, x, x_nodes) for i in range(n)]))


def plotting(x, n, i, axs, x_ch_nodes):
    y = runge(x)
    y_cheb = runge(x_ch_nodes)
    print(len(x_ch_nodes))
    x_default = np.linspace(-1, 1, 100)
    axs[0].plot(x, y, 'ro', markersize=5)
    axs[0].plot(x_default, runge(x_default), label=r'$f (x)$')
    axs[0].set_ylabel('Y')
    axs[0].set_xlabel('X')
    axs[0].grid(True)
    axs[0].plot(x_default, [lagrange_interpolant(i, x, y) for i in x_default], 'm-', linewidth=2, label = r'$L (x)$')
    axs[1].plot(x_default, [lagrange_interpolant(i, x_ch_nodes, y_cheb) for i in x_default], 'm-', linewidth = 2, label = r'$L (x)$')
    axs[1].plot(x_ch_nodes, y_cheb, 'ro', markersize=5)
    axs[1].plot(x_default, runge(x_default), label=r'$f (x)$')
    axs[1].grid(True)
    axs[1].set_ylabel('Y')
    axs[1].set_xlabel('X')
    axs[0].legend()
    axs[1].legend()

if __name__ == '__main__':
    n = [5, 8, 11, 14, 17, 20, 23]
    fig, axs = plt.subplots(1, 2, figsize=(12, 3))
    axs[0].set_title('Равномерно распределенные узлы')
    axs[1].set_title('Чебышевские узлы')
    for i in range(len(n)):
        x_ch_nodes = np.array(chebish(n[i]))
        lag_x_nodes = np.linspace(-1, 1, n[i])
        plotting(lag_x_nodes, n[i], i, axs, x_ch_nodes)
        plt.show()
        fig, axs = plt.subplots(1, 2, figsize=(12, 3))
        axs[0].set_title('Равномерно распределенные узлы')




