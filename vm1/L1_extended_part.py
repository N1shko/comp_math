import numpy as np
import matplotlib.pyplot as plt
import time
def pade(n, m, a, b, i, x):
    temp = 0
    temp1 = 1
    for j in range(m[i]):
        temp += a[i][j] * x**j
    for j in range(n[i]):
        temp1 += b[i][j] * x**j
    return temp/temp1


def chebish(n) -> list:
    return([np.cos((2 * i - 1) / (2 * n) * np.pi) for i in range(1, n+1)])


def l_i(i, x, x_nodes):
    iter = np.delete(x_nodes, i)
    return np.prod((x - iter) / (x_nodes[i] - iter))


def lagrange_interpolant(x, x_nodes, y_nodes):
    n = len(x_nodes)
    return np.sum(y_nodes * np.array([l_i(i, x, x_nodes) for i in range(n)]))


def plotting(x, n, m, a, b, k, x_cheb, dot):
    x_default = np.linspace(-1, 1, dot)
    y_default = pade(n, m, a, b, k, x_default)
    y_linear = pade(n, m, a, b, k, x)  # n, m, a, b, i, x
    y_linear_interpolar = [lagrange_interpolant(i, x,  y_linear) for i in x_default]
    y_cheb = pade(n, m, a, b, k,  x_cheb)
    y_cheb_interpolar = [lagrange_interpolant(i, x_cheb, y_cheb) for i in x_default]
    ax.plot(x_default,  y_default, '-', color='#aaa', linewidth=1, label=r'$f_{n, m}(x)$')  # о
    ax.plot(x_default, y_linear_interpolar, 'g-', linewidth=1, label="L(x) (равномерно распределенные узлы)")
    ax.grid()
    ax.set_xlabel('X', fontsize=20)
    ax.set_ylabel('Y', fontsize=20)
    ax.plot(x_default, y_cheb_interpolar, 'm-', linewidth=1, label="L(x) (Чебышевские узлы))")
    ax.legend(prop={'size': 15})


def distance(n, m, a, b, fun_number, dot):
    max_distance = [[], []]
    max_distance[0] = [0.] * 30
    max_distance[1] = [0.] * 30
    var_default = []
    for i in range(1, 31):
        var_default.append(i)
    for i in range(1, 31):
        x_ch_nodes = np.array(chebish(i))
        lag_x_nodes = np.linspace(-1, 1, i)
        x_default = np.linspace(-1, 1, dot)
        y_default = pade(n, m, a, b, fun_number, x_default)
        y_linear = pade(n, m, a, b, fun_number, lag_x_nodes)  # n, m, a, b, i, x
        y_linear_interpolar = [lagrange_interpolant(i, lag_x_nodes, y_linear) for i in x_default]
        y_cheb = pade(n, m, a, b, fun_number, x_ch_nodes)
        y_cheb_interpolar = [lagrange_interpolant(i, x_ch_nodes, y_cheb) for i in x_default]
        max_distance[0][i - 1] = max(np.absolute(y_default[iterator] - y_linear_interpolar[iterator]) for iterator in range(dot))
        max_distance[1][i - 1] = max(np.absolute(y_default[iterator] - y_cheb_interpolar[iterator]) for iterator in range(dot))
    ax[0].semilogy(var_default, max_distance[0], color='#aaa')
    ax[1].semilogy(var_default, max_distance[1], 'g-')
    ax[0].set_xlabel('N', fontsize=15)
    ax[1].set_ylabel(r'$||f_{n, m}(x) - L(x)||$', fontsize=13)
    ax[1].set_xlabel('N', fontsize=15)
    ax[0].set_ylabel(r'$||f_{n, m}(x) - L(x)||$', fontsize=13)
    ax[1].set_title('Расстояние по чебышевским узлам', fontsize=13)
    ax[0].set_title('Расстояние по равномерно распределенным узлам', fontsize=13)
    ax[0].grid()
    ax[1].grid()


if __name__ == '__main__':
    #np.random.seed(1000)
    #start_time = time.time()
    dot = 200
    set = [24, 45, 67, 88]
    n_node = 30
    n = np.random.randint(7, 16, 100)
    m = np.random.randint(7, 16, 100)
    a = [[]] * 100
    b = [[]] * 100
    for i in range(100):
        a[i] = np.random.rand(m[i])
        b[i] = np.random.rand(n[i])
    x_ch_nodes = np.array(chebish(n_node))
    lag_x_nodes = np.linspace(-1, 1, n_node)
    for k in range(len(set)):
        fig, ax = plt.subplots(1, 1, figsize=(12, 6))
        plotting(lag_x_nodes, n, m, a, b, k, x_ch_nodes, dot)
        plt.show()
        fig, ax = plt.subplots(1, 2, figsize=(12, 4))
        distance(n, m, a, b, set[k], dot)
        plt.show()
    #print("--- %s seconds ---" % (time.time() - start_time))




