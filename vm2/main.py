import numpy as np
import matplotlib.pyplot as plt


def diff2(x_0, h, f):
    return (f(x_0 + h) - f(x_0 - h)) / (2 * h)


def diff4(x_0, h, f):
    return (8 * f(x_0 + h) - 8 * f(x_0 - h) + f(x_0 - 2 * h) - f(x_0 + 2 * h)) / (12 * h)


def compos_avg(a, b, n, f):
    h = (b - a) / (n - 1)
    x_nodes = np.linspace(a, b, n)
    result_sum = 0
    result_sum += f(x_nodes[0])
    result_sum += f(x_nodes[n - 1])
    for i in range( n ):
            result_sum +=  f(x_nodes[i])
    result_sum *= h
    return result_sum


def compos_trapec(a, b, n, f):
    h = (b - a) / (n - 1)
    x_nodes = np.linspace(a, b, n)
    result_sum = 0
    result_sum += f(x_nodes[0])
    result_sum += f(x_nodes[n - 1])
    for i in range(2, n - 1):
            result_sum += 2 * f(x_nodes[i])
    result_sum *= (h / 2)
    return result_sum


def composite_simpson(a, b, n, f):
    h = (b - a)/(n - 1)
    x_nodes = np.linspace(a, b, n)
    result_sum = 0
    result_sum += f(x_nodes[0])
    result_sum += f(x_nodes[n - 1])
    for i in range(1, n - 1):
        if i % 2 == 0:
            result_sum += 2 * f(x_nodes[i])
        else:
            result_sum += 4 * f(x_nodes[i])
    result_sum *= (h / 3)
    return result_sum


def gauss_quad5(f):
    left_lim = 0
    right_lim = 2
    det = (right_lim - left_lim) / 2
    det1 = (right_lim + left_lim) / 2
    x_1 = -np.sqrt(3/5)
    x_2 = 0
    x_3 = np.sqrt(3/5)
    t_1 = det1 + x_1 * det
    t_2 = det1 + x_2 * det
    t_3 = det1 + x_3 * det
    c_1 = 5/9
    c_2 = 8/9
    c_3 = 5/9
    return det * (c_1 * f(t_1) + c_2 * f(t_2) + c_3 * f(t_3))


def gauss_quad(f):
    x_1 = -np.sqrt(3/5)
    x_2 = 0
    x_3 = np.sqrt(3/5)
    c_1 = 5/9
    c_2 = 8/9
    c_3 = 5/9
    return  (c_1 * f(x_1) + c_2 * f(x_2) + c_3 * f(x_3))


def real_integral(a, x, x1):
    integral1 = 0
    integral2 = 0
    for i in range(7):
        integral1 += (a[i] * x ** (i + 1) / (i + 1))
        integral2 += (a[i] * x1 ** (i + 1) / (i + 1))
    return integral2 - integral1


def base_int():
    g2 = lambda x: x * x * np.sin(3 * x)
    print (g2)
    a = 0
    b = np.pi
    n = []
    for i in range (3, 1000, 10):
        n.append(i)
    for i in range(999, 10000, 100):
        n.append(i)
    for i in range(10001, 30000, 1000):
        n.append(i)
    g2_int_real = (- b**2 * np.cos(3 * b) + a**2 * np.cos(3 * a)) / 3 + \
                (2 * b * np.sin(3 * b) - 2 * a * np.sin(3 *a )) / 9\
                +(2 * np.cos(3 * b) - 2 * np.cos(3 * a)) / 27
    print('Real integral')
    print(g2_int_real)
    err_int = []
    h = []
    x_for_plotting = np.logspace(-4, 0, 50)
    for el in n:
        h.append((b - a) / el)
        g2_int = composite_simpson(a, b, el, g2)
        err_int.append(abs(g2_int - g2_int_real))
    fig, ax = plt.subplots(1, 1, figsize=(12, 6))
    ax.grid()
    ax.set_xlabel(r'$h$', fontsize = 18)
    ax.set_ylabel(r'$\Delta$', fontsize = 18)
    ax.loglog(h, err_int, 'ro', markersize=6, label = '$\Delta$')
    ax.loglog(x_for_plotting, 10**(0) * x_for_plotting**4, label='$O(h^4)$')
    #ax.loglog(x_for_plotting, x_for_plotting**2, label='$O(h^2)$')
    #ax.loglog(x_for_plotting, 0.1 * x_for_plotting ** 3, label='$O(h^3)$')
    ax.legend()
    plt.show()


def base_diff():
    g1 = lambda x: x * np.exp(x)
    print(g1)
    x_0 = 2
    h = np.logspace(-16, 0, 200)
    g1_d_real = np.exp(x_0) + x_0 * np.exp(x_0)
    err_dif = []
    err_dif2 = []
    for el in h:
        g1_d = diff4(x_0, el, g1)
        err_dif.append(abs(g1_d - g1_d_real))
        g1_d2 = diff2(x_0, el, g1)
        err_dif2.append(abs(g1_d2 - g1_d_real))
    fig, ax = plt.subplots(1, 1, figsize=(12, 6))
    ax.grid()
    ax.set_xlabel('$h$', fontsize = 18)
    ax.set_ylabel('$\Delta$', fontsize = 18)
    ax.loglog(h, err_dif, 'o', label='$\Delta(g_4\'(x_0))$')
    ax.loglog(h, err_dif2, 'o', label='$\Delta(g_2\'(x_0))$')
    ax.loglog(h[-20:], [el ** 4 for el in h[-20:]], label='$O(h^4)$', color = 'red' )
    ax.loglog(h[-30:], [el ** 2 for el in h[-30:]], label='$O(h^2)$')
    ax.legend()
    plt.show()


def gauss_int():
    left_lim = 0
    right_lim = 2
    a = []
    for i in range(7):
        a.append([0] * 7)
        for j in range(i + 1):
            a[i][j] = np.random.randn()
        print(a[i])
    print(type(a[0][0]))
    poli_list = []
    real_int = []
    for i in range(7):
        poli_list.append(lambda x: a[i][0] * x ** 0 +
                                   a[i][1] * x ** 1 +
                                   a[i][2] * x ** 2 +
                                   a[i][3] * x ** 3 +
                                   a[i][4] * x ** 4 +
                                   a[i][5] * x ** 5 +
                                   a[i][6] * x ** 6)
        real_int.append(real_integral(a[i], 0, 2))
    err_dif = []
    for i in range(7):
        err_dif.append(abs(gauss_quad5(poli_list[i]) - real_int[i]))
        print(gauss_quad5(poli_list[i]), real_int[i])
    print(err_dif)
    print("Степень полинома      Аналитическое значение    Значение по квадратуре     Абсолютная погрешность")
    for i in range(7):
        print('%12d           %18.18f      %18.18f      %18.18f' % (i, real_int[i], gauss_quad5(poli_list[i]), err_dif[i]))


if (__name__ == '__main__'):
    base_diff()
    base_int()
    np.random.seed(101)
    gauss_int()
#