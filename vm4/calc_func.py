import numpy as np



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
            res += y[i*(p+1)+j] * lagrange_basis([i, j], cur_values, value0_bounds, p, m)
    return res