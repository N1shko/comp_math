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


def lagrange_basis(i_offset, cur_values, value0_bonds, p=2, m=2, debug=False):
    res=1
    for i in range(p+1):
        for j in range(m):
            if i_offset[j] != i:
                res *= (p*(cur_values[j]-value0_bonds[j][0])/(value0_bonds[j][1]-value0_bonds[j][0]) - i) / (i_offset[j] - i)
                if debug:
                    print("basis, 0: ",value0_bonds, p)
                    print("basis, cur_values:", cur_values)
                    print("basis: ", i, j, i_offset[j])
                    print("chisl:", (p*(cur_values[j]-value0_bonds[j][0])))
                    print("znam:", i_offset[j] - i)
                    print("basis: ", (p*(cur_values[j]-value0_bonds[j][0])/(value0_bonds[j][1]-value0_bonds[j][0]) - i) / (i_offset[j] - i), "\n")
            # else:
            #     print("i am here", i, j, i_offset[j], cur_values)
    return res


def lagrange_interpolant(y, cur_values, value0_bounds, p=2, m=2, debug=False):
    res = 0
    for i in range(p+1):
        for j in range(p+1):
            if debug:
                print(y[i*(p+1)+j] * lagrange_basis([i, j], cur_values, value0_bounds, p, m, debug=True))
                print(y, cur_values, y[i*(p+1)+j], [i, j])
            res += y[i*(p+1)+j] * lagrange_basis([j, i], cur_values, value0_bounds, p, m)
    if debug:
        print(res)
    return res