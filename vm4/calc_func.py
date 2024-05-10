import numpy as np
from functools import lru_cache


class CalcTools:
    def __init__ (self):
        self.root=None
        
    def rk4(self, x_0, t_n, f, h, koeffs=[]):
        t = np.arange(0, t_n, h)
        n = t.size
        if len(koeffs):
            for i in range(n):
                k1 = h * f(t, x_0, koeffs)
                k2 = h * f(t + h / 2, x_0 + k1 / 2, koeffs)
                k3 = h * f(t + h / 2, x_0 + k2 / 2, koeffs)
                k4 = h * f(t + h, x_0 + k3, koeffs)
                koef = (k1 + 2 * k2 + 2 * k3 + k4) / 6
                x_0=x_0 + koef        
        else:
            for i in range(n):
                k1 = h * f(t, x_0)
                k2 = h * f(t + h / 2, x_0 + k1 / 2)
                k3 = h * f(t + h / 2, x_0 + k2 / 2)
                k4 = h * f(t + h, x_0 + k3)
                koef = (k1 + 2 * k2 + 2 * k3 + k4) / 6
                x_0=x_0 + koef
        return x_0

    def pure_rk4(self, x_0, t_n, f, h, koeffs=[]):
        t = np.arange(0, t_n, h)
        n = t.size - 1
        res=[]
        res.append(x_0)
        if len(koeffs):
            for i in range(n):
                k1 = h * f(t, x_0, koeffs)
                k2 = h * f(t + h / 2, x_0 + k1 / 2, koeffs)
                k3 = h * f(t + h / 2, x_0 + k2 / 2, koeffs)
                k4 = h * f(t + h, x_0 + k3, koeffs)
                koef = (k1 + 2 * k2 + 2 * k3 + k4) / 6
                x_0=x_0 + koef
                res.append(x_0)      
        else:
            for i in range(n):
                k1 = h * f(t, x_0)
                k2 = h * f(t + h / 2, x_0 + k1 / 2)
                k3 = h * f(t + h / 2, x_0 + k2 / 2)
                k4 = h * f(t + h, x_0 + k3)
                koef = (k1 + 2 * k2 + 2 * k3 + k4) / 6
                x_0=x_0 + koef
                res.append(x_0)      
        return res
    # @lru_cache
    def rk4_iteration(self, x_0, t, f, h, koeffs=[]):
        if len(koeffs):
            k1 = h * f(t, x_0, koeffs)
            k2 = h * f(t + h / 2, x_0 + k1 / 2, koeffs)
            k3 = h * f(t + h / 2, x_0 + k2 / 2, koeffs)
            k4 = h * f(t + h, x_0 + k3, koeffs)
            koef = (k1 + 2 * k2 + 2 * k3 + k4) / 6
        else:
            k1 = h * f(t, x_0)
            k2 = h * f(t + h / 2, x_0 + k1 / 2)
            k3 = h * f(t + h / 2, x_0 + k2 / 2)
            k4 = h * f(t + h, x_0 + k3)
            koef = (k1 + 2 * k2 + 2 * k3 + k4) / 6
        return x_0 + koef

    def lagrange_basis(self, i_offset, cur_values, value0_bonds, p=2, m=2, debug=False):
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
        return res


    def lagrange_interpolant_2d(self, y, cur_values, value0_bounds, p=2, m=2, debug=False):
        res = 0
        for i in range(p+1):
            for j in range(p+1):
                if debug:
                    print(y[i*(p+1)+j] * self.lagrange_basis([i, j], cur_values, value0_bounds, p, m, debug=True))
                    print(y, cur_values, y[i*(p+1)+j], [i, j])
                a = self.lagrange_basis([j, i], cur_values, value0_bounds, p, m)
                res += y[i*(p+1)+j] * a
        if debug:
            print(res)
        return res

    
    def lagrange_interpolant_1d(self, y, cur_values, value0_bounds, p=2, m=2, debug=False):
        res = 0
        for i in range(p+1):
                res += y[i] * self.lagrange_basis([i], cur_values, value0_bounds, p, m)
        return res