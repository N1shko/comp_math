import numpy as np


def optimal_diff_step_2(eps, M):
    return (3.*eps /M)**(1./3)


def optimal_diff_step_4(eps, M):
    return (45.*eps /(4.*M))**(1./5)


x = 2
M_3 = (3 + x) * np.exp(x)
M_5 = (5 + x) * np.exp(x)
h_3_opt = optimal_diff_step_2(np.finfo(np.float64).eps, M_3)
h_5_opt = optimal_diff_step_4(np.finfo(np.float64).eps, M_5)
print(h_3_opt)
print(h_5_opt)
e_3 = (np.finfo(np.float64).eps) / h_3_opt + (h_3_opt ** 2) * M_3 / 6
e_5 = (3 * (np.finfo(np.float64).eps)) / (12 * h_5_opt) + (h_5_opt ** 4) * M_5 / 30
print(e_3)
print(e_5)