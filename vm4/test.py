# import matplotlib.pyplot as plt
# import numpy as np

# x = np.linspace(0,4)
# y = np.linspace(0,1)

# def f(x, y):
#     return y * np.sin(x) 

# X, Y = np.meshgrid(x,y)
# Z = np.zeros((50,50))

# for i in range(50):
#    for j in range(50):
#        Z[i,j] = f(X[i,j],Y[i,j])

# plt.pcolor(X, Y, Z)
# plt.show()

import pandas as pd
import numpy as np
a = [ [1,2], [2,9], [3,7] ]
na = np.array(a)
print(na[:,1])
# array([1, 2, 3])

import numpy as np
# interval = [(0, 1), (0, 1)]
interval = [(0, 1), (0, 1), (0, 1)]
# interval = [(0, 1), (0, 1), (0, 1), (0, 1)]
# interval = [(0, 1), (0, 1), (0, 1), (0, 1), (0, 1)]
p = 3
points = [np.linspace(interval[i][0], interval[i][1], p) for i in range(len(interval))]
grid = np.meshgrid(*points)
dots = np.vstack([grid[i].ravel() for i in range(len(interval))]).T

print(dots)
for i in range(len(dots)):
    print(dots[i], "\n")




    # for i in range(len(koef1)):
    #     temp = []
    #     gor_grid_x.append([])
    #     gor_grid_y.append([])
    #     for j in range(len(koef2)):
    #         temp1 = koef1[i]
    #         temp2 = koef2[i]
    #         t_temp_x = []
    #         t_temp_y = []
    #         for t_ in t:
    #             res1 = rk4_iteration([temp1, temp2], t_, f, h)
    #             temp1 = res1[0]
    #             temp2 = res1[1]
    #             t_temp_x.append(temp1)
    #             t_temp_y.append(temp2)
    #         temp_x.append(t_temp_x)
    #         temp_y.append(t_temp_y)
    #         # res1 = rk4([koef1[i], koef2[j]], t_n, f, h)
    #         #     #ax.plot(res1[0, :], res1[1, :], '-', color='steelblue')
    #         # temp_x.append(res1[0])
    #         # temp_y.append(res1[1])
    #         # ax[0].plot(time, res1[0, :], color='cornflowerblue')
    #         # #ax[0].plot(time, res1[1, :], color='indigo')
    #         # ax[1].scatter(res1[0][-1], res1[1][-1], color='red')
    #         # #ax[1].plot()
    #         # #print(res1[0][1])
    #     gor_grid_x[i].append(temp_x[-1])
    #     gor_grid_y[i].append(temp_y[-1])