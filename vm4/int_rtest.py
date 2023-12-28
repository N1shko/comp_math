import intvalpy as ip
ip.precision.extendedPrecisionQ = False
import numpy as np


iplt = ip.IPlot(figsize=(15, 15))
fig, ax = iplt.subplots(nrows=2, ncols=2)


#########################################################################
A, b = ip.Shary(2)
shary_uni = ip.IntLinIncR2(A, b, show=False)
shary_tol = ip.IntLinIncR2(A, b, consistency='tol', show=False)

axindex = (0, 0)
ax[axindex].set_title('United and tolerable solution sets for the Shary interval system')
ax[axindex].title.set_size(15)
iplt.IntLinIncR2(shary_uni, color='gray', alpha=0.5, s=0, axindex=axindex)
iplt.IntLinIncR2(shary_tol, color='blue', alpha=0.3, s=10, axindex=axindex)

#########################################################################
A = ip.Interval([
    [[-1, 1], [-1, 1]],
    [[-1, -1], [-1, 1]]
])
b = ip.Interval([[1, 1], [-2, 2]])
unconstrained_set = ip.IntLinIncR2(A, b, show=False)

axindex = (0, 1)
ax[axindex].set_title('Unbounded set')
ax[axindex].title.set_size(15)
iplt.IntLinIncR2(unconstrained_set, color='darkolivegreen', alpha=0.3, s=10, axindex=axindex)

#########################################################################
A = -np.array([[-3, -1],
              [-2, -2],
              [-1, -3],
              [1, -3],
              [2, -2],
              [3, -1],
              [3, 1],
              [2, 2],
              [1, 3],
              [-1, 3],
              [-2, 2],
              [-3, 1]])
b = -np.array([18,16,18,18,16,18,18,16,18,18,16,18])
duodecagon = ip.lineqs(A, b, show=False)

axindex = (1, 0)
ax[axindex].set_title('Duodecagon')
ax[axindex].title.set_size(15)
iplt.lineqs(duodecagon, color='peru', alpha=0.3, s=10, axindex=axindex)

#########################################################################
x = ip.Interval([[1, 1.2], [1.9, 2.7], [1.7, 1.95], [3.5, 3.5],
                 [4.5, 5.5], [6, 6], [6.5, 7.5], [7, 7.8]])
y = ip.Interval([[4, 4.3], [4.5, 5.3], [4.6, 4.8], [5.1, 6],
                 [6, 6.5], [7, 7], [6.7, 7.4], [6.8, 8]])

axindex = (1, 1)
ax[axindex].set_title('Interval scatterplot')
ax[axindex].title.set_size(15)
iplt.scatter(x, y, color='gray', alpha=0.7, s=10, axindex=axindex)