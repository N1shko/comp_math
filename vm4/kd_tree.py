# import numpy as np
# import matplotlib.pyplot as plt

# #k1 = [10, 50]
# #k2 = [10, 50]

# k1 = np.arange(10, 50, 1)
# k2 = np.arange(10, 50, 1)
# d = 0
# i = 1
# print(k1)
# import numpy as np
# from sklearn.neighbors import KDTree
# rng = np.random.RandomState(0)
# X = rng.random_sample((10, 3))  # 10 points in 3 dimensions
# print(X)
# tree = KDTree(X, leaf_size=2)
# print(tree.get_arrays())     
# print(tree.query_radius(X[:1], r=0.3, count_only=True))
# ind = tree.query_radius(X[:1], r=0.3)  
# print(ind)  # indices of neighbors within distance 0.3

import numpy as np
from scipy.spatial import KDTree

data = np.array([[40, 20], 
                [30, 60], 
                [50, 20], 
                [50, 70],
                [20, 80],
                [10, 70],
                [70, 60],
                [30, 30],
                [20, 20],
                [70, 30]])

a = KDTree(data)
