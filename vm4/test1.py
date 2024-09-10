import numpy as np
borders = [[0.0, 1.0], [1.0, 2.0], [2.0, 3.0]]
p = 3
points = [np.linspace(borders[i][0], borders[i][1], p+1) for i in range(len(borders))]
grid = np.meshgrid(*points)
dots = np.vstack([grid[i].ravel() for i in range(len(borders))]).T

print(dots)


random_dots = np.zeros((p, len(borders)))
for i in range(len(borders)):
    random_dots[:, i] = np.random.uniform(borders[i][0], borders[i][1], p)

print(random_dots)

i = 4
m = 3
index = []
for _ in range(m):
    i, remainder = divmod(i, p+1)
    index.append(remainder)
print(index[::-1])