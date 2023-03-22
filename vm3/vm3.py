import numpy as np
import matplotlib.pyplot as plt
# import scipy.linalg

A1 = np.array([[1, 1, 0, 3],
               [2, 1, -1, 1],
               [3, -1, -1, 2],
               [-1, 2, 3, -1]], dtype=np.float16)

b1 = np.array([[4],
               [1],
               [-3],
               [4]], dtype=np.float16)

A2 = np.array([[3, 1, -3],
               [6, 2, 5],
               [1, 4, -3]], dtype=np.float64)

b2 = np.array([[-16],
               [12],
               [-39]], dtype=np.float64)


def lu_p(A: list, permute: bool):
    n = len(A)
    C = np.array(A.copy())
    P = np.array([np.array([0**(abs(i-j)) for j in range(n)]) for i in range(n)], dtype=np.float64)
    for i in range(n):
        max_abs = 0
        max_row = -1
        for j in range(i, n):
            if(abs(C[j][i]) > max_abs):
                max_abs = abs(C[j][i])
                max_row = j
        if(max_abs!=0):
            if(permute):
                P[[max_row, i]] = P[[i, max_row]]
                C[[max_row, i]] = C[[i, max_row]]
            for j in range(i+1, n):
               C[j][i] /= C[i][i]
               for k in range(i+1, n):
                   C[j][k] -= C[j][i] * C[i][k]
    U = np.triu(C)
    L = np.tril(C, -1)
    return P, L, U

def lu(A):
    n = len(A)
    C = np.array(A.copy())
    for i in range(n):
        for j in range(i+1, n):
            C[j][i] /= C[i][i]
            for k in range(i+1, n):
                C[j][k] -= C[j][i] * C[i][k]
                print(C)
    U = np.triu(C)
    L = np.tril(C, -1)
    return L, U

def solve(L, U, P, b):
    n = len(b)
    x = np.zeros(n, dtype=np.float64)
    y = np.zeros(n, dtype=np.float64)
    b = np.matmul(P, b)
    for i in range(n):
        y[i] = b[i] - sum([L[i][k] * y[k] for k in range(i)])
    for i in range(n - 1, -1, -1):
        x[i] = (y[i] - sum([U[i][k] * x[k] for k in range(i, n)])) / U[i][i]
    return x

def sol(L, U, b):
    n = len(b)
    X = np.zeros(n, dtype=np.float64)
    Y = np.zeros(n, dtype=np.float64)
    for i in range(n):
        Y[i] = b[i] - sum([L[i][k] * Y[k] for k in range(i)])
    for i in range(n - 1, -1, -1):
        X[i] = (Y[i] - sum([U[i][k] * X[k] for k in range(i, n)])) / U[i][i]
    return X



def norm(A):
    k = len(A)
    summary = 0
    for i in range(k):
        summary += A[i]**2
    return(summary**(1/2))


if __name__ == '__main__':
    P, L, U = lu_p(A2, True)
    X = solve(L, U, P, b2)
    # L1, U1 = lu(A1)
    print(X)
    # X1 = sol(L1, U1, b1)
    # print(X1)
    temp_A = A2[0][0]
    temp_b = b2[0][0]
    result = []
    result1 = []
    y_nodes = []
    y_nodes1 = []
    x_nodes = []
    # for p in range(13):
    #     A2[0][0] = temp_A + 10**((-1) * p)
    #     b2[0][0] = temp_b + 10**((-1) * p)
    #     P, L, U = lu_wiki(A2, True)
    #     # print(b2)
    #     # print(A2)
    #     result.append(solve(L, U, P, b2))
    #     # print(p, end=" & ")
    #     # print(result[p][0], " & ", result[p][1], " & ", result[p][2], "\\", "\\")
    #     print("%d & %.16f & %.16f & %.16f \\\\\\hline" % (p, result[p][0], result[p][1], result[p][2]))
    #     y_nodes.append(abs(norm(X) - norm(result[p])) / norm(X))
    #     x_nodes.append(p)
    p_ax = np.linspace(0, 12, 101)
    for i in range(101):
        A2[0][0] = temp_A + 10**(-p_ax[i])
        b2[0][0] = temp_b + 10**(-p_ax[i])
        P, L, U = lu_p(A2, False)
        result.append(solve(L, U, P, b2))
        P, L, U = lu_p(A2, True)
        result1.append(solve(L, U, P, b2))
        print(p_ax[i], " ", result[i])
        y_nodes.append(abs(norm(X) - norm(result[i])) / norm(X))
        y_nodes1.append(abs(norm(X) - norm(result1[i])) / norm(X))
        x_nodes.append(p_ax[i])
    for i in range(13):
        print(y_nodes[i])
    fig, ax = plt.subplots(1, 1, figsize=(12, 6))
    ax.grid()
    ax.set_xlabel('$p$', size = 18)
    ax.set_ylabel('$E$', size = 18)
    ax.semilogy(p_ax, y_nodes, 'ro', markersize=6, label="Без частичного выбора главного элемента")
    ax.semilogy(p_ax, y_nodes1, 'bo', markersize=6, label="С частичным выбором главного элемента")
    ax.legend(prop={'size': 15})

    plt.show()
