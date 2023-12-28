import numpy as np



def conjugate_gradients(A, b):
    x_k = np.array([0, 0])
    d_k = np.array([0, 0])
    g_k = b
    g_k_next = np.matmul(A, x_k.T)- b
    d_k_next = -g_k_next + (np.dot(g_k_next.T, g_k_next) / np.dot(g_k.T, g_k)) * d_k
    s_k_next = np.dot(d_k_next, g_k_next) / np.matmul(np.matmul(d_k_next.T, A), d_k_next)

    x_k_next = np.array([round(i, 2) for i in x_k - s_k_next * d_k_next])

    while np.linalg.norm(x_k_next - x_k, ord=np.inf) >= 10 ** (-8):
        x_k = x_k_next
        d_k = d_k_next
        g_k = g_k_next
        g_k_next = np.matmul(A, x_k.T) - b
        d_k_next = -g_k_next + (np.dot(g_k_next.T, g_k_next) / np.dot(g_k.T, g_k)) * d_k
        s_k_next = np.dot(d_k_next, g_k_next) / np.matmul(np.matmul(d_k_next.T, A), d_k_next)
        x_k_next = np.array([round(i, 2) for i in x_k - s_k_next * d_k_next])
    return x_k_next


if __name__ == "__main__":
    A = np.array([[1.0, 0.5], [0.5, 0.33]])
    b = np.array([0.24, 0.13])
    print(conjugate_gradients(A, b))