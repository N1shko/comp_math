import numpy as np
import matplotlib.pyplot as plt
import time

#функция вычиления базисного полинома Лагранжа
def l_i(i, x, x_nodes):
    if i >= len(x_nodes):
        exit(1)
    l_i = 1
    for j in range(len(x_nodes)):
        if j!=i:
             l_i *= (x - x_nodes[j])/(x_nodes[i] - x_nodes[j])
    return l_i


#функция вычисления полинома Лагранжа
def L(x, x_nodes, y_nodes):
    L = 0
    for i in range(len(x_nodes)):
        L += y_nodes[i] * l_i(i, x, x_nodes)
    return L


#функция вычисления чебышевских узлов
def cheb(n) -> list:
    return(np.array([np.cos((2*i - 1) / (2 * n) * np.pi) for i in range(1, n + 1)]))


#класс параметров функции Паде
class Params:
    def __init__(self, n, m, a_j, b_k):
        self.n = n
        self.m = m
        self.a_j = a_j
        self.b_k = b_k


#генерируем параметры для функции Паде
def generate_params():
    n = np.random.randint(7, 15)
    m = np.random.randint(7, 15)
    a_j = np.random.sample(m, )
    b_k = np.random.sample(n, )
    return Params(n, m, a_j, b_k)


#функция вычисления значения функции Паде
def pade(x, *params_list):
    sum_a_j = 0
    sum_b_k = 0
    for obj in params_list:
        for j in range(obj.m):
            sum_a_j += obj.a_j[j] * x ** j
        for k in range(1, obj.n):
            sum_b_k += obj.b_k[k] * x ** k
    return sum_a_j / (1 + sum_b_k)


#функция вычисления и вывода функции Паде и интеполяционных полиномов
# def plotting(x, l, x_cheb, kter):
#     x_def = np.linspace(-1, 1, 200)
#     y_def = pade(x_def, params_list[kter])
#     y_lin = pade(x, params_list[kter])
#     y_lin_interpolar = [L(i, x, y_lin) for i in x_def]
#     y_cheb = pade(x_cheb, params_list[kter])
#     y_cheb_interpolar = [L(i, x_cheb, y_cheb) for i in x_def]
#     ax[l, 0].plot(x, y_lin, 'ro', markersize = 1) #nodes
#     ax[l, 0].plot(x_def, y_def, 'g-', linewidth=1, label='Pade')
#     ax[l, 0].plot(x_def, y_lin_interpolar, 'b-', linewidth=0.5, label='line ')
#     ax[l, 0].plot(x_cheb, y_cheb, 'go', linewidth=1, label='cheb ')
#     ax[l, 0].plot(x_def, y_cheb_interpolar, '-', color = "#aaa",  markersize=1, label='cheb nodes')
#     ax[l, 0].grid()
#     ax[l, 0].legend()


#функция вычисления расстояния и построения графиков
# def distance(numb, fun_i):
#     var = []
#     max_dist = [[0.] * 30] * 2 #В пространстве 𝐿∞ нормой является равномерная норма,
#                                 # которую можно определить как ||𝑔(𝑥)||∞ = max𝑥∈[𝑎;𝑏]
#     for i in range(1, 31):
#         var.append(i)
#     for i in range(1, 31):
#         x_def = np.linspace(-1, 1, 200)
#         x_ch = np.array(cheb(i))  # генерируем список чебышевских узлов
#         x_eq = np.linspace(-1, 1, i)  # генерируем список равномерно расположеных узлов
#         y_def = pade(x_def, params_list[numb])
#         y_cheb = pade(x_ch_nodes, params_list[numb])
#         y_eq = pade(x_eq, params_list[numb])
#         y_eq_interpolar = [L(it, x_eq_nodes, y_eq) for it in x_def]
#         y_cheb_interpolar = [L(it, x_ch_nodes, y_cheb) for it in x_def]
#         max_dist[0][i - 1] = max(np.abs(y_def[iter] - y_eq_interpolar[iter]) for iter in range(200))
#         max_dist[1][i - 1] = max(np.abs(y_def[iter] - y_cheb_interpolar[iter]) for iter in range(200))
#     ax[fun_i, 1].semilogy(var, max_dist[0], color = "aaa")
#     ax[fun_i, 1].grid()
#     ax[fun_i, 2].semilogy(var, max_dist[1], 'r-')
#     ax[fun_i, 2].grid()


if __name__ == '__main__':
    start_time = time.time()
    #np.random.seed(1000)
    params_list = [] #список параметров для 100 функций Паде
    N = 7 #должно быть не меньше 5
    var = [] #массив значений N = {1, ... , 30}
    max_dist = [0.] * 30   # В пространстве 𝐿∞ нормой является равномерная норма,
    max_dist1 = [0.] * 30                      # которую можно определить как ||𝑔(𝑥)||∞ = max𝑥∈[𝑎;𝑏]
    for i in range(1, 31):
        var.append(i)

    for i in range(100): #В цикле генерируем параметры для 100 функций
        params_list.append(generate_params())
    x_ch_nodes = np.array(cheb(N)) #генерируем список чебышевских узлов
    x_eq_nodes = np.linspace(-1, 1, N) #генерируем список равномерно расположеных узлов
    fig, ax = plt.subplots(4, 3, figsize=(24, 24))
    for kter in range(100): #проходимся по 100 функциям
        if kter % 25 == 0: #будем выводить каждую 25-ую функцию
            l = kter // 25
            #plotting(x_eq_nodes, l, x_ch_nodes, kter)
            x_def = np.linspace(-1, 1, 200)
            y_def = [pade(x, params_list[kter]) for x in x_def]
            y_eq = [pade(x, params_list[kter]) for x in x_eq_nodes]
            y_eq_interp = [L(i, x_eq_nodes, y_eq) for i in x_def]
            y_cheb = [pade(x, params_list[kter]) for x in x_ch_nodes]
            y_cheb_interp = [L(i, x_ch_nodes, y_cheb) for i in x_def]
            ax[l][0].plot(x_eq_nodes, y_eq, 'o', color = "#9400D3", markersize = 3)  # nodes
            ax[l][0].plot(x_def, y_def, '-', color = "#40E0D0", label = 'Pade approximation')
            ax[l][0].plot(x_def, y_eq_interp, '-', color = "#BA55D3", label = 'Interpolar eq nodes (purple nodes)')
            ax[l][0].plot(x_def, y_cheb_interp, '-', color = "#4169E1", label = 'Interpolar cheb nodes (grey nodes)')
            ax[l][0].plot(x_ch_nodes, y_cheb, 'o', color = "#aaa", markersize = 3)
            ax[l][0].grid()
            ax[l][0].legend()
            #plt.show()
            for i in range(1, 31):
                x_def = np.linspace(-1, 1, 200)
                x_ch = np.array(cheb(i))  # генерируем список чебышевских узлов
                x_eq = np.linspace(-1, 1, i)  # генерируем список равномерно расположеных узлов
                y_def = [pade(x, params_list[kter]) for x in x_def]
                y_cheb = [pade(x, params_list[kter]) for x in x_ch]
                y_eq = [pade(x, params_list[kter]) for x in x_eq]
                y_eq_interpolar = [L(it, x_eq, y_eq) for it in x_def]
                y_cheb_interpolar = [L(it, x_ch, y_cheb) for it in x_def]
                max_dist[i - 1] = max([np.abs(y_def[iter] - y_eq_interpolar[iter]) for iter in range(200)])
                max_dist1[i - 1] = max([np.abs(y_def[iter] - y_cheb_interpolar[iter]) for iter in range(200)])
            ax[l, 1].semilogy(var, max_dist, '-', color = "#FF00FF", label = "The distance on\n eq nodes")
            ax[l, 1].grid()
            ax[l, 1].legend()
            ax[l, 2].semilogy(var, max_dist1, '-', color = "#1E90FF", label = "The distance on\n cheb nodes")
            ax[l, 2].grid()
            ax[l, 2].legend()
    plt.show()
    print("--- %s seconds ---" % (time.time() - start_time))