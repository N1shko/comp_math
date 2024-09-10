import numpy as np
import matplotlib.pyplot as plt
import intvalpy as ip
from self_made_kd_tree import KDTree
from calc_func import CalcTools
plt.rcParams.update({'font.size': 20})
EPS = 1e-3

k = 1 
b = 0.3
c = 0.7

fl = 1
random_dots_amount = 10
def func(x, koeffs):
    #return np.array([-k * x[0], 2*k*x[0]])
    #return np.array([-koeffs[0] * x[0] + b * x[0] * x[1], c * x[1] - koeffs[1] * x[0] * x[1]], dtype = np.float32)
    return np.array([koeffs[0] * x[1], 
                     -koeffs[0] * x[1], 
                     0.5 * koeffs[0] * x[1]], 
                     dtype = np.float32)



if __name__ == '__main__':
    tool = CalcTools()
    p = 4
    m = 1
    f = lambda t, x, koeffs: func(x, koeffs)
    x_0 = [0.0, 150.0, 0.0]
    
    t_n = 10.5
    h = 0.5
    t_h = h
    t = np.arange(0, t_n, t_h, dtype = np.float32)
    n = t.size

    min_val = np.zeros((n, 3))
    max_val = np.zeros((n, 3))
    tree = KDTree([[0.0, 5.0]], m, p, f, h, len(x_0))

    separation_candidates = []
    for k in range(len(t)):
        is_set = False
        global_flag = False
        for leaf in tree.leaf_list:
            tree.nodes[leaf].calculation_needed = True
        while not global_flag:
            local_flag = True
            # if t[k] == t[-1] or fl:
            #     fig, ax = plt.subplots(1, 1, figsize=(10, 10))
            #     ax.grid()
            #     ax.set_xlabel(f"t={t[k]}")
            # ax.set_xlim(0,3)
            # ax.set_ylim(0,3)
            for i in range(len(separation_candidates)):
                tree.separate_node(separation_candidates[i])
            separation_candidates = []
            for leaf in tree.leaf_list:
                # print(tree.nodes[leaf].dots)
                gor_grid_x = np.zeros(len(tree.nodes[leaf].dots), dtype=np.float32)
                gor_grid_y = np.zeros(len(tree.nodes[leaf].dots), dtype=np.float32)
                gor_grid_z = np.zeros(len(tree.nodes[leaf].dots), dtype=np.float32)
                if tree.nodes[leaf].is_new or tree.nodes[leaf].just_made:
                    #print("overlooking: ", tree.nodes[leaf].dots)
                    for i in range(len(tree.nodes[leaf].dots)):
                        dot = tree.nodes[leaf].dots[i]
                        #print(dot)
                        res = tool.rk4(x_0, t[k], f, t_h, koeffs=dot)
                        # print("here", res, dot)
                        tree.nodes[leaf].dots_val[i] = res.copy()
                        gor_grid_x[i] = res[0]
                        gor_grid_y[i] = res[1]
                        gor_grid_z[i] = res[2]
                        if not is_set:
                            min_val[k][0] = res[0]
                            min_val[k][1] = res[1]
                            min_val[k][2] = res[2]
                            max_val[k][0] = res[0]
                            max_val[k][1] = res[1]
                            max_val[k][2] = res[2]
                            is_set = True
                        else:
                            min_val[k][0] = min(min_val[k][0], res[0])
                            min_val[k][1] = min(min_val[k][1], res[1])
                            min_val[k][2] = min(min_val[k][2], res[2])
                            max_val[k][0] = max(max_val[k][0], res[0])
                            max_val[k][1] = max(max_val[k][1], res[1]) 
                            max_val[k][2] = max(max_val[k][2], res[2]) 
                    # if t[k] == t[-1] or fl:
                    #     ax.scatter(gor_grid_x, gor_grid_y, c="yellow")
                    tree.nodes[leaf].is_new = False
                    tree.nodes[leaf].swap_needed = False
                else:
                    tree.nodes[leaf].swap_needed = True
                    for ik in range(len(tree.nodes[leaf].dots_val)):
                        dot = tree.nodes[leaf].dots[ik]
                        res = tool.rk4_iteration(tree.nodes[leaf].dots_val[ik], t[k], f, t_h, koeffs=dot)
                        # print(res)
                        gor_grid_x[ik] = res[0]
                        gor_grid_y[ik] = res[1]
                        gor_grid_z[ik] = res[2]
                        tree.nodes[leaf].temp_val[ik] = res.copy()
                        if not is_set:
                            min_val[k][0] = res[0]
                            min_val[k][1] = res[1]
                            min_val[k][2] = res[2]
                            max_val[k][0] = res[0]
                            max_val[k][1] = res[1]
                            max_val[k][2] = res[2]
                            is_set = True
                        else:
                            min_val[k][0] = min(min_val[k][0], res[0])
                            min_val[k][1] = min(min_val[k][1], res[1])
                            min_val[k][2] = min(min_val[k][2], res[2])
                            max_val[k][0] = max(max_val[k][0], res[0])
                            max_val[k][1] = max(max_val[k][1], res[1]) 
                            max_val[k][2] = max(max_val[k][2], res[2]) 

                    # if t[k] == t[-1] or fl:
                    #     #print("wtf", gor_grid_x[-1], gor_grid_y[-1])
                    #     ax.scatter(gor_grid_x, gor_grid_y, c="black")                
                x_for_plotting = np.zeros(len(tree.nodes[leaf].plot_dots), dtype=np.float32)
                y_for_plotting = np.zeros(len(tree.nodes[leaf].plot_dots), dtype=np.float32)
                # if t[k] == t[-1] or fl:
                #     for i in range(len(tree.nodes[leaf].plot_dots)):
                #         dot = tree.nodes[leaf].plot_dots[i]
                #         x_for_plotting[i] = tool.lagrange_interpolant_1d(gor_grid_x, dot, tree.nodes[leaf].borders, p, m)
                #         y_for_plotting[i] = tool.lagrange_interpolant_1d(gor_grid_y, dot, tree.nodes[leaf].borders, p, m)
                #     # ax.scatter(x_for_plotting, y_for_plotting, c="green", s=1)
                #     num = tree.nodes[leaf].plot_dots_num
                #     h_plot = num // p
                #     for i in range(0, num+1, h_plot):
                #         ax.plot(x_for_plotting[i*(num+1):(i+1)*(num+1)], y_for_plotting[i*(num+1):(i+1)*(num+1)])
                #         ax.plot(x_for_plotting[i::num+1], y_for_plotting[i::num+1])

                x_test_real = np.zeros(len(tree.nodes[leaf].random_dots), dtype=np.float32)
                y_test_real = np.zeros(len(tree.nodes[leaf].random_dots), dtype=np.float32)
                z_test_real = np.zeros(len(tree.nodes[leaf].random_dots), dtype=np.float32)

                x_test_interpolation = np.zeros(len(tree.nodes[leaf].random_dots), dtype=np.float32)
                y_test_interpolation = np.zeros(len(tree.nodes[leaf].random_dots), dtype=np.float32)
                z_test_interpolation = np.zeros(len(tree.nodes[leaf].random_dots), dtype=np.float32)
                #tree.nodes[leaf].generate_random_dots(random_dots_amount, t[k])
                for i in range(len(tree.nodes[leaf].random_dots)):
                    dot = tree.nodes[leaf].random_dots[i]
                    res = tool.rk4(x_0, t[k], f, t_h, koeffs=dot)
                    x_test_real[i] = res[0]
                    y_test_real[i] = res[1]
                    z_test_real[i] = res[2]
                    
                    x_test_interpolation[i] = tool.lagrange_interpolant_nd(gor_grid_x, dot, tree.nodes[leaf].borders, p, m)
                    y_test_interpolation[i] = tool.lagrange_interpolant_nd(gor_grid_y, dot, tree.nodes[leaf].borders, p, m)
                    z_test_interpolation[i] = tool.lagrange_interpolant_nd(gor_grid_z, dot, tree.nodes[leaf].borders, p, m)
                # print(gor_grid_x)
                # if k>2: 
                #     exit()
                # if t[k] == t[-1] or fl:
                #     ax.text(gor_grid_x[0], gor_grid_y[0], f"{leaf}", fontsize=12)

                #     ax.scatter(x_test_real, y_test_real, c="blue")

                #     #if t[k] == t[-1]:
                #     ax.scatter(x_test_interpolation, y_test_interpolation, c="red")
                #     for i in range(len(x_test_real)):
                #         ax.plot([x_test_interpolation[i], x_test_real[i]] ,[y_test_interpolation[i], y_test_real[i]])
                error_x = 0
                error_y = 0
                error_z = 0
                for i in range(len(x_test_real)):
                    error_x += np.abs(x_test_real[i] - x_test_interpolation[i])
                    error_y += np.abs(y_test_real[i] - y_test_interpolation[i])
                    error_z += np.abs(z_test_real[i] - z_test_interpolation[i])
                avg_error_x = error_x / len(x_test_real)
                avg_error_y = error_y / len(y_test_real)
                avg_error_z = error_z / len(z_test_real)
                if ((avg_error_x + avg_error_y + avg_error_z)/3)>EPS:
                    local_flag = False
                    #tree.separate_node(leaf)
                    separation_candidates.append(leaf)
                # print(avg_error_y, avg_error_x, t[k], local_flag, global_flag)
                # if leaf == 9:
                #     print(tree.nodes[9].dots)   
                #     print(tree.nodes[9].dots_val)
                #     print(tree.nodes[9].temp_val)
            if local_flag == True:
                global_flag = True
            # if t[k] == t[-1] or fl:
            #     plt.show()
        for leaf in tree.leaf_list:
            tree.nodes[leaf].just_made = False
            if tree.nodes[leaf].swap_needed:
                tree.nodes[leaf].dots_val = tree.nodes[leaf].temp_val.copy()
        

    tree.plot_grid()


    fig, ax = plt.subplots(3, 1, figsize=(10, 10))
    ax[0].grid()
    ax[1].grid()
    ax[2].grid()
    ax[2].set_xlabel("Время, сек")
    ax[1].set_ylabel("Концентрация, мг/л")
    ax[0].plot(t, [x[0] for x in min_val], label=r"$H_2O_{min}$")
    ax[0].plot(t, [x[0] for x in max_val], label=r"$H_2O_{max}$")

    ax[1].plot(t, [x[1] for x in min_val], label=r"$H_2O_{2min}$")
    ax[1].plot(t, [x[1] for x in max_val], label=r"$H_2O_{2max}$")

    ax[2].plot(t, [x[2] for x in min_val], label=r"$O_{2min}$")
    ax[2].plot(t, [x[2] for x in max_val], label=r"$O_{2max}$")
    ax[0].legend()
    ax[1].legend()
    ax[2].legend()
    plt.show()