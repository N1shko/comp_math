import numpy as np
import matplotlib.pyplot as plt
from self_made_kd_tree import KDTree
from calc_func import CalcTools
plt.rcParams.update({'font.size': 20})
EPS = 1e-2

k = 1 
a = -0.9
b = 0.5
c = 0.7
d = -0.8




fl = 0
random_dots_amount = 15
def func(x):
    #return np.array([-k * x[0], 2*k*x[0]])
    return np.array([a * x[0] + b * x[0] * x[1], c * x[1] + d * x[0] * x[1]], dtype = np.float64)

def calc_interpolation():
    tree_separate_counter = 0
    total_rk4_full = 0
    total_rk4_part = 0
    p = 7
    m = 2
    f = lambda t, x: func(x)
    
    tool = CalcTools()
    t_n = 15.5
    h = 0.5
    t_h = h
    t = np.arange(0, t_n, t_h, dtype = np.float64)
    print(t)
    n = t.size
    print("s", len(t))
    print("s", t.size)
    
    min_val = np.zeros((n, m))
    max_val = np.zeros((n, m))
    tree = KDTree([[0.5, 2.5], [0.5, 2.5]], m, p, f, h, 2)
    separation_candidates = []
    
    for k in range(len(t)):
        is_set = False
        print("3")
        global_flag = False
        # for leaf in tree.leaf_list:
        #     tree.nodes[leaf].calculation_needed = True
        while not global_flag:
            local_flag = True
            if t[k] == t[-1] or fl:
                fig, ax = plt.subplots(1, 1, figsize=(10, 10))
                ax.grid()
                ax.set_xlabel(f"t = {t[k]}")
            for i in range(len(separation_candidates)):
                tree_separate_counter+=1
                tree.separate_node(separation_candidates[i])
            separation_candidates = []
            print("tree:", tree.leaf_list, t[k])
            for leaf in tree.leaf_list:
                gor_grid_x = np.zeros(len(tree.nodes[leaf].dots), dtype=np.float64)
                gor_grid_y = np.zeros(len(tree.nodes[leaf].dots), dtype=np.float64)
                if tree.nodes[leaf].is_new or tree.nodes[leaf].just_made:
                    for i in range(len(tree.nodes[leaf].dots)):
                        dot = tree.nodes[leaf].dots[i]
                        res = tool.rk4([dot[0], dot[1]], t[k], f, t_h)
                        tree.nodes[leaf].dots_val[i] = res.copy()
                        gor_grid_x[i] = res[0]
                        gor_grid_y[i] = res[1]
                        total_rk4_full += 1

                    print("first")
                    # if t[k] == t[-1] or fl:
                    #     ax.scatter(gor_grid_x, gor_grid_y, c="yellow")
                    tree.nodes[leaf].is_new = False
                    tree.nodes[leaf].swap_needed = False
                else:
                    tree.nodes[leaf].swap_needed = True
                    print("second")
                    # if tree.nodes[leaf].calculation_needed:
                    for ik in range(len(tree.nodes[leaf].dots_val)):
                        res = tool.rk4_iteration([tree.nodes[leaf].dots_val[ik][0], tree.nodes[leaf].dots_val[ik][1]], t[k], f, t_h)
                        gor_grid_x[ik] = res[0]
                        gor_grid_y[ik] = res[1]
                        tree.nodes[leaf].temp_val[ik] = res.copy()
                        total_rk4_full += 1
                        if not is_set:
                            min_val[k][0] = res[0]
                            min_val[k][1] = res[1]
                            max_val[k][0] = res[0]
                            max_val[k][1] = res[1]
                            is_set = True
                        else:
                            min_val[k][0] = min(min_val[k][0], res[0])
                            min_val[k][1] = min(min_val[k][1], res[1])
                            max_val[k][0] = max(max_val[k][0], res[0])
                            max_val[k][1] = max(max_val[k][1], res[1])  
                        #tree.nodes[leaf].calculation_needed = False
                    # if t[k] == t[-1] or fl and tree.nodes[leaf].cal:
                    #     ax.scatter(gor_grid_x, gor_grid_y, c="black")                
                x_for_plotting = np.zeros(len(tree.nodes[leaf].plot_dots), dtype=np.float64)
                y_for_plotting = np.zeros(len(tree.nodes[leaf].plot_dots), dtype=np.float64)
                if t[k] == t[-1] or fl:
                    for i in range(len(tree.nodes[leaf].plot_dots)):
                        dot = tree.nodes[leaf].plot_dots[i]
                        x_for_plotting[i] = tool.lagrange_interpolant_2d(gor_grid_x, dot, tree.nodes[leaf].borders, p, m)
                        y_for_plotting[i] = tool.lagrange_interpolant_2d(gor_grid_y, dot, tree.nodes[leaf].borders, p, m)
                    #ax.scatter(x_for_plotting, y_for_plotting, c="green", s=1)
                    num = tree.nodes[leaf].plot_dots_num
                    h_plot = num // (p)
                    for i in range(0, num+1, h_plot):
                        ax.plot(x_for_plotting[i*(num+1):(i+1)*(num+1)], y_for_plotting[i*(num+1):(i+1)*(num+1)])
                        ax.plot(x_for_plotting[i::num+1], y_for_plotting[i::num+1])

                x_test_real = np.zeros(len(tree.nodes[leaf].random_dots), dtype=np.float64)
                y_test_real = np.zeros(len(tree.nodes[leaf].random_dots), dtype=np.float64)

                x_test_interpolation = np.zeros(len(tree.nodes[leaf].random_dots), dtype=np.float64)
                y_test_interpolation = np.zeros(len(tree.nodes[leaf].random_dots), dtype=np.float64)
                for i in range(len(tree.nodes[leaf].random_dots)):
                    dot = tree.nodes[leaf].random_dots[i]
                    res = tool.rk4([dot[0], dot[1]], t[k], f, t_h)
                    total_rk4_part += 1
                    x_test_real[i] = res[0]
                    y_test_real[i] = res[1]
                    x_test_interpolation[i] = tool.lagrange_interpolant_2d(gor_grid_x, dot, tree.nodes[leaf].borders, p, m)
                    y_test_interpolation[i] = tool.lagrange_interpolant_2d(gor_grid_y, dot, tree.nodes[leaf].borders, p, m)
                # if t[k] == t[-1] or fl:
                #     # ax.text(gor_grid_x[0], gor_grid_y[0], f"{leaf}", fontsize=12)
                #     # ax.scatter(x_test_real, y_test_real, c="blue")
                #     # ax.scatter(x_test_interpolation, y_test_interpolation, c="red")
                #     for i in range(len(x_test_real)):
                #         ax.plot([x_test_interpolation[i], x_test_real[i]] ,[y_test_interpolation[i], y_test_real[i]])
                error_x = 0
                error_y = 0
                for i in range(len(x_test_real)):
                    error_x += np.abs(x_test_real[i] - x_test_interpolation[i])
                    error_y += np.abs(y_test_real[i] - y_test_interpolation[i])
                avg_error_x = error_x / len(x_test_real)
                avg_error_y = error_y / len(y_test_real)
                if (avg_error_x + avg_error_y/2)>EPS:
                    local_flag = False
                    separation_candidates.append(leaf)
                print(avg_error_y, avg_error_x, t[k], local_flag, global_flag)
            if local_flag == True:
                print("i am here")
                global_flag = True
            if t[k] == t[-1] or fl:
                plt.show()
        for leaf in tree.leaf_list:
            tree.nodes[leaf].just_made = False
            if tree.nodes[leaf].swap_needed:
                tree.nodes[leaf].dots_val = tree.nodes[leaf].temp_val.copy()
    tree.plot_grid()
    return t, min_val, max_val

if __name__ == '__main__':
    t, min_val, max_val = calc_interpolation()
    fig, ax = plt.subplots(2, 1, figsize=(10, 10))
    ax[0].grid()
    ax[1].grid()
    ax[0].set_xlabel("t")
    ax[1].set_xlabel("t")
    ax[0].set_ylabel("x")
    ax[1].set_ylabel("y")
    ax[0].plot(t, [x[0] for x in min_val], label="x_min")
    ax[0].plot(t, [x[0] for x in max_val], label="x_max")
    ax[1].plot(t, [x[1] for x in min_val], label="y_min")
    ax[1].plot(t, [x[1] for x in max_val], label="y_max")



    ax[0].legend()
    ax[1].legend()
    plt.show()


    # print("tree_separate_counter: ", tree_separate_counter)
    # print("total_rk4_full: ", total_rk4_full)
    # print("total_rk4_part: ", total_rk4_part)