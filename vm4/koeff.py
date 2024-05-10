import numpy as np
import matplotlib.pyplot as plt
from self_made_kd_tree import KDTree
from calc_func import CalcTools
plt.rcParams.update({'font.size': 20})
EPS = 1e-3

k = 1 
b = 0.3
c = 0.7

fl = 0
random_dots_amount = 10
def func(x, koeffs):
    #return np.array([-k * x[0], 2*k*x[0]])
    return np.array([-koeffs[0] * x[0] + b * x[0] * x[1], c * x[1] - koeffs[1] * x[0] * x[1]], dtype = np.float32)


if __name__ == '__main__':
    tool = CalcTools()
    p = 4
    m = 2
    f = lambda t, x, koeffs: func(x, koeffs)
    x_0 = 0.8
    y_0 = 0.6
    
    t_n = 15.5
    h = 0.5
    t_h = h
    t = np.arange(0, t_n, t_h, dtype = np.float32)
    n = t.size


    tree = KDTree([[0.9, 1.3], [0.6, 1.0]], m, p, f, h, 2)

    separation_candidates = []
    for k in range(len(t)):
        print("3")
        global_flag = False
        for leaf in tree.leaf_list:
            tree.nodes[leaf].calculation_needed = True
        while not global_flag:
            local_flag = True
            if t[k] == t[-1] or fl:
                fig, ax = plt.subplots(1, 1, figsize=(10, 10))
                ax.grid()
                ax.set_xlabel(f"t={t[k]}")
            # ax.set_xlim(0,3)
            # ax.set_ylim(0,3)
            for i in range(len(separation_candidates)):
                tree.separate_node(separation_candidates[i])
            separation_candidates = []
            print("tree:", tree.leaf_list, t[k])
            for leaf in tree.leaf_list:
                gor_grid_x = np.zeros(len(tree.nodes[leaf].dots), dtype=np.float32)
                gor_grid_y = np.zeros(len(tree.nodes[leaf].dots), dtype=np.float32)
                if tree.nodes[leaf].is_new or tree.nodes[leaf].just_made:
                    #print("overlooking: ", tree.nodes[leaf].dots)
                    for i in range(len(tree.nodes[leaf].dots)):
                        dot = tree.nodes[leaf].dots[i]
                        #print(dot)
                        res = tool.rk4([x_0, y_0], t[k], f, t_h, koeffs=dot)
                        tree.nodes[leaf].dots_val[i] = res.copy()
                        gor_grid_x[i] = res[0]
                        gor_grid_y[i] = res[1]

                    # if t[k] == t[-1] or fl:
                    #     ax.scatter(gor_grid_x, gor_grid_y, c="yellow")
                    tree.nodes[leaf].is_new = False
                    tree.nodes[leaf].swap_needed = False
                else:
                    tree.nodes[leaf].swap_needed = True
                    for ik in range(len(tree.nodes[leaf].dots_val)):
                        dot = tree.nodes[leaf].dots[ik]
                        res = tool.rk4_iteration([tree.nodes[leaf].dots_val[ik][0], tree.nodes[leaf].dots_val[ik][1]], t[k], f, t_h, koeffs=dot)
                        gor_grid_x[ik] = res[0]
                        gor_grid_y[ik] = res[1]
                        tree.nodes[leaf].temp_val[ik] = res.copy()

                    # if t[k] == t[-1] or fl:
                    #     #print("wtf", gor_grid_x[-1], gor_grid_y[-1])
                    #     ax.scatter(gor_grid_x, gor_grid_y, c="black")                
                x_for_plotting = np.zeros(len(tree.nodes[leaf].plot_dots), dtype=np.float32)
                y_for_plotting = np.zeros(len(tree.nodes[leaf].plot_dots), dtype=np.float32)
                if t[k] == t[-1] or fl:
                    for i in range(len(tree.nodes[leaf].plot_dots)):
                        dot = tree.nodes[leaf].plot_dots[i]
                        x_for_plotting[i] = tool.lagrange_interpolant_2d(gor_grid_x, dot, tree.nodes[leaf].borders, p, m)
                        y_for_plotting[i] = tool.lagrange_interpolant_2d(gor_grid_y, dot, tree.nodes[leaf].borders, p, m)
                    # ax.scatter(x_for_plotting, y_for_plotting, c="green", s=1)
                    num = tree.nodes[leaf].plot_dots_num
                    h_plot = num // p
                    for i in range(0, num+1, h_plot):
                        ax.plot(x_for_plotting[i*(num+1):(i+1)*(num+1)], y_for_plotting[i*(num+1):(i+1)*(num+1)])
                        ax.plot(x_for_plotting[i::num+1], y_for_plotting[i::num+1])

                x_test_real = np.zeros(len(tree.nodes[leaf].random_dots), dtype=np.float32)
                y_test_real = np.zeros(len(tree.nodes[leaf].random_dots), dtype=np.float32)

                x_test_interpolation = np.zeros(len(tree.nodes[leaf].random_dots), dtype=np.float32)
                y_test_interpolation = np.zeros(len(tree.nodes[leaf].random_dots), dtype=np.float32)
                #tree.nodes[leaf].generate_random_dots(random_dots_amount, t[k])
                for i in range(len(tree.nodes[leaf].random_dots)):
                    dot = tree.nodes[leaf].random_dots[i]
                    res = tool.rk4([x_0, y_0], t[k], f, t_h, koeffs=dot)
                    x_test_real[i] = res[0]
                    y_test_real[i] = res[1]
                    x_test_interpolation[i] = tool.lagrange_interpolant_2d(gor_grid_x, dot, tree.nodes[leaf].borders, p, m)
                    y_test_interpolation[i] = tool.lagrange_interpolant_2d(gor_grid_y, dot, tree.nodes[leaf].borders, p, m)
                # if t[k] == t[-1] or fl:
                #     ax.text(gor_grid_x[0], gor_grid_y[0], f"{leaf}", fontsize=12)

                #     ax.scatter(x_test_real, y_test_real, c="blue")

                #     #if t[k] == t[-1]:
                #     ax.scatter(x_test_interpolation, y_test_interpolation, c="red")
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
                    #tree.separate_node(leaf)
                    separation_candidates.append(leaf)
                print(avg_error_y, avg_error_x, t[k], local_flag, global_flag)
                if leaf == 9:
                    print(tree.nodes[9].dots)   
                    print(tree.nodes[9].dots_val)
                    print(tree.nodes[9].temp_val)
            if local_flag == True:
                global_flag = True
            if t[k] == t[-1] or fl:
                plt.show()
        for leaf in tree.leaf_list:
            tree.nodes[leaf].just_made = False
            if tree.nodes[leaf].swap_needed:
                tree.nodes[leaf].dots_val = tree.nodes[leaf].temp_val.copy()
    tree.plot_grid()