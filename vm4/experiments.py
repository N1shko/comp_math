import numpy as np
import matplotlib.pyplot as plt
from self_made_kd_tree import KDTree
from calc_func import CalcTools
import random
plt.rcParams.update({'font.size': 20})
EPS = 1e-3

k = 1 
a = -0.9
b = 0.5
c = 0.7
d = -0.8


t_n = 5.5
h = 0.25
t_h = h
p = 6
m = 2

fl = 0
fl2 = 0
random_dots_amount = 15

random.seed(0)
max_dispersion = 0.1

intervals = [[1.5, 1.6], [1.5, 1.6]]

def func(x):
    #return np.array([-k * x[0], 2*k*x[0]])
    return np.array([a * x[0] + b * x[0] * x[1], c * x[1] + d * x[0] * x[1]], dtype = np.float64)

def calc_interpolation(inter):
    f = lambda t, x: func(x)
    
    tool = CalcTools()

    t = np.arange(0, t_n, t_h, dtype = np.float64)
    n = t.size
    min_val = np.zeros((n, m))
    max_val = np.zeros((n, m))
    tree = KDTree(inter, m, p, f, h, 2)
    separation_candidates = []
    
    for k in range(len(t)):
        is_set = False
        global_flag = False
        while not global_flag:
            local_flag = True
            if (t[k] == t[-1] or fl) and fl2:
                fig, ax = plt.subplots(1, 1, figsize=(10, 10))
                ax.grid()
                ax.set_xlabel(f"t = {t[k]}")
            for i in range(len(separation_candidates)):
                tree.separate_node(separation_candidates[i])
            separation_candidates = []
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
                    tree.nodes[leaf].is_new = False
                    tree.nodes[leaf].swap_needed = False
                else:
                    tree.nodes[leaf].swap_needed = True
                    for ik in range(len(tree.nodes[leaf].dots_val)):
                        res = tool.rk4_iteration([tree.nodes[leaf].dots_val[ik][0], tree.nodes[leaf].dots_val[ik][1]], t[k], f, t_h)
                        gor_grid_x[ik] = res[0]
                        gor_grid_y[ik] = res[1]
                        tree.nodes[leaf].temp_val[ik] = res.copy()
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
                x_for_plotting = np.zeros(len(tree.nodes[leaf].plot_dots), dtype=np.float64)
                y_for_plotting = np.zeros(len(tree.nodes[leaf].plot_dots), dtype=np.float64)
                if (t[k] == t[-1] or fl) and fl2:
                    for i in range(len(tree.nodes[leaf].plot_dots)):
                        dot = tree.nodes[leaf].plot_dots[i]
                        x_for_plotting[i] = tool.lagrange_interpolant_2d(gor_grid_x, dot, tree.nodes[leaf].borders, p, m)
                        y_for_plotting[i] = tool.lagrange_interpolant_2d(gor_grid_y, dot, tree.nodes[leaf].borders, p, m)
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
                    x_test_real[i] = res[0]
                    y_test_real[i] = res[1]
                    x_test_interpolation[i] = tool.lagrange_interpolant_2d(gor_grid_x, dot, tree.nodes[leaf].borders, p, m)
                    y_test_interpolation[i] = tool.lagrange_interpolant_2d(gor_grid_y, dot, tree.nodes[leaf].borders, p, m)
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
            if local_flag == True:
                global_flag = True
            if (t[k] == t[-1] or fl) and fl2:
                plt.show()
        for leaf in tree.leaf_list:
            tree.nodes[leaf].just_made = False
            if tree.nodes[leaf].swap_needed:
                tree.nodes[leaf].dots_val = tree.nodes[leaf].temp_val.copy()
    # tree.plot_grid()
    return t, min_val, max_val


if __name__ == '__main__':
    t, min_val, max_val = calc_interpolation(intervals)
    fixed_0 = (np.array([1.5, 1.5]))
    tool = CalcTools()
    f = lambda t, x: func(x)
    fixed_result = tool.pure_rk4(fixed_0, t_n, f, h)
    for i in range(len(fixed_result)-1):
        while True:
            check = fixed_result[i][0] + random.uniform(-max_dispersion, max_dispersion) 
            if check >= min_val[i][0] and check <= max_val[i][0]:
                fixed_result[i][0] = check
                break
        while True:
            check = fixed_result[i][1] + random.uniform(-max_dispersion, max_dispersion) 
            if check >= min_val[i][1] and check <= max_val[i][1]:
                fixed_result[i][1] = check
                break
    # print(fixed_result[0])
    print("fixed:", fixed_result[0])
    dot_count = 10
    start_sum = sum(x[1]-x[0] for x in intervals)
    for i in range(1):
        # fig, ax = plt.subplots(2, 1, figsize=(10, 10))
        # ax[0].plot(t, [x[0] for x in fixed_result[:-1]], color='cornflowerblue', label=r'$x(t)$')
        # ax[1].plot(t, [x[1] for x in fixed_result[:-1]], color='indigo', label=r'$y(t)$')
        # ax[0].grid()
        # ax[1].grid()
        # ax[0].set_xlabel("t")
        # ax[1].set_xlabel("t")
        # ax[0].set_ylabel("x")
        # ax[1].set_ylabel("y")
        # ax[0].plot(t, [x[0] for x in min_val], label="x_min")
        # ax[0].plot(t, [x[0] for x in max_val], label="x_max")
        # ax[1].plot(t, [x[1] for x in min_val], label="y_min")
        # ax[1].plot(t, [x[1] for x in max_val], label="y_max")
        new_indices = np.linspace(start=0, stop=len(t), num=i+1).astype(int)
        new_linspace = t[new_indices]
        print("new", new_linspace, new_indices)
        # for j in range(i):
        #     ax[0].scatter(new_linspace[j], fixed_result[new_indices[j]][0])
        #     ax[1].scatter(new_linspace[j], fixed_result[new_indices[j]][1])
        # ax[0].legend()
        # ax[1].legend()
        plt.show()
        local_intervals = intervals
        result_intervals = intervals
        temp_sum = sum(x[1]-x[0] for x in local_intervals)
        # print(temp_sum)
        # num_random = 100
        # for experiments in range(num_random):
        #     new_interval = []
        #     for interval in local_intervals:
        #         new_a = random.uniform(interval[0], interval[1])
        #         new_b = random.uniform(interval[0], interval[1])
        #         # print("New", new_a, new_b)
        #         random_subinterval = [min(new_a, new_b), max(new_a, new_b)]
        #         new_interval.append(random_subinterval)
        #     # print("New interval", new_interval)
        #     t, new_min_val, new_max_val = calc_interpolation(new_interval)
        #     # fig, ax = plt.subplots(2, 1, figsize=(10, 10))
        #     # ax[0].plot(t, [x[0] for x in fixed_result[:-1]], color='cornflowerblue', label=r'$x(t)$')
        #     # ax[1].plot(t, [x[1] for x in fixed_result[:-1]], color='indigo', label=r'$y(t)$')
        #     # ax[0].grid()
        #     # ax[1].grid()
        #     # ax[0].set_xlabel("t")
        #     # ax[1].set_xlabel("t")
        #     # ax[0].set_ylabel("x")
        #     # ax[1].set_ylabel("y")
        #     # ax[0].plot(t, [x[0] for x in new_min_val], label="x_min")
        #     # ax[0].plot(t, [x[0] for x in new_max_val], label="x_max")
        #     # ax[1].plot(t, [x[1] for x in new_min_val], label="y_min")
        #     # ax[1].plot(t, [x[1] for x in new_max_val], label="y_max")
        #     # for j in range(i):
        #     #     ax[0].scatter(new_linspace[j], fixed_result[new_indices[j]][0])
        #     #     ax[1].scatter(new_linspace[j], fixed_result[new_indices[j]][1])
        #     # ax[0].legend()
        #     # ax[1].legend()
        #     # plt.show()
        #     for checkpoint in range(len(new_linspace)):
        #         # print(new_min_val[new_indices[checkpoint]][0], new_max_val[new_indices[checkpoint]][0], new_min_val[new_indices[checkpoint]][1], new_max_val[new_indices[checkpoint]][1])
        #         if fixed_result[new_indices[checkpoint]][0] < new_min_val[new_indices[checkpoint]][0] or \
        #             fixed_result[new_indices[checkpoint]][0] > new_max_val[new_indices[checkpoint]][0] or \
        #             fixed_result[new_indices[checkpoint]][1] < new_min_val[new_indices[checkpoint]][1] or \
        #             fixed_result[new_indices[checkpoint]][1] > new_max_val[new_indices[checkpoint]][1]:
        #             # print("am i supposed to be here?", fixed_result[new_indices[checkpoint]][0], fixed_result[new_indices[checkpoint]][1])
        #             break
        #         new_temp_sum = sum(x[1]-x[0] for x in new_interval)
        #         # print("new sum", new_temp_sum)
        #         if new_temp_sum < temp_sum:
        #             temp_sum = new_temp_sum
        #             result_intervals = new_interval
        new_sum = sum(x[1]-x[0] for x in result_intervals)
        prev_sum = new_sum
        steps = [[0.01, 0.01], [0.01, 0.01]]

        num_determine = 100
        print("Первичное приближение: ", result_intervals)
        stop_eps = 1e-5
        stop_flag = 0
        for iterations in range(num_determine):
            prev_sum = temp_sum
            prev_sum_no_success = temp_sum
            for it in range(len(result_intervals)):
                for k in range(2):
                    new_interval = result_intervals
                    if k == 0:
                        while new_interval[it][k] + steps[it][k] > new_interval[it][1]:
                            steps[it][k] /= 2 
                        new_interval[it][k] += steps[it][k]
                    else:
                        while new_interval[it][k] - steps[it][k] < new_interval[it][0]:
                            steps[it][k] /= 2 
                        new_interval[it][k] -= steps[it][k]
                    t, new_min_val, new_max_val = calc_interpolation(result_intervals)
                    # fig, ax = plt.subplots(2, 1, figsize=(10, 10))
                    # ax[0].plot(t, [x[0] for x in fixed_result], color='cornflowerblue', label=r'$x(t)$')
                    # ax[1].plot(t, [x[1] for x in fixed_result], color='indigo', label=r'$y(t)$')
                    # ax[0].grid()
                    # ax[1].grid()
                    # ax[0].set_xlabel("t")
                    # ax[1].set_xlabel("t")
                    # ax[0].set_ylabel("x")
                    # ax[1].set_ylabel("y")
                    # ax[0].plot(t, [x[0] for x in new_min_val], label="x_min")
                    # ax[0].plot(t, [x[0] for x in new_max_val], label="x_max")
                    # ax[1].plot(t, [x[1] for x in new_min_val], label="y_min")
                    # ax[1].plot(t, [x[1] for x in new_max_val], label="y_max")
                    # for j in range(i):
                    #     ax[0].scatter(new_linspace[j], fixed_result[new_indices[j]][0])
                    #     ax[1].scatter(new_linspace[j], fixed_result[new_indices[j]][1])
                    # ax[0].legend()
                    # ax[1].legend()
                    # plt.show()
                    allow=True
                    for checkpoint in range(len(new_linspace)):
                        print("check:", fixed_result[new_indices[checkpoint]])
                        if fixed_result[new_indices[checkpoint]][0] < new_min_val[new_indices[checkpoint]][0] or \
                            fixed_result[new_indices[checkpoint]][0] > new_max_val[new_indices[checkpoint]][0] or \
                            fixed_result[new_indices[checkpoint]][1] < new_min_val[new_indices[checkpoint]][1] or \
                            fixed_result[new_indices[checkpoint]][1] > new_max_val[new_indices[checkpoint]][1]:
                            steps[i][k] /= 2
                            allow=False
                            break
                    if allow:  
                        new_temp_sum = sum(x[1]-x[0] for x in new_interval)
                        if new_temp_sum < temp_sum:
                            temp_sum = new_temp_sum
                            result_intervals = new_interval
                            steps[i][k] *= 2
                            if abs(new_temp_sum - prev_sum) < stop_eps:
                                stop_flag = 1
                    else:
                        new_temp_sum_no_success = sum(x[1]-x[0] for x in new_interval)
                        # print("l", new_temp_sum_no_success, prev_sum_no_success)
                        if abs(new_temp_sum_no_success - prev_sum_no_success) < stop_eps or abs(new_temp_sum_no_success - prev_sum) < stop_eps:
                            print("here")
                            stop_flag = 1
                    print(result_intervals)
                
            if stop_flag:
                break        
            
        print(f"result for {i} dots: {result_intervals}, efficency: {new_sum/start_sum}")


            

