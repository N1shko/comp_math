import numpy as np
import matplotlib.pyplot as plt
from self_made_kd_tree import KDTree
from calc_func import CalcTools
import random
plt.rcParams.update({'font.size': 20})

k = 1 
b = 0.3
c = 0.7
x_0 = 0.8
y_0 = 0.6
    
fl = 0
fl2 = 0
def func(x, koeffs):
    #return np.array([-k * x[0], 2*k*x[0]])
    return np.array([-koeffs[0] * x[0] + b * x[0] * x[1], c * x[1] - koeffs[1] * x[0] * x[1]], dtype = np.float32)


def calc_interpolation(inter, function, time_integr, h_integr, integr_eps, m, p_degree):
    
    tool = CalcTools()

    t = np.arange(0, time_integr, h_integr, dtype = np.float64)
    n = t.size
    min_val = np.zeros((n, m))
    max_val = np.zeros((n, m))
    tree = KDTree(inter, m, p_degree, function, h_integr, 2)
    separation_candidates = []
    for k in range(len(t)):
        print(k)
        is_set = False
        global_flag = False
        while not global_flag:
            print("here")
            local_flag = True
            if (t[k] == t[-1] or fl) and fl2:
                fig, ax = plt.subplots(1, 1, figsize=(10, 10))
                ax.grid()
                ax.set_xlabel(f"t = {t[k]}")    
            for i in range(len(separation_candidates)):
                tree.separate_node(separation_candidates[i])
            separation_candidates = []
            for leaf in tree.leaf_list:
                # print("ok not here")
                gor_grid_x = np.zeros(len(tree.nodes[leaf].dots), dtype=np.float64)
                gor_grid_y = np.zeros(len(tree.nodes[leaf].dots), dtype=np.float64)
                if tree.nodes[leaf].is_new or tree.nodes[leaf].just_made:
                    for i in range(len(tree.nodes[leaf].dots)):
                        dot = tree.nodes[leaf].dots[i]
                        print(dot)
                        res = tool.rk4([x_0, y_0], t[k], function, h_integr, koeffs=dot)
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
                        dot = tree.nodes[leaf].dots[ik]
                        res = tool.rk4_iteration([tree.nodes[leaf].dots_val[ik][0], tree.nodes[leaf].dots_val[ik][1]], t[k], function, h_integr, koeffs=dot)
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
                        x_for_plotting[i] = tool.lagrange_interpolant_2d(gor_grid_x, dot, tree.nodes[leaf].borders, p_degree, m)
                        y_for_plotting[i] = tool.lagrange_interpolant_2d(gor_grid_y, dot, tree.nodes[leaf].borders, p_degree, m)
                    num = tree.nodes[leaf].plot_dots_num
                    h_plot = num // (p_degree)
                    for i in range(0, num+1, h_plot):
                        ax.plot(x_for_plotting[i*(num+1):(i+1)*(num+1)], y_for_plotting[i*(num+1):(i+1)*(num+1)])
                        ax.plot(x_for_plotting[i::num+1], y_for_plotting[i::num+1])

                x_test_real = np.zeros(len(tree.nodes[leaf].random_dots), dtype=np.float64)
                y_test_real = np.zeros(len(tree.nodes[leaf].random_dots), dtype=np.float64)

                x_test_interpolation = np.zeros(len(tree.nodes[leaf].random_dots), dtype=np.float64)
                y_test_interpolation = np.zeros(len(tree.nodes[leaf].random_dots), dtype=np.float64)
                for i in range(len(tree.nodes[leaf].random_dots)):
                    dot = tree.nodes[leaf].random_dots[i]
                    res = tool.rk4([x_0, y_0], t[k], function, h_integr, koeffs=dot)
                    x_test_real[i] = res[0]
                    y_test_real[i] = res[1]
                    x_test_interpolation[i] = tool.lagrange_interpolant_2d(gor_grid_x, dot, tree.nodes[leaf].borders, p_degree, m)
                    y_test_interpolation[i] = tool.lagrange_interpolant_2d(gor_grid_y, dot, tree.nodes[leaf].borders, p_degree, m)
                error_x = 0
                error_y = 0
                for i in range(len(x_test_real)):
                    error_x += np.abs(x_test_real[i] - x_test_interpolation[i])
                    error_y += np.abs(y_test_real[i] - y_test_interpolation[i])
                avg_error_x = error_x / len(x_test_real)
                avg_error_y = error_y / len(y_test_real)
                if (avg_error_x + avg_error_y/2)>integr_eps:
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
        exit()
        
    # tree.plot_grid()
    return t, min_val, max_val


def optimize_intervals(start_intervals, function, time_integr, h_integr, integr_eps, steps=None, p_degree=3):
    m = len(start_intervals)
    print("am i even here")
    t, min_val, max_val = calc_interpolation(start_intervals, function, time_integr, h_integr, integr_eps, m, p_degree)
    print("am i even here")
    fig, ax = plt.subplots(m, 1, figsize=(10, 10))
    for i in range(m):
        ax[i].grid()
        ax[i].set_xlabel("t")
        ax[i].set_ylabel("x")
        ax[i].plot(t, [x[0] for x in min_val], label="x_min")
        ax[i].plot(t, [x[0] for x in max_val], label="x_max")
        ax[i].legend()
    plt.show()
    fixed_0 = (np.array([1.0, 0.7]))
    tool = CalcTools()
    
    fixed_result = tool.pure_rk4([x_0, y_0], time_integr, function, h_integr, koeffs=fixed_0)
    print("lel")
    print(min_val, max_val, fixed_result)
    for i in range(len(fixed_result)-1):
        print(i)
        max_dispersion = (max_val[i][0]- min_val[i][0]) * 0.05
        while True:
            check = fixed_result[i][0] + random.uniform(-max_dispersion, max_dispersion) 
            # print(check)
            if check >= min_val[i][0] and check <= max_val[i][0]:
                fixed_result[i][0] = check
                break
        max_dispersion = (max_val[i][1]- min_val[i][1]) * 0.05
        while True:
            check = fixed_result[i][1] + random.uniform(-max_dispersion, max_dispersion) 
            if check >= min_val[i][1] and check <= max_val[i][1]:
                fixed_result[i][1] = check
                break
    # print(fixed_result[0])
    print("fixed:", fixed_result[0])
    start_sum = sum(x[1]-x[0] for x in start_intervals)
    for i in range(5,6):
        print("am i even here")

        new_indices = np.linspace(start=0, stop=len(t)-1, num=i+1).astype(int)
        new_linspace = t[new_indices]
        plt.show()
        local_intervals = np.copy(start_intervals)
        result_intervals = np.copy(start_intervals)
        temp_sum = sum(x[1]-x[0] for x in local_intervals)
        new_sum = sum(x[1]-x[0] for x in result_intervals)
        prev_sum = new_sum
        

        num_determine = 20
        stop_eps = 1e-5
        stop_flag = 0
        for iterations in range(num_determine):
            prev_sum = temp_sum
            prev_sum_no_success = temp_sum
            for it in range(m):
                for k in range(2):
                    steps_temp = np.copy(steps)
                    new_interval = np.copy(result_intervals)
                    if k == 0:
                        while new_interval[it][0] + steps[it][0] >= new_interval[it][1]:
                            steps[it][0] /= 2
                        new_interval[it][0] += steps[it][0]
                    else:
                        while new_interval[it][1] - steps[it][1] <= new_interval[it][0]:
                            steps[it][1] /= 2
                        new_interval[it][1] -= steps[it][1]
                    # print("steps before", steps)
                    t, new_min_val, new_max_val = calc_interpolation(new_interval, function, time_integr, h_integr, integr_eps, m, p_degree)
                    if iterations == num_determine-1 or iterations == 0:
                        fig, ax = plt.subplots(2, 1, figsize=(10, 10))
                        # ax[0].plot(t, [x[0] for x in fixed_result], color='cornflowerblue', label=r'$x(t)$')
                        # ax[1].plot(t, [x[1] for x in fixed_result], color='indigo', label=r'$y(t)$')
                        ax[0].grid()
                        ax[1].grid()
                        ax[0].set_xlabel("t")
                        ax[1].set_xlabel("t")
                        ax[0].set_ylabel("x")
                        ax[1].set_ylabel("y")
                        ax[0].plot(t, [x[0] for x in new_min_val], label="x_min")
                        ax[0].plot(t, [x[0] for x in new_max_val], label="x_max")
                        ax[1].plot(t, [x[1] for x in new_min_val], label="y_min")
                        ax[1].plot(t, [x[1] for x in new_max_val], label="y_max")
                        for j in range(i + 1):
                            if not j:
                                ax[0].scatter(new_linspace[j], fixed_result[new_indices[j]][0], color='cornflowerblue', label=r'$x^*$')
                                ax[1].scatter(new_linspace[j], fixed_result[new_indices[j]][1], color='red', label=r"$y^*$")
                            else:
                                ax[0].scatter(new_linspace[j], fixed_result[new_indices[j]][0], color='cornflowerblue')
                                ax[1].scatter(new_linspace[j], fixed_result[new_indices[j]][1], color='red')                               
                        ax[0].legend()
                        ax[1].legend()
                        plt.show()
                    allow=True
                    for checkpoint in range(len(new_linspace)):
                        for iterate_through_m in range(m):
                            if fixed_result[new_indices[checkpoint]][iterate_through_m] < new_min_val[new_indices[checkpoint]][iterate_through_m] or \
                                fixed_result[new_indices[checkpoint]][iterate_through_m] > new_max_val[new_indices[checkpoint]][iterate_through_m]:
                                allow=False
                                break
                        if not allow:
                            break
                    if allow:  
                        new_temp_sum = sum(x[1]-x[0] for x in new_interval)
                        if new_temp_sum < temp_sum:
                            temp_sum = new_temp_sum
                            result_intervals = np.copy(new_interval)
                            steps[it][k] *= 2
                            if abs(new_temp_sum - prev_sum) < stop_eps:
                                stop_flag = 1
                    else:
                        steps = np.copy(steps_temp)
                        steps[it][k] /= 2
                        new_temp_sum_no_success = sum(x[1]-x[0] for x in new_interval)
                        if abs(new_temp_sum_no_success - prev_sum_no_success) < stop_eps or abs(new_temp_sum_no_success - prev_sum) < stop_eps:
                            stop_flag = 1
            if stop_flag:
                fig, ax = plt.subplots(2, 1, figsize=(10, 10))
                ax[0].plot(t, [x[0] for x in fixed_result], color='cornflowerblue', label=r'$x(t)$')
                ax[1].plot(t, [x[1] for x in fixed_result], color='indigo', label=r'$y(t)$')
                ax[0].grid()
                ax[1].grid()
                ax[0].set_xlabel("t")
                ax[1].set_xlabel("t")
                ax[0].set_ylabel("x")
                ax[1].set_ylabel("y")
                ax[0].plot(t, [x[0] for x in new_min_val], label="x_min")
                ax[0].plot(t, [x[0] for x in new_max_val], label="x_max")
                ax[1].plot(t, [x[1] for x in new_min_val], label="y_min")
                ax[1].plot(t, [x[1] for x in new_max_val], label="y_max")
                for j in range(i + 1):
                    ax[0].scatter(new_linspace[j], fixed_result[new_indices[j]][0])
                    ax[1].scatter(new_linspace[j], fixed_result[new_indices[j]][1])
                ax[0].legend()
                ax[1].legend()
                plt.show()
                print(iterations)
                break     

        print(f"result for {i+1} dots: {result_intervals}, efficency: {new_temp_sum/start_sum}")


if __name__ == '__main__':
    start_intervals = np.array([[0.9, 1.3], [0.6, 1.0]])
    steps = np.array([[0.01, 0.01], [0.01, 0.01]], dtype=np.float64)
    random.seed(0)
    max_dispersion = 0.05
    f = lambda t, x, koeffs: func(x, koeffs)
    optimize_intervals(start_intervals, function=f, time_integr=5.5, h_integr=0.5, integr_eps=1e-3, steps=steps, p_degree=4)