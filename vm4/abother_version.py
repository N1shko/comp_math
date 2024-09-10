import numpy as np
import matplotlib.pyplot as plt
from chut_menee_govnoderevo import KDTree
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


def benzilamin(x, koeffs):
    o1 = koeffs[0] * x[0] * x[1]
    o2 = koeffs[1] * x[2]
    o3 = koeffs[2] * x[4] * x[0]
    o4 = koeffs[3] * x[7] * x[5]
    return np.array([
                    -o1 - o3,
                    -o1,
                    o1 - o2,
                    o1,
                    o2-o3,
                    -o2-o4,
                    o3,
                    o3-o4,
                    o4
            ], dtype=np.float32)

def func(x, koeffs):
    #return np.array([-k * x[0], 2*k*x[0]])
    return np.array([-koeffs[0] * x[0] + b * x[0] * x[1], 
                     c * x[1] - koeffs[1] * x[0] * x[1]], dtype=np.float32)


def chem(x, koeffs):
    #return np.array([-k * x[0], 2*k*x[0]])
    #return np.array([-koeffs[0] * x[0] + b * x[0] * x[1], c * x[1] - koeffs[1] * x[0] * x[1]], dtype = np.float32)
    return np.array([koeffs[0] * x[1], 
                     -koeffs[0] * x[1], 
                     0.5 * koeffs[0] * x[1]], 
                     dtype = np.float32)



def calc_interpolation(inter, function, time_integr, h_integr, integr_eps, m, p_degree, eq_amount, initial_descrete):
    tool = CalcTools()

    t = np.arange(0, time_integr, h_integr, dtype = np.float64)
    k = t.size
    min_val = np.zeros((k, eq_amount))
    max_val = np.zeros((k, eq_amount))
    tree = KDTree(inter, m, p_degree, function, h_integr, eq_amount)
    separation_candidates = []
    for k in range(len(t)):
        # print(k)
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
                # print("ok not here")
                gor_grid = np.zeros((eq_amount, len(tree.nodes[leaf].dots)), dtype=np.float64)
                if tree.nodes[leaf].is_new or tree.nodes[leaf].just_made:
                    for i in range(len(tree.nodes[leaf].dots)):
                        dot = tree.nodes[leaf].dots[i]
                        res = tool.rk4(initial_descrete, t[k], function, h_integr, koeffs=dot)
                        # print(eq_amount, len(res), len(gor_grid), len(tree.nodes[leaf].dots_val[i]))
                        tree.nodes[leaf].dots_val[i] = res.copy()
                        for r in range(eq_amount):
                            gor_grid[r][i] = res[r]
                        if not is_set:
                            for r in range(len(res)):    
                                min_val[k][r] = res[r]
                                max_val[k][r] = res[r]
                            is_set = True
                        else:
                            for r in range(len(res)):   
                                min_val[k][r] = min(min_val[k][r], res[r])
                                max_val[k][r] = max(max_val[k][r], res[r]) 
                    tree.nodes[leaf].is_new = False
                    tree.nodes[leaf].swap_needed = False
                else:
                    tree.nodes[leaf].swap_needed = True
                    for ik in range(len(tree.nodes[leaf].dots_val)):
                        dot = tree.nodes[leaf].dots[ik]
                        res = tool.rk4_iteration(tree.nodes[leaf].dots_val[ik], t[k], function, h_integr, koeffs=dot)
                        for r in range(eq_amount):
                            gor_grid[r][ik] = res[r]
                        tree.nodes[leaf].temp_val[ik] = res.copy()
                        if not is_set:
                            for r in range(len(res)):    
                                min_val[k][r] = res[r]
                                max_val[k][r] = res[r]
                            is_set = True
                        else:
                            for r in range(len(res)):   
                                min_val[k][r] = min(min_val[k][r], res[r])
                                max_val[k][r] = max(max_val[k][r], res[r])            
                for_plotting = np.zeros((eq_amount, len(tree.nodes[leaf].plot_dots)), dtype=np.float64)
                # if (t[k] == t[-1] or fl) and fl2:
                #     for i in range(len(tree.nodes[leaf].plot_dots)):
                #         dot = tree.nodes[leaf].plot_dots[i]
                #         for r in range(eq_amount):
                #             for_plotting[r][i] = tool.lagrange_interpolant_nd(gor_grid[r], dot, tree.nodes[leaf].borders, p_degree, m)
                #             for_plotting[r][i] = tool.lagrange_interpolant_nd(gor_grid[r], dot, tree.nodes[leaf].borders, p_degree, m)
                #     num = tree.nodes[leaf].plot_dots_num
                #     h_plot = num // (p_degree)
                #     for i in range(0, num+1, h_plot):
                #         ax.plot(for_plotting[0][i*(num+1):(i+1)*(num+1)], for_plotting[1][i*(num+1):(i+1)*(num+1)])
                #         ax.plot(for_plotting[0][i::num+1], for_plotting[1][i::num+1])
                # print(gor_grid)
                # a = input()
                test_real = np.zeros((eq_amount, len(tree.nodes[leaf].random_dots)), dtype=np.float64)

                test_interpolation = np.zeros((eq_amount, len(tree.nodes[leaf].random_dots)), dtype=np.float64)
                for i in range(len(tree.nodes[leaf].random_dots)):
                    dot = tree.nodes[leaf].random_dots[i]
                    res = tool.rk4(initial_descrete, t[k], function, h_integr, koeffs=dot)
                    for r in range(eq_amount):
                        test_real[r][i] = res[r]
                        test_interpolation[r][i] = tool.lagrange_interpolant_nd(gor_grid[r], dot, tree.nodes[leaf].borders, p_degree, m)
                error = np.zeros(eq_amount)
                for i in range(len(test_real[0])):
                    # print(test_real[r], test_interpolation[r])
                    # a = input()
                    for r in range(eq_amount):
                        error[r] += np.abs(test_real[r][i] - test_interpolation[r][i])
                avg_error = np.zeros(eq_amount)
                for r in range(eq_amount):
                    avg_error[r] = error[r] / len(test_real[r])
                    avg_error[r] = error[r] / len(test_real[r])
                # print(np.sum(avg_error))
                if (np.sum(avg_error))>integr_eps:
                    local_flag = False
                    separation_candidates.append(leaf)
            if local_flag == True:
                global_flag = True
            # if (t[k] == t[-1] or fl) and fl2:
            #     plt.show()
        for leaf in tree.leaf_list:
            tree.nodes[leaf].just_made = False
            if tree.nodes[leaf].swap_needed:
                tree.nodes[leaf].dots_val = tree.nodes[leaf].temp_val.copy()
    # tree.print_nodes()
    return t, min_val, max_val


def optimize_intervals(function, eq_amount, time_integr, h_integr, integr_eps, steps=None, p_degree=3, koeff_intervals=[], initial_intervals=[], initial_discrete=None):
    m = len(koeff_intervals) + len(initial_intervals)
    t, min_val, max_val = calc_interpolation(koeff_intervals, function, time_integr, h_integr, integr_eps, m, p_degree, eq_amount, initial_discrete)
    fixed_0 = (np.array([0.2025, 4.85, 14.75, 2.2]))
    tool = CalcTools()
    if eq_amount <= 3:
        fig, ax = plt.subplots(eq_amount, 1, figsize=(10, 10))
        for graphic in range(eq_amount):
            ax[graphic].grid()
            ax[graphic].set_xlabel("t")
            ax[graphic].set_ylabel("x")
            ax[graphic].plot(t, [x[graphic] for x in min_val], label="x_min")
            ax[graphic].plot(t, [x[graphic] for x in max_val], label="x_max")
            ax[graphic].legend()
    else:
        fig, ax = plt.subplots(eq_amount//3, eq_amount//3, figsize=(10, 10))
        for graphic in range(eq_amount):
            ax[graphic//3][graphic%3].grid()
            ax[graphic//3][graphic%3].set_xlabel("t")
            ax[graphic//3][graphic%3].set_ylabel("x")
            ax[graphic//3][graphic%3].plot(t, [x[graphic] for x in min_val], label="x_min")
            ax[graphic//3][graphic%3].plot(t, [x[graphic] for x in max_val], label="x_max")
            ax[graphic//3][graphic%3].legend()
    plt.show()  
    fixed_result = tool.pure_rk4(initial_discrete, time_integr, function, h_integr, koeffs=fixed_0)
    # exit()
    for i in range(len(fixed_result)-1):
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
    print(fixed_result[0])
    start_sum = sum(x[1]-x[0] for x in koeff_intervals)
    for i in range(5,6):

        new_indices = np.linspace(start=0, stop=len(t)-1, num=i+1).astype(int)
        new_linspace = t[new_indices]
        plt.show()
        local_intervals = np.copy(koeff_intervals)
        result_intervals = np.copy(koeff_intervals)
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
                    t, new_min_val, new_max_val = calc_interpolation(new_interval, function, time_integr, h_integr, integr_eps, m, p_degree, eq_amount, initial_discrete)
                    if iterations == num_determine-1:
                        if eq_amount <= 3:
                            fig, ax = plt.subplots(eq_amount, 1, figsize=(10, 10))
                            for graphic in range(eq_amount):
                                ax[graphic].grid()
                                ax[graphic].set_xlabel("t")
                                ax[graphic].set_ylabel("x")
                                ax[graphic].plot(t, [x[graphic] for x in new_min_val], label="x_min")
                                ax[graphic].plot(t, [x[graphic] for x in new_max_val], label="x_max")
                                for j in range(i + 1):
                                    if not j:
                                        ax[graphic].scatter(new_linspace[j], fixed_result[new_indices[j]][graphic], color='cornflowerblue', label=r'$x^*$')
                                    else:
                                        ax[graphic].scatter(new_linspace[j], fixed_result[new_indices[j]][graphic], color='cornflowerblue')
                                ax[graphic].legend()
                        else:
                            
                            for a in range(3):
                                fig, ax = plt.subplots(3, 1, figsize=(10, 10))
                                for graphic in range(3):
                                    ax[graphic].grid()
                                    ax[graphic].set_xlabel("t")
                                    ax[graphic].set_ylabel(f"x{graphic + 3*a + 1}")
                                    ax[graphic].plot(t, [x[graphic + 3*a] for x in new_min_val], label=f"x{graphic + 3*a + 1}_min")
                                    ax[graphic].plot(t, [x[graphic + 3*a] for x in new_max_val], label=f"x{graphic + 3*a + 1}_max")
                                    for j in range(i + 1):
                                        if not j:
                                            ax[graphic].scatter(new_linspace[j], fixed_result[new_indices[j]][graphic + 3*a], color='cornflowerblue', label=r'$x^*$')
                                        else:
                                            ax[graphic].scatter(new_linspace[j], fixed_result[new_indices[j]][graphic + 3*a], color='cornflowerblue')
                                    ax[graphic].legend()
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
                if eq_amount <= 3:
                    fig, ax = plt.subplots(eq_amount, 1, figsize=(10, 10))
                    for graphic in range(eq_amount):
                        ax[graphic].grid()
                        ax[graphic].set_xlabel("t")
                        ax[graphic].set_ylabel("x")
                        ax[graphic].plot(t, [x[graphic] for x in new_min_val], label="x_min")
                        ax[graphic].plot(t, [x[graphic] for x in new_max_val], label="x_max")
                        for j in range(i + 1):
                            if not j:
                                ax[graphic].scatter(new_linspace[j], fixed_result[new_indices[j]][graphic], color='cornflowerblue', label=r'$x^*$')
                            else:
                                ax[graphic].scatter(new_linspace[j], fixed_result[new_indices[j]][graphic], color='cornflowerblue')
                        ax[graphic].legend()
                else:
                    for a in range(3):
                        fig, ax = plt.subplots(3, 1, figsize=(10, 10))
                        for graphic in range(3):
                            ax[graphic].grid()
                            ax[graphic].set_xlabel("t")
                            ax[graphic].set_ylabel(f"x{graphic + 3*a + 1}")
                            ax[graphic].plot(t, [x[graphic + 3*a] for x in new_min_val], label=f"x{graphic + 3*a + 1}_min")
                            ax[graphic].plot(t, [x[graphic + 3*a] for x in new_max_val], label=f"x{graphic + 3*a + 1}_max")
                            for j in range(i + 1):
                                if not j:
                                    ax[graphic].scatter(new_linspace[j], fixed_result[new_indices[j]][graphic + 3*a], color='cornflowerblue', label=r'$x^*$')
                                else:
                                    ax[graphic].scatter(new_linspace[j], fixed_result[new_indices[j]][graphic + 3*a], color='cornflowerblue')
                            ax[graphic].legend()
                        plt.show()   
                plt.show() 
                break     

        print(f"result for {i+1} dots: {result_intervals}, efficency: {new_temp_sum/start_sum}")


if __name__ == '__main__':
    random.seed(0)

    koeff_intervals = np.array([[0.1525, 0.3206], [4.772, 5], [14.49, 15], [2.0, 2.34]])
    steps = np.array([[0.005, 0.005], [0.005, 0.005], [0.005, 0.005], [0.005, 0.005]], dtype=np.float64)
    max_dispersion = 0.05
    f = lambda t, x, koeffs: benzilamin(x, koeffs)
    optimize_intervals(koeff_intervals=koeff_intervals, eq_amount=9, function=f, time_integr=5.5, h_integr=0.5, integr_eps=1, steps=steps, p_degree=3,
                       initial_discrete=np.array([0.3738, 0.6262, 0, 0, 0, 0, 0, 0, 0]))
       
    # koeff_intervals = np.array([[0, 5.0]])
    # steps = np.array([[0.005, 0.005]], dtype=np.float64)
    # max_dispersion = 0.05
    # f = lambda t, x, koeffs: chem(x, koeffs)
    # optimize_intervals(koeff_intervals=koeff_intervals, eq_amount=3, function=f, time_integr=15.5, h_integr=0.5, integr_eps=1e-2, steps=steps, p_degree=2,
    #                    initial_discrete=np.array([0, 150, 0])) 

    # koeff_intervals = np.array([[0.9, 1.3], [0.6, 1.0]])
    # steps = np.array([[0.01, 0.01], [0.01, 0.01]], dtype=np.float64)
    # max_dispersion = 0.05
    # f = lambda t, x, koeffs: func(x, koeffs)
    # optimize_intervals(koeff_intervals=koeff_intervals, eq_amount=2, function=f, time_integr=5.5, h_integr=0.5, integr_eps=1e-3, steps=steps, p_degree=4,
    #                    initial_discrete=np.array([0.8, 0.6]))
