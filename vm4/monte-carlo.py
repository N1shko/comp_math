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
        #             result_intervals = new_intervals