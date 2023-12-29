import numpy as np
import matplotlib.pyplot as plt


class Node:
    
    
    def __init__(self, number, borders, p, k, d, ancestor=None, multiplier=10):
        self.number = number
        self.borders = borders
        self.k = k
        self.d = d
        self.ancestor = ancestor
        self.successors = []
        self.is_leaf = True
        self.set_grid(p)
        self.plot_dots_num = p * multiplier
        self.set_plot_dots(self.plot_dots_num)
    

    def set_grid(self, p):
        points = [np.linspace(self.borders[i][0], self.borders[i][1], p+1) for i in range(len(self.borders))]
        grid = np.meshgrid(*points)
        self.dots = np.vstack([grid[i].ravel() for i in range(len(self.borders))]).T


    def set_plot_dots(self, p):
        points = [np.linspace(self.borders[i][0], self.borders[i][1], p+1) for i in range(len(self.borders))]
        grid = np.meshgrid(*points)
        self.plot_dots = np.vstack([grid[i].ravel() for i in range(len(self.borders))]).T
    


class KDTree:
    plot_var = []
    nodes = []
    node_count = 0


    def __init__(self, init_boundaries, m, p):
        self.nodes.append(Node(0, init_boundaries, p, 0, 0))
        self.node_count+=1
        self.m = m
        self.p = p
        if m==2:
            self.plot_var.append([[init_boundaries[0][0], init_boundaries[0][1]], [init_boundaries[1][0], init_boundaries[1][0]]])
            self.plot_var.append([[init_boundaries[0][0], init_boundaries[0][1]], [init_boundaries[1][1], init_boundaries[1][1]]])
            self.plot_var.append([[init_boundaries[0][0], init_boundaries[0][0]], [init_boundaries[1][0], init_boundaries[1][1]]])
            self.plot_var.append([[init_boundaries[0][1], init_boundaries[0][1]], [init_boundaries[1][0], init_boundaries[1][1]]])


    def separate_node(self, i):
        cur_d = self.nodes[i].d
        self.nodes[i].is_leaf = False
        border_to_change = self.nodes[i].borders[cur_d].copy()
        middle = (border_to_change[0] + border_to_change[1]) / 2
        new_temp_1 = [border_to_change[0], middle]
        new_temp_2 = [middle, border_to_change[1]]
        
        new_border_1 = self.nodes[i].borders.copy()
        new_border_2 = self.nodes[i].borders.copy()
        
        new_border_1[cur_d] = new_temp_1.copy()
        new_border_2[cur_d] = new_temp_2.copy()

        self.nodes.append(Node(self.node_count, new_border_1, self.p, 0, (cur_d + 1) % self.m, i))
        self.nodes[i].successors.append(self.node_count)
        self.node_count += 1

        self.nodes.append(Node(self.node_count, new_border_2, self.p, 0, (cur_d + 1) % self.m, i))
        self.nodes[i].successors.append(self.node_count)
        self.node_count += 1

        
        if self.m == 2:
            if cur_d == 0:
                self.plot_var.append([[middle, middle], [new_border_1[1][0], new_border_1[1][1]]])
            else:
                self.plot_var.append([[new_border_1[0][0], new_border_1[0][1]], [middle, middle]])



    def print_nodes(self):
        for i in range(self.node_count):
            print(f"{i} {self.nodes[i].borders} {self.nodes[i].ancestor} {self.nodes[i].successors}")


    def plot_grid(self):
        print(self.m)
        if self.m == 2:
            fig, ax = plt.subplots(1, 1, figsize=(10, 10))
            for i in range(len(self.plot_var)):
                ax.plot(self.plot_var[i][0], self.plot_var[i][1], linewidth=3)
            for i in range(self.node_count):
                for j in range(len(self.nodes[i].dots)):
                    ax.scatter(self.nodes[i].dots[j][0], self.nodes[i].dots[j][1])
            ax.grid()
            plt.show()



if __name__ == "__main__":
    tree = KDTree([[0.5, 2.5], [0.5, 2.5], [0.5, 2.5]], 2, 2)
    tree.separate_node(0)
    print(tree.nodes[0].dots)
    tree.separate_node(1)
    tree.separate_node(2)
    tree.separate_node(4)
    tree.separate_node(5)
    tree.separate_node(7)
    tree.separate_node(12)
    tree.separate_node(13)
    tree.separate_node(16)
    tree.print_nodes()
    tree.plot_grid()