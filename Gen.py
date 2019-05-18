from DistributionsUtil import distribute_st, get_row_col_num, generate_samples, within
from matplotlib import pyplot as plt

class Gen:
    SAMPLE_SIZE = 100000
    ITER_SIZE = 50

    def __init__(self, param_dist, split=10):
        self.all_dists = distribute_st(param_dist)
        self.n_split = split
        self.all_dists_borders = {}
        self.border_info = {}

        for d in self.all_dists.keys():
            curr_d_param_borders = {}
            for p in self.all_dists[d][0].keys():
                p_len = len(self.all_dists[d][0][p])
                p_border_inds = [(p_len-1)*i//self.n_split for i in range(0, self.n_split+1)]
                curr_d_param_borders[p] = [self.all_dists[d][0][p][ind] for ind in p_border_inds]

            curr_d_st_borders = distribute_st({d : curr_d_param_borders})[d][1]

            self.all_dists_borders[d] = (curr_d_param_borders, curr_d_st_borders)



        # centers map here
        for d in self.all_dists.keys():
            curr_d_param_borders = self.all_dists_borders[d][0]
            curr_d_st_borders = self.all_dists_borders[d][1]

            p1_centers = []
            p2_centers = []
            # len=2 case
            curr_d_param_labels = sorted(curr_d_param_borders.keys())
            for i in range(1, len(curr_d_param_borders[curr_d_param_labels[0]])):
                p1_centers.append((curr_d_param_borders[curr_d_param_labels[0]][i] +
                                   curr_d_param_borders[curr_d_param_labels[0]][i - 1]) / 2.)

            for i in range(1, len(curr_d_param_borders[curr_d_param_labels[1]])):
                p2_centers.append((curr_d_param_borders[curr_d_param_labels[1]][i] +
                                   curr_d_param_borders[curr_d_param_labels[1]][i - 1]) / 2.)

            for p1 in p1_centers:
                for p2 in p2_centers:
                    # plt.plot(p1, p2, 'go')
                    for k in range(self.ITER_SIZE):
                        sstt = generate_samples(d, [p1, p2], self.SAMPLE_SIZE)
                        curr_el = (d, p1, p2)

                        for in_d in sorted(self.all_dists_borders.keys()):
                            st_borders_check_list = self.all_dists_borders[in_d][1]

                            for i in range(self.n_split):
                                for j in range(self.n_split):
                                    st_borders_check = (st_borders_check_list[self.n_split * i + j],
                                                        st_borders_check_list[self.n_split * i + j + 1],
                                                        st_borders_check_list[self.n_split * (i + 1) + j],
                                                        st_borders_check_list[self.n_split * (i + 1) + j + 1])
                                    if within(sstt, st_borders_check):
                                        if st_borders_check not in self.border_info:
                                            self.border_info[st_borders_check] = {}
                                        else:
                                            if curr_el not in self.border_info[st_borders_check]:
                                                self.border_info[st_borders_check][curr_el] = 1
                                            else:
                                                self.border_info[st_borders_check][curr_el] += 1
                        # print(k, sstt)
                    # print()

            # plt.show()
            print(self.border_info)



    def plot_params(self, borders):
        for i, d in enumerate(self.all_dists.keys(), 1):
            plt.figure(i)

            currd_params = self.all_dists[d][0]
            currd_param_labels = sorted(currd_params.keys())

            if len(currd_params) == 1:
                p_list = currd_params[currd_param_labels[0]]
                pass

            if len(currd_params) == 2:
                p1_list = currd_params[currd_param_labels[0]]
                p2_list = currd_params[currd_param_labels[1]]

                for p1 in p1_list:
                    for p2 in p2_list:
                        plt.plot(p1, p2, 'bo')

                plt.xlabel(currd_param_labels[0])
                plt.ylabel(currd_param_labels[1])

                if borders:
                    borders_arr = self.all_dists_borders[d][0]
                    p1_border_list = borders_arr[currd_param_labels[0]]
                    p2_border_list = borders_arr[currd_param_labels[1]]

                    # a bit of an overkill, can be achieved with less operations
                    for j, p1 in enumerate(p1_border_list):
                        for k, p2 in enumerate(p2_border_list):
                            plt.plot([p1, p1], [p2, p2_border_list[(k+1)%self.n_split]], 'k')
                            plt.plot([p1, p1_border_list[(k+1)%self.n_split]], [p2, p2], 'k')

        plt.show()


    def plot_st_map(self, true_borders, line_borders):
        for d in self.all_dists:
            st_arr = self.all_dists[d][1]

            plt.plot([s for s, t in st_arr], [t for s, t in st_arr], 'bo')

        if line_borders:
            st_borders_arr = self.all_dists_borders[d][1]
            st_borders_rows, st_borders_cols = get_row_col_num(st_borders_arr)

            for j in range(st_borders_cols):
                plt.plot([s for s, t in st_borders_arr[j*st_borders_rows:(j+1)*st_borders_rows]],
                         [t for s, t in st_borders_arr[j * st_borders_rows:(j + 1) * st_borders_rows]], 'y')

                plt.plot([s for s, t in [st_borders_arr[j + st_borders_rows * l] for l in range(st_borders_rows)]],
                         [t for s, t in [st_borders_arr[j + st_borders_rows * l] for l in range(st_borders_rows)]], 'y')

        if true_borders:

            pass

        plt.xlabel('S')
        plt.ylabel('T')

        plt.show()

    def identify_dist(self, st_coord):
        answer = (None, None)

        return answer
