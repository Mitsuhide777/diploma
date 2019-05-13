from DistributionsUtil import distribute_st
from matplotlib import pyplot as plt

class Gen:
    def __init__(self, param_dist, split=10):
        self.all_dists = distribute_st(param_dist)
        self.n_split = split
        self.all_dists_borders = {}

        for d in self.all_dists.keys():
            curr_d_param_borders = {}
            for p in self.all_dists[d][0].keys():
                p_len = len(self.all_dists[d][0][p])
                p_border_inds = [(p_len-1)*i//self.n_split for i in range(0, self.n_split+1)]
                curr_d_param_borders[p] = [self.all_dists[d][0][p][ind] for ind in p_border_inds]

            curr_d_st_borders = distribute_st({d : curr_d_param_borders})[d][1]

            self.all_dists_borders[d] = (curr_d_param_borders, curr_d_st_borders)

        for i, d in enumerate(self.all_dists.keys(), 1):
            param_list = self.all_dists[d][0]
            st_arr = self.all_dists[d][1]
            st_borders_arr = self.all_dists_borders[d][1]

            plt.figure(i)
            plt.plot([s for s, t in st_arr], [t for s, t in st_arr], 'o')
            plt.plot([s for s, t in st_borders_arr], [t for s, t in st_borders_arr]) # st as matrix?

        plt.show()