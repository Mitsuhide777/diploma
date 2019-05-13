from scipy.stats import beta, betaprime
OCTILES = [0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875]

def calculate_st(sample):
    n = len(sample)
    sample_sorted = sorted(sample)

    curr_e = [sample_sorted[int(j * n)] for j in OCTILES]
    curr_s = (curr_e[5] - 2 * curr_e[3] + curr_e[1]) / (curr_e[5] - curr_e[1])
    curr_t = (curr_e[6] - curr_e[4] + curr_e[2] - curr_e[0]) / (curr_e[5] - curr_e[1])

    return (curr_s, curr_t)

def distribute_st(param_dict):
    all_dists = {}
    for d in param_dict.keys():
        curr_dist_st = []
        if d == 'beta':
            all_dists['beta'] = (param_dict['beta'], generate_beta(param_dict['beta']))
        elif d == 'betaprime':
            all_dists['betaprime'] = (param_dict['betaprime'], generate_beta(param_dict['betaprime']))
        else:
            continue

    return all_dists

def generate_beta(param_list):
    a_arr = param_list['a']
    b_arr = param_list['b']

    st_arr = []
    for curr_a in a_arr:
        for curr_b in b_arr:
            curr_e = [beta.ppf(o, curr_a, curr_b) for o in OCTILES]
            curr_s = (curr_e[5] - 2 * curr_e[3] + curr_e[1]) / (curr_e[5] - curr_e[1])
            curr_t = (curr_e[6] - curr_e[4] + curr_e[2] - curr_e[0]) / (curr_e[5] - curr_e[1])

            # ab_arr.append((curr_a, curr_b))
            st_arr.append((curr_s, curr_t))

    return st_arr


def generate_betaprime(param_list):
    a_arr = param_list['a']
    b_arr = param_list['b']

    st_arr = []
    for curr_a in a_arr:
        for curr_b in b_arr:
            curr_e = [betaprime.ppf(o, curr_a, curr_b) for o in OCTILES]
            curr_s = (curr_e[5] - 2 * curr_e[3] + curr_e[1]) / (curr_e[5] - curr_e[1])
            curr_t = (curr_e[6] - curr_e[4] + curr_e[2] - curr_e[0]) / (curr_e[5] - curr_e[1])

            # ab_arr.append((curr_a, curr_b))
            st_arr.append((curr_s, curr_t))

    return st_arr