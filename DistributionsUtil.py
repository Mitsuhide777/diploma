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
            all_dists['betaprime'] = (param_dict['betaprime'], generate_betaprime(param_dict['betaprime']))
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

def sgn(x):
    if x > 0:
        return 1
    elif x < 0:
        return -1
    else:
        return 0

def get_row_col_num(st_arr):
    if len(st_arr) == 0:
        return (0, 0)
    if len(st_arr) == 1:
        return (1, 0)

    row_num = 2

    prev_diff = sgn(st_arr[1][1] - st_arr[0][1])
    for i in range(2, len(st_arr)):
        curr_diff = sgn(st_arr[i][1] - st_arr[i-1][1])
        if curr_diff == prev_diff:
            prev_diff = curr_diff
            row_num += 1
        else:
            break

    return (row_num, len(st_arr) // row_num)

def within(st, st_borders_list):
    ratio02 = (st_borders_list[2][1]-st_borders_list[0][1]) / (st_borders_list[2][0]-st_borders_list[0][0])
    bias02 = st_borders_list[0][1] - ratio02*st_borders_list[0][0]
    ratio23 = (st_borders_list[3][1] - st_borders_list[2][1]) / (st_borders_list[3][0] - st_borders_list[2][0])
    bias23 = st_borders_list[2][1] - ratio23 * st_borders_list[2][0]
    ratio31 = (st_borders_list[1][1] - st_borders_list[3][1]) / (st_borders_list[1][0] - st_borders_list[3][0])
    bias31 = st_borders_list[1][1] - ratio31 * st_borders_list[1][0]
    ratio10 = (st_borders_list[0][1] - st_borders_list[1][1]) / (st_borders_list[0][0] - st_borders_list[1][0])
    bias10 = st_borders_list[1][1] - ratio10 * st_borders_list[1][0]

    cond1 = st[1] >= ratio02*st[0] + bias02
    cond2 = st[1] <= ratio23*st[0] + bias23
    cond3 = st[1] <= ratio31*st[0] + bias31
    cond4 = st[1] >= ratio10*st[0] + bias10
    return cond1 and cond2 and cond3 and cond4

def generate_samples(dist, curr_params, n):
    if dist == 'beta':
        a, b = curr_params
        curr_sample = beta.rvs(a, b, size=n)
    elif dist == 'betaprime':
        a, b = curr_params
        curr_sample = betaprime.rvs(a, b, size=n)
    else:
        return (None, None)

    curr_sample.sort()
    curr_e = [curr_sample[int(j * n)] for j in OCTILES]
    curr_s = (curr_e[5] - 2 * curr_e[3] + curr_e[1]) / (curr_e[5] - curr_e[1])
    curr_t = (curr_e[6] - curr_e[4] + curr_e[2] - curr_e[0]) / (curr_e[5] - curr_e[1])

    return (curr_s, curr_t)
