from scipy.stats import beta, betaprime
OCTILES = [0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875]


def calculate_st(sample):
    """
    This method calculates quantile skewness and quantile kurtosis of a sample
    via sorting the sample and getting its empiric octiles which are used in formulas
    for quantile skewness and kurtosis
    :param sample: list[float]
        A sample of float numbers
    :return: tuple[float. float]
        A pair (quantile skewness, quantile kurtosis) of the sample
    """
    n = len(sample)
    sample_sorted = sorted(sample)

    curr_e = [sample_sorted[int(j * n)] for j in OCTILES]
    curr_s = (curr_e[5] - 2 * curr_e[3] + curr_e[1]) / (curr_e[5] - curr_e[1])
    curr_t = (curr_e[6] - curr_e[4] + curr_e[2] - curr_e[0]) / (curr_e[5] - curr_e[1])

    return (curr_s, curr_t)


def distribute_st(param_dict):
    """
    This method generates a list of pairs (quantile skewness, quantile, kurtosis)
    for given distributions with their lists of parameters
    :param param_dict: dict[str, dict[str, list[float]]]
        A dictionary with keys which correspond to various distributions and values
        that store specific distribution's parameter ranges
    :return: dict[str, tuple[dict[str, list[float]], list[float]]
        A dictionary with keys which correspond to various distributions and values
        that store specific distribution's parameter ranges and pairs (quantile skewness,
        quantile kurtosis) for each specific combination of parameters
    """
    all_dists = {}
    for d in param_dict.keys():
        if d == 'beta':
            all_dists['beta'] = (param_dict['beta'], generate_beta(param_dict['beta']))
        elif d == 'betaprime':
            all_dists['betaprime'] = (param_dict['betaprime'], generate_betaprime(param_dict['betaprime']))
        else:
            continue

    return all_dists


def generate_beta(param_list):
    """
    This method generates a list of pairs (quantile skewness, quantile kurtosis)
    for beta distribution with given parameter ranges
    :param param_list: dict[str, list[float]]
        A dictionary with keys which correspond to the names of parameters of a
        given distribution ('a' and 'b' in case of beta distribution) and values
        that store specific parameter ranges
    :return: list[tuple[float, float]]
        A list of pairs (quantile skewness, quantile kurtosis) which can be viewed
        as a point and act as a representation of beta distribution with specific
        parameters on a Pearson System
    """
    a_arr = param_list['a']
    b_arr = param_list['b']

    st_arr = []
    for curr_a in a_arr:
        for curr_b in b_arr:
            curr_e = [beta.ppf(o, curr_a, curr_b) for o in OCTILES]
            curr_s = (curr_e[5] - 2 * curr_e[3] + curr_e[1]) / (curr_e[5] - curr_e[1])
            curr_t = (curr_e[6] - curr_e[4] + curr_e[2] - curr_e[0]) / (curr_e[5] - curr_e[1])

            st_arr.append((curr_s, curr_t))

    return st_arr


def generate_betaprime(param_list):
    """
    This method generates a list of pairs (quantile skewness, quantile kurtosis)
    for inverse beta distribution with given parameters
    :param param_list: dict[str, list[float]]
        A dictionary with keys which correspond to the names of parameters of a
        given distribution ('a' and 'b' in case of inverse beta distribution) and
        values that store specific parameter ranges
    :return: list[tuple[float, float]]
        A list of pairs (quantile skewness, quantile kurtosis) which can be viewed
        as a point and act as a representation of inverse beta distribution with
        specific parameters on a Pearson System
    """
    a_arr = param_list['a']
    b_arr = param_list['b']

    st_arr = []
    for curr_a in a_arr:
        for curr_b in b_arr:
            curr_e = [betaprime.ppf(o, curr_a, curr_b) for o in OCTILES]
            curr_s = (curr_e[5] - 2 * curr_e[3] + curr_e[1]) / (curr_e[5] - curr_e[1])
            curr_t = (curr_e[6] - curr_e[4] + curr_e[2] - curr_e[0]) / (curr_e[5] - curr_e[1])

            st_arr.append((curr_s, curr_t))

    return st_arr


def sgn(x):
    """
    This method returns an int which specifies a number's sign
    :param x: float
        A number
    :return: int
        1, if x is positive
        0, if x equals to zero
        -1, if x is negative
    """
    if x > 0:
        return 1
    elif x < 0:
        return -1
    else:
        return 0


def get_row_col_num(st_arr):
    """
    This method calculates of how many rows and columns of points consists an area
     representing a distribution on a Pearson System
    :param st_arr: list[tuple[float, float]]
        A list of pairs (quantile skewness, quantile kurtosis) which can be viewed
        as a point and act as a representation of a distribution with specific
        parameters on a Pearson System
    :return: tuple[int, int]
        A pair (number of rows, number of columns)
    """
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
    """
    This method determines whether a point is inside a curvilinear rectangle
    :param st: tuple[float, float]
        A pair (quantile skewness, quantile kurtosis) which can be viewed
        as a point and act as a representation of a distribution with specific
        parameters on a Pearson System
    :param st_borders_list: list[tuple[float, float]]
        A list of four pairs (quantile skewness, quantile kurtosis) which act as points
        on a Pearson System and as vertices for curvilinear rectangles
    :return: bool
        true, if st is inside a curvilinear rectangle created by points in st_borders_list
        false, otherwise
    """
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
    """
    This method generates a sample from a given distribution with
    given parameters
    :param dist: str
        Type of distribution from which to draw a sample
    :param curr_params: tuple[float, float]
        Parameters for a distribution
    :param n: int
        Length of generated sample
    :return: tuple[float, float]
        A pair (quantile skewness, quantile kurtosis) which can be viewed
        as a point and act as a representation of the generated sample on a Pearson System
    """
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
