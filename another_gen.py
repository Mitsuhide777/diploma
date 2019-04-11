from numpy import arange
from scipy.stats import beta
from math import sqrt

N_SPLIT = 10
SAMPLE_SIZE = 300000
ITER_SIZE = 1000
OCTILES = [0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875]
# A_MIN = 21
# A_MAX = 61
# B_MIN = 31
# B_MAX = 71
# SCALE = 10.
#
# a_arr = [i / SCALE for i in range(A_MIN, A_MAX+1)]
# b_arr = [i / SCALE for i in range(B_MIN, B_MAX+1)]

A_MIN = 2.1
A_MAX = 6.1
A_STEP = 0.1
B_MIN = 3.1
B_MAX = 7.1
B_STEP = 0.1

a_arr = arange(A_MIN, A_MAX+A_STEP, A_STEP)
b_arr = arange(B_MIN, B_MAX+B_STEP, B_STEP)

a_len = len(a_arr)
b_len = len(b_arr)
ab_len = a_len * b_len

ab_arr = []
st_arr = []

# calculate S_T for beta distribution with given a_b
for a in a_arr:
    for b in b_arr:
        curr_e = [beta.ppf(o, a, b) for o in OCTILES]
        curr_s = (curr_e[5] - 2 * curr_e[3] + curr_e[1]) / (curr_e[5] - curr_e[1])
        curr_t = (curr_e[6] - curr_e[4] + curr_e[2] - curr_e[0]) / (curr_e[5] - curr_e[1])

        ab_arr.append((a, b))
        st_arr.append((curr_s, curr_t))
        # print()

# calculate border indices for a_b
a_border_inds = [(a_len-1)*i//N_SPLIT for i in range(0, N_SPLIT+1)]
b_border_inds = [(b_len-1)*i//N_SPLIT for i in range(0, N_SPLIT+1)]

# and take border values for a_b
a_borders = [a_arr[ind] for ind in a_border_inds]
b_borders = [b_arr[ind] for ind in b_border_inds]
# s_borders = [st_arr[a*b_len + b][0] for a, b in zip(a_border_inds, b_border_inds)]
# t_borders = [st_arr[a*b_len +b][1] for a, b in zip(a_border_inds, b_border_inds)]
st_borders = [[st_arr[a*b_len +b] for b in b_border_inds] for a in a_border_inds]

# get a_b centers
ab_centers = []
# for i in range(5, N_SPLIT-3):
#     for j in range(5, N_SPLIT-3):
#         ab_centers.append(((a_borders[i] + a_borders[i-1]) / 2., (b_borders[j] + b_borders[j-1]) / 2.))
ab_centers.append((a_borders[N_SPLIT//2], b_borders[N_SPLIT//2]))

def dist(x, y):
    x1, x2 = x
    y1, y2 = y
    return sqrt((x1-y1)**2. + (x2-y2)**2.)

# def st_dist(st, st_borders_list):
#     norm = dist(st_borders_list[0], st_borders_list[1]) + dist(st_borders_list[0], st_borders_list[2])\
#            + dist(st_borders_list[0], st_borders_list[3])
#     curr_dist = dist(st, st_borders_list[0]) + dist(st, st_borders_list[1])\
#            + dist(st, st_borders_list[2]) + dist(st, st_borders_list[3])
#     return curr_dist <= norm

def st_dist(st, st_borders_list):
    center = ((st_borders_list[0][0]+st_borders_list[3][0])/2.,(st_borders_list[1][1]+st_borders_list[2][1])/2.)
    # norm = dist(st, st_borders_list[0]) + dist(st, st_borders_list[1])\
    #        + dist(st, st_borders_list[2]) + dist(st, st_borders_list[3])
    # return (dist(st, center) <= max(dist(center, st_borders_list[0]), dist(center, st_borders_list[1])), dist(st, center))
    return

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
    # cond1 = (st[1] >= st_borders_list[0][1])
    # cond2 = (st[0] >= st_borders_list[2][0])
    # cond3 = (st[1] <= st_borders_list[3][1])
    # cond4 = (st[0] <= st_borders_list[1][0])
    # cond5 = (st[0] <= st_borders_list[0][0]) and (st[1] <= st_borders_list[2][1])
    # cond6 = (st[1] >= st_borders_list[2][1]) and (st[0] <= st_borders_list[3][0])
    # cond7 = (st[0] >= st_borders_list[3][0]) and (st[1] >= st_borders_list[1][1])
    # cond8 = (st[1] <= st_borders_list[1][1]) and (st[0] >= st_borders_list[0][0])

    # return cond1 & cond2 & cond3 & cond4 # square area
    return cond1 and cond2 and cond3 and cond4
# generate samples
sample_s_arr = []
sample_t_arr = []
# junk (for checking purposes)
true_s_arr = []
true_t_arr = []
# junk
# print(st_arr[0], st_arr[1])
hits = []
for n in range(10):
    for a, b in ab_centers:
        # junk (for checking purposes)
        curr_e = [beta.ppf(o, a, b) for o in OCTILES]
        curr_s = (curr_e[5] - 2 * curr_e[3] + curr_e[1]) / (curr_e[5] - curr_e[1])
        curr_t = (curr_e[6] - curr_e[4] + curr_e[2] - curr_e[0]) / (curr_e[5] - curr_e[1])

        ab_arr.append((a, b))
        true_s_arr.append(curr_s)
        true_t_arr.append(curr_t)
        # junk

        curr_sample = beta.rvs(a, b, size=SAMPLE_SIZE)
        curr_sample.sort()
        curr_e = [curr_sample[int(j * SAMPLE_SIZE)] for j in OCTILES]
        curr_s = (curr_e[5] - 2 * curr_e[3] + curr_e[1]) / (curr_e[5] - curr_e[1])
        curr_t = (curr_e[6] - curr_e[4] + curr_e[2] - curr_e[0]) / (curr_e[5] - curr_e[1])

        sample_s_arr.append(curr_s)
        sample_t_arr.append(curr_t)

        min_dist = 1000.
        hit = -1
        for i in range(N_SPLIT):
            for j in range(N_SPLIT):
                curr_dist = st_dist((curr_s, curr_t),
                                      [st_borders[i][j], st_borders[i][j+1], st_borders[i+1][j], st_borders[i+1][j+1]])
                # if curr_dist < min_dist:
                #     min_dist = curr_dist
                # if curr_dist[0] & (curr_dist[1] < min_dist):
                #     min_dist = curr_dist[1]
                if within((curr_s, curr_t),
                                      [st_borders[i][j], st_borders[i][j+1], st_borders[i+1][j], st_borders[i+1][j+1]]):
                    hit = i*N_SPLIT + j
        hits.append(hit)

# for i in range(len(true_s_arr)):
#     print((true_s_arr[i], true_t_arr[i]))
#     print((sample_s_arr[i], sample_t_arr[i]))
print(hits)

from matplotlib import pyplot as plt

# #plot a_b
# plt.plot([el[0] for el in ab_arr], [el[1] for el in ab_arr], 'o')
#
# # draw a_b borders
# for a_border, b_border in zip(a_borders, b_borders):
#         plt.plot([a_borders[0], a_borders[-1]], [b_border]*2, 'k')
#         plt.plot([a_border]*2, [b_borders[0], b_borders[-1]], 'k')
#
# # draw a_b centers
# plt.plot([el[0] for el in ab_centers], [el[1] for el in ab_centers], 'ro')


# plot s_t
plt.plot([s for s, t in st_arr], [t for s, t in st_arr], 'o')

# draw s_t borders for a lines
for a_ind in a_border_inds:
    plt.plot([s for s, t in st_arr][a_ind*b_len:(a_ind+1)*b_len], [t for s, t in st_arr][a_ind*b_len:(a_ind+1)*b_len], 'k')

# draw s_t borders for b lines
for b_ind in b_border_inds:
    plt.plot([[s for s, t in st_arr][b_ind+i*b_len] for i in range(a_len)],
             [[t for s, t in st_arr][b_ind+i*b_len] for i in range(a_len)], 'k')

# draw sample s_t
plt.plot(sample_s_arr, sample_t_arr, 'ro')
# plt.plot([[s for s, t in st_row] for st_row in st_borders], [[t for s, t in st_row] for st_row in st_borders], 'ro')

plt.show()