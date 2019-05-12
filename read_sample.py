from DistributionsUtil import calculate_st

def read_sample(filename):
    sample = []
    with open(filename, 'r') as file:
        lines = file.readlines()
        ll = lines[10000].strip('\n')
        sample = [float(l.strip('\n')) for l in lines if l != '\n']
    return calculate_st(sample)