from distributions_util import calculate_st

def read_sample(filename):
    """
    This method reads a sample from a file and calculates its
    quantile skewness and kurtosis
    :param filename: str
        A path to a file from which to read a sample
    :return: tuple[float, float]
        A pair (quantile skewness, quantile kurtosis) of the sample
        read from a file
    """
    with open(filename, 'r') as file:
        lines = file.readlines()
        sample = [float(line.strip('\n')) for line in lines if line != '\n']
    return calculate_st(sample)
