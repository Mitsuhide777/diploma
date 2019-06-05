from .distributions_util import calculate_st
from os.path import exists
from os.path import isfile
from .exceptions import FileParseException

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
    if exists(filename):
        if isfile(filename):
            with open(filename, 'r') as file:
                try:
                    lines = file.readlines()
                    sample = [float(line.strip('\n')) for line in lines if line != '\n']
                    print('Sample read from file.')
                except Exception:
                    raise FileParseException('Something went wrong when reading a sample. Check the contents of the sample file.')
        else:
            raise FileNotFoundError(filename + ' is not a file.')
    else:
        raise FileExistsError('File ' + filename + ' does not exist.')

    return calculate_st(sample)
