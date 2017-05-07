"""
Collection of utility methods related to wavelet analysis
"""
from tools.utils import *
import mlpy.wavelet as wavelet

# ------------------------------------------


def haar(arr):
    """
    Wrapper function: Refer to http://mlpy.sourceforge.net/docs/3.5/wavelet.html
    :param arr: 1-dimensional numeric numpy array
    :return: 1-dimensional numeric numpy array consisting of Haar wavelet transformed data
    THROWS EXCEPTION IF length of arr not a power of 2
    """
    try:
        check_if_numeric_array(arr)
    except AssertionError:
        raise Exception("Problem with input arguments: One dimensional numeric Numpy array expected")

    try:
        wt_array = wavelet.dwt(x=arr, wf=b'h', k=2)
    except:  # length of arr should be a power of 2, most common cause of exception here
        raise Exception("Length of input array should be a power of 2")
    return wt_array


# ------------------------------------------


def get_detail_coeffs(wt_array):
    """
    wt_array is the result of invoking Discrete Wavelet Transform via mlpy.wavelet.dwt
    Reference: http://mlpy.sourceforge.net/docs/3.5/wavelet.html
    Keyword arguments:
    wt_array -- a one dimensional numpy array storing all the detail coefficients
             -- it's first element is the smoothing coefficient
    Output: a tuple:
        First element is an integer = number of levels in wavelet analysis (denoted by J in reference documentation)
        Second object is a dictionary with level number as key and corresponding detail coeffs as numpy array
        Note that level count starts from zero (ranges from 0 .. J-1)
    """

    try:
        check_if_numeric_array(wt_array)
    except AssertionError:
        raise Exception("Problem with input arguments: One dimensional numeric Numpy array expected")

    #todo - check length of wt_array should be a power of 2. Not an issue if wt_array is the output from haar()

    n = len(wt_array)
    number_of_levels = int(np.log2(n))
    detail_coefficients = {}

    detail_coefficient_count_array = np.zeros(number_of_levels)
    for level in range(number_of_levels):
        detail_coefficient_count = 2**level
        detail_coefficient_count_array[level] = detail_coefficient_count

    # Now for each level, extract the corresponding coefficients

    for level in range(number_of_levels):
        end_pos = int(sum(detail_coefficient_count_array[0:level+1]))
        detail_coefficients[level] = wt_array[end_pos-int(2**level)+1:end_pos+1]

    return number_of_levels, detail_coefficients


# ------------------------------------------


def get_cutoffs(detail_coefficients):
    """
    Keyword arguments:
    detail_coefficients -- a dictionary (one of the return objects from get_detail_coeffs function)
    Output: a dictionary with level number as key and corresponding cutoffs as value
    """

    level_cutoffs = {}
    constant_factor = (1.0/0.6745)*np.sqrt(2.0*len(detail_coefficients))
    # Note: len(detail_coefficients) is the number of levels
    for level in detail_coefficients:
        detail_coeff_array = detail_coeffs[level]
        median_absolute_deviation = np.median(np.absolute(detail_coeff_array - np.median(detail_coeff_array)))
        cutoff = constant_factor*median_absolute_deviation
        level_cutoffs[level] = cutoff
    return level_cutoffs


# ------------------------------------------


if __name__ == "__main__":
    wt = np.array([-1, 0, 10, 11, 20, 21, 22, 23, 30, 31, 32, 33, 34, 35, 36, 37])
    # The above array is completely arbitrary, the idea is that level 0 should have one coefficient,
    # level 1 should have 2 (10, 11), level 2 should have 4 (20, 21, 22, 23), etc.
    level_count, detail_coeffs = get_detail_coeffs(wt)
    print(level_count)  # wt has 16 numbers =  2^4, so level_count = 4
    print(detail_coeffs)
    print(get_cutoffs(detail_coeffs))
