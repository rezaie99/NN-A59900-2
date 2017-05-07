import numpy as np
from scipy import stats
import matplotlib.pyplot as pl
#
# todo: Check that length of both arrays are the same for relevant functions
#

# -----------------------------------------


def check_if_numeric_array(arr):
    """
    :param arr: Hopefully a  1-dimensional real numpy array, but that's what we are checking for
    :return: throws exception if not so
    """
    assert (isinstance(arr, np.ndarray)), "Numpy array expected"
    assert (1 == len(arr.shape)), "Expecting a one dimensional array"
    bool_array = np.isreal(arr)
    assert (np.all(bool_array) == True), "Real valued array expected"

# -----------------------------------------


def MAE(true_hrs, predicted_hrs,plot=False):
    """
    :param true_hrs: 1-dimensional numeric numpy array
    :param predicted_hrs: 1-dimensional numeric numpy array
    :return: mean absolute error
    """
    try:
        check_if_numeric_array(true_hrs)
        check_if_numeric_array(predicted_hrs)
    except AssertionError:
        raise Exception("Problem with input arguments: One dimensional numeric Numpy array expected")
    error = np.abs(true_hrs - predicted_hrs)

    if plot==True:
        bins=np.arange(0,100,1)
        h,b = np.histogram(error,bins=bins,density=True)
        pl.figure()
        pl.plot(b[:-1],100*h,'-+')
        pl.plot(b[:-1],100*np.cumsum(h),'-+')
        pl.grid()
        pl.xlabel('error (BPM)')
        pl.ylabel('% occurence')
        pl.legend(['PDF','CDF'])

    return round(np.mean(error), 2),round(np.std(100*(error/true_hrs)), 2)

# -----------------------------------------


def MAPE(true_hrs, predicted_hrs,plot=False, fileNameTitle=None):
    """
    :param true_hrs: 1-dimensional numeric numpy array
    :param predicted_hrs: 1-dimensional numeric numpy array
    :return: mean absolute error percentage
    """
    try:
        check_if_numeric_array(true_hrs)
        check_if_numeric_array(predicted_hrs)
    except AssertionError:
        raise Exception("Problem with input arguments: One dimensional numeric Numpy array expected")
    error = np.abs(true_hrs - predicted_hrs)

    if plot==True:
        bins=np.arange(0,101,1)
        h,b = np.histogram(100*(error/true_hrs),bins=bins,density=True)
        pl.figure()
        pl.plot(b[:-1],100*h,'-+')
        pl.plot(b[:-1],100*np.cumsum(h),'-+')
        pl.grid()
        pl.xlabel('% error')
        pl.ylabel('% occurence')
        pl.legend(['PDF','CDF'])

        # save pdf and cdf to a file
        if fileNameTitle is None:
            fileNameTitle = 'pdfCdf.txt'
        else:
            fileNameTitle = fileNameTitle + '.pdfCdf.txt'
        z = np.hstack((np.array(b[:-1])[np.newaxis].T, np.array(h)[np.newaxis].T, np.array(np.cumsum(h))[np.newaxis].T))
        np.savetxt(fileNameTitle, z, fmt=['%i','%f','%f'],delimiter=',')

    return round(np.mean(100*(error/true_hrs)), 2)

# -----------------------------------------


def scatter2D(true_hrs, predicted_hrs,plot=False):
    """
    :param true_hrs: 1-dimensional numeric numpy array
    :param predicted_hrs: 1-dimensional numeric numpy array
    :return: mean absolute error percentage
    """
    try:
        check_if_numeric_array(true_hrs)
        check_if_numeric_array(predicted_hrs)
    except AssertionError:
        raise Exception("Problem with input arguments: One dimensional numeric Numpy array expected")
    error = np.abs(true_hrs - predicted_hrs)

    if plot==True:
        bins=np.arange(0,200,2)

        # note flipping of axes, due to how histogram2d handles the axes
        hist, xedges, yedges = np.histogram2d(predicted_hrs, true_hrs, bins=(bins, bins))

        scatterFig = pl.figure()
        pl.title('HR Error 2D Histogram')
        pl.xlabel('Chest Strap')
        pl.ylabel('Whoop')
        im = pl.imshow(hist, interpolation='nearest', origin='low', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
        return scatterFig
    return None

# -----------------------------------------


def RMSE(true_hrs, predicted_hrs):
    """
    :param true_hrs: 1-dimensional numeric numpy array
    :param predicted_hrs: 1-dimensional numeric numpy array
    :return: square root of the mean of the square of the error
    """
    try:
        check_if_numeric_array(true_hrs)
        check_if_numeric_array(predicted_hrs)
    except AssertionError:
        raise Exception("Problem with input arguments: One dimensional numeric Numpy array expected")
    error_sq = (true_hrs - predicted_hrs)**2
    return round(np.sqrt(np.mean(error_sq)), 2)

# -----------------------------------------


def get_error_percentiles(true_hrs, predicted_hrs):
    """
    :param true_hrs: 1-dimensional numeric numpy array
    :param predicted_hrs: 1-dimensional numeric numpy array
    :return: 25th, 50th, 75th, 80th, 85th, 90th and 95th percentiles of the absolute error as a python list
    """
    try:
        check_if_numeric_array(true_hrs)
        check_if_numeric_array(predicted_hrs)
    except AssertionError:
        raise Exception("Problem with input arguments: One dimensional numeric Numpy array expected")
    error = np.abs(true_hrs - predicted_hrs)
    return [np.percentile(error, 25), np.percentile(error, 50), np.percentile(error, 75), np.percentile(error, 80),
            np.percentile(error, 85), np.percentile(error, 90), np.percentile(error, 95)]

# -----------------------------------------


def get_percentage_error_percentiles(true_hrs, predicted_hrs):
    """
    :param true_hrs: 1-dimensional numeric numpy array
    :param predicted_hrs: 1-dimensional numeric numpy array
    :return: 25th, 50th, 75th, 80th, 85th, 90th and 95th percentiles of the % error as a python list
    """
    try:
        check_if_numeric_array(true_hrs)
        check_if_numeric_array(predicted_hrs)
    except AssertionError:
        raise Exception("Problem with input arguments: One dimensional numeric Numpy array expected")
    error = np.abs(true_hrs - predicted_hrs)
    percentage_error = 100*(error/true_hrs)
    return [np.percentile(percentage_error, 25), np.percentile(percentage_error, 50), np.percentile(percentage_error, 75), np.percentile(percentage_error, 80),
            np.percentile(percentage_error, 85), np.percentile(percentage_error, 90), np.percentile(percentage_error, 95)]

# -----------------------------------------


def score_at_percentile(arr, percentile=50):
    """
    :param arr: 1-dimensional numeric numpy array
    :param precentile: numeric between 0 and 100, default value set for median
    :return: The numeric corresponding to the percentile specified. Wrapper around the scipy function
    """
    try:
        check_if_numeric_array(arr)
    except AssertionError:
        raise Exception("Problem with input argument: One dimensional numeric Numpy array expected")
    try:
        score = stats.scoreatpercentile(arr, percentile)
    except:
        raise Exception("Percentile parameter needs to a number between 0 and 100")
    return score

# -----------------------------------------


def get_ST_metric(true_hrs, predicted_hrs, percentage=10.0):
    """
    :param true_hrs: 1-dimensional numeric numpy array
    :param predicted_hrs: 1-dimensional numeric numpy array
    :param percentage: Numeric between 0 and 100.0
    :return: % of time the error is within 10% of the original
    """
    try:
        check_if_numeric_array(true_hrs)
        check_if_numeric_array(predicted_hrs)
    except AssertionError:
        raise Exception("Problem with input arguments: One dimensional numeric Numpy array expected")

    bool_var = 0.0 <= percentage <= 100.0
    if not bool_var:
        raise Exception("Percentage has to lie in [0, 100] range")

    error = np.abs(true_hrs - predicted_hrs)
    percentage_error = 100*(error/true_hrs)
    candidates = [elem for elem in percentage_error if elem <= percentage]
    return 100*len(candidates)/len(percentage_error)

