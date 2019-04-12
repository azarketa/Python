"""
Created on Fri Apr 12 01:35:10 2019

@author: azarketa
This file is intended to act as a side-module to scripts that aim at processing TDMS files. It includes the mathematical tools that are
necessary when processing TDMS files that store MGEP's wind tunnel data.
"""

######################################################################################################################
########################################################PACKAGES######################################################
######################################################################################################################

# numpy is intended to perform numeric calculations.
import numpy as np

######################################################################################################################
######################################################FUNCTIONS#######################################################
######################################################################################################################

# get_cdf() function.
def get_cdf(data):
    '''Obtains the normalized cummulative distribution function of the input.
    
    - **parameters**, **return**, **return types**::

        :param data: data from which to extract the cdf.         
        :return: list containing the cdf, array containing the linear space on the interval [0, 1] with a number of points equal
        to the length of the input argument 'data'.
        :rtype: list, array
    '''
    
    # Return statement computing the normalized cdf (by considering the minimum and maximum values of data through the 'sorted()'
    # function), and the linear space on the interval [0, 1].
    return [(b - sorted(data)[0])/(sorted(data)[-1] - sorted(data)[0]) for b in sorted(data)], np.linspace(0, 1, len(data))