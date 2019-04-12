"""
Created on Fri Nov 2 15:22:35 2018

@author: azarketa

This file is intended to act as a side-module to scripts that aim at processing TDMS files. It includes the necessary custom enum classes
to provide input information to the classes that process the TDMS files that store MGEP's wind tunnel data.
"""

######################################################################################################################
########################################################PACKAGES######################################################
######################################################################################################################

# enum package necessary to declare custom enumerator classes.
import enum

######################################################################################################################
########################################################CLASSES#######################################################
######################################################################################################################

# Public class tdmsFilereadMode.
class tdmsFilereadMode(enum.Enum):
    """This class provides options for allowing different reading modes of TDMS files."""

    standard = 1
    means = 2
    means_only = 3
    projected = 4
    projected_means_only = 5
    wake_rake = 6

# Public class deviceAxis.
class deviceAxis(enum.Enum):
    """This class provides options for defining axis information about devices."""

    x = 1
    y = 2
    z = 3
    angle = 4

# Public class runsPerParamMode.
class runsPerParamMode(enum.Enum):
    """This class provides options for setting the information pertaining the group to be added when performing a
    classification of a TDMS file data according to the positional information."""

    kistler = 1
    wake_rake = 2

######################################################################################################################
#####################################################END OF FILE######################################################
######################################################################################################################
