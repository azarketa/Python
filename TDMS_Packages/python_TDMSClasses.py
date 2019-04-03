# PACKAGES
######################################################################################################################
######################################################################################################################
######################################################################################################################

# numpy and scipy are intended to perform numeric calculations.
import numpy as np
import scipy as sc
import scipy.integrate

# The entities classmethod and staticmethod are imported from the builtins package as function or class decorators.
# The classmethod decorator provides an akin functionality to the C# overloading of a function.
# The staticmethod decorator provides an akin functionality to the C# static method, which is to say that the object
# containing the static method need not be instantiated in order to use that method; what's more, a static method does
# not provide access to the parent object nor to other internals.
from builtins import staticmethod

# nptdms is a package that provides basic functionalities for reading a TDMS file (https://pypi.org/project/npTDMS/).
import nptdms

# TDMS_Packages.python_TDMSEnums is a custom package containing customized Enum objects used on the workflow herein.
from TDMS_Packages.python_TDMSEnums import *

# The library dataclasses implements the functionality of the dataclass decorator, which is used to create C#-like
# structs.
from dataclasses import *

# The library attr implements the functionality of attribute decorators
from attr import *

# The os library is intended to obtain information and perform operations with operative-system-related objects.
import os

# Serialization for saving files in binary format.
import pickle


######################################################################################################################
######################################################################################################################
######################################################################################################################

# CLASS HIERARCHY
# The following hierarchy intends to reproduce the structure of a generic TDMS file. A generic TDMS file structure
# consists of a root entity and an entity containing a number of groups, with each group containing a number of
# channels:
#
# ROOT
#  |
#  |
#  --GROUPS
#      |
#      |
#      --GROUP_1
#          |
#          |
#          --CHANNEL_1
#          --CHANNEL_2
#               .
#               .
#               .
#          --CHANNEL_N
#      --GROUP_2
#          |
#          |
#          --CHANNEL_1
#          --CHANNEL_2
#               .
#               .
#               .
#          --CHANNEL_M
#          .
#          .
#          .
#      --GROUP_L
#          |
#          |
#          --CHANNEL_1
#          --CHANNEL_2
#               .
#               .
#               .
#          --CHANNEL_K
#
# The structure of a TDMS file is reproduced on a class-based fashion; the root entity (tdmsFileRoot) is the only 
# 'public' class (on a C# terminology) accessible by external code. The other classes (groups-containing, 
# __tdmsFileGroups__; group-representing, __tdmsFileGroup__; channels-containing, __TdmsGroupChannels__; channel-
# representing, __TdmsGroupChannel__), are private internal classes (on C# terminology). As Python does not own an
# analogous construction for defining such classes, they are defined by notation convention (the double underscores
# before and aft the class declarations are meant, by convention, to warn potential coders and/or users that those
# classes are private by definition).
######################################################################################################################
######################################################################################################################
######################################################################################################################

# Public class tdmsFileRoot.
class TdmsFileRoot:
    """Public class tdmsFileRoot representing the core entity of a TDMS file."""

    # __init__ (constructor) of the public class tdmsFileRoot.
    def __init__(self, file_path, tdms_file_read_mode, file_alias, ref_mag):
        """Initializes an instance of the tdmsFileRoot object.

        The tdmsFileRoot object is the basic object that contains the data of a TDMS file.        

        - **parameters**, **return**, **return types**::

        :param file_path: path of the file to be read.
        :param tdms_file_read_mode: custom enum object specifying the reading mode:
            -Standard: reads and damps the TDMS structure, as a whole, to runtime variables.
            -Means: reads and damps the TDMS structure to runtime variables, adding to the structure itself
            the average values of the read data.
            -MeansOnly: reads and damps the TDMS structure to runtime variables, considering solely the
            average values of the read data.
            -Projected: to be used when force measurements are included; it projects the forces according to
            the measurement angle to provide the loads on wind tunnel axes.
            -ProjectedMeansOnly: to be used when force measurements are included; it projects the force
            averages according to the measurement angle to provide the average loads on wind tunnel axes.
        :param file_alias: string by which the data pertaining to a specific file is recognized on further charts
        or tables.
        :param ref_mag: an instantiation of the dataclass RefMagnitudes for non-dimensionalizing purposes.
        :return: tdmsFileRoot object structured in accordance to the provided reading mode.
        :rtype: TdmsFileRoot
        
        """

        # Conditional structure selecting in accordance to reading mode.
        # Standard reading mode.
        if tdms_file_read_mode == tdmsFilereadMode.Standard:
            self.__standard_init__(file_path, file_alias, ref_mag)
        # Means reading mode.
        elif tdms_file_read_mode == tdmsFilereadMode.Means:
            self.__means_init__(file_path, file_alias, ref_mag)
        # MeansOnly reading mode.
        elif tdms_file_read_mode == tdmsFilereadMode.MeansOnly:
            self.__means_only_init__(file_path, file_alias, ref_mag)
        # Projected reading mode.
        elif tdms_file_read_mode == tdmsFilereadMode.Projected:
            self.__projected_init__(file_path, file_alias, ref_mag)
        # ProjectedMeansOnly reading mode.
        elif tdms_file_read_mode == tdmsFilereadMode.ProjectedMeansOnly:
            self.__projected_means_only_init__(file_path, file_alias, ref_mag)
        # WakeRake reading mode.
        elif tdms_file_read_mode == tdmsFilereadMode.WakeRake:
            self.__wake_rake_init__(file_path, file_alias, ref_mag)

    # standard reading mode initializer.
    def __standard_init__(self, file_path, file_alias, ref_mag):
        """Initializes an instance of the tdmsFileRoot object in "Standard" reading mode.

        The instantiation of a tdmsFileRoot object in "Standard" reading mode leads to the damping of the
        structure of the provided TDMS file, as a whole, to runtime variables.

        - **parameters**::
        
        :param file_path: path of the file to be read.
        :param file_alias: string by which the data pertaining to a specific file is recognized on further charts
        or tables.

        """

        # Sets the 'file_alias' parameter as the new attribute '__Alias__' of the tdmsFileRoot object being
        # instantiated.
        self.__setattr__('__alias__', file_alias)
        # Sets the 'file_path' parameter as the new attribute '__path__' of the tdmsFileroot object being instantiated.
        self.__setattr__('__path__', file_path)
        # Sets the 'refMag' parameter as the new attribute '__refMag__' of the tdmsFileRoot object being instantiated.
        self.__setattr__('__ref_mag__', ref_mag)
        # Calls the internal __setRootProperties__ function to set the root properties.
        self.__set_root_properties__(file_path)
        # Instantiates a new internal __tdmsFileGroups__ object that stores the group structure of the TDMS file, and
        # sets that group into the groups_original variable of the tdmsFileRoot object being instantiated.
        self.groups_original = __TdmsFileGroups__(self.__file__, self.__file__.groups(), generate_isolated_groups=False)
        self.__non_dimensionalize__(self.__file__, ref_mag)

    # means reading mode initializer.
    def __means_init__(self, file_path, file_alias, ref_mag):
        """Initializes an instance of the tdmsFileRoot object in "Means" reading mode.

        The instantiation of a tdmsFileRoot object in "Means" reading mode leads to the damping of the
        structure of the provided TDMS file, in addition to the averages of the read data, to runtime variables.

        - **parameters**::

        :param file_path: path of the file to be read.
        :param file_alias: string by which the data pertaining to a specific file is recognized on further charts
        or tables.

        """

        # Calls the internal __standardInit__ method to set the structure of the TDMS file, as a whole, in the runtime
        # variable.
        self.__standard_init__(file_path, file_alias, ref_mag)
        # Calls the internal __setMeansAttributes__ method to set the necessary structure for dealing with mean values
        # of data.
        # self.__setMeansAttributes__()
        # Calls the internal __setMeans__ for setting the channel properties and their
        # correspondent values into variables contained within the means_channels variable.
        self.__set_means__()

    # means_only reading mode initializer.
    def __means_only_init__(self, file_path, file_alias, ref_mag):
        """Initializes an instance of the tdmsFileRoot object in "Means" reading mode.

        The instantiation of a tdmsFileRoot object in "MeansOnly" reading mode leads to the damping of the structure
        of the provided TDMS file, considering solely the averaged values of the read data.

        - **parameters**::

        :param file_path: path of the file to be read.
        :param file_alias: string by which the data pertaining to a specific file is recognized on further charts
        or tables.

        """

        # Calls the internal __meansInit__ method to set the structure of the TDMS file, as a whole, in the runtime
        # variable.
        self.__means_init__(file_path, file_alias, ref_mag)
        # Deleting the groups_original variable, with groups_added being the one remaining.
        if hasattr(self, 'groups_original'):
            self.__setattr__('__groups_original__', self.__getattribute__('groups_original'))
            self.__delattr__('groups_original')

    # projected reading mode initializer.
    def __projected_init__(self, file_path, file_alias, ref_mag):
        """Initializes an instance of the tdmsFileRoot object in "Projected" reading mode.

        The instantiation of a tdmsFileRoot object in "Projected" reading mode is to be used when force measurements
        are included; it projects the forces according to the measurement angle to provide the loads on wind tunnel
        axes.

        - **parameters**::

        :param file_path: path of the file to be read.
        :param file_alias: string by which the data pertaining to a specific file is recognized on further charts
        or tables.

        """

        # This initializer undertakes a process similar to that of the standard method, which is why it relies on the
        # previously defined __standardInit__ method for performing the data damping.
        self.__standard_init__(file_path, file_alias, ref_mag)
        # Calls the internal __projectForces__ method for performing the projecting operation on the loads.
        self.__project_forces__(self.__file__, False)

    # projected_means_only reading mode initializer.
    def __projected_means_only_init__(self, file_path, file_alias, ref_mag):
        """Initializes an instance of the tdmsFileRoot object in "ProjectedMeansOnly" reading mode.

        The instantiation of a tdmsFileRoot object in "ProjectedMeansOnly" reading mode is to be used when force
        measurements are included; it projects the averaged forces according to the measurement angle to provide the
        loads on wind tunnel axes.

        - **parameters**::

        :param file_path: path of the file to be read.
        :param file_alias: string by which the data pertaining to a specific file is recognized on further charts
        or tables.

        """

        # As this initializer performs a task similar to that of the __meansOnlyInit__ method, it relies on it to
        # undertake the data loading part.
        self.__means_only_init__(file_path, file_alias, ref_mag)
        # Calls the internal __projectForces__ method for performing the projecting operation on the loads.
        self.__project_forces__(self.__file__, True)

    # wake_rake reading mode initializer.
    def __wake_rake_init__(self, file_path, file_alias, ref_mag):
        """Initializes an instance of the tdmsFileRoot object in "WakeRake" reading mode.

        The instantiation of a tdmsFileRoot object in "WakeRake" reading mode is to be used when wake rake surveys are
        undertaken; it builds the attribute structure corresponding to the wake rake surveying mode.

        -**parameters**::

        :param file_path: path of the file to be read.
        :param file_alias: string by which the data pertaining to a specific file is recognized on further charts
        or tables.

        """

        # As this initializer performs a task similar to that of the __meansOnlyInit__ method, it relies on it to
        # undertake the data loading part.
        self.__means_only_init__(file_path, file_alias, ref_mag)
        # Calls the internal __setWakeRake__ method for performing the attribute structuring task.
        self.__runs_counter__(runsPerParamMode.WakeRake)

    # Internal __set_root_properties__ method.
    def __set_root_properties__(self, file_path):
        """Programmatically declares runtime variables within a tdmsFileRoot object.

        The properties found at the root level of a TDMS file are programmatically declared when instantiating a
        tdmsFileRoot object. The declared name matches the one coming from the file, and the values are set as "get"
        properties of the declared variables.

        - **properties**::

        :param file_path: path of the file to be read.

        """

        # An auxiliary variable (file) is declared as an instantiation of the object TdmsFile coming from the nptdms
        # package. The purpose of this variable is to provide access to the entities underlying the TDMS file to be
        # read and to allow performing the subsequent allocation of properties, both at the root level and at the
        # group/channel levels.
        #
        # Notice that the instantiation of the TDMS file is done, exclusively, by this method;
        # i.e. all the initializers are meant to own, within their methods, either the implementation of this method
        # or a method that internally calls this one.
        #
        # Correspondingly, the variable "__file__" is meant to be destroyed at the finalization stage of the
        # initializers.
        self.__file__ = nptdms.TdmsFile(file_path)

        # At this level, the following "for" loop runs over the properties of the root entity of the TDMS file and
        # retrieves their names and values for the declaration and allocation of runtime variables.
        #
        # The string replacements are meant to avoid having non-allowed characters on the declared variables, such as
        # white spaces or special characters (!"·$%&/()=?¿\).
        for prop in list(self.__file__.object().properties.items()):
            self.__setattr__('__' + prop[0].replace(' ', '_').replace('(', '').replace(')', '') + '__', prop[1])

    # Internal __add_isolated_groups__ method.
    def __add_isolated_groups__(self, isolated_groups_name):
        """Programmatically declares an entity intended to store the groups coming from a TDMS file.

        The information found at the root level is extended by adding an entity intended to store the groups coming
        from the TDMS file being read. The declared name for this variable is a settable parameter of the method.

        - **parameters**, **return**, **return types**::

        :param isolated_groups_name: name of the programmatically declared variable.
        :return: __tdmsFileGroups__ object.
        :rtype: __TdmsFileGroups__

        """

        # Declaring an instantiation of the __tdmsFileGroups__ object and setting it into the "isolatedGroupsName"
        # variable.
        self.__setattr__(isolated_groups_name, __TdmsFileGroups__(self.__file__, None, generate_isolated_groups=True))

    # Internal __set_standard_data__ method.
    def __set_standard_data__(self):
        """Method that loads standard TDMS data in case the initialization has not done so."""

        if hasattr(self, '__groups_original__'):
            self.__setattr__('groups_original', self.__getattribute__('__groups_original__'))
            self.__delattr__('__groups_original__')
        else:
            # Conditional that assigns the attribute '__file__' to an instantiation of the nptmds.TdmsFile object in order
            # to accomplish the task of further attribute assignments.
            if hasattr(self, '__file__'):
                pass
            else:
                self.__setattr__('__file__', nptdms.TdmsFile(self.__path__))
            # Conditional that checks whether the "groups_original" attribute is already assigned to the object. If so,
            # the function does not perform further actions. Otherwise, the subsequent attribute sturcture is set.
            if hasattr(self, "groups_original"):
                return
            else:
                self.__setattr__("groups_original",
                                 __TdmsFileGroups__(self.__file__, self.__file__.groups(), generate_isolated_groups=False))

        # Calls the internal __nondimensionalize__ method with the correspondent refMag object (referential magnitudes)
        # and the nonDimensaionlizeGroups.original flag indicating that just the original group is to be
        # non-dimensionalized.
        self.__non_dimensionalize__(file_obj=self.__file__, ref_mag=self.__ref_mag__)

    # Internal __set_means_attributes__ method.
    def __set_mean_attributes__(self):
        """Sets the attributes necessary to deal with mean values of data.

        -**returns, return types**::

        :return: boolean indicating whether the "means_group" attribute was already set before the function call.
        :rtype: boolean.

        """

        # Conditional that assigns the attribute '__file__' to an instantiation of the nptmds.TdmsFile object in order
        # to accomplish the task of further attribute assignments.
        if hasattr(self, '__file__'):
            pass
        else:
            self.__setattr__('__file__', nptdms.TdmsFile(self.__path__))
        # Calls the internal __addIsolatedGroups__ method for setting an additional internal object __tdmsFileGroups__
        # into the groups_added variable of the tdmsFileRoot object being instantiated. The groups_added variable
        # reproduces the data structure of the groups_original variable, but with derived data (means, projections,
        # Fourier transformations...).
        if hasattr(self, "groups_added"):
            pass
        else:
            self.__add_isolated_groups__("groups_added")
        # Calls the internal __addIsolatedGroup__ method for setting an additional internal object __tdmsFileGroup__
        # into the means_group variable of the groups_added variable of the tdmsFileRoot object being instantiated. The
        # means_group variable is intended to store the averaged values of the TDMS file data. Instead of storing the
        # data of each group separately, averaged data of the channels of all groups are stored in channels
        # corresponding to each of the variables.
        if hasattr(self.groups_added, "means_group"):
            return False
        else:
            self.groups_added.__add_isolated_group__("means_group")
            # Calls the internal __addIsolatedChannels__ method for setting an additional internal object
            # __TdmsGroupChannels__ into the means_group variable; this class intends to store the channels that
            # correspond
            # to the data categories of each of the groups of the TDMS file.
            self.groups_added.means_group.__add_isolated_channels__("means_channels")

        # "True" statement return in case the structure was not set before the function call.
        return True

    # Internal __set_means__ method.
    def __set_means__(self):
        """Rearranges group properties and channel values for their allocation on the groups_added variable.

        When employing reading modes that allocate average values of data on variables, group and channel properties
        need to be rearranged so that all the information pertaining a group (both at its level and at its channels'
        level) is included on variables within the groups_added variable.
        As an illustrating example: the positional information of a measurement (X_pos, Y_pos, Z_pos and Angle_pos),
        which is a group level information, is the same for all the child channels of a particular group. If, say, the
        averaged channels' data is to be matched with the positional information, both the channel data and the
        positional information are required to be on a same hierarchical level. This level is provided by the
        groups_added variable, which stores positional information and averaged channel data on one-to-one corresponding
        lists.

        """

        # Calling the internal method __setMeansAttributes__() in case it is necessary to set the attribute structure,
        # and storing the return boolean in the boolean_flag variable.
        boolean_flag = self.__set_mean_attributes__()

        if boolean_flag:
            # A groups level "for" loop that runs over the groups on the TDMS file (excluding "Initial drift" and "Final
            # drift" groups that are not meant to provide data-analysis information).
            for group in [group for group in self.__file__.groups() if group not in ["Initial drift", "Final drift"]]:
                # A group level "for" loop that runs over the properties on that group and looks for the properties
                # containing the positional information (X_pos, Y_pos, Z_pos and Angle_pos).
                for pos_prop in [prop for prop in dir(self.groups_original.__getattribute__(group)) if
                                 (prop.__repr__() in ["'X_pos'", "'Y_pos'", "'Z_pos'", "'Angle_pos'"])]:
                    # Conditional that checks whether the groups_added variable already contains a child variable
                    # containing the positional information found. If so, that variable's list is appended with the
                    # newly found positional information piece; otherwise, a new variable is declared with its name
                    # matching that of the positional information, and its value being a list with a single element
                    # matching the value of the positional information.
                    if hasattr(self.groups_added.means_group.means_channels, pos_prop):
                        self.groups_added.means_group.means_channels.__getattribute__(pos_prop).append(
                            self.groups_original.__getattribute__(group).__getattribute__(pos_prop))
                    else:
                        self.groups_added.means_group.means_channels.__setattr__(pos_prop, [
                            self.groups_original.__getattribute__(group).__getattribute__(pos_prop)])
                # A channel level "for" loop that runs over the channels of the group.
                # for channel in self.__file__.group_channels(group):
                for channel in [not_dunder_property for not_dunder_property in dir(
                        self.groups_original.__getattribute__(group).channels) if
                                ('__' not in not_dunder_property.__repr__())]:
                    # Conditional that checks whether the groups_added variable already contains a child variable
                    # containing the channel information found. If so, that variable's list is appended with the newly
                    # found channel information piece; otherwise, a new variable is declared with its name matching that
                    # of the channel information, and its value being a list with a single element matching the value of
                    # the averaged channel information.
                    if hasattr(self.groups_added.means_group.means_channels, channel.replace(" ", "_")):
                        self.groups_added.means_group.means_channels.__getattribute__(channel.replace(" ", "_")).append(
                            np.average(
                                self.groups_original.__getattribute__(group).channels.__getattribute__(channel).data))
                    else:
                        self.groups_added.means_group.means_channels.__setattr__(channel.replace(" ", "_"), [np.average(
                            self.groups_original.__getattribute__(group).channels.__getattribute__(channel).data)])

    # TODO: DOCUMENT METHOD CORRECTLY.
    # Internal __set_wake_rake__ method.
    def __set_wake_rake__(self, eval_string, runs_info, wake_rake_port_ordering=tuple((np.arange(1, 19))),
                          wake_rake_cross_axis=deviceAxis.Y):
        """ """

        actual_param = eval_string.split('__getattribute__')[-1].replace('(', '').replace(')', '').replace("'", "")
        actual_param_index = [runInfo[0] for runInfo in runs_info if runInfo[2] == actual_param][0]
        parent_eval_string = ".".join(eval_string.split('.')[:-1])
        parent_param = parent_eval_string.split('__getattribute__')[-1].replace('(', '').replace(')', '').replace(
            "'", "")
        parent_param_index = [runInfo[1] for runInfo in runs_info if runInfo[2] == parent_param][0]
        wake_rake_value_list = []

        for i in range(0, len(
                self.groups_added.means_group.means_channels.__getattribute__(parent_param + '_pos')[
                ::parent_param_index])):
            for _ in np.arange(parent_param_index * i, parent_param_index * (i + 1)):
                wake_rake_pos = []
                wake_rake_values = []
                for actualIndex in range(0, actual_param_index):
                    specific_index = np.arange(parent_param_index * i, parent_param_index * (i + 1))[actualIndex]
                    for port in wake_rake_port_ordering:
                        wake_rake_pos.append(
                            eval(eval_string + ".__getattribute__(" + (actual_param + "_pos").__repr__() + ")")[
                                actualIndex] + (-wake_rake_port_ordering.index(port) + 9) * 2.505 / (
                                    1000 * self.__ref_mag__.length) - 1.25 / (1000 * self.__ref_mag__.length))
                        if wake_rake_port_ordering.index(port) == 12:
                            if 0.5 * (self.groups_added.means_group.means_channels.__getattribute__(
                                    "Port_" + str(wake_rake_port_ordering[11]) + "_signal")[specific_index] +
                                      self.groups_added.means_group.means_channels.__getattribute__(
                                          "Port_" + str(wake_rake_port_ordering[13]) + "_signal")[
                                          specific_index]) > 1.0:
                                wake_rake_values.append(1.0)
                            else:
                                wake_rake_values.append(0.5 * (
                                        self.groups_added.means_group.means_channels.__getattribute__(
                                            "Port_" + str(wake_rake_port_ordering[11]) + "_signal")[specific_index] +
                                        self.groups_added.means_group.means_channels.__getattribute__(
                                            "Port_" + str(wake_rake_port_ordering[13]) + "_signal")[specific_index]))
                        else:
                            if self.groups_added.means_group.means_channels.__getattribute__(
                                    "Port_" + str(port) + "_signal")[specific_index] > 1.0:
                                wake_rake_values.append(1.0)
                            else:
                                wake_rake_values.append(self.groups_added.means_group.means_channels.__getattribute__(
                                    "Port_" + str(port) + "_signal")[specific_index])
                    eval(eval_string + ".__setattr__(" + (
                            wake_rake_cross_axis.name.upper() + "_rake_pos").__repr__() + ", wake_rake_pos)")
            wake_rake_value_list.append(wake_rake_values)
        eval(eval_string + ".__setattr__(" + (
                wake_rake_cross_axis.name.upper() + "_rake_values").__repr__() + ", wake_rake_value_list)")

        integrate_integrate_x = [int_int_x for int_int_x in eval(
            eval_string + ".__getattribute__(" + (wake_rake_cross_axis.name.upper() + "_rake_pos").__repr__() + ")")]
        integration_pre_process = [eval_string + ".__getattribute__(" + (
                wake_rake_cross_axis.name.upper() + "_rake_values").__repr__() + ")[" + str(i) + "]" for i in list(
            np.arange(0, len(eval(eval_string + ".__getattribute__(" + (
                    wake_rake_cross_axis.name.upper() + "_rake_values").__repr__() + ")"))))]
        cd = []
        integrates_y = []

        parent_for_average_list = eval(
            parent_eval_string +
            '.__getattribute__(' +
            (parent_param + '_pos').__repr__() +
            ')')

        parent_for_average_indices = [parent_for_average_list.index(x) for x in
                                      filter(lambda x: np.abs(x) <= 5, parent_for_average_list)]

        actual_for_average_list = eval(
            eval_string +
            '.__getattribute__(' +
            (wake_rake_cross_axis.name.upper() + "_rake_values").__repr__() +
            ')')

        actual_average_list = []

        for parent_for_average_index in parent_for_average_indices:
            actual_average_list.append(np.average(actual_for_average_list[parent_for_average_index][:10]))
            actual_average_list.append(np.average(actual_for_average_list[parent_for_average_index][-10:]))

        actual_average = np.average(actual_average_list)

        for integrate_pre_process in integration_pre_process:
            # post-process=0-5Mean
            int_list = [
                np.sqrt(int_y) * (1.0 - np.sqrt(int_y)) -
                np.sqrt(actual_average) * (1.0 - np.sqrt(actual_average))
                for int_y in eval(integrate_pre_process)]
            x, y = zip(*sorted(zip(integrate_integrate_x, int_list)))
            integrate_y = [int_y for int_y in y]

            def get_cdf(data):
              return [(b - sorted(data)[0])/(sorted(data)[-1] - sorted(data)[0]) for b in sorted(data)], np.linspace(0, 1, len(data))

            cdf = get_cdf(integrate_y)
            ave = np.average([(list(np.diff(cdf[0]))[i]**2 + list(np.diff(cdf[1]))[i]**2)**0.5 for i in range(0, len(np.diff(cdf[0])))])
            modulus_list = [(list(np.diff(cdf[0]))[i]**2 + list(np.diff(cdf[1]))[i]**2)**0.5 for i in range(0, len(np.diff(cdf[0])))]
            cdf_lower_limit_index = list(np.where(modulus_list[2:] > ave)[0])[0] + 2
            cdf_lower_limit = cdf[0][cdf_lower_limit_index]
            for i in range(0, len(integrate_y)):
                if (integrate_y[i] - np.min(integrate_y))/(np.max(integrate_y) - np.min(integrate_y)) < cdf_lower_limit:
                    integrate_y[i] = 0

            integrates_y.append(integrate_y)
            cd.append(2.0 * sc.integrate.simps(integrate_y, x=x))

        eval(eval_string + ".__setattr__(" + (wake_rake_cross_axis.name.upper() + "_rake_pos").__repr__() + ", x)")
        eval(eval_string + ".__setattr__(" + (
                wake_rake_cross_axis.name.upper() + "_rake_int_values").__repr__() + ", integrates_y)")
        eval(parent_eval_string + ".__setattr__(" + (parent_param + "_cd").__repr__() + ", cd)")

    # TODO: DOCUMENT METHOD CORRECTLY.
    # Internal __runs_counter__ method.
    def __runs_counter__(self, runs_per_param_mode, device_axis=deviceAxis.Y):
        """ """

        if not hasattr(self, "groups_added"):
            self.__set_means__()
        elif not hasattr(self.groups_added, "means_group"):
            self.__set_means__()

        poss_info = [axis.name for axis in deviceAxis]
        runs_info = []

        for pos_info in poss_info:
            runs_info.append(
                # True value of in-line conditional for a 4-member tuple: (Member 1, Member 2, Member 3, Member 4).
                (
                    # Member 1
                    len(self.groups_added.means_group.means_channels.__getattribute__(pos_info + "_pos")) //
                    min(
                        [self.groups_added.means_group.means_channels.__getattribute__(pos_info + "_pos").count(
                            run_value)
                            for run_value in
                            self.groups_added.means_group.means_channels.__getattribute__(pos_info + "_pos")]
                    ),
                    # Member 2
                    [i != self.groups_added.means_group.means_channels.__getattribute__(pos_info + "_pos")[0]
                     for i in
                     self.groups_added.means_group.means_channels.__getattribute__(pos_info + "_pos")].index(True),
                    # Member 3
                    pos_info,
                    # Member 4
                    "consecutive runs per" + pos_info
                )
                # Condition.
                if
                min(
                    [self.groups_added.means_group.means_channels.__getattribute__(pos_info + "_pos").count(run_value)
                     for run_value in
                     self.groups_added.means_group.means_channels.__getattribute__(pos_info + "_pos")]
                ) !=
                len(
                    self.groups_added.means_group.means_channels.__getattribute__(pos_info + "_pos"))
                # False value of in-line conditional for a 4-member tuple: (Member 1, Member 2, Member 3, Member 4).
                else
                # Member 1.
                (np.NaN if (pos_info != device_axis.name) else 1,
                 # Member 2.
                 np.NaN if (pos_info != device_axis.name) else 1,
                 # Member 3.
                 "" if (pos_info != device_axis.name) else pos_info,
                 # Member 4.
                 "" if (pos_info != device_axis.name) else "consecutive runs per" + pos_info)
            )

        sorted_runs_info = sorted([run_info for run_info in runs_info if not np.isnan(run_info[1])],
                                  key=lambda x: x[1], reverse=True)

        to_assign_param_names = [run_info[2] for run_info in sorted_runs_info]
        to_assign_con_values = [run_info[1] for run_info in sorted_runs_info]
        to_assign_values = [runInfo[0] for runInfo in sorted_runs_info]

        added_group = runs_per_param_mode.name[0].lower() + runs_per_param_mode.name[1:] + "Group"
        runs_per_param_group = __RunsPerParam__(added_group)
        self.groups_added.__setattr__(added_group, runs_per_param_group)
        self.groups_added.__getattribute__(added_group).__delattr__(added_group)

        print(self.__path__.split("/")[-2])
        if (np.nanprod([runs_per_param[0] for runs_per_param in sorted_runs_info]) ==
                len(self.groups_added.means_group.means_channels.Z_pos)):
            print(
                "Sequential runs: check correct --> Ordering: ("
                + ",".join([str(sorted_run_info[2]) for sorted_run_info in sorted_runs_info]) +
                ") --> Numbering: (" + ",".join([str(sorted_run_info[0]) for sorted_run_info in sorted_runs_info]) +
                ")"
            )
            for param_name in to_assign_param_names:
                if to_assign_param_names.index(param_name) == 0:
                    runs_per_param = __RunsPerParam__(
                        param_name + "_pos",
                        self.groups_added.means_group.means_channels.__getattribute__(param_name + "_pos")
                        [::to_assign_con_values[to_assign_param_names.index(param_name)]]
                    )
                    self.groups_added.__getattribute__(added_group).__setattr__(param_name, runs_per_param)
                else:
                    eval_string = "self.groups_added.__getattribute__(" + added_group.__repr__() + ")."
                    for previous_param_name in to_assign_param_names[:to_assign_param_names.index(param_name)]:
                        eval_string += "__getattribute__(" + previous_param_name.__repr__() + ")."
                    eval(
                        eval_string +
                        "__setattr__(" +
                        param_name.__repr__() +
                        ", __RunsPerParam__(" +
                        (param_name + "_pos").__repr__() +
                        ", self.groups_added.means_group.means_channels.__getattribute__(" +
                        (param_name + "_pos").__repr__() +
                        ")[:to_assign_values[to_assign_param_names.index(" +
                        param_name.__repr__() + ")]:to_assign_con_values[to_assign_param_names.index(" +
                        param_name.__repr__() + ")]]))"
                    )
                    if to_assign_param_names.index(param_name) == len(to_assign_param_names) - 1:
                        eval_string = \
                            eval_string.split("__setattr__")[0] + \
                            "__getattribute__(" + param_name.__repr__() + ")"
                        if runs_per_param_mode == runs_per_param_mode.WakeRake:
                            self.__set_wake_rake__(eval_string, runs_info)
        else:
            print("Check not correct. File not liable to be ordered coherently by levelled runs.")

    # Internal __project_forces__ method.
    def __project_forces__(self, file_obj, means_only=False):
        """Projects force-related data on the TDMS file into wind tunnel axes.

        The projection of force-related data on the TDMS file into wind tunnel axes is done according to the positional
        information on each of the measurement groups. Mind that the load balance is assumed to be aligned with wind
        tunnel axes for a 0º configuration, with x- and y-axis being flow-wise and traverse directions, respectively.
        Clockwise angles are assumed positive.

        - **parameters**::

        :param file_obj: a nptdms.TdmsFile object containing non_drift_group/channel information to be projected.
        :param means_only: boolean flag signaling whether averaged values are meant to be considered exclusively.

        """

        # A groups level "for" loop that runs over the groups on the TDMS file (excluding "Initial drift" and "Final
        # drift" groups that are not meant to provide data-analysis information).
        for non_drift_group in [group for group in file_obj.groups() if group not in ["Initial drift", "Final drift"]]:
            # A group_obj looping variable is declared that instantiates a group object of the fileObj variable.
            group_obj = file_obj.object(non_drift_group)
            # A channel_list looping variable is declared that contains a list of the channels on the group
            # "non_drift_group" of the file_Obj variable.
            channel_list = [channel.channel for channel in file_obj.group_channels(non_drift_group)]
            # Conditional for checking whether the non_drift_group "non_drift_group" contains, in its channels,
            # drift-compensated load measurements.
            if all(elem in channel_list for elem in ["Kistler signal time", "Corrected Kistler Fx signal"]):
                # Getting the angle at which the measurement was taken from the positional information, and converting
                # it to radian units.
                angle = group_obj.properties.get("Angle pos") * np.pi / 180
                # Computing drag and lift projections according to:
                # d = Fx·cos(a) - Fy·sin(a)
                # length = Fx·sin(a) + Fy·cos(a)
                drag = np.cos(angle) * file_obj.object(non_drift_group, "Corrected Kistler Fx signal").data - \
                    np.sin(angle) * file_obj.object(non_drift_group, "Corrected Kistler Fy signal").data
                lift = np.sin(angle) * file_obj.object(non_drift_group, "Corrected Kistler Fx signal").data + \
                    np.cos(angle) * file_obj.object(non_drift_group, "Corrected Kistler Fy signal").data
                # Conditional that checks on the boolean flag; if false (temporal information is meant to be damped),
                # additional attributes of "drag" and "lift" are declared at the channel level of the groups_original
                # variable, and their values set to the computed values of drag and lift.
                if not means_only:
                    self.groups_original.__getattribute__(non_drift_group).__setattr__("drag", drag)
                    self.groups_original.__getattribute__(non_drift_group).__setattr__("lift", lift)
                # Conditional that checks on the boolean flag; if true (averaged information is meant to be damped),
                # a further conditional is encountered which performs the damping operation of the averaged information
                # on the groups_added variable.
                if means_only:
                    # Conditional that checks whether the groups_added variable already contains a child variable
                    # containing the projected information. If so, that variable's list is appended with the newly
                    # found projected information piece; otherwise, a new variable is declared with its name matching
                    # that of the projected information, and its value being a list with a single element matching the
                    # value of the projected information.
                    if hasattr(self.groups_added.means_group.means_channels, "drag"):
                        self.groups_added.means_group.means_channels.drag.append(np.average(drag))
                        self.groups_added.means_group.means_channels.lift.append(np.average(lift))
                    else:
                        self.groups_added.means_group.means_channels.__setattr__("drag", [np.average(drag)])
                        self.groups_added.means_group.means_channels.__setattr__("lift", [np.average(lift)])

    # Internal __non_dimensionalize__ method.
    def __non_dimensionalize__(self, file_obj, ref_mag):
        """Non-dimensionalizes data coming from the TDMS file.

        - **parameters**::

        :param file_obj: a nptdms.TdmsFile object containing group/channel information to be projected.
        :param ref_mag: an instantiation of the dataclass RefMagnitudes for non-dimensionalizing purposes.

        """

        # Conditional for checking whether the non-dimensionalization of pressure is necessary. This code block is
        # intended to be hit when loading the TDMS file the first time.
        if "Atmospheric_pressure_Pa" in dir(self):
            ref_mag.p = self.__getattribute__("Atmospheric_pressure_Pa")

        # Conditional for checking whether a non-dimensionalization is needed on the group of data coming from the
        # original TDMS file. If it happens that the data has already been non-dimensionalized, this code block is
        # skipped (i.e. when calculating means from non-dimensionalized original data).
        for item in file_obj.objects.items():
            item_name_list = item[0].replace("/", "").replace("''", '"').replace("'", "").split('"')
            if (len(item_name_list) == 1) and (item_name_list[0] != ""):
                if hasattr(self.groups_original.__getattribute__(
                        item_name_list[0].replace(" ", "_")).__getattribute__(
                        "channels"), "Thermocouple_signal"):
                    ref_mag.temp = np.average(
                        self.groups_original.__getattribute__(
                            item_name_list[0].replace(" ", "_")).__getattribute__(
                            "channels").__getattribute__("Thermocouple_signal").data) + 273
                if hasattr(self.groups_original.__getattribute__(
                        item_name_list[0].replace(" ", "_")).__getattribute__(
                        "channels"), "Fix_probe_signal"):
                    ref_mag.u = np.average(
                        self.groups_original.__getattribute__(
                            item_name_list[0].replace(" ", "_")).__getattribute__(
                            "channels").__getattribute__(
                            "Fix_probe_signal").data)
                self.groups_original.__getattribute__(item_name_list[0].replace(" ", "_")).X_pos /=\
                    (1000 * ref_mag.length)
                self.groups_original.__getattribute__(item_name_list[0].replace(" ", "_")).Y_pos /=\
                    (1000 * ref_mag.length)
                self.groups_original.__getattribute__(item_name_list[0].replace(" ", "_")).Z_pos /=\
                    (1000 * ref_mag.length)
            if "unit_string" in item[1].properties:
                if item[1].property("unit_string") == "Pa":
                    self.groups_original.__getattribute__(
                        item_name_list[0].replace(" ", "_")).channels.__getattribute__(
                        item_name_list[1].replace(" ", "_")).data /= ref_mag.q
                elif item[1].property("unit_string") == "N":
                    self.groups_original.__getattribute__(
                        item_name_list[0].replace(" ", "_")).channels.__getattribute__(
                        item_name_list[1].replace(" ", "_")).data /= (ref_mag.q * ref_mag.length * ref_mag.span)
                elif item[1].property("unit_string") == "Nm":
                    self.groups_original.__getattribute__(
                        item_name_list[0].replace(" ", "_")).channels.__getattribute__(
                        item_name_list[1].replace(" ", "_")).data /=\
                        (ref_mag.q * (ref_mag.length ** 2.0) * ref_mag.span)
                elif item[1].property("unit_string") == "m/s":
                    self.groups_original.__getattribute__(
                        item_name_list[0].replace(" ", "_")).channels.__getattribute__(
                        item_name_list[1].replace(" ", "_")).data /= ref_mag.u
                elif item[1].property("unit_string") == "ºC":
                    self.groups_original.__getattribute__(
                        item_name_list[0].replace(" ", "_")).channels.__getattribute__(
                        item_name_list[1].replace(" ", "_")).data += 273
                    self.groups_original.__getattribute__(
                        item_name_list[0].replace(" ", "_")).channels.__getattribute__(
                        item_name_list[1].replace(" ", "_")).data /= ref_mag.temp
                else:
                    pass

    # TODO: Implement corrections to variables according to bibliography.
    # Internal __corrector__ method.
    def __corrector__(self, file_obj, ref_mag):
        raise NotImplementedError
        #        ######################################################################################################
        #        ######################################################################################################
        #        # Correction parameters coming from "Wind Tunnel Testing Airfoils at Low Reynolds Numbers" (Selig).
        #        # Maybe it should go on a separate method.
        #        # Solid blockage, sb = (k1 · Mv ) / (Ats)^(3/2)
        #        k1 = 0.52
        #        Mv = 0.00291296
        #        Ats = 0.9*0.75
        #        sb = (k1*Mv)/(Ats**1.5)
        #        # Wake blockage, wb = (c/(2hts))*ndDrag
        #        hts = 0.9
        #        wb = (c/(2*hts))*ndDrag
        #        # Blockage, bl = sb + wb
        #        bl = sb + wb
        #        # Sigma, sigma = (pi/48)*(c/hts)**2
        #        sigma = (np.pi/48)*(c/hts)**2
        #        # ndLift corrected, ndLift = ndLiftUnc * (1 - sigma)/(1 - blockage)**2
        #        # ndDrag corrected, ndDrag = ndDragUnc * (1 - solid blockage)/(1 - blockage)**2
        #        # ndPitch corrected, ndMz = (ndMzUnc - ndLiftUnc*sigma*(1 - sigma)/4)/(1 - blockage)**2
        #        ndMzCorr = (ndMz - ndLift*sigma*(1 - sigma)/4)/((1 - bl)**2)
        #        ndLiftCorr = ndLift*(1 - sigma)/((1 - bl)**2)
        #        ndDragCorr = ndDrag*(1 - sb)/((1 - bl)**2)

# Internal class __TdmsFileGroups__, which contains the group objects of the tdmsFileRoot object.
class __TdmsFileGroups__:
    """Internal class __tdmsFileGroups__ representing the entity of groups of a TDMS file."""

    # __init__ (constructor) of the internal class __tdmsFileGroups__.
    def __init__(self, file_obj, group_objs=None, generate_isolated_groups=False):
        """Initializes an instance of the __tdmsFileGroups__ object.

        The __tdmsFileGroups__ object is the entity object that contains the groups of a TDMS file.

        - **parameters**, **return**, **return types**::

        :param file_obj: an instantiation of the nptdms.TdmsFile object.
        :param group_objs: a list containing strings of group names to be instantiated as nptdms.GroupObject objects

        :param generate_isolated_groups: boolean flag signaling whether an isolated group object is to be instantiated.
        :return: __tdmsFileGroups__ object structured in accordance to the groups contained on the TDMS file.
        :rtype: __tdmsFileGroups__.
        """

        # Conditional for checking whether an isolated groups is to be instantiated. The code hits this block if the
        # boolean flag is false.
        if group_objs is None:
            group_objs = []
        if not generate_isolated_groups:
            # Declaring a group_index variable for running over the group objects passed as parameters.
            group_index = 0
            # For loop running over the group objects. For each group object, a __tdmsFileGroup__ object is
            # instantiated and set as an attribute on the parent __tdmsFileGroups__ object.
            for groupObj in group_objs:
                self.__setattr__(groupObj.replace(" ", "_"),
                                 __TdmsFileGroup__(nptdms.GroupObject(groupObj), group_index, file_obj, False))
                group_index += 1
        else:
            pass

    # Internal __addIsolatedGroup__ method.
    def __add_isolated_group__(self, isolated_group_name):
        """Instantiates an isolated __TdmsFileGroup__ object.

        The method is a side call to the nptdms library built-in method nptdms.GroupObject(), which is the proper
        instantiator of a group object. This is done so to allow declaring isolated group instantiations from the
        parent group __TdmsFileGroups__. Declaring isolated groups turns necessary when running the root initializer
        with certain reading modes.

        - **parameters**

        :param isolated_group_name: name given to the isolated group to be instantiated.

        """

        # Call to the static internal method __TdmsFileGroup__.__generateIsolatedGroup__() which is, itself, a side
        # call to the proper group instantiator nptdms.GroupObject().
        self.__setattr__(isolated_group_name, __TdmsFileGroup__.__generate_isolated_group__(isolated_group_name))


# Internal class __TdmsFileGroup__, which contains individual groups that populate a __tdmsFileGroups__ object.
class __TdmsFileGroup__:
    """Internal class __TdmsFileGroup__ representing a group entity of a TDMS file."""

    # __init__ (constructor) of the internal class __tdmsFileGroup__.
    def __init__(self, group_obj, group_index, file_obj, generate_isolated_group=False):
        """Initializes an instance of the __TdmsFileGroup__ object.

        The __TdmsFileGroup__ object is the entity object representing a group of a TDMS file.

        - **parameters**, **return**, **return types**::

        :param group_obj: a string literal containing the declared name of the group object to be instantiated.
        :param group_index: index for identifying the specific group within the groups of the fileObj variable and thus
        retrieving the properties of the group for its proper instantiation.
        :param file_obj: an instantiation of the nptdms.TdmsFile object representing a TDMS file.
        :param generate_isolated_group: boolean flag signaling whether an isolated group is to be instantiated.
        :return: __TdmsFileGroup__ object structured in accordance to the properties and channels contained within
        that group.
        :rtype: __TdmsFileGroup__.

        """

        # Conditional for checking whether an isolated group is to be instantiated. The code hits this block if the
        # boolean flag is false.
        if not generate_isolated_group:
            # Conditional for allowing the instantiation of nptdms.GroupObject objects outside the scope of a specific
            # nptmds.TdmsFile object.
            if file_obj is None:
                pass
            else:
                # Declaring the attributes of "group" and "group_name".
                self.group = group_obj
                self.group_name = group_obj.group
                # Declaring the list of properties "prop_list".
                prop_list = list((file_obj.object(file_obj.groups()[group_index])).properties.items())
                # A "for" loop running over the properties list and setting the correspondent attributes.
                for prop in prop_list:
                    self.__setattr__(prop[0].replace(" ", "_"), prop[1])
                # Declaring the attribute "channels" as an instantiation of the __TdmsGroupChannels__ object, which
                # represents the collection of channels contained within a group.
                self.channels = __TdmsGroupChannels__(self.group, file_obj.group_channels(self.group_name),
                                                      generate_isolated_channels=False)
        else:
            pass

    # Static internal __generateIsolatedGroup__ method.
    @staticmethod
    def __generate_isolated_group__(isolated_group_name):
        """Static internal method __generateIsolatedGroup__ for generating isolated groups.

        The groups are instantiations of the nptdms.GroupObject objects.

        - **parameters**, **returns**, **return types**

        :param isolated_group_name: string literal providing the declaration name of the instantiated object.
        :return: __tdmsFileGroup__ object reflecting a nptdms.GroupObject object.
        :rtype: __TdmsFileGroup__

        """

        # Return statement instantiating a __tdmsFileGroup__ object.
        return __TdmsFileGroup__(nptdms.GroupObject(isolated_group_name), None, None, False)

        # Internal __addIsolatedChannels__ method.

    def __add_isolated_channels__(self, isolated_channels_name):
        """Programmatically declares an entity intended to store the channels coming from a TDMS group.

        The information found at the group level is extended by adding an entity intended to store the channels coming
        from the TDMS group being read. The declared name for this variable is a settable parameter of the method.

        - **parameters**, **return**, **return types**::

        :param isolated_channels_name: name of the programmatically declared variable.
        :return: __tdmsFileChannels__ object.
        :rtype: __tdmsFileChannels__

        """

        # Declaring an instantiation of the __tdmsFileChannels__ object and setting it into the "isolatedChannelsName"
        # variable.
        self.__setattr__(isolated_channels_name, __TdmsGroupChannels__(None, None, generate_isolated_channels=True))


# Internal class __TdmsGroupChannels__, which contains the channel objects of a __tdmsFileGroup__ object.
class __TdmsGroupChannels__:
    """Internal class __TdmsGroupChannels__ representing the entity of channels of a TDMS group."""

    # __init__ (constructor) of the internal class __TdmsGroupChannels__.
    def __init__(self, group_obj, channel_objs, generate_isolated_channels=False):
        """Initializes an instance of the __TdmsGroupChannels__ object.

        The __TdmsGroupChannels__ object is the entity containing the channels of a TDMS group.

        - **parameters**, **return**, **return types**::

        :param group_obj: an instantiation of the nptdms.GroupObject object representing a TDMS group for which the
        channels are to be determined.
        :param channel_objs: a list of string literals containing the names of the channels to be declared as
        instantiations of the __TdmsGroupChannel__ object.
        :param generate_isolated_channels: boolean flag signaling whether an isolated channel is to be instantiated.
        :return: __TdmsGroupChannels__ object structured in accordance to the channels contained on a TDMS group.
        :rtype: __TdmsGroupChannels__.

        """

        # Conditional for checking whether an isolated channel is to be instantiated. This block of code is hit solely
        # when the boolean flag is false.
        if not generate_isolated_channels:
            # Declaring the channel_obj looping variable for running over the channel names and declaring them as
            # instantiations of the __TdmsGroupChannel__ object.
            for channel_obj in channel_objs:
                self.__setattr__(channel_obj.channel.replace(" ", "_"),
                                 __TdmsGroupChannel__(channel_obj.channel.replace(" ", "_"),
                                                      nptdms.ChannelObject(channel_obj, group_obj, channel_obj.data)))
        else:
            pass

    # Internal __add_isolated_channel__ method.
    def __add_isolated_channel__(self, isolated_channel_name):
        """Instantiates an isolated __TdmsGroupChannel__ object.

        The method is a call to the initializer of the __TdmsGroupChannel__ object.

        - **parameters**

        :param isolated_channel_name: name given to the isolated channel to be instantiated.

        """

        # Call to the static internal method __tdmsGroupchannel__.__generateIsolatedchannel__(), which declares an
        # instantiation of the __TdmsGroupChannel__ object.
        self.__setattr__(
            isolated_channel_name, __TdmsGroupChannel__.__generate_isolated_channel__(isolated_channel_name))


# Internal class __tdmsFileChannel__, which contains individual channels that populate a __tdmsFileChannels__ object.
class __TdmsGroupChannel__:
    """Internal class __TdmsGroupChannel__ representing a channel entity of a TDMS group."""

    # __init__ (constructor) of the internal class __TdmsGroupChannel__.
    def __init__(self, channel_name, channel_obj, generate_isolated_channel=False):
        """Initializes an instance of the __TdmsGroupChannel__ object.

        The __TdmsGroupChannel__ object is the entity representing a channel of a TDMS group.

        - **parameters**, **return**, **return types**::

        :param channel_name: declared name of the channel to be instantiated.
        :param channel_obj: an instantiation of the nptdms.group_channel object containing a TDMS channel information.
        :param generate_isolated_channel: boolean flag signaling whether an isolated channel is to be instantiated.
        :return: __TdmsGroupChannel__ object structured in accordance to the properties contained within itself.
        :rtype: __TdmsGroupChannel__.

        """

        # Conditional for checking whether an isolated channel is to be instantiated.
        if not generate_isolated_channel:
            # Declaring channel_name, data and dataMean variables in case the boolean flag owns a false value.
            self.channelName = channel_name
            self.data = channel_obj.data
        else:
            # Declaring channel_name, data (empty array) and dataMean (NaN) in case the boolean flag owns a true value.
            self.channelName = channel_name
            self.data = np.array(())

    # Static internal __generate_isolated_channel__ method.
    @staticmethod
    def __generate_isolated_channel__(isolated_channel_name):
        """Static internal method __generate_isolated_channel__ for generating isolated groups.

        The channels are instantiations of this same object (__tdmsGroupChannel__). This approach is sensible insofar
        it turns necessary to auto-declare certain channels.

        - **parameters**, **returns**, **return types**

        :param isolated_channel_name: string literal providing the declaration name of the instantiated object.
        :return: __tdmsGroupChannel__ object reflecting an instantiation of this same object.
        :rtype: __TdmsGroupChannel__

        """

        # Return statement instantiating a __TdmsGroupChannel__ object.
        return __TdmsGroupChannel__(isolated_channel_name, None, True)


# Internal class __wakeRake__, which contains the basic individual structure for constructing the wake rake processor.
class __RunsPerParam__:
    """Internal class __RunsPerParam__, which contains the basic individual structure for constructing the
    attribute-based structure that classifies the data of a file according to its positional information."""

    # __init__ (constructor) of the internal class __wakeRake__.
    def __init__(self, name, value=None):
        """Initializes an instance of the __RunsPerParam__ object."""

        # Setting an attribute with the given name to the given value.
        if value is None:
            value = []
        self.__setattr__(name, value)


# TODO: DOCUMENT CLASS CORRECTLY.
# Public class XFoilRoot, which contains the basic structure for constructing a data entity coming from an XFoil file.
class XFoilRoot:
    """ """

    #
    angles = None
    #
    lift = None
    #
    drag = None
    #
    eff = None
    #
    pitch = None

    def __init__(self, angles, lift, drag, eff, pitch):
        self.angles = angles
        self.lift = lift
        self.drag = drag
        self.eff = eff
        self.pitch = pitch

# TODO: DOCUMENT CLASS CORRECTLY.
# Public class TdmsFileData, which provides the necessary metadata to read a TDMS file.
class TdmsFileData:
    """ """

    #
    file_path = ""
    #
    file_alias = ""
    #
    file_read_mode = ""
    #
    file_title = ""
    #
    file_ref_mag = ""

    def __init__(self, path, alias, title=None, read_mode=tdmsFilereadMode.Standard, ref_mag=""):
        self.file_path = path
        self.file_alias = alias
        self.file_title = title
        self.file_read_mode = read_mode
        self.file_ref_mag = ref_mag

# TODO: DOCUMENT CODE CORRECTLY.
# Public class XFoilFileData, which provides the necessary metadata to read a XFoil file.
class XFoilFileData:
    """ """

    #
    file_path = ""
    #
    file_alias = ""
    #
    file_title = ""

    def __init__(self, path, alias, title):
        self.file_path = path
        self.file_alias = alias
        self.file_title = title

    # Public class RefMagnitudes, which is a C#-like struct class for non-dimensionalizing purposes.


@dataclass
class RefMagnitudes:
    """Public class for setting referential or characteristic magnitudes for non-dimensionalizing purposes."""

    __length__: float = 1  # Internal member characteristic length (chordwise)
    __span__: float = 1  # Internal member characteristic length (spanwise)
    __thick__: float = 1  # Internal member characteristic length (thickness)
    __u__: float = 1  # Internal member characteristic velocity
    __p__: float = 1  # Internal member characteristic pressure
    __temp__: float = 1  # Internal member characteristic temperature
    __r_const__: float = 1  # Internal member gas constant
    __rho__: float = field(init=False)  # Internal member density
    __mu__: float = field(init=False)  # Internal member viscosity
    __q__: float = field(init=False)  # Internal member dynamic pressure

    # Internal __post_init__ method for computing derived magnitudes such as density or dynamic pressure.
    def __post_init__(self):
        self.__rho__ = self.p / (self.r_const * self.temp)  # Density computation by perfect gas law.
        self.__mu__ = (1.456e-6 * self.temp ** (3. / 2.)) / (self.temp + 110.4)  # Viscosity computation by
        # Sutherland's law.
        self.__q__ = 0.5 * self.rho * self.u * self.u  # Dynamic pressure computation.

    # Declaration of get property length for characteristic length (chordwise).
    @property
    def length(self) -> float:
        return self.__length__

    # Declaration of set property length for characteristic length (chordwise).
    @length.setter
    def length(self, value: float) -> None:
        self.__length__ = value

    # Declaration of get property span for characteristic length (spanwise).
    @property
    def span(self) -> float:
        return self.__span__

    # Declaration of set property span for characteristic length (spanwise).
    @span.setter
    def span(self, value: float) -> None:
        self.__span__ = value

    # Declaration of get property thick for characteristic length (thickness).
    @property
    def thick(self) -> float:
        return self.__thick__

    # Declaration of set property thick for characteristic length (thickness).
    @thick.setter
    def thick(self, value: float) -> None:
        self.__thick__ = value

    # Declaration of get property u for characteristic velocity.
    @property
    def u(self) -> float:
        return self.__u__

    # Declaration of set property u for characteristic velocity.
    @u.setter
    def u(self, value: float) -> None:
        self.__u__ = value
        self.__post_init__()

    # Declaration of get property p for characteristic pressure.
    @property
    def p(self) -> float:
        return self.__p__

    # Declaration of set property p for characteristic pressure.
    @p.setter
    def p(self, value: float) -> None:
        self.__p__ = value
        self.__post_init__()

    # Declaration of get property temp for characteristic temperature.
    @property
    def temp(self) -> float:
        return self.__temp__

    # Declaration of set property temp for characteristic temperature.
    @temp.setter
    def temp(self, value: float) -> None:
        self.__temp__ = value
        self.__post_init__()

    # Declaration of get property r_const for gas constant.
    @property
    def r_const(self) -> float:
        return self.__r_const__

    # Declaration of set property for gas constant.
    @r_const.setter
    def r_const(self, value: float) -> None:
        self.__r_const__ = value
        self.__post_init__()

    # Declaration of get property rho for derived density.
    @property
    def rho(self) -> float:
        return self.__rho__

    # Declaration of get property mu for derived viscosity.
    @property
    def mu(self) -> float:
        return self.__mu__

    # Declaration of get property q for dervied dynamic pressure.
    @property
    def q(self) -> float:
        return self.__q__

    def __init__(self, length=1., span=1., thick=1., vel=0., press=101325., temp=273., r_const=287.):
        """Initializer of data class RefMagnitudes.

        By design philosophy, an initializer should not be necessary when employing the dataclass decorator. However,
        the drawback
        of such a decorator is to lack of a proper signature display on intellisense when instantiating the object.
        The purpose
        of the custom initializer is to provide users with such an aid when using the object.

        - **parameters**, **returns**, **return types**

        :param length: characteristic length (chordwise, m).
        :param span: characteristic length (spanwise, m).
        :param thick: characteristic length (thickness, m).
        :param vel: characteristic velocity (m/s).
        :param press: characteristic pressure (Pa).
        :param temp: characteristic temperature (ºK).
        :param r_const: characteristic gas constant (J/(kg·ºK)).
        :return: an instantiation of the dataclass refMagnitude object.
        :rtype: refMagnitude.

        """

        # Parameter assignment.
        self.__length__ = length
        self.__span__ = span
        self.__thick__ = thick
        self.__u__ = vel
        self.__p__ = press
        self.__temp__ = temp
        self.__r_const__ = r_const
        # __post_init__ method call for obtaining derived magnitudes.
        self.__post_init__()


######################################################################################################################
######################################################################################################################
######################################################################################################################

# OUT-OF-HIERARCHY METHODS
# Out-of-hierarchy methods are intended to provide general functionalities on data manipulation, such as automated
# file opening, programmatic variable declaration or fast plotting.
######################################################################################################################
######################################################################################################################
######################################################################################################################
# TODO: CHECK IF WORKS WITH global_dict set to eval("globals()").
# Public method open_list_of_files.
def open_list_of_files(list_of_files, global_dict):
    """Opens a list of TDMS files in a provided reading mode.

    A number of TDMS files are meant to be passed as a list-like argument, together with a custom enum variable
    determining the reading mode.

    - **parameters**, **returns**, **return types**::

    :param list_of_files: list of instances of the fileData object providing TDMS file metadata.
    :param global_dict: dictionary-typed variable within which the instantiated TDMS files are to be programmatically
    allocated. This variable is intended to match with the globals() dictionary statically declared at runtime.
    :returns: a list of tdmsFileRoot objects representing the instantiations of the provided TDMS files.
    :rtypes: tdmsFileRoot[].

    """

    # Deleting any previously declared variable of the tdmsFileRoot type.
    if "list_of_file_variables" in global_dict.keys():
        globals_deleter([glob for glob in global_dict if (glob in global_dict["list_of_file_variables"])], global_dict)
    # Deleting the list_of_file_variables to avoid having extra files to open due to previously introduced file paths on
    # that variable.
    globals_deleter(["list_of_file_variables"], global_dict)

    # Re-declaring the list_of_file_variables variable.
    list_of_file_variables = []

    # A "for" loop running over the files on the list_of_file_variables variable and programmatically declaring them
    # with sequential names ranging from "file1" to "fileN", where N stands for the number of files in list_of_file
    # variables. The declarations are done as instantiations of the tdmsFileRoot object.
    for file in list_of_files:
        if type(file) == TdmsFileData:
            path = file.file_path
            file_alias_no_whitespace = file.file_alias.replace(" ", "_")
            parent_path = "/".join(path.split("/")[:-1])
            child_path = path.split("/")[-1].split(".")[0]
            parent_path_files = os.listdir(parent_path)
            parent_path_file_extensions = [file.split(".")[-1] for file in parent_path_files]
            if 'p' in parent_path_file_extensions:
                pickle_file = parent_path_files[parent_path_file_extensions.index('p')]
                pickle_file_path = parent_path + "/" + pickle_file
                print("Pickle file found: {0}. Loading data from serialized object.".format(pickle_file))
                global_dict[file_alias_no_whitespace] = pickle.load(open(pickle_file_path, "rb"))
            else:
                global_dict[file_alias_no_whitespace] = TdmsFileRoot(
                    file.file_path, file.file_read_mode, file.file_alias, file.file_ref_mag)
                pickle_file = child_path + '.p'
                pickle_file_path = parent_path + "/" + pickle_file
                print("Saving serialized object in Pickle file: {0}".format(pickle_file))
                pickle.dump(global_dict[file_alias_no_whitespace], open(pickle_file_path, "wb"))
            # Printing statement for monitoring the programmatic declaration/instantiation process.
            print("TDMS file " +
                  str(list_of_files.index(file) + 1) +
                  " of " +
                  str(len(list_of_files)) +
                  " opened in " +
                  file.file_read_mode.name +
                  " mode; stored in " +
                  file_alias_no_whitespace +
                  " global variable.")
            list_of_file_variables.append(file_alias_no_whitespace)
        elif type(file) == XFoilFileData:
            angles, lift, drag, eff, pitch = np.loadtxt(file.filePath,
                                                        delimiter='\t', usecols=(0, 1, 2, 3, 4), unpack=True)
            global_dict[file.file_alias.replace(" ", "_")] = XFoilRoot(angles, lift, drag, eff, pitch)
            list_of_file_variables.append(file.file_alias.replace(" ", "_"))
            # Printing statement for monitoring the programmatic declaration/instantiation process.
            print("XFoil file " +
                  str(list_of_files.index(file) + 1) +
                  " of " +
                  str(len(list_of_files)) +
                  " opened; stored in " +
                  file.fileAlias.replace(" ", "_") +
                  " global variable.")

    # Damping the list_of_file_variables variable to the global_dict dictionary.
    global_dict["list_of_file_variables"] = list_of_file_variables

# TODO: CHECK IF WORKS WITH global_dict set to eval("globals()").
# Public method globals_deleter.
def globals_deleter(look_up_list, global_dict):
    """Deletes entries from a dictionary according to a given list.
    
    This method is intended to update the globals() dictionary statically declared at runtime.

    - **parameters**::

    :param look_up_list: list from which the entries to be deleted are retrieved.
    :param global_dict: dictionary-typed variable within which the instantiated TDMS files are to be programmatically
    allocated. This variable is intended to match with the globals() dictionary statically declared at runtime.
    
    """

    # A "for" loop that runs over the string literals on the look_up_list variables and deletes them from the global
    # dict dictionary, if found.
    for glob in look_up_list:
        if glob in global_dict.keys():
            # Printing statement for monitoring the deleting process.
            print(glob + " global variable has been deleted.")
            del global_dict[glob]


# Public method setStandardData.
def process_standard_data(TdmsFileRoot_obj):
    """Public method that externally sets the attribute structure pertaining a TDMS file values in case the
    initialization has not done so.
   
    -**parameters**::

    :param TdmsFileRoot_obj: an instantiation of the TdmsFileRoot object upon which the attribute structure is to be
    set.
   
    """

    # Calling the internal __setStandardData__() method to set the averaged values.
    TdmsFileRoot_obj.__set_standard_data__()


# Public method meansProcessing.
def process_means(TdmsFileRoot_obj):
    """Public method that externally sets the attribute structure pertaining the averaged values.
    
    -**parameters**::

    :param TdmsFileRoot_obj: an instantiation of the TdmsFileRoot object upon which the attribute structure is to be
    set.

    """

    # Calling the internal __setMeans__() method to set the averaged values.
    TdmsFileRoot_obj.__set_means__()


# Public method wakeRakeProcessing.
def process_wake_rake(TdmsFileRoot_obj, wake_rake_cross_axis=deviceAxis.Y):
    """Public method that externally sets the attribute structure pertaining wake rake values.
    
    -**parameters**::

    :param TdmsFileRoot_obj: an instantiation of the TdmsFileRoot object upon which the attribute structure is to be
    set.
    :param wake_rake_cross_axis: an instance of the custom deviceAxis class specifying the traverse direction, relative
    to wind tunnel axes, of the rake.
    
    """

    # Calling the internal __setWakeRake__ method for performing the task of attribute structuring.
    TdmsFileRoot_obj.__runs_counter__(runsPerParamMode.WakeRake)

######################################################################################################################
######################################################################################################################
######################################################################################################################

def a():
    pass
