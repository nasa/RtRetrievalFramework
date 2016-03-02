from builtins import object
# Base class which defines some common behaviour for plot routines

import re
import types
from copy import copy
import inspect

import numpy

from matplotlib import pyplot
from matplotlib import axes

class RoutinesBase(object):
    def __init__(self, **kwargs):
        self.analysis_env = kwargs.pop("analysis_env", None)

class PlotRoutinesBase(RoutinesBase):
    def new_axis(self, axis=None, **kwargs):
        if axis != None:
            if hasattr(axis, "pop"):
                return axis.pop()
            elif isinstance(axis, axes.Axes):
                return axis
            else:
                raise TypeError("Unknown object passed as existing axis to use: %s" % axis)
        else:
            fig = pyplot.figure()
            new_axis = fig.add_subplot(111)
            return new_axis

    def _add_diff_stat_text(self, ax, dataset_diff):
        valid_idx = numpy.where(numpy.isnan(dataset_diff) == False)
        diff_mean = numpy.mean(dataset_diff[valid_idx])
        diff_std = numpy.std(dataset_diff[valid_idx])
        
        ax.text(0.02, 0.92, "Mean diff: %7.4f\nDiff std dev: %7.4f" % (diff_mean, diff_std),
                transform=ax.transAxes)

 
class PlotMaker(object):
    """Defines a method of automatically creating plot routines given a dataset name
    by using_create_dataset_routines. This class should be inherited and the __init__
    populated with the dersired calls to _create_dataset_routines"""
                
    def _create_dataset_routines(self, routine_prefix, data_name, value_name, translate_func=None, source_datasets=None, **static_kwargs):
        """Creates a named wrapper  for __dataset__ routines that can be picked up
        by analyze_l2.py and passed the correct dataset values."""
        
        if translate_func == None:
            translate_func = lambda x: x

        if source_datasets == None:
            source_datasets = [ data_name ]

        # Inspects the methods of this class with the pattern __dataset__ in their name
        # then creates wrapper method for the instance that will call these routines but with
        # new names and arguments matching those to be exposed to analyze_l2.py
        created_routines = []
        for method_name, method_obj in inspect.getmembers(self, inspect.ismethod):
            if re.search("__dataset__", method_name):
                # Skip this dataset if it is a _diff_ dataset and we know
                # only one set of data is going to be present
                if (re.search("_diff_", method_name) or re.search("_multi", method_name)) and hasattr(self, "analysis_env") and len(self.analysis_env.obj_names) < 2:
                    continue
                
                # Set up name for new routine and the name for the dataset argument to the method
                routine_name = re.sub("_dataset_", routine_prefix, method_name)
                routine_name = re.sub("^_", "", routine_name)

                # Look for optional dataset in routines and add those so they will be
                # Passed as keyword arguments
                opt_ds_args = []
                if hasattr(method_obj, "optional_datasets"):
                    for ds_name in method_obj.optional_datasets:
                        if self.analysis_env.data_objs[0].get(ds_name, None):
                            opt_ds_args.append( re.sub("/", "__", ds_name.strip("/")) )
                opt_ds_arg_spec = ", ".join(opt_ds_args)
                opt_ds_arg_call = ", ".join([ nm + "=" + nm for nm in opt_ds_args ])
                if len(opt_ds_arg_spec) > 0:
                    opt_ds_arg_spec += "," 
                    opt_ds_arg_call += ","

                datasets_arg_spec = ", ".join( [ re.sub("/", "__", dn.strip("/")) for dn in source_datasets ] )

                # Define a function to perform the operations needed to call the method
                # we really want to call, which is one of the __dataset__ routines
                # Note that method_obj is an argument here because if we counted on
                # closures it would be bound to the last iteration of this loop
                def call_function(method_obj, obj_names, sounding_ids, *dataset_values, **kwargs):
                    data_translated = translate_func(*dataset_values)
                    pass_kwargs = copy(static_kwargs)
                    pass_kwargs.update(kwargs)
                    return method_obj(data_name, value_name, obj_names, sounding_ids, data_translated, **pass_kwargs)
                # Create a lambda function that calls our helper function that we can bind to the class
                # and that we can define the argument names more freely
                method_lambda = "lambda self, obj_names, sounding_id, {datasets_arg_spec}, {opt_ds_arg_spec} **kwargs: call_function(method_obj, obj_names, sounding_id, {datasets_arg_spec}, {opt_ds_arg_call} **kwargs)".format(datasets_arg_spec=datasets_arg_spec, opt_ds_arg_spec=opt_ds_arg_spec, opt_ds_arg_call=opt_ds_arg_call)
                                               
                # Bind the lambda as a class method
                method_decl = eval(method_lambda, {'call_function':call_function,'method_obj':method_obj})
                method_decl.__doc__ = method_obj.__doc__
                method_decl.__name__ = routine_name
                created_routines.append(method_decl)
                setattr(self, method_decl.__name__, types.MethodType(method_decl, self.__class__))

        return created_routines
