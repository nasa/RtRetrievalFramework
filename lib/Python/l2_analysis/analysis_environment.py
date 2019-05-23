from builtins import zip
from builtins import map
from builtins import str
from builtins import range
from builtins import object
import os
import re
import sys
import inspect
import textwrap
import subprocess
from types import GeneratorType
from functools import wraps
import collections

import numpy
import h5py

from matplotlib.pyplot import *
from matplotlib import pylab

import full_physics.acos_file as acos_file
from functools import reduce
import six

def get_term_size():
    try:
        proc = subprocess.Popen(["stty", "size"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return [ int(val) for val in proc.stdout.read().split() ]
    except OSError:
        return None
   
class AnalysisEnvironment(object):
    """Creates an environment where multiple acos_file.L2 files can be analyzed
    using analysis routines where the data passed to the routines are matched
    by sounding id.

    The routines can be run individually or an interactive shell created
    for interactive inspection of the data.
    """
    
    def __init__(self, filenames, routine_objs=None, obj_names=None, addl_files=None):
        # Load files 
        self.data_objs = [ acos_file.L2(fn) for fn in filenames ]

        # Load additional files
        self.addl_objs = []
        if addl_files != None:
            for curr_addl in addl_files:
                # Try adding file as one where we need to find
                # the sounding id dataset, otherwise assume the
                # file has the sounding id dimension as first one
                # in each file ie ECMWF file
                try:
                    curr_obj = acos_file.IdDatasetFinder(curr_addl)
                except (LookupError, AttributeError):
                    curr_obj = acos_file.SoundingFirstFile(curr_addl)
                self.addl_objs.append(curr_obj)

        if len(self.addl_objs) > 0 and len(self.data_objs) != len(self.addl_objs):
            raise Exception("If additional files are supplied, there must be one for each data object")
                
        # Names to use when plotting
        if obj_names != None:
            if len(obj_names) != len(filenames):
                raise ValueError("Number of object names passed must be same as number of filenames")
            self.obj_names = obj_names
        else:
            self.obj_names = self._obj_names_from_filenames(filenames)
            
        # Make sure we can iterate one or more objects
        if not hasattr(routine_objs, "__iter__"):
            routine_objs = [ routine_objs ]

        # Get routines present in the objects
        self.routines = {}
        if routine_objs != None:
            self.add_routines(routine_objs)

        # Correlate the soundings in the supplied files
        self.__dict__.update(self.correlate_soundings())

        # Initialize the global name space
        self._init_global_namespace()        

        # Init local namespace
        self._init_local_namespace()
                          
    def _obj_names_from_filenames(self, filenames):
        # Using filenames, filter out common parts to files
        # This algorithm could come into trouble if a path part occurs
        # more than once

        obj_names = []
        if len(filenames) > 1:
            part_occur = {}
            for fn in filenames:
                for file_part in fn.split(os.path.sep):
                    part_occur[file_part] = part_occur.get(file_part, 0) + 1

            common_parts = []
            for part, count in list(part_occur.items()):
                if count == len(filenames):
                    common_parts.append(part)

            for fn in filenames:
                for part in common_parts:
                    fn = re.sub("/?"+part+"/?", "", fn)
                obj_names.append(fn)
        else:
            obj_names = [ os.path.basename(fn) for fn in filenames ]

        return obj_names

    def _init_local_namespace(self):
        # Set up namespaces used for executing interactively or by scripts
        self._local_namespace = {}

        # Add analysis (self) object to locals
        self._local_namespace['analysis_env'] = self

        # Add a short cut to analysis_env.get_object_data
        self._local_namespace['get_data'] = self.get_object_data

    def _init_global_namespace(self):
        # Import pylab, numpy, h5py into global namespace
        self._global_namespace = {}
        self._global_namespace.update( __import__("matplotlib.pylab", {}, {}, "*").__dict__ )
        self._global_namespace.update( __import__("numpy", {}, {}, "*").__dict__ )
        self._global_namespace["h5py"] = h5py

    def add_routines(self, routine_objs):
        if not hasattr(routine_objs, "__iter__"):
            routine_objs = (routine_objs,)

        for obj in routine_objs:
            if inspect.ismethod(obj) or inspect.isfunction(obj):
                self.routines[obj.__name__] = obj
                self._add_helper_function(obj.__name__, obj)
            elif hasattr(obj, "__class__"):
                for routine_name, routine_method in inspect.getmembers(obj, inspect.ismethod):
                    # Ignore methods starting with _
                    if routine_name[0] != "_":
                        self.routines[routine_name] = routine_method
                        self._add_helper_function(routine_name, routine_method)
            else:
                raise ValueError("The object %s is can not be used as a routine" % obj)

    def _indexes_from_lists(self, comparison_ids, sounding_lists):
        # For each sounding list make an index list of where the ids appear in their files
        # From the result of nonzero get the value using [0] so we can turn it into a tuple
        # and use for indexing
        indexes = []
        for snd_list in sounding_lists:
            where_same = reduce(lambda x,y: x+y, [snd_id == snd_list for snd_id in comparison_ids ])
            if len(where_same.shape) > 1:
                dim_indexes = [ dim for dim in zip(*numpy.nonzero(where_same)) ]
            else:
                dim_indexes = numpy.nonzero(where_same)[0]
            indexes.append(tuple(dim_indexes))
        return tuple(indexes)
    
    def correlate_soundings(self, in_sounding_lists=None, l2_sounding_lists=None, addl_id_lists=None):
        # Ensure we are working on a container of containser
        if l2_sounding_lists == None:
            l2_sounding_lists = [ obj.get_sounding_ids()[:] for obj in self.data_objs ]
        
        if in_sounding_lists == None:
            in_sounding_lists = l2_sounding_lists
        elif not isinstance(in_sounding_lists, collections.Iterable) or isinstance(in_sounding_lists, str):
            in_sounding_lists = [ (in_sounding_lists,)]
        elif isinstance(in_sounding_lists, collections.Container) and not isinstance(in_sounding_lists[0], collections.Container):
            in_sounding_lists = [in_sounding_lists]
       
        # Get a list of the ids common to all sets
        comparison_ids = list(reduce(set.intersection, list(map(set, in_sounding_lists))))

        # Sort this list since we used a set and bad things can happen when this list is not sorted
        comparison_ids.sort()

        # Get indexes of those soundings in the L2 sounding lists
        if len(comparison_ids) > 0:
            data_id_indexes = self._indexes_from_lists(comparison_ids, l2_sounding_lists)
        else:
            data_id_indexes = None

        # Indexes into additional objects if loaded
        addl_id_indexes = None
        if len(comparison_ids) > 0 and len(self.addl_objs) > 0:
            if addl_id_lists == None:
                addl_id_lists = [ obj.get_sounding_ids()[:] for obj in self.addl_objs ]
            addl_id_indexes = self._indexes_from_lists(comparison_ids, addl_id_lists)
            
        return {'comparison_ids': comparison_ids,
                'data_id_indexes': data_id_indexes,
                'addl_id_indexes': addl_id_indexes,
                }

    def routine_names(self):
        return list(self.routines.keys())

    def formatted_routine_names(self, group_len=2, part_sep="_", name_indent="  "):
        # Wrap routine names at the screen width if we can obtain it
        term_size = get_term_size()
        if term_size != None and len(term_size) >= 2:
            term_width = term_size[1]
        else:
            term_width = 75
        wrapper = textwrap.TextWrapper(width=term_width, subsequent_indent=name_indent)
        
        # Group routine names by common part(s) which is a string
        # with the parts are seperated by the part_sep character
        # up to group_len
        name_groups = {}
        for r_name in self.routine_names():
            g_name = part_sep.join(r_name.split(part_sep)[:group_len])
            name_groups[g_name] = name_groups.get(g_name, [])
            name_groups[g_name].append(r_name)

        # Return a string containing the routine names formatted and sorted
        routine_name_str = ""
        for g_name in sorted(name_groups.keys()):
            group_str = ", ".join(name_groups[g_name]) + "\n"
            group_str = wrapper.fill(group_str)
            routine_name_str += group_str + "\n"
                
        return routine_name_str

    def get_object_data(self, data_name, ids=None, **kwargs):
        # If a filter was supplied use it for getting data indexes
        if ids == None:
            use_data_indexes = self.data_id_indexes
            use_addl_indexes = self.addl_id_indexes
        else:
            correlate_res = self.correlate_soundings(ids)
            use_data_indexes = correlate_res['data_id_indexes']
            use_addl_indexes = correlate_res['addl_id_indexes']

        # Get the data to be passed to a routine for each object given a specific 
        # data name
        obj_vals = []
        for d_idx, d_obj, d_indexes in zip(list(range(len(self.data_objs))), self.data_objs, use_data_indexes):
            if len(d_indexes) == 0:
                raise ValueError("No sounding ids avaliable from file: %s for data name: %s" % (d_obj.filename, data_name))

            try:
                obj_vals.append( d_obj.get_sounding_data(data_name, indexes=d_indexes, preserve_shape=True, **kwargs) )
            except ValueError as err:
                # Try to find dataset in an additional file
                if d_idx < len(self.addl_objs):
                    try:
                        a_obj = self.addl_objs[d_idx]
                        obj_vals.append( a_obj.get_sounding_data(data_name, indexes=use_addl_indexes[d_idx], **kwargs) )
                    except ValueError as err:
                        raise ValueError("Error getting data named %s from either %s or %s : %s" % (data_name, d_obj.filename, a_obj.filename, str(err)))
                else:
                    raise ValueError("Error getting data named %s from file %s : %s" % (data_name, d_obj.filename, str(err)))

        return obj_vals

    def call_analysis_routine(self, routine_name, *arg_list, **kwargs):
        """Calls a desired routine with the appropriate data values passed to it"""

        # Find the requested routine  and output a helpful message if not found
        try:
            routine_method = self.routines[routine_name]
        except KeyError:
            raise KeyError("Could not find routine named: %s" % routine_name)

        # Get list of arguments from routine
        argspec = inspect.getargspec(routine_method)
        routine_args = list(argspec.args)

        if argspec.defaults != None:
            defaults = list(argspec.defaults)
        else:
            defaults = None

        # Remove arguments that have default values but not a value in kwargs
        while(defaults != None and len(defaults) > 0):
            defaults.pop()
            if not routine_args[-1] in list(kwargs.keys()):
                routine_args.pop()

        # For each argument not called self, search for data
        # that matches its name and extract it using the
        # indexes for that file of the common sounding ids
        arg_list = list(arg_list) # Change tuple var arg into list
        for arg_name in routine_args:
            if hasattr(self, arg_name):
                # Take argument from analyis environment argument
                arg_list.append(getattr(self, arg_name))
            elif arg_name in list(kwargs.keys()):
                # Take argument from kwargs, remove so we do not
                # get duplicates should we pass it as **kwargs
                # below
                arg_list.append(kwargs.pop(arg_name))
            elif arg_name == "sounding_id":
                obj_vals = []
                for d_obj in self.data_objs:
                    passed_ids = kwargs.get('ids', None)
                    if passed_ids != None:
                        # passed_ids should be a list of lists to match
                        # how sounding ids are passed around internally
                        if not hasattr(passed_ids, "__iter__") or not len(passed_ids) == 1:
                            raise TypeError("sounding ids passed through the 'ids' keyword should be an iterable with one element containing the list of ids")
                        obj_vals.append(passed_ids[0])
                    else:
                        obj_vals.append( self.comparison_ids )
                arg_list.append(obj_vals)
            elif arg_name != 'self' and arg_name[0] != "_":
                # Look up argument in hdf files
                arg_list.append(self.get_object_data(arg_name, **kwargs))

        # Call routine with constructed data, a list of data from each file
        # for each argument
        try:
            if argspec.keywords != None:
                return routine_method(*arg_list, **kwargs)
            else:
                return routine_method(*arg_list)
        except TypeError as e:
            # Raise error preserving traceback
            type, value, traceback = sys.exc_info()
            err_msg = "%s when calling %s with arguments named: %s with value list of size: %d" % (e, routine_name, routine_args, len(arg_list))
            six.reraise(TypeError, err_msg, traceback)

    def _add_helper_function(self, routine_name, routine_obj):
        @wraps(routine_obj)
        def routine_call_helper(*varargs, **kwargs):
            result = self.call_analysis_routine(routine_name, *varargs, **kwargs)

            # Evaluate the generator results so plots show up
            # immediately in an interactive environment
            if type(result) is GeneratorType:
                return tuple(result)
            else:
                return result
        
        # Add a helper for a given routine as it is added
        self._local_namespace[routine_name] = routine_call_helper

    def launch_shell(self, argv=[]):
        # Configure prompts and messages
        in_template  = 'L2A: In <\\#>: '
        in_template2 = '   .\\D.:'
        out_template = 'L2A: Out<\\#>: '
        banner = '*** Launcing L2 Analysis Shell ***\nAvailable analysis routines:\n%s' % self.formatted_routine_names()
        exit_msg = '*** Bye ***'
        
        # Set pylab environment as interactive
        pylab.interactive(True)

        # Try using an older version of IPython first 
        try:
            argv += [ '-pi1', in_template, '-pi2', in_template2, '-po', out_template ]

            from IPython.Shell import IPShellEmbed

            ipshell = IPShellEmbed(argv,banner=banner,exit_msg=exit_msg)
            ipshell(local_ns=self._local_namespace, global_ns=self._global_namespace)

        except ImportError as imp_err:
            # Newer version of IPython, 0.11 onward use this interface
            from IPython.config.loader import Config

            cfg = Config()
            prompt_config = cfg.PromptManager
            prompt_config.in_template = in_template
            prompt_config.in2_template = in_template2 
            prompt_config.out_template = out_template 

            from IPython.frontend.terminal.embed import InteractiveShellEmbed

            ipshell = InteractiveShellEmbed(config=cfg, banner1=banner, exit_msg=exit_msg)

            # There is no access to global namespace in this version of IPython
            # put everything into the local namespace
            namespace = {}
            namespace.update(self._local_namespace)
            namespace.update(self._global_namespace)
            ipshell(local_ns=namespace)


    def run_script(self, script_file):
        """Reads in script file and wraps it in a function so that closures work and the 
        execution of the script works the same way as %run under the interactive shell"""

        indent = "\t"
        script_code = [ "def main():" ]
        with open(script_file,"r") as f_obj:
            for script_line in f_obj.readlines():
                script_code.append( indent + script_line.rstrip() )
        script_code.append("main()")
        code = compile("\n".join(script_code), script_file, 'exec')
        exec(code, self._local_namespace, self._global_namespace)
