from builtins import range
import itertools
from functools import wraps, update_wrapper
import inspect

import numpy

# Free functions for PlotRoutines
def call_data_combinations(call_method, pair_arg_idxs, *varargs, **kwargs):

    def pair_combinations(arg):
        return [v for v in list(itertools.combinations(arg, 2)) if id(v[0]) == id(arg[0])]

    num_pairs = len(pair_combinations(varargs[pair_arg_idxs[0]]))

    call_args = []
    for arg_idx, arg_val in enumerate(varargs):
        if arg_idx in pair_arg_idxs:
            call_args.append( pair_combinations(arg_val) )
        else:
            call_args.append( [ arg_val for x in range(num_pairs) ] )

    res = []
    for curr_args in zip(*call_args):
        res.append( call_method(*curr_args, **kwargs) )
    return res

# This was used at one time but due to shuffling code it got orphaned...
def check_same_ids(sounding_ids):
    last_ids = sounding_ids[0]

    for curr_ids in sounding_ids[1:]:
        for chk_id in curr_ids:
            if not any(chk_id == last_ids):
                return False
        last_ids = curr_ids
    return True

def masq_arguments(wrapped):
    "Wraps a function but also masquerades as having the same argument specification"
    argspec = inspect.getargspec(wrapped)
    arg_names = argspec.args
    if argspec.varargs != None: arg_names.append( "*%s" % argspec.varargs )
    if argspec.keywords != None: arg_names.append( "**%s" % argspec.keywords )
    def masq_helper(call_function):
        wrapper_txt = "lambda {arg_spec}: call_function({arg_spec})".format(arg_spec=", ".join(arg_names))
        wrapper = eval(wrapper_txt, {"call_function": call_function})
        update_wrapper(wrapper, wrapped)
        return wrapper
    return masq_helper

# decorator for __dataset__ routines
def call_data_pairs(pair_arg_idxs):
    """Decorator to check the number of sets of data passed and if greater than 2 then
    call the method for each combination of the sets of data"""
    def wrapper_outer(method):
        @masq_arguments(method)
        def wrapper_inner(*varargs, **kwargs):
            if numpy.any(numpy.array(pair_arg_idxs) > len(varargs)-1):
                raise Exception("Pair argument indexes %s exceeds number of arguments for function: %d" % (pair_arg_idxs, len(varargs)))
            test_arg = pair_arg_idxs[0]
            if len(varargs[test_arg]) == 1:
                raise Exception("Can not execute %s with only one set of data" % method.__name__)
            elif len(varargs[test_arg]) != 2:
                return call_data_combinations(method, pair_arg_idxs, *varargs, **kwargs)
            else:
                return method(*varargs, **kwargs)
        return wrapper_inner

    return wrapper_outer

# Add attributes specifying optional dataset names to method
def add_optional_datasets(*opt_datasets):
    def apply_wrapped(method):
        method.optional_datasets = opt_datasets
        return method
    return apply_wrapped
