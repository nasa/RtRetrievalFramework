# This file was automatically generated by SWIG (http://www.swig.org).
# Version 3.0.7
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.





from sys import version_info
if version_info >= (3, 0, 0):
    new_instancemethod = lambda func, inst, cls: _cost_minimizer_gsl.SWIG_PyInstanceMethod_New(func)
else:
    from new import instancemethod as new_instancemethod
if version_info >= (2, 6, 0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_cost_minimizer_gsl', [dirname(__file__)])
        except ImportError:
            import _cost_minimizer_gsl
            return _cost_minimizer_gsl
        if fp is not None:
            try:
                _mod = imp.load_module('_cost_minimizer_gsl', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _cost_minimizer_gsl = swig_import_helper()
    del swig_import_helper
else:
    import _cost_minimizer_gsl
del version_info
try:
    _swig_property = property
except NameError:
    pass  # Python < 2.2 doesn't have 'property'.


def _swig_setattr_nondynamic(self, class_type, name, value, static=1):
    if (name == "thisown"):
        return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'SwigPyObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name, None)
    if method:
        return method(self, value)
    if (not static):
        object.__setattr__(self, name, value)
    else:
        raise AttributeError("You cannot add attributes to %s" % self)


def _swig_setattr(self, class_type, name, value):
    return _swig_setattr_nondynamic(self, class_type, name, value, 0)


def _swig_getattr_nondynamic(self, class_type, name, static=1):
    if (name == "thisown"):
        return self.this.own()
    method = class_type.__swig_getmethods__.get(name, None)
    if method:
        return method(self)
    if (not static):
        return object.__getattr__(self, name)
    else:
        raise AttributeError(name)

def _swig_getattr(self, class_type, name):
    return _swig_getattr_nondynamic(self, class_type, name, 0)


def _swig_repr(self):
    try:
        strthis = "proxy of " + self.this.__repr__()
    except:
        strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

try:
    _object = object
    _newclass = 1
except AttributeError:
    class _object:
        pass
    _newclass = 0



def _swig_setattr_nondynamic_method(set):
    def set_attr(self, name, value):
        if (name == "thisown"):
            return self.this.own(value)
        if hasattr(self, name) or (name == "this"):
            set(self, name, value)
        else:
            raise AttributeError("You cannot add attributes to %s" % self)
    return set_attr


try:
    import weakref
    weakref_proxy = weakref.proxy
except:
    weakref_proxy = lambda x: x



_cost_minimizer_gsl.SHARED_PTR_DISOWN_swigconstant(_cost_minimizer_gsl)
SHARED_PTR_DISOWN = _cost_minimizer_gsl.SHARED_PTR_DISOWN

def _new_from_init(cls, version, *args):
    '''For use with pickle, covers common case where we just store the
    arguments needed to create an object. See for example HdfFile'''
    if(cls.pickle_format_version() != version):
      raise RuntimeException("Class is expecting a pickled object with version number %d, but we found %d" % (cls.pickle_format_version(), version))
    inst = cls.__new__(cls)
    inst.__init__(*args)
    return inst

def _new_from_set(cls, version, *args):
    '''For use with pickle, covers common case where we use a set function 
    to assign the value'''
    if(cls.pickle_format_version() != version):
      raise RuntimeException("Class is expecting a pickled object with version number %d, but we found %d" % (cls.pickle_format_version(), version))
    inst = cls.__new__(cls)
    inst.__init__()
    inst.set(*args)
    return inst

import full_physics_swig.cost_minimizer
import full_physics_swig.iterative_solver
import full_physics_swig.generic_object
import full_physics_swig.problem_state
class CostMinimizerGSL(full_physics_swig.cost_minimizer.CostMinimizer):
    """

    C++ includes: cost_minimizer_gsl.h

    """
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr

    def __init__(self, max_cost_function_calls, dx_tol_abs, dx_tol_rel, size_tol, p, init_step_size, vrbs=False):
        """

        CostMinimizerGSL::CostMinimizerGSL(int max_cost_function_calls, double dx_tol_abs, double dx_tol_rel,
        double size_tol, const boost::shared_ptr< CostFunc > &p, const
        blitz::Array< double, 1 > &init_step_size, bool vrbs=false)
        Initializes the minimizer.

        Parameters:
        -----------

        max_cost_function_calls:  Input value

        dx_tol_abs:  Input value

        dx_tol_rel:  Input value

        size_tol:

        p:  Input value

        init_step_size:  The initial step stize

        vrbs:  Input value 
        """
        _cost_minimizer_gsl.CostMinimizerGSL_swiginit(self, _cost_minimizer_gsl.new_CostMinimizerGSL(max_cost_function_calls, dx_tol_abs, dx_tol_rel, size_tol, p, init_step_size, vrbs))
    __swig_destroy__ = _cost_minimizer_gsl.delete_CostMinimizerGSL
CostMinimizerGSL_swigregister = _cost_minimizer_gsl.CostMinimizerGSL_swigregister
CostMinimizerGSL_swigregister(CostMinimizerGSL)



