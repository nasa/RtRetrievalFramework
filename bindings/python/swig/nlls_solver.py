# This file was automatically generated by SWIG (http://www.swig.org).
# Version 3.0.12
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.

from sys import version_info as _swig_python_version_info
if _swig_python_version_info >= (3, 0, 0):
    new_instancemethod = lambda func, inst, cls: _nlls_solver.SWIG_PyInstanceMethod_New(func)
else:
    from new import instancemethod as new_instancemethod
if _swig_python_version_info >= (2, 7, 0):
    def swig_import_helper():
        import importlib
        pkg = __name__.rpartition('.')[0]
        mname = '.'.join((pkg, '_nlls_solver')).lstrip('.')
        try:
            return importlib.import_module(mname)
        except ImportError:
            return importlib.import_module('_nlls_solver')
    _nlls_solver = swig_import_helper()
    del swig_import_helper
elif _swig_python_version_info >= (2, 6, 0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_nlls_solver', [dirname(__file__)])
        except ImportError:
            import _nlls_solver
            return _nlls_solver
        try:
            _mod = imp.load_module('_nlls_solver', fp, pathname, description)
        finally:
            if fp is not None:
                fp.close()
        return _mod
    _nlls_solver = swig_import_helper()
    del swig_import_helper
else:
    import _nlls_solver
del _swig_python_version_info

try:
    _swig_property = property
except NameError:
    pass  # Python < 2.2 doesn't have 'property'.

try:
    import builtins as __builtin__
except ImportError:
    import __builtin__

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


def _swig_getattr(self, class_type, name):
    if (name == "thisown"):
        return self.this.own()
    method = class_type.__swig_getmethods__.get(name, None)
    if method:
        return method(self)
    raise AttributeError("'%s' object has no attribute '%s'" % (class_type.__name__, name))


def _swig_repr(self):
    try:
        strthis = "proxy of " + self.this.__repr__()
    except __builtin__.Exception:
        strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)


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
except __builtin__.Exception:
    weakref_proxy = lambda x: x


SHARED_PTR_DISOWN = _nlls_solver.SHARED_PTR_DISOWN

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

import full_physics_swig.iterative_solver_der
import full_physics_swig.iterative_solver
import full_physics_swig.generic_object
import full_physics_swig.cost_func_diff
import full_physics_swig.cost_func
import full_physics_swig.problem_state
class NLLSSolver(full_physics_swig.iterative_solver_der.IterativeSolverDer):
    """

    The base class for the solvers of the Nonlinear-Least-Squares Problem.

    This is the base class for all Nonlinear-Least-Squares solvers.

    This class is associated with a problem (NLLSProblem) because the
    problem interface is determined: provide a point in the parameter
    space

    evaluate the residual function (a vector function) at the point

    evaluate the Jacobian (a matrix function) of the residual at the point

    C++ includes: nlls_solver.h 
    """

    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')

    def __init__(self, *args, **kwargs):
        raise AttributeError("No constructor defined - class is abstract")
    __repr__ = _swig_repr
    __swig_destroy__ = _nlls_solver.delete_NLLSSolver

    def _v_nlls_problem(self):
        """

        const boost::shared_ptr<NLLSProblem>& FullPhysics::NLLSSolver::nlls_problem() const
        Returns the Nonlinear Least Squares problem.

        This method returns the Nonlinear Least Squares problem that is passed
        to the constructor of the solver.

        Nonlinear least squares problem 
        """
        return _nlls_solver.NLLSSolver__v_nlls_problem(self)


    @property
    def nlls_problem(self):
        return self._v_nlls_problem()

NLLSSolver._v_nlls_problem = new_instancemethod(_nlls_solver.NLLSSolver__v_nlls_problem, None, NLLSSolver)
NLLSSolver_swigregister = _nlls_solver.NLLSSolver_swigregister
NLLSSolver_swigregister(NLLSSolver)



