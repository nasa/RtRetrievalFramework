# This file was automatically generated by SWIG (http://www.swig.org).
# Version 3.0.12
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.

from sys import version_info as _swig_python_version_info
if _swig_python_version_info >= (3, 0, 0):
    new_instancemethod = lambda func, inst, cls: _nlls_problem_scaled.SWIG_PyInstanceMethod_New(func)
else:
    from new import instancemethod as new_instancemethod
if _swig_python_version_info >= (2, 7, 0):
    def swig_import_helper():
        import importlib
        pkg = __name__.rpartition('.')[0]
        mname = '.'.join((pkg, '_nlls_problem_scaled')).lstrip('.')
        try:
            return importlib.import_module(mname)
        except ImportError:
            return importlib.import_module('_nlls_problem_scaled')
    _nlls_problem_scaled = swig_import_helper()
    del swig_import_helper
elif _swig_python_version_info >= (2, 6, 0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_nlls_problem_scaled', [dirname(__file__)])
        except ImportError:
            import _nlls_problem_scaled
            return _nlls_problem_scaled
        try:
            _mod = imp.load_module('_nlls_problem_scaled', fp, pathname, description)
        finally:
            if fp is not None:
                fp.close()
        return _mod
    _nlls_problem_scaled = swig_import_helper()
    del swig_import_helper
else:
    import _nlls_problem_scaled
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


SHARED_PTR_DISOWN = _nlls_problem_scaled.SHARED_PTR_DISOWN

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

import full_physics_swig.nlls_problem
import full_physics_swig.cost_func_diff
import full_physics_swig.cost_func
import full_physics_swig.problem_state
import full_physics_swig.generic_object
class NLLSProblemScaled(full_physics_swig.nlls_problem.NLLSProblem):
    """

    C++ includes: nlls_problem_scaled.h

    """

    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr

    def __init__(self, s, p):
        """

        NLLSProblemScaled::NLLSProblemScaled(const blitz::Array< double, 1 > &s, const boost::shared_ptr<
        NLLSProblem > &p)
        Default Constructor. 
        """
        _nlls_problem_scaled.NLLSProblemScaled_swiginit(self, _nlls_problem_scaled.new_NLLSProblemScaled(s, p))
    __swig_destroy__ = _nlls_problem_scaled.delete_NLLSProblemScaled

    def _v_residual(self):
        """

        Array< double, 1 > NLLSProblemScaled::residual()
        Return the residual of the scaled NLLS problem at the current set
        point.

        Residual 
        """
        return _nlls_problem_scaled.NLLSProblemScaled__v_residual(self)


    @property
    def residual(self):
        return self._v_residual()


    def _v_jacobian(self):
        """

        Array< double, 2 > NLLSProblemScaled::jacobian()
        Return the Jacobian of the residual of the scaled NLLS problem at the
        current set point.

        The Jacobian of the residual function. 
        """
        return _nlls_problem_scaled.NLLSProblemScaled__v_jacobian(self)


    @property
    def jacobian(self):
        return self._v_jacobian()


    @property
    def residual_size(self):
        return self._v_residual_size()


    @property
    def expected_parameter_size(self):
        return self._v_expected_parameter_size()


    def _v_parameters(self, *args):
        """

        virtual blitz::Array<double, 1> FullPhysics::NLLSProblemScaled::parameters() const
        Just returns the current values of parameters.

        This method is redefined here (see the root base class) because of a
        compiler bug; otherwise, there should be no need for its redefinition.

        Current parameter values 
        """
        return _nlls_problem_scaled.NLLSProblemScaled__v_parameters(self, *args)


    @property
    def parameters(self):
        return self._v_parameters()

    @parameters.setter
    def parameters(self, value):
      self._v_parameters(value)


    def _v_nlls_problem(self):
        """

        boost::shared_ptr<NLLSProblem> FullPhysics::NLLSProblemScaled::nlls_problem()

        """
        return _nlls_problem_scaled.NLLSProblemScaled__v_nlls_problem(self)


    @property
    def nlls_problem(self):
        return self._v_nlls_problem()


    def scale_parameters(self, x):
        """

        Array< double, 1 > NLLSProblemScaled::scale_parameters(const blitz::Array< double, 1 > &x) const
        If x is the input to the NLLS problem that this class is trying to
        scale, then this method scales the input to be used by this class,
        i.e.

        this->parameters(x). The reason for scaling x outside of
        this->parameters(x) is that we can also scale an already scaled NLLS
        problem.

        In summary, the input x is a correct input directly to the NLLS
        problem being scaled. The returned value is correctly scaled to be
        used as input to this scaled NLLS problem. 
        """
        return _nlls_problem_scaled.NLLSProblemScaled_scale_parameters(self, x)


    def unscale_parameters(self, x):
        """

        Array< double, 1 > NLLSProblemScaled::unscale_parameters(const blitz::Array< double, 1 > &x) const
        The input is correctly scaled to be used as input to this scaled NLLS
        problem.

        The returned value is a correct direct input to the NLLS problem being
        scaled. 
        """
        return _nlls_problem_scaled.NLLSProblemScaled_unscale_parameters(self, x)

NLLSProblemScaled._v_residual = new_instancemethod(_nlls_problem_scaled.NLLSProblemScaled__v_residual, None, NLLSProblemScaled)
NLLSProblemScaled._v_jacobian = new_instancemethod(_nlls_problem_scaled.NLLSProblemScaled__v_jacobian, None, NLLSProblemScaled)
NLLSProblemScaled._v_parameters = new_instancemethod(_nlls_problem_scaled.NLLSProblemScaled__v_parameters, None, NLLSProblemScaled)
NLLSProblemScaled._v_nlls_problem = new_instancemethod(_nlls_problem_scaled.NLLSProblemScaled__v_nlls_problem, None, NLLSProblemScaled)
NLLSProblemScaled.scale_parameters = new_instancemethod(_nlls_problem_scaled.NLLSProblemScaled_scale_parameters, None, NLLSProblemScaled)
NLLSProblemScaled.unscale_parameters = new_instancemethod(_nlls_problem_scaled.NLLSProblemScaled_unscale_parameters, None, NLLSProblemScaled)
NLLSProblemScaled_swigregister = _nlls_problem_scaled.NLLSProblemScaled_swigregister
NLLSProblemScaled_swigregister(NLLSProblemScaled)



