# This file was automatically generated by SWIG (http://www.swig.org).
# Version 3.0.12
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.

from sys import version_info as _swig_python_version_info
if _swig_python_version_info >= (3, 0, 0):
    new_instancemethod = lambda func, inst, cls: _composite_perturbation.SWIG_PyInstanceMethod_New(func)
else:
    from new import instancemethod as new_instancemethod
if _swig_python_version_info >= (2, 7, 0):
    def swig_import_helper():
        import importlib
        pkg = __name__.rpartition('.')[0]
        mname = '.'.join((pkg, '_composite_perturbation')).lstrip('.')
        try:
            return importlib.import_module(mname)
        except ImportError:
            return importlib.import_module('_composite_perturbation')
    _composite_perturbation = swig_import_helper()
    del swig_import_helper
elif _swig_python_version_info >= (2, 6, 0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_composite_perturbation', [dirname(__file__)])
        except ImportError:
            import _composite_perturbation
            return _composite_perturbation
        try:
            _mod = imp.load_module('_composite_perturbation', fp, pathname, description)
        finally:
            if fp is not None:
                fp.close()
        return _mod
    _composite_perturbation = swig_import_helper()
    del swig_import_helper
else:
    import _composite_perturbation
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


SHARED_PTR_DISOWN = _composite_perturbation.SHARED_PTR_DISOWN

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

import full_physics_swig.perturbation
import full_physics_swig.generic_object
class PerturbationBuilder(object):
    """

    Class that builds a perturbation to use for a finite difference
    Jacobian.

    We use a std::vector here rather than a blitz::Array just because it
    is easier to add something to the end of std::vector than
    blitz::Array. CompositePerturbation converts this to a blitz::Array
    before finishing the construction of the initial guess.

    C++ includes: composite_perturbation.h 
    """

    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')

    def __init__(self, *args, **kwargs):
        raise AttributeError("No constructor defined - class is abstract")
    __repr__ = _swig_repr
    __swig_destroy__ = _composite_perturbation.delete_PerturbationBuilder

    def _v_number_element(self):
        """

        virtual int FullPhysics::PerturbationBuilder::number_element() const =0
        Number of elements we will be adding to the perturbation.

        0 is a legal value, if we are changing elements but not adding any. 
        """
        return _composite_perturbation.PerturbationBuilder__v_number_element(self)


    @property
    def number_element(self):
        return self._v_number_element()


    def build_perturbation(self, v, index):
        """

        virtual void FullPhysics::PerturbationBuilder::build_perturbation(blitz::Array< double, 1 > &v, int index) const =0
        Called when we need this class to do its part in setting up the
        perturbation array.

        Parameters:
        -----------

        v:  Perturbation vector that should be updated in place.

        index:  Since we are often adding to the end of the state vector,
        index is passed in. This is the sum of the number_elements() of all
        the PerturbationBuilder that appear before this object in the list. 
        """
        return _composite_perturbation.PerturbationBuilder_build_perturbation(self, v, index)

PerturbationBuilder._v_number_element = new_instancemethod(_composite_perturbation.PerturbationBuilder__v_number_element, None, PerturbationBuilder)
PerturbationBuilder.build_perturbation = new_instancemethod(_composite_perturbation.PerturbationBuilder_build_perturbation, None, PerturbationBuilder)
PerturbationBuilder_swigregister = _composite_perturbation.PerturbationBuilder_swigregister
PerturbationBuilder_swigregister(PerturbationBuilder)

class CompositePerturbation(full_physics_swig.perturbation.Perturbation):
    """

    A common way to create a perturbation is to have other classes
    responsible for portions of the state vector (e.g., an Atmosphere
    class creates the portion of the initial guess that handles the
    description of the atmosphere layers).

    This class implements this division.

    This is an example of the "Builder" design pattern. This class is
    what is commonly called the "Director", and the PerturbationBuilder
    classes are the "Builder" classes.

    Note that the PerturbationBuilder objects are called in the order they
    are added to the CompositePerturbation object. A common
    PerturbationBuilder adds additional values to the end state vector, so
    the order is important.

    C++ includes: composite_perturbation.h 
    """

    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr

    @property
    def perturbation(self):
        return self._v_perturbation()


    def add_builder(self, B):
        """

        void FullPhysics::CompositePerturbation::add_builder(const boost::shared_ptr< PerturbationBuilder > &B)
        Add a builder to the build list. 
        """
        return _composite_perturbation.CompositePerturbation_add_builder(self, B)


    def remove_builder(self, B):
        """

        void FullPhysics::CompositePerturbation::remove_builder(const boost::shared_ptr< PerturbationBuilder > &B)
        Remove a builder to the build list. 
        """
        return _composite_perturbation.CompositePerturbation_remove_builder(self, B)


    def __init__(self):
        _composite_perturbation.CompositePerturbation_swiginit(self, _composite_perturbation.new_CompositePerturbation())
    __swig_destroy__ = _composite_perturbation.delete_CompositePerturbation
CompositePerturbation.add_builder = new_instancemethod(_composite_perturbation.CompositePerturbation_add_builder, None, CompositePerturbation)
CompositePerturbation.remove_builder = new_instancemethod(_composite_perturbation.CompositePerturbation_remove_builder, None, CompositePerturbation)
CompositePerturbation_swigregister = _composite_perturbation.CompositePerturbation_swigregister
CompositePerturbation_swigregister(CompositePerturbation)



