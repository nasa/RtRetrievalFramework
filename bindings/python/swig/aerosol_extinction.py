# This file was automatically generated by SWIG (http://www.swig.org).
# Version 3.0.12
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.

from sys import version_info as _swig_python_version_info
if _swig_python_version_info >= (3, 0, 0):
    new_instancemethod = lambda func, inst, cls: _aerosol_extinction.SWIG_PyInstanceMethod_New(func)
else:
    from new import instancemethod as new_instancemethod
if _swig_python_version_info >= (2, 7, 0):
    def swig_import_helper():
        import importlib
        pkg = __name__.rpartition('.')[0]
        mname = '.'.join((pkg, '_aerosol_extinction')).lstrip('.')
        try:
            return importlib.import_module(mname)
        except ImportError:
            return importlib.import_module('_aerosol_extinction')
    _aerosol_extinction = swig_import_helper()
    del swig_import_helper
elif _swig_python_version_info >= (2, 6, 0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_aerosol_extinction', [dirname(__file__)])
        except ImportError:
            import _aerosol_extinction
            return _aerosol_extinction
        try:
            _mod = imp.load_module('_aerosol_extinction', fp, pathname, description)
        finally:
            if fp is not None:
                fp.close()
        return _mod
    _aerosol_extinction = swig_import_helper()
    del swig_import_helper
else:
    import _aerosol_extinction
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


class SwigPyIterator(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')

    def __init__(self, *args, **kwargs):
        raise AttributeError("No constructor defined - class is abstract")
    __repr__ = _swig_repr
    __swig_destroy__ = _aerosol_extinction.delete_SwigPyIterator
    def __iter__(self):
        return self
SwigPyIterator.value = new_instancemethod(_aerosol_extinction.SwigPyIterator_value, None, SwigPyIterator)
SwigPyIterator.incr = new_instancemethod(_aerosol_extinction.SwigPyIterator_incr, None, SwigPyIterator)
SwigPyIterator.decr = new_instancemethod(_aerosol_extinction.SwigPyIterator_decr, None, SwigPyIterator)
SwigPyIterator.distance = new_instancemethod(_aerosol_extinction.SwigPyIterator_distance, None, SwigPyIterator)
SwigPyIterator.equal = new_instancemethod(_aerosol_extinction.SwigPyIterator_equal, None, SwigPyIterator)
SwigPyIterator.copy = new_instancemethod(_aerosol_extinction.SwigPyIterator_copy, None, SwigPyIterator)
SwigPyIterator.next = new_instancemethod(_aerosol_extinction.SwigPyIterator_next, None, SwigPyIterator)
SwigPyIterator.__next__ = new_instancemethod(_aerosol_extinction.SwigPyIterator___next__, None, SwigPyIterator)
SwigPyIterator.previous = new_instancemethod(_aerosol_extinction.SwigPyIterator_previous, None, SwigPyIterator)
SwigPyIterator.advance = new_instancemethod(_aerosol_extinction.SwigPyIterator_advance, None, SwigPyIterator)
SwigPyIterator.__eq__ = new_instancemethod(_aerosol_extinction.SwigPyIterator___eq__, None, SwigPyIterator)
SwigPyIterator.__ne__ = new_instancemethod(_aerosol_extinction.SwigPyIterator___ne__, None, SwigPyIterator)
SwigPyIterator.__iadd__ = new_instancemethod(_aerosol_extinction.SwigPyIterator___iadd__, None, SwigPyIterator)
SwigPyIterator.__isub__ = new_instancemethod(_aerosol_extinction.SwigPyIterator___isub__, None, SwigPyIterator)
SwigPyIterator.__add__ = new_instancemethod(_aerosol_extinction.SwigPyIterator___add__, None, SwigPyIterator)
SwigPyIterator.__sub__ = new_instancemethod(_aerosol_extinction.SwigPyIterator___sub__, None, SwigPyIterator)
SwigPyIterator_swigregister = _aerosol_extinction.SwigPyIterator_swigregister
SwigPyIterator_swigregister(SwigPyIterator)

SHARED_PTR_DISOWN = _aerosol_extinction.SHARED_PTR_DISOWN

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

import full_physics_swig.observer
import full_physics_swig.generic_object
import full_physics_swig.state_vector
class ObservableAerosolExtinction(full_physics_swig.generic_object.GenericObject):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')

    def __init__(self, *args, **kwargs):
        raise AttributeError("No constructor defined - class is abstract")
    __repr__ = _swig_repr
    __swig_destroy__ = _aerosol_extinction.delete_ObservableAerosolExtinction
ObservableAerosolExtinction.add_observer_and_keep_reference = new_instancemethod(_aerosol_extinction.ObservableAerosolExtinction_add_observer_and_keep_reference, None, ObservableAerosolExtinction)
ObservableAerosolExtinction.add_observer = new_instancemethod(_aerosol_extinction.ObservableAerosolExtinction_add_observer, None, ObservableAerosolExtinction)
ObservableAerosolExtinction.remove_observer = new_instancemethod(_aerosol_extinction.ObservableAerosolExtinction_remove_observer, None, ObservableAerosolExtinction)
ObservableAerosolExtinction_swigregister = _aerosol_extinction.ObservableAerosolExtinction_swigregister
ObservableAerosolExtinction_swigregister(ObservableAerosolExtinction)

class ObserverAerosolExtinction(full_physics_swig.generic_object.GenericObject):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr

    def __init__(self):
        _aerosol_extinction.ObserverAerosolExtinction_swiginit(self, _aerosol_extinction.new_ObserverAerosolExtinction())
    __swig_destroy__ = _aerosol_extinction.delete_ObserverAerosolExtinction
ObserverAerosolExtinction.notify_update = new_instancemethod(_aerosol_extinction.ObserverAerosolExtinction_notify_update, None, ObserverAerosolExtinction)
ObserverAerosolExtinction.notify_add = new_instancemethod(_aerosol_extinction.ObserverAerosolExtinction_notify_add, None, ObserverAerosolExtinction)
ObserverAerosolExtinction.notify_remove = new_instancemethod(_aerosol_extinction.ObserverAerosolExtinction_notify_remove, None, ObserverAerosolExtinction)
ObserverAerosolExtinction_swigregister = _aerosol_extinction.ObserverAerosolExtinction_swigregister
ObserverAerosolExtinction_swigregister(ObserverAerosolExtinction)

class AerosolExtinction(full_physics_swig.state_vector.StateVectorObserver, ObservableAerosolExtinction):
    """

    This class maps the state vector to the aerosol extinction on each
    level.

    Other objects may depend on the AerosolExtinction, and should be
    updated when the AerosolExtinction is updated. To facilitate that,
    this class in an Oberverable, and objects can add themselves as
    Observers to be notified when the AerosolExtinction is updated.

    When implementing a new class, you almost always will want to derive
    from AerosolExtinctionImpBase rather than from this class. See that
    class for a description.

    C++ includes: aerosol_extinction.h 
    """

    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')

    def __init__(self, *args, **kwargs):
        raise AttributeError("No constructor defined - class is abstract")
    __repr__ = _swig_repr
    __swig_destroy__ = _aerosol_extinction.delete_AerosolExtinction

    def clone(self, *args):
        """

        virtual boost::shared_ptr<AerosolExtinction> FullPhysics::AerosolExtinction::clone(const boost::shared_ptr< Pressure > &Press) const =0
        This version of clone takes a pressure to use.

        The intent is that the pressure has been cloned from the original
        pressure (although this class has no way to verify this). This allows
        sets of objects to be cloned using a common Pressure clone, e.g.
        Atmosphere. 
        """
        return _aerosol_extinction.AerosolExtinction_clone(self, *args)


    def extinction_for_layer(self, i):
        """

        virtual AutoDerivative<double> FullPhysics::AerosolExtinction::extinction_for_layer(int i) const =0
        Extinction for given layer. 
        """
        return _aerosol_extinction.AerosolExtinction_extinction_for_layer(self, i)


    def _v_aerosol_name(self):
        """

        virtual std::string FullPhysics::AerosolExtinction::aerosol_name() const =0
        Name of aerosol. 
        """
        return _aerosol_extinction.AerosolExtinction__v_aerosol_name(self)


    @property
    def aerosol_name(self):
        return self._v_aerosol_name()

AerosolExtinction.clone = new_instancemethod(_aerosol_extinction.AerosolExtinction_clone, None, AerosolExtinction)
AerosolExtinction.extinction_for_layer = new_instancemethod(_aerosol_extinction.AerosolExtinction_extinction_for_layer, None, AerosolExtinction)
AerosolExtinction._v_aerosol_name = new_instancemethod(_aerosol_extinction.AerosolExtinction__v_aerosol_name, None, AerosolExtinction)
AerosolExtinction.__str__ = new_instancemethod(_aerosol_extinction.AerosolExtinction___str__, None, AerosolExtinction)
AerosolExtinction.print_desc = new_instancemethod(_aerosol_extinction.AerosolExtinction_print_desc, None, AerosolExtinction)
AerosolExtinction_swigregister = _aerosol_extinction.AerosolExtinction_swigregister
AerosolExtinction_swigregister(AerosolExtinction)

class vector_aerosol_extinction(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __iter__(self):
        return self.iterator()

    def __init__(self, *args):
        _aerosol_extinction.vector_aerosol_extinction_swiginit(self, _aerosol_extinction.new_vector_aerosol_extinction(*args))
    __swig_destroy__ = _aerosol_extinction.delete_vector_aerosol_extinction
vector_aerosol_extinction.iterator = new_instancemethod(_aerosol_extinction.vector_aerosol_extinction_iterator, None, vector_aerosol_extinction)
vector_aerosol_extinction.__nonzero__ = new_instancemethod(_aerosol_extinction.vector_aerosol_extinction___nonzero__, None, vector_aerosol_extinction)
vector_aerosol_extinction.__bool__ = new_instancemethod(_aerosol_extinction.vector_aerosol_extinction___bool__, None, vector_aerosol_extinction)
vector_aerosol_extinction.__len__ = new_instancemethod(_aerosol_extinction.vector_aerosol_extinction___len__, None, vector_aerosol_extinction)
vector_aerosol_extinction.__getslice__ = new_instancemethod(_aerosol_extinction.vector_aerosol_extinction___getslice__, None, vector_aerosol_extinction)
vector_aerosol_extinction.__setslice__ = new_instancemethod(_aerosol_extinction.vector_aerosol_extinction___setslice__, None, vector_aerosol_extinction)
vector_aerosol_extinction.__delslice__ = new_instancemethod(_aerosol_extinction.vector_aerosol_extinction___delslice__, None, vector_aerosol_extinction)
vector_aerosol_extinction.__delitem__ = new_instancemethod(_aerosol_extinction.vector_aerosol_extinction___delitem__, None, vector_aerosol_extinction)
vector_aerosol_extinction.__getitem__ = new_instancemethod(_aerosol_extinction.vector_aerosol_extinction___getitem__, None, vector_aerosol_extinction)
vector_aerosol_extinction.__setitem__ = new_instancemethod(_aerosol_extinction.vector_aerosol_extinction___setitem__, None, vector_aerosol_extinction)
vector_aerosol_extinction.pop = new_instancemethod(_aerosol_extinction.vector_aerosol_extinction_pop, None, vector_aerosol_extinction)
vector_aerosol_extinction.append = new_instancemethod(_aerosol_extinction.vector_aerosol_extinction_append, None, vector_aerosol_extinction)
vector_aerosol_extinction.empty = new_instancemethod(_aerosol_extinction.vector_aerosol_extinction_empty, None, vector_aerosol_extinction)
vector_aerosol_extinction.size = new_instancemethod(_aerosol_extinction.vector_aerosol_extinction_size, None, vector_aerosol_extinction)
vector_aerosol_extinction.swap = new_instancemethod(_aerosol_extinction.vector_aerosol_extinction_swap, None, vector_aerosol_extinction)
vector_aerosol_extinction.begin = new_instancemethod(_aerosol_extinction.vector_aerosol_extinction_begin, None, vector_aerosol_extinction)
vector_aerosol_extinction.end = new_instancemethod(_aerosol_extinction.vector_aerosol_extinction_end, None, vector_aerosol_extinction)
vector_aerosol_extinction.rbegin = new_instancemethod(_aerosol_extinction.vector_aerosol_extinction_rbegin, None, vector_aerosol_extinction)
vector_aerosol_extinction.rend = new_instancemethod(_aerosol_extinction.vector_aerosol_extinction_rend, None, vector_aerosol_extinction)
vector_aerosol_extinction.clear = new_instancemethod(_aerosol_extinction.vector_aerosol_extinction_clear, None, vector_aerosol_extinction)
vector_aerosol_extinction.get_allocator = new_instancemethod(_aerosol_extinction.vector_aerosol_extinction_get_allocator, None, vector_aerosol_extinction)
vector_aerosol_extinction.pop_back = new_instancemethod(_aerosol_extinction.vector_aerosol_extinction_pop_back, None, vector_aerosol_extinction)
vector_aerosol_extinction.erase = new_instancemethod(_aerosol_extinction.vector_aerosol_extinction_erase, None, vector_aerosol_extinction)
vector_aerosol_extinction.push_back = new_instancemethod(_aerosol_extinction.vector_aerosol_extinction_push_back, None, vector_aerosol_extinction)
vector_aerosol_extinction.front = new_instancemethod(_aerosol_extinction.vector_aerosol_extinction_front, None, vector_aerosol_extinction)
vector_aerosol_extinction.back = new_instancemethod(_aerosol_extinction.vector_aerosol_extinction_back, None, vector_aerosol_extinction)
vector_aerosol_extinction.assign = new_instancemethod(_aerosol_extinction.vector_aerosol_extinction_assign, None, vector_aerosol_extinction)
vector_aerosol_extinction.resize = new_instancemethod(_aerosol_extinction.vector_aerosol_extinction_resize, None, vector_aerosol_extinction)
vector_aerosol_extinction.insert = new_instancemethod(_aerosol_extinction.vector_aerosol_extinction_insert, None, vector_aerosol_extinction)
vector_aerosol_extinction.reserve = new_instancemethod(_aerosol_extinction.vector_aerosol_extinction_reserve, None, vector_aerosol_extinction)
vector_aerosol_extinction.capacity = new_instancemethod(_aerosol_extinction.vector_aerosol_extinction_capacity, None, vector_aerosol_extinction)
vector_aerosol_extinction_swigregister = _aerosol_extinction.vector_aerosol_extinction_swigregister
vector_aerosol_extinction_swigregister(vector_aerosol_extinction)



