# This file was automatically generated by SWIG (http://www.swig.org).
# Version 2.0.9
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.



from sys import version_info
if version_info >= (3,0,0):
    new_instancemethod = lambda func, inst, cls: _named_spectrum.SWIG_PyInstanceMethod_New(func)
else:
    from new import instancemethod as new_instancemethod
if version_info >= (2,6,0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_named_spectrum', [dirname(__file__)])
        except ImportError:
            import _named_spectrum
            return _named_spectrum
        if fp is not None:
            try:
                _mod = imp.load_module('_named_spectrum', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _named_spectrum = swig_import_helper()
    del swig_import_helper
else:
    import _named_spectrum
del version_info
try:
    _swig_property = property
except NameError:
    pass # Python < 2.2 doesn't have 'property'.
def _swig_setattr_nondynamic(self,class_type,name,value,static=1):
    if (name == "thisown"): return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'SwigPyObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    if (not static):
        self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)

def _swig_setattr(self,class_type,name,value):
    return _swig_setattr_nondynamic(self,class_type,name,value,0)

def _swig_getattr(self,class_type,name):
    if (name == "thisown"): return self.this.own()
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError(name)

def _swig_repr(self):
    try: strthis = "proxy of " + self.this.__repr__()
    except: strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

try:
    _object = object
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0


def _swig_setattr_nondynamic_method(set):
    def set_attr(self,name,value):
        if (name == "thisown"): return self.this.own(value)
        if hasattr(self,name) or (name == "this"):
            set(self,name,value)
        else:
            raise AttributeError("You cannot add attributes to %s" % self)
    return set_attr


try:
    import weakref
    weakref_proxy = weakref.proxy
except:
    weakref_proxy = lambda x: x


class SwigPyIterator(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    def __init__(self, *args, **kwargs): raise AttributeError("No constructor defined - class is abstract")
    __repr__ = _swig_repr
    __swig_destroy__ = _named_spectrum.delete_SwigPyIterator
    def __iter__(self): return self
SwigPyIterator.value = new_instancemethod(_named_spectrum.SwigPyIterator_value,None,SwigPyIterator)
SwigPyIterator.incr = new_instancemethod(_named_spectrum.SwigPyIterator_incr,None,SwigPyIterator)
SwigPyIterator.decr = new_instancemethod(_named_spectrum.SwigPyIterator_decr,None,SwigPyIterator)
SwigPyIterator.distance = new_instancemethod(_named_spectrum.SwigPyIterator_distance,None,SwigPyIterator)
SwigPyIterator.equal = new_instancemethod(_named_spectrum.SwigPyIterator_equal,None,SwigPyIterator)
SwigPyIterator.copy = new_instancemethod(_named_spectrum.SwigPyIterator_copy,None,SwigPyIterator)
SwigPyIterator.next = new_instancemethod(_named_spectrum.SwigPyIterator_next,None,SwigPyIterator)
SwigPyIterator.__next__ = new_instancemethod(_named_spectrum.SwigPyIterator___next__,None,SwigPyIterator)
SwigPyIterator.previous = new_instancemethod(_named_spectrum.SwigPyIterator_previous,None,SwigPyIterator)
SwigPyIterator.advance = new_instancemethod(_named_spectrum.SwigPyIterator_advance,None,SwigPyIterator)
SwigPyIterator.__eq__ = new_instancemethod(_named_spectrum.SwigPyIterator___eq__,None,SwigPyIterator)
SwigPyIterator.__ne__ = new_instancemethod(_named_spectrum.SwigPyIterator___ne__,None,SwigPyIterator)
SwigPyIterator.__iadd__ = new_instancemethod(_named_spectrum.SwigPyIterator___iadd__,None,SwigPyIterator)
SwigPyIterator.__isub__ = new_instancemethod(_named_spectrum.SwigPyIterator___isub__,None,SwigPyIterator)
SwigPyIterator.__add__ = new_instancemethod(_named_spectrum.SwigPyIterator___add__,None,SwigPyIterator)
SwigPyIterator.__sub__ = new_instancemethod(_named_spectrum.SwigPyIterator___sub__,None,SwigPyIterator)
SwigPyIterator_swigregister = _named_spectrum.SwigPyIterator_swigregister
SwigPyIterator_swigregister(SwigPyIterator)

SHARED_PTR_DISOWN = _named_spectrum.SHARED_PTR_DISOWN
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

import full_physics_swig.spectrum
import full_physics_swig.generic_object
class ObservableNamedSpectrum(full_physics_swig.generic_object.GenericObject):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    def __init__(self, *args, **kwargs): raise AttributeError("No constructor defined - class is abstract")
    __repr__ = _swig_repr
    __swig_destroy__ = _named_spectrum.delete_ObservableNamedSpectrum
ObservableNamedSpectrum.add_observer_and_keep_reference = new_instancemethod(_named_spectrum.ObservableNamedSpectrum_add_observer_and_keep_reference,None,ObservableNamedSpectrum)
ObservableNamedSpectrum.add_observer = new_instancemethod(_named_spectrum.ObservableNamedSpectrum_add_observer,None,ObservableNamedSpectrum)
ObservableNamedSpectrum.remove_observer = new_instancemethod(_named_spectrum.ObservableNamedSpectrum_remove_observer,None,ObservableNamedSpectrum)
ObservableNamedSpectrum_swigregister = _named_spectrum.ObservableNamedSpectrum_swigregister
ObservableNamedSpectrum_swigregister(ObservableNamedSpectrum)

class ObservablePtrNamedSpectrum(full_physics_swig.generic_object.GenericObject):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    def __init__(self, *args, **kwargs): raise AttributeError("No constructor defined - class is abstract")
    __repr__ = _swig_repr
    __swig_destroy__ = _named_spectrum.delete_ObservablePtrNamedSpectrum
ObservablePtrNamedSpectrum.add_observer_and_keep_reference = new_instancemethod(_named_spectrum.ObservablePtrNamedSpectrum_add_observer_and_keep_reference,None,ObservablePtrNamedSpectrum)
ObservablePtrNamedSpectrum.add_observer = new_instancemethod(_named_spectrum.ObservablePtrNamedSpectrum_add_observer,None,ObservablePtrNamedSpectrum)
ObservablePtrNamedSpectrum.remove_observer = new_instancemethod(_named_spectrum.ObservablePtrNamedSpectrum_remove_observer,None,ObservablePtrNamedSpectrum)
ObservablePtrNamedSpectrum_swigregister = _named_spectrum.ObservablePtrNamedSpectrum_swigregister
ObservablePtrNamedSpectrum_swigregister(ObservablePtrNamedSpectrum)

class ObserverNamedSpectrum(full_physics_swig.generic_object.GenericObject):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self): 
        _named_spectrum.ObserverNamedSpectrum_swiginit(self,_named_spectrum.new_ObserverNamedSpectrum())
    __swig_destroy__ = _named_spectrum.delete_ObserverNamedSpectrum
ObserverNamedSpectrum.notify_update = new_instancemethod(_named_spectrum.ObserverNamedSpectrum_notify_update,None,ObserverNamedSpectrum)
ObserverNamedSpectrum.notify_add = new_instancemethod(_named_spectrum.ObserverNamedSpectrum_notify_add,None,ObserverNamedSpectrum)
ObserverNamedSpectrum.notify_remove = new_instancemethod(_named_spectrum.ObserverNamedSpectrum_notify_remove,None,ObserverNamedSpectrum)
ObserverNamedSpectrum_swigregister = _named_spectrum.ObserverNamedSpectrum_swigregister
ObserverNamedSpectrum_swigregister(ObserverNamedSpectrum)

class NamedSpectrum(full_physics_swig.spectrum.Spectrum):
    """
    Adds name and spec index fields to a Spectrum.

    Useful for sending Spectrum files to output files.

    Note that there are a few closely related classes, with similar
    sounding names. See spectrum_doxygen for a description of each of
    these.

    C++ includes: named_spectrum.h 
    """
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        """
        FullPhysics::NamedSpectrum::NamedSpectrum()
        Default constructor needed for SWIG. 
        """
        _named_spectrum.NamedSpectrum_swiginit(self,_named_spectrum.new_NamedSpectrum(*args))
    def _v_name(self):
        """
        virtual const std::string& FullPhysics::NamedSpectrum::name() const
        Name that makes this a named spectrum. 
        """
        return _named_spectrum.NamedSpectrum__v_name(self)

    @property
    def name(self):
        return self._v_name()

    def _v_index(self):
        """
        virtual int FullPhysics::NamedSpectrum::index() const
        An reference index for the spectrum, ie a spectrometer index. 
        """
        return _named_spectrum.NamedSpectrum__v_index(self)

    @property
    def index(self):
        return self._v_index()

    __swig_destroy__ = _named_spectrum.delete_NamedSpectrum
NamedSpectrum._v_name = new_instancemethod(_named_spectrum.NamedSpectrum__v_name,None,NamedSpectrum)
NamedSpectrum._v_index = new_instancemethod(_named_spectrum.NamedSpectrum__v_index,None,NamedSpectrum)
NamedSpectrum.__str__ = new_instancemethod(_named_spectrum.NamedSpectrum___str__,None,NamedSpectrum)
NamedSpectrum_swigregister = _named_spectrum.NamedSpectrum_swigregister
NamedSpectrum_swigregister(NamedSpectrum)

class vector_named_spectrum(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __iter__(self): return self.iterator()
    def __init__(self, *args): 
        _named_spectrum.vector_named_spectrum_swiginit(self,_named_spectrum.new_vector_named_spectrum(*args))
    __swig_destroy__ = _named_spectrum.delete_vector_named_spectrum
vector_named_spectrum.iterator = new_instancemethod(_named_spectrum.vector_named_spectrum_iterator,None,vector_named_spectrum)
vector_named_spectrum.__nonzero__ = new_instancemethod(_named_spectrum.vector_named_spectrum___nonzero__,None,vector_named_spectrum)
vector_named_spectrum.__bool__ = new_instancemethod(_named_spectrum.vector_named_spectrum___bool__,None,vector_named_spectrum)
vector_named_spectrum.__len__ = new_instancemethod(_named_spectrum.vector_named_spectrum___len__,None,vector_named_spectrum)
vector_named_spectrum.pop = new_instancemethod(_named_spectrum.vector_named_spectrum_pop,None,vector_named_spectrum)
vector_named_spectrum.__getslice__ = new_instancemethod(_named_spectrum.vector_named_spectrum___getslice__,None,vector_named_spectrum)
vector_named_spectrum.__setslice__ = new_instancemethod(_named_spectrum.vector_named_spectrum___setslice__,None,vector_named_spectrum)
vector_named_spectrum.__delslice__ = new_instancemethod(_named_spectrum.vector_named_spectrum___delslice__,None,vector_named_spectrum)
vector_named_spectrum.__delitem__ = new_instancemethod(_named_spectrum.vector_named_spectrum___delitem__,None,vector_named_spectrum)
vector_named_spectrum.__getitem__ = new_instancemethod(_named_spectrum.vector_named_spectrum___getitem__,None,vector_named_spectrum)
vector_named_spectrum.__setitem__ = new_instancemethod(_named_spectrum.vector_named_spectrum___setitem__,None,vector_named_spectrum)
vector_named_spectrum.append = new_instancemethod(_named_spectrum.vector_named_spectrum_append,None,vector_named_spectrum)
vector_named_spectrum.empty = new_instancemethod(_named_spectrum.vector_named_spectrum_empty,None,vector_named_spectrum)
vector_named_spectrum.size = new_instancemethod(_named_spectrum.vector_named_spectrum_size,None,vector_named_spectrum)
vector_named_spectrum.clear = new_instancemethod(_named_spectrum.vector_named_spectrum_clear,None,vector_named_spectrum)
vector_named_spectrum.swap = new_instancemethod(_named_spectrum.vector_named_spectrum_swap,None,vector_named_spectrum)
vector_named_spectrum.get_allocator = new_instancemethod(_named_spectrum.vector_named_spectrum_get_allocator,None,vector_named_spectrum)
vector_named_spectrum.begin = new_instancemethod(_named_spectrum.vector_named_spectrum_begin,None,vector_named_spectrum)
vector_named_spectrum.end = new_instancemethod(_named_spectrum.vector_named_spectrum_end,None,vector_named_spectrum)
vector_named_spectrum.rbegin = new_instancemethod(_named_spectrum.vector_named_spectrum_rbegin,None,vector_named_spectrum)
vector_named_spectrum.rend = new_instancemethod(_named_spectrum.vector_named_spectrum_rend,None,vector_named_spectrum)
vector_named_spectrum.pop_back = new_instancemethod(_named_spectrum.vector_named_spectrum_pop_back,None,vector_named_spectrum)
vector_named_spectrum.erase = new_instancemethod(_named_spectrum.vector_named_spectrum_erase,None,vector_named_spectrum)
vector_named_spectrum.push_back = new_instancemethod(_named_spectrum.vector_named_spectrum_push_back,None,vector_named_spectrum)
vector_named_spectrum.front = new_instancemethod(_named_spectrum.vector_named_spectrum_front,None,vector_named_spectrum)
vector_named_spectrum.back = new_instancemethod(_named_spectrum.vector_named_spectrum_back,None,vector_named_spectrum)
vector_named_spectrum.assign = new_instancemethod(_named_spectrum.vector_named_spectrum_assign,None,vector_named_spectrum)
vector_named_spectrum.resize = new_instancemethod(_named_spectrum.vector_named_spectrum_resize,None,vector_named_spectrum)
vector_named_spectrum.insert = new_instancemethod(_named_spectrum.vector_named_spectrum_insert,None,vector_named_spectrum)
vector_named_spectrum.reserve = new_instancemethod(_named_spectrum.vector_named_spectrum_reserve,None,vector_named_spectrum)
vector_named_spectrum.capacity = new_instancemethod(_named_spectrum.vector_named_spectrum_capacity,None,vector_named_spectrum)
vector_named_spectrum_swigregister = _named_spectrum.vector_named_spectrum_swigregister
vector_named_spectrum_swigregister(vector_named_spectrum)

class ObservableStokesUpdate(full_physics_swig.generic_object.GenericObject):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    def __init__(self, *args, **kwargs): raise AttributeError("No constructor defined - class is abstract")
    __repr__ = _swig_repr
    __swig_destroy__ = _named_spectrum.delete_ObservableStokesUpdate
ObservableStokesUpdate.add_observer_and_keep_reference = new_instancemethod(_named_spectrum.ObservableStokesUpdate_add_observer_and_keep_reference,None,ObservableStokesUpdate)
ObservableStokesUpdate.add_observer = new_instancemethod(_named_spectrum.ObservableStokesUpdate_add_observer,None,ObservableStokesUpdate)
ObservableStokesUpdate.remove_observer = new_instancemethod(_named_spectrum.ObservableStokesUpdate_remove_observer,None,ObservableStokesUpdate)
ObservableStokesUpdate_swigregister = _named_spectrum.ObservableStokesUpdate_swigregister
ObservableStokesUpdate_swigregister(ObservableStokesUpdate)

class ObserverStokesUpdate(full_physics_swig.generic_object.GenericObject):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self): 
        if self.__class__ == ObserverStokesUpdate:
            _self = None
        else:
            _self = self
        _named_spectrum.ObserverStokesUpdate_swiginit(self,_named_spectrum.new_ObserverStokesUpdate(_self, ))
    __swig_destroy__ = _named_spectrum.delete_ObserverStokesUpdate
    def __disown__(self):
        self.this.disown()
        _named_spectrum.disown_ObserverStokesUpdate(self)
        return weakref_proxy(self)
ObserverStokesUpdate.notify_update = new_instancemethod(_named_spectrum.ObserverStokesUpdate_notify_update,None,ObserverStokesUpdate)
ObserverStokesUpdate.notify_add = new_instancemethod(_named_spectrum.ObserverStokesUpdate_notify_add,None,ObserverStokesUpdate)
ObserverStokesUpdate.notify_remove = new_instancemethod(_named_spectrum.ObserverStokesUpdate_notify_remove,None,ObserverStokesUpdate)
ObserverStokesUpdate_swigregister = _named_spectrum.ObserverStokesUpdate_swigregister
ObserverStokesUpdate_swigregister(ObserverStokesUpdate)



