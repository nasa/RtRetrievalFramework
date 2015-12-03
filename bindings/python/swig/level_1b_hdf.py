# This file was automatically generated by SWIG (http://www.swig.org).
# Version 3.0.7
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.





from sys import version_info
if version_info >= (3, 0, 0):
    new_instancemethod = lambda func, inst, cls: _level_1b_hdf.SWIG_PyInstanceMethod_New(func)
else:
    from new import instancemethod as new_instancemethod
if version_info >= (2, 6, 0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_level_1b_hdf', [dirname(__file__)])
        except ImportError:
            import _level_1b_hdf
            return _level_1b_hdf
        if fp is not None:
            try:
                _mod = imp.load_module('_level_1b_hdf', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _level_1b_hdf = swig_import_helper()
    del swig_import_helper
else:
    import _level_1b_hdf
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



_level_1b_hdf.SHARED_PTR_DISOWN_swigconstant(_level_1b_hdf)
SHARED_PTR_DISOWN = _level_1b_hdf.SHARED_PTR_DISOWN

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

import full_physics_swig.level_1b
import full_physics_swig.generic_object
class Level1bHdf(full_physics_swig.level_1b.Level1b):
    """

    This is an abstract interface class which mainly serves to reduce the
    amount of code needed for HDF related Level 1B classes by containing
    the common values needed.

    With the HDF files, we typically generate the uncertainty through a
    NoiseModel. Because it is convenient, we include this functionality in
    this class. This is separate from the handling of common values, so if
    needed we can pull this out as a separate mixin. But for now, there
    doesn't seem to be any point to do this.

    C++ includes: level_1b_hdf.h 
    """
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')

    def __init__(self, *args, **kwargs):
        raise AttributeError("No constructor defined - class is abstract")
    __repr__ = _swig_repr
    __swig_destroy__ = _level_1b_hdf.delete_Level1bHdf

    def _v_noise_model(self, *args):
        """

        void FullPhysics::Level1bHdf::noise_model(const boost::shared_ptr< NoiseModel > &Noise_model)
        Set the noise model. 
        """
        return _level_1b_hdf.Level1bHdf__v_noise_model(self, *args)


    @property
    def noise_model(self):
        return self._v_noise_model()

    @noise_model.setter
    def noise_model(self, value):
      self._v_noise_model(value)

Level1bHdf._v_noise_model = new_instancemethod(_level_1b_hdf.Level1bHdf__v_noise_model, None, Level1bHdf)
Level1bHdf_swigregister = _level_1b_hdf.Level1bHdf_swigregister
Level1bHdf_swigregister(Level1bHdf)



