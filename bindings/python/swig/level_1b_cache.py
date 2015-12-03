# This file was automatically generated by SWIG (http://www.swig.org).
# Version 3.0.7
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.





from sys import version_info
if version_info >= (3, 0, 0):
    new_instancemethod = lambda func, inst, cls: _level_1b_cache.SWIG_PyInstanceMethod_New(func)
else:
    from new import instancemethod as new_instancemethod
if version_info >= (2, 6, 0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_level_1b_cache', [dirname(__file__)])
        except ImportError:
            import _level_1b_cache
            return _level_1b_cache
        if fp is not None:
            try:
                _mod = imp.load_module('_level_1b_cache', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _level_1b_cache = swig_import_helper()
    del swig_import_helper
else:
    import _level_1b_cache
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



_level_1b_cache.SHARED_PTR_DISOWN_swigconstant(_level_1b_cache)
SHARED_PTR_DISOWN = _level_1b_cache.SHARED_PTR_DISOWN

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
class Level1bCache(full_physics_swig.level_1b.Level1b):
    """

    This is a Level1b implementation that just saves the data read from
    another Level1b object.

    We then allow these values to be changed if desired. This can be
    useful when setting up special run in Python, among other uses.

    C++ includes: level_1b_cache.h 
    """
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr

    def __init__(self, L1_in):
        """

        Level1bCache::Level1bCache(const Level1b &L1_in)
        Constructor. 
        """
        _level_1b_cache.Level1bCache_swiginit(self, _level_1b_cache.new_Level1bCache(L1_in))

    def set_latitude(self, i, V):
        """

        void FullPhysics::Level1bCache::set_latitude(int i, const DoubleWithUnit &V)
        Change value. 
        """
        return _level_1b_cache.Level1bCache_set_latitude(self, i, V)


    def set_longitude(self, i, V):
        """

        void FullPhysics::Level1bCache::set_longitude(int i, const DoubleWithUnit &V)
        Change value. 
        """
        return _level_1b_cache.Level1bCache_set_longitude(self, i, V)


    def set_sounding_zenith(self, i, V):
        """

        void FullPhysics::Level1bCache::set_sounding_zenith(int i, const DoubleWithUnit &V)
        Change value. 
        """
        return _level_1b_cache.Level1bCache_set_sounding_zenith(self, i, V)


    def set_sounding_azimuth(self, i, V):
        """

        void FullPhysics::Level1bCache::set_sounding_azimuth(int i, const DoubleWithUnit &V)
        Change value. 
        """
        return _level_1b_cache.Level1bCache_set_sounding_azimuth(self, i, V)


    def set_stokes_coefficient(self, i, V):
        """

        void FullPhysics::Level1bCache::set_stokes_coefficient(int i, const blitz::Array< double, 1 > &V)
        Change value. 
        """
        return _level_1b_cache.Level1bCache_set_stokes_coefficient(self, i, V)


    def set_solar_zenith(self, i, V):
        """

        void FullPhysics::Level1bCache::set_solar_zenith(int i, const DoubleWithUnit &V)
        Change value. 
        """
        return _level_1b_cache.Level1bCache_set_solar_zenith(self, i, V)


    def set_solar_azimuth(self, i, V):
        """

        void FullPhysics::Level1bCache::set_solar_azimuth(int i, const DoubleWithUnit &V)
        Change value. 
        """
        return _level_1b_cache.Level1bCache_set_solar_azimuth(self, i, V)


    def set_altitude(self, i, V):
        """

        void FullPhysics::Level1bCache::set_altitude(int i, const DoubleWithUnit &V)
        Change value. 
        """
        return _level_1b_cache.Level1bCache_set_altitude(self, i, V)


    def set_relative_velocity(self, i, V):
        """

        void FullPhysics::Level1bCache::set_relative_velocity(int i, const DoubleWithUnit &V)
        Change value. 
        """
        return _level_1b_cache.Level1bCache_set_relative_velocity(self, i, V)


    def set_spectral_coefficient(self, i, V):
        """

        void FullPhysics::Level1bCache::set_spectral_coefficient(int i, const ArrayWithUnit< double, 1 > &V)
        Change value. 
        """
        return _level_1b_cache.Level1bCache_set_spectral_coefficient(self, i, V)


    def set_time(self, i, V):
        """

        void FullPhysics::Level1bCache::set_time(int i, const Time &V)
        Change value. 
        """
        return _level_1b_cache.Level1bCache_set_time(self, i, V)


    def radiance(self, Spec_index):
        """

        virtual SpectralRange FullPhysics::Level1bCache::radiance(int i) const

        """
        return _level_1b_cache.Level1bCache_radiance(self, Spec_index)


    def set_radiance(self, *args):
        """

        void FullPhysics::Level1bCache::set_radiance(int i, const SpectralRange &V, const std::vector< int > &Plist)
        Change value, but only for a subset of pixels.

        This might come from the ForwardModelSpectralGrid for example. 
        """
        return _level_1b_cache.Level1bCache_set_radiance(self, *args)


    def set_sounding_id(self, V):
        """

        void FullPhysics::Level1bCache::set_sounding_id(int64_t V)
        Change value. 
        """
        return _level_1b_cache.Level1bCache_set_sounding_id(self, V)


    def set_exposure_index(self, V):
        """

        void FullPhysics::Level1bCache::set_exposure_index(int V)
        Change value. 
        """
        return _level_1b_cache.Level1bCache_set_exposure_index(self, V)

    __swig_destroy__ = _level_1b_cache.delete_Level1bCache
Level1bCache.set_latitude = new_instancemethod(_level_1b_cache.Level1bCache_set_latitude, None, Level1bCache)
Level1bCache.set_longitude = new_instancemethod(_level_1b_cache.Level1bCache_set_longitude, None, Level1bCache)
Level1bCache.set_sounding_zenith = new_instancemethod(_level_1b_cache.Level1bCache_set_sounding_zenith, None, Level1bCache)
Level1bCache.set_sounding_azimuth = new_instancemethod(_level_1b_cache.Level1bCache_set_sounding_azimuth, None, Level1bCache)
Level1bCache.set_stokes_coefficient = new_instancemethod(_level_1b_cache.Level1bCache_set_stokes_coefficient, None, Level1bCache)
Level1bCache.set_solar_zenith = new_instancemethod(_level_1b_cache.Level1bCache_set_solar_zenith, None, Level1bCache)
Level1bCache.set_solar_azimuth = new_instancemethod(_level_1b_cache.Level1bCache_set_solar_azimuth, None, Level1bCache)
Level1bCache.set_altitude = new_instancemethod(_level_1b_cache.Level1bCache_set_altitude, None, Level1bCache)
Level1bCache.set_relative_velocity = new_instancemethod(_level_1b_cache.Level1bCache_set_relative_velocity, None, Level1bCache)
Level1bCache.set_spectral_coefficient = new_instancemethod(_level_1b_cache.Level1bCache_set_spectral_coefficient, None, Level1bCache)
Level1bCache.set_time = new_instancemethod(_level_1b_cache.Level1bCache_set_time, None, Level1bCache)
Level1bCache.radiance = new_instancemethod(_level_1b_cache.Level1bCache_radiance, None, Level1bCache)
Level1bCache.set_radiance = new_instancemethod(_level_1b_cache.Level1bCache_set_radiance, None, Level1bCache)
Level1bCache.set_sounding_id = new_instancemethod(_level_1b_cache.Level1bCache_set_sounding_id, None, Level1bCache)
Level1bCache.set_exposure_index = new_instancemethod(_level_1b_cache.Level1bCache_set_exposure_index, None, Level1bCache)
Level1bCache_swigregister = _level_1b_cache.Level1bCache_swigregister
Level1bCache_swigregister(Level1bCache)



