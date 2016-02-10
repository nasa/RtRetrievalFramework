# This file was automatically generated by SWIG (http://www.swig.org).
# Version 3.0.7
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.





from sys import version_info
if version_info >= (3, 0, 0):
    new_instancemethod = lambda func, inst, cls: _spectral_domain.SWIG_PyInstanceMethod_New(func)
else:
    from new import instancemethod as new_instancemethod
if version_info >= (2, 6, 0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_spectral_domain', [dirname(__file__)])
        except ImportError:
            import _spectral_domain
            return _spectral_domain
        if fp is not None:
            try:
                _mod = imp.load_module('_spectral_domain', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _spectral_domain = swig_import_helper()
    del swig_import_helper
else:
    import _spectral_domain
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



_spectral_domain.SHARED_PTR_DISOWN_swigconstant(_spectral_domain)
SHARED_PTR_DISOWN = _spectral_domain.SHARED_PTR_DISOWN

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

import full_physics_swig.generic_object
class SpectralDomain(full_physics_swig.generic_object.GenericObject):
    """

    For different instruments, it is more natural to either work with
    wavenumbers (e.g., GOSAT) or wavelength (e.g., OCO).

    Most of our code doesn't care if we are using wavenumber or
    wavelength, so we have this one class that can be either. For code
    that needs one or the other, we supply conversion routines to present
    the data as either wavenumber or wavelength.

    As far as I can determine, there isn't any commonly used name that
    means "either wavelength or wavenumber". We've named this
    "SpectralDomain", where "Domain" is used like "Domain and Range"
    of a function, i.e., this is the X-axis of a spectral plot. Perhaps a
    better name will arise and we can rename this class.

    This class is essentially just a blitz::array with units, and the
    additional functionality to convert to wavenumber or wavelength.

    Note that there are a few closely related classes, with similar
    sounding names. See spectrum_doxygen for a description of each of
    these.

    C++ includes: spectral_domain.h 
    """
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    PREFER_WAVENUMBER = _spectral_domain.SpectralDomain_PREFER_WAVENUMBER
    PREFER_WAVELENGTH = _spectral_domain.SpectralDomain_PREFER_WAVELENGTH

    def __init__(self, *args):
        """

        FullPhysics::SpectralDomain::SpectralDomain()
        Default constructor needed for SWIG. 
        """
        _spectral_domain.SpectralDomain_swiginit(self, _spectral_domain.new_SpectralDomain(*args))

    def _v_data(self):
        """

        const blitz::Array<double, 1>& FullPhysics::SpectralDomain::data() const
        Return data.

        This is either wavenumber or wavelength. This member function is
        intended for classes that don't care which one we are using. If you do
        care, then you should call either wavenumber or wavelength.

        Note that this is a reference to the actual data, so if you intend on
        modifying this you should make a deep copy. 
        """
        return _spectral_domain.SpectralDomain__v_data(self)


    @property
    def data(self):
        return self._v_data()


    def _v_sample_index(self):
        """

        const blitz::Array<int, 1>& FullPhysics::SpectralDomain::sample_index() const
        Return sample index.

        This may be empty if we aren't dealing with the lower resolution grid
        that maps to a sample index. If present, this gives the sample index
        for each of the data() points.

        Note that by convention the sample index is 1 based, so the first
        sample index is 1 rather than zero. 
        """
        return _spectral_domain.SpectralDomain__v_sample_index(self)


    @property
    def sample_index(self):
        return self._v_sample_index()


    def _v_units(self):
        """

        const Unit FullPhysics::SpectralDomain::units() const
        Units that go with data() 
        """
        return _spectral_domain.SpectralDomain__v_units(self)


    @property
    def units(self):
        return self._v_units()


    def _v_type_preference(self):
        """

        TypePreference FullPhysics::SpectralDomain::type_preference() const
        Indicate if this class prefers wavelength or wavenumber.

        This is what data() is. 
        """
        return _spectral_domain.SpectralDomain__v_type_preference(self)


    @property
    def type_preference(self):
        return self._v_type_preference()


    def convert_wave(self, *args):
        """

        Array< double, 1 > SpectralDomain::convert_wave(const Unit &Units) const
        Return data as the supplied the units. 
        """
        return _spectral_domain.SpectralDomain_convert_wave(self, *args)


    def wavenumber(self, *args):
        """

        Array< double, 1 > SpectralDomain::wavenumber(const Unit &Units=units::inv_cm) const
        Return data as wavenumbers.

        You can optionally supply the units to use. Throws an error if the the
        optionally supplied units are not commensurate with cm^-1 
        """
        return _spectral_domain.SpectralDomain_wavenumber(self, *args)


    def wavelength(self, *args):
        """

        Array< double, 1 > SpectralDomain::wavelength(const Unit &Units=units::micron) const
        Return data as wavelengths You can optionally supply the units to use.

        Throws an error if the the optionally supplied units are not
        commensurate with microns 
        """
        return _spectral_domain.SpectralDomain_wavelength(self, *args)


    def photon_to_radiance_factor(self):
        """

        ArrayWithUnit< double, 1 > SpectralDomain::photon_to_radiance_factor() const
        We may want to convert from photon number per second to radiance
        units.

        This gives the factor to use in converting. 
        """
        return _spectral_domain.SpectralDomain_photon_to_radiance_factor(self)

    __swig_destroy__ = _spectral_domain.delete_SpectralDomain
SpectralDomain._v_data = new_instancemethod(_spectral_domain.SpectralDomain__v_data, None, SpectralDomain)
SpectralDomain._v_sample_index = new_instancemethod(_spectral_domain.SpectralDomain__v_sample_index, None, SpectralDomain)
SpectralDomain._v_units = new_instancemethod(_spectral_domain.SpectralDomain__v_units, None, SpectralDomain)
SpectralDomain._v_type_preference = new_instancemethod(_spectral_domain.SpectralDomain__v_type_preference, None, SpectralDomain)
SpectralDomain.convert_wave = new_instancemethod(_spectral_domain.SpectralDomain_convert_wave, None, SpectralDomain)
SpectralDomain.wavenumber = new_instancemethod(_spectral_domain.SpectralDomain_wavenumber, None, SpectralDomain)
SpectralDomain.wavelength = new_instancemethod(_spectral_domain.SpectralDomain_wavelength, None, SpectralDomain)
SpectralDomain.photon_to_radiance_factor = new_instancemethod(_spectral_domain.SpectralDomain_photon_to_radiance_factor, None, SpectralDomain)
SpectralDomain.__str__ = new_instancemethod(_spectral_domain.SpectralDomain___str__, None, SpectralDomain)
SpectralDomain_swigregister = _spectral_domain.SpectralDomain_swigregister
SpectralDomain_swigregister(SpectralDomain)



