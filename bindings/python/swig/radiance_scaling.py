# This file was automatically generated by SWIG (http://www.swig.org).
# Version 2.0.9
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.



from sys import version_info
if version_info >= (3,0,0):
    new_instancemethod = lambda func, inst, cls: _radiance_scaling.SWIG_PyInstanceMethod_New(func)
else:
    from new import instancemethod as new_instancemethod
if version_info >= (2,6,0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_radiance_scaling', [dirname(__file__)])
        except ImportError:
            import _radiance_scaling
            return _radiance_scaling
        if fp is not None:
            try:
                _mod = imp.load_module('_radiance_scaling', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _radiance_scaling = swig_import_helper()
    del swig_import_helper
else:
    import _radiance_scaling
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


SHARED_PTR_DISOWN = _radiance_scaling.SHARED_PTR_DISOWN
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

import full_physics_swig.instrument_correction
import full_physics_swig.state_vector
import full_physics_swig.generic_object
class RadianceScaling(full_physics_swig.instrument_correction.InstrumentCorrection):
    """
    This abstract class provides the generic capabilities for applying a
    radiance scaling to a Radiance.

    The radiance scaling slope is referenced to a reference
    wavelength/wavenumber

    This class can support both a scale factor polynomial and a single
    radiance wide offset.

    C++ includes: radiance_scaling.h 
    """
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    def __init__(self, *args, **kwargs): raise AttributeError("No constructor defined - class is abstract")
    __repr__ = _swig_repr
    __swig_destroy__ = _radiance_scaling.delete_RadianceScaling
    def print_desc(self, *args):
        """
        void RadianceScaling::print(std::ostream &Os) const

        """
        return _radiance_scaling.RadianceScaling_print_desc(self, *args)

    def apply_scaling(self, *args):
        """
        void RadianceScaling::apply_scaling(const SpectralDomain &Grid, SpectralRange &Radiance) const
        Apply scaling and offset coefficients to Radiance. 
        """
        return _radiance_scaling.RadianceScaling_apply_scaling(self, *args)

    def _v_radiance_scaling_coeff(self):
        """
        virtual blitz::Array<double, 1> FullPhysics::RadianceScaling::radiance_scaling_coeff() const
        Return radiance scaling coefficients for the output file. 
        """
        return _radiance_scaling.RadianceScaling__v_radiance_scaling_coeff(self)

    @property
    def radiance_scaling_coeff(self):
        return self._v_radiance_scaling_coeff()

    def _v_radiance_scaling_coeff_uncertainty(self):
        """
        virtual blitz::Array<double, 1> FullPhysics::RadianceScaling::radiance_scaling_coeff_uncertainty() const =0
        Return radiance scaling coefficients uncertainty for the output file.

        """
        return _radiance_scaling.RadianceScaling__v_radiance_scaling_coeff_uncertainty(self)

    @property
    def radiance_scaling_coeff_uncertainty(self):
        return self._v_radiance_scaling_coeff_uncertainty()

    def _v_radiance_offset(self):
        """
        virtual double FullPhysics::RadianceScaling::radiance_offset() const
        Return radiance scaling offset for the output file. 
        """
        return _radiance_scaling.RadianceScaling__v_radiance_offset(self)

    @property
    def radiance_offset(self):
        return self._v_radiance_offset()

RadianceScaling.print_desc = new_instancemethod(_radiance_scaling.RadianceScaling_print_desc,None,RadianceScaling)
RadianceScaling.apply_scaling = new_instancemethod(_radiance_scaling.RadianceScaling_apply_scaling,None,RadianceScaling)
RadianceScaling._v_radiance_scaling_coeff = new_instancemethod(_radiance_scaling.RadianceScaling__v_radiance_scaling_coeff,None,RadianceScaling)
RadianceScaling._v_radiance_scaling_coeff_uncertainty = new_instancemethod(_radiance_scaling.RadianceScaling__v_radiance_scaling_coeff_uncertainty,None,RadianceScaling)
RadianceScaling._v_radiance_offset = new_instancemethod(_radiance_scaling.RadianceScaling__v_radiance_offset,None,RadianceScaling)
RadianceScaling_swigregister = _radiance_scaling.RadianceScaling_swigregister
RadianceScaling_swigregister(RadianceScaling)



