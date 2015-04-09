# This file was automatically generated by SWIG (http://www.swig.org).
# Version 2.0.9
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.



from sys import version_info
if version_info >= (3,0,0):
    new_instancemethod = lambda func, inst, cls: _l_rad_rt.SWIG_PyInstanceMethod_New(func)
else:
    from new import instancemethod as new_instancemethod
if version_info >= (2,6,0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_l_rad_rt', [dirname(__file__)])
        except ImportError:
            import _l_rad_rt
            return _l_rad_rt
        if fp is not None:
            try:
                _mod = imp.load_module('_l_rad_rt', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _l_rad_rt = swig_import_helper()
    del swig_import_helper
else:
    import _l_rad_rt
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


SHARED_PTR_DISOWN = _l_rad_rt.SHARED_PTR_DISOWN
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

import full_physics_swig.radiative_transfer_single_wn
import full_physics_swig.radiative_transfer_fixed_stokes_coefficient
import full_physics_swig.radiative_transfer
import full_physics_swig.generic_object
import full_physics_swig.observer
import full_physics_swig.named_spectrum
import full_physics_swig.state_vector
class LRadRt(full_physics_swig.radiative_transfer_single_wn.RadiativeTransferSingleWn):
    """
    This class drives the LRAD code, which gives a polarization correction
    to scalar intensity and jacobians.

    This can also be used on its own to provide a single scatter
    approximation to the RadiativeTransfer (i.e., without also running
    LIDORT).

    C++ includes: l_rad_rt.h 
    """
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    def __init__(self, *args, **kwargs): raise AttributeError("No constructor defined - class is abstract")
    __repr__ = _swig_repr
    @property
    def number_stokes(self):
        return self._v_number_stokes()

    @property
    def number_stream(self):
        return self._v_number_stream()

    def _v_surface_type(self):
        """
        virtual int FullPhysics::LRadRt::surface_type() const
        Returns an integer with l_rad's representation of surface type. 
        """
        return _l_rad_rt.LRadRt__v_surface_type(self)

    @property
    def surface_type(self):
        return self._v_surface_type()

    def _v_radiative_transfer(self):
        """
        const boost::shared_ptr<RadiativeTransferSingleWn>& FullPhysics::LRadRt::radiative_transfer() const

        """
        return _l_rad_rt.LRadRt__v_radiative_transfer(self)

    @property
    def radiative_transfer(self):
        return self._v_radiative_transfer()

    __swig_destroy__ = _l_rad_rt.delete_LRadRt
LRadRt._v_surface_type = new_instancemethod(_l_rad_rt.LRadRt__v_surface_type,None,LRadRt)
LRadRt._v_radiative_transfer = new_instancemethod(_l_rad_rt.LRadRt__v_radiative_transfer,None,LRadRt)
LRadRt_swigregister = _l_rad_rt.LRadRt_swigregister
LRadRt_swigregister(LRadRt)


