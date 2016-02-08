# This file was automatically generated by SWIG (http://www.swig.org).
# Version 3.0.7
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.





from sys import version_info
if version_info >= (3, 0, 0):
    new_instancemethod = lambda func, inst, cls: _merra_aerosol.SWIG_PyInstanceMethod_New(func)
else:
    from new import instancemethod as new_instancemethod
if version_info >= (2, 6, 0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_merra_aerosol', [dirname(__file__)])
        except ImportError:
            import _merra_aerosol
            return _merra_aerosol
        if fp is not None:
            try:
                _mod = imp.load_module('_merra_aerosol', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _merra_aerosol = swig_import_helper()
    del swig_import_helper
else:
    import _merra_aerosol
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



_merra_aerosol.SHARED_PTR_DISOWN_swigconstant(_merra_aerosol)
SHARED_PTR_DISOWN = _merra_aerosol.SHARED_PTR_DISOWN

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
import full_physics_swig.state_vector
import full_physics_swig.initial_guess
class MerraAerosol(full_physics_swig.generic_object.GenericObject):
    """

    This class is used to create the Aerosol from a Merra climatology
    file.

    C++ includes: merra_aerosol.h 
    """
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr

    def __init__(self, Merra_climatology, Aerosol_property, Latitude, Longitude, Press, Rh, Aerosol_cov, Max_aod=0.2, Exp_aod=0.8, Min_types=2, Max_types=4, Linear_aod=False, Max_residual=0.005, Reference_wn=1):
        """

        MerraAerosol::MerraAerosol(const HdfFile &Merra_climatology, const HdfFile &Aerosol_property,
        DoubleWithUnit Latitude, DoubleWithUnit Longitude, const
        boost::shared_ptr< Pressure > &Press, const boost::shared_ptr<
        RelativeHumidity > &Rh, const blitz::Array< double, 2 > &Aerosol_cov,
        double Max_aod=0.2, double Exp_aod=0.8, int Min_types=2, int
        Max_types=4, bool Linear_aod=false, double Max_residual=0.005, double
        Reference_wn=1e4/0.755)
        Constructor.

        Parameters:
        -----------

        Merra_climatology:  The Merra climatology file

        Aerosol_property:  The Aerosol property file

        Latitude:  The latitude of the ground point

        Longitude:  The longitude of the ground point

        Press:  The Pressure object that gives the pressure grid.

        Rh:  The RelativeHumidity object that gives the relative humidity.

        Aerosol_cov:  The covariance matrix to use for each Aerosol.

        Max_aod:  Maximum AOD cap for each composite type

        Exp_aod:  Threshold for explained fraction of AOD

        Min_types:  Minimum number of types to be selected

        Max_types:  Maximum number of types to be selected

        Linear_aod:  If true, then use linear aod rather than logarithmic aod.

        Max_residual:  Maximum resdidual in total AOD (either > threshold of
        fraction or < max residual will suffice

        Reference_wn:  The wavenumber that Aext is given for. This is
        optional, the default value matches the reference band given in the
        ATB. 
        """
        _merra_aerosol.MerraAerosol_swiginit(self, _merra_aerosol.new_MerraAerosol(Merra_climatology, Aerosol_property, Latitude, Longitude, Press, Rh, Aerosol_cov, Max_aod, Exp_aod, Min_types, Max_types, Linear_aod, Max_residual, Reference_wn))

    def _v_aerosol(self):
        """

        boost::shared_ptr< Aerosol > MerraAerosol::aerosol() const
        Return the aerosol setup generated from this class. 
        """
        return _merra_aerosol.MerraAerosol__v_aerosol(self)


    @property
    def aerosol(self):
        return self._v_aerosol()


    def _v_initial_guess(self):
        """

        boost::shared_ptr<InitialGuessBuilder> FullPhysics::MerraAerosol::initial_guess() const
        Return IntitialGuessBuilder for the Aerosol. 
        """
        return _merra_aerosol.MerraAerosol__v_initial_guess(self)


    @property
    def initial_guess(self):
        return self._v_initial_guess()


    def _v_number_merra_particle(self):
        """

        int FullPhysics::MerraAerosol::number_merra_particle() const
        Number of merra particles. 
        """
        return _merra_aerosol.MerraAerosol__v_number_merra_particle(self)


    @property
    def number_merra_particle(self):
        return self._v_number_merra_particle()

    __swig_destroy__ = _merra_aerosol.delete_MerraAerosol
MerraAerosol._v_aerosol = new_instancemethod(_merra_aerosol.MerraAerosol__v_aerosol, None, MerraAerosol)
MerraAerosol._v_initial_guess = new_instancemethod(_merra_aerosol.MerraAerosol__v_initial_guess, None, MerraAerosol)
MerraAerosol._v_number_merra_particle = new_instancemethod(_merra_aerosol.MerraAerosol__v_number_merra_particle, None, MerraAerosol)
MerraAerosol_swigregister = _merra_aerosol.MerraAerosol_swigregister
MerraAerosol_swigregister(MerraAerosol)



