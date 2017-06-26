# This file was automatically generated by SWIG (http://www.swig.org).
# Version 3.0.7
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.





from sys import version_info
if version_info >= (3, 0, 0):
    new_instancemethod = lambda func, inst, cls: _l_rad_driver.SWIG_PyInstanceMethod_New(func)
else:
    from new import instancemethod as new_instancemethod
if version_info >= (2, 6, 0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_l_rad_driver', [dirname(__file__)])
        except ImportError:
            import _l_rad_driver
            return _l_rad_driver
        if fp is not None:
            try:
                _mod = imp.load_module('_l_rad_driver', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _l_rad_driver = swig_import_helper()
    del swig_import_helper
else:
    import _l_rad_driver
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



_l_rad_driver.SHARED_PTR_DISOWN_swigconstant(_l_rad_driver)
SHARED_PTR_DISOWN = _l_rad_driver.SHARED_PTR_DISOWN

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
class LRadDriver(object):
    """

    This class drives the LRAD code, which gives a polarization correction
    to scalar intensity and jacobians.

    The correction used here is described in the paper "A fast linearized
    pseudo-spherical two orders of scattering model to account for
    polarization in vertically inhomogeneous scattering-absorbing media"
    by Vijah Natrah and Robert Spurr in Journal of Quantitative
    Spectroscopy & Radiative Transfer 107 (2007) 263-293 (a copy of this
    paper can be found in the source tree at 'doc/LRAD Paper.pdf').

    C++ includes: l_rad_driver.h 
    """
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    REGULAR = _l_rad_driver.LRadDriver_REGULAR
    ENHANCED = _l_rad_driver.LRadDriver_ENHANCED
    PLANE_PARALLEL = _l_rad_driver.LRadDriver_PLANE_PARALLEL
    DETECT = _l_rad_driver.LRadDriver_DETECT

    def __init__(self, *args):
        """

        LRadDriver::LRadDriver(int Number_stream, int Number_stokes, int surface_type, bool
        Tms_Correction=false, bool Pure_nadir=false, const PsMode
        ps_mode=DETECT)

        """
        _l_rad_driver.LRadDriver_swiginit(self, _l_rad_driver.new_LRadDriver(*args))

    def _v_number_stokes(self):
        """

        virtual int FullPhysics::LRadDriver::number_stokes() const

        """
        return _l_rad_driver.LRadDriver__v_number_stokes(self)


    @property
    def number_stokes(self):
        return self._v_number_stokes()


    def _v_number_stream(self):
        """

        virtual int FullPhysics::LRadDriver::number_stream() const

        """
        return _l_rad_driver.LRadDriver__v_number_stream(self)


    @property
    def number_stream(self):
        return self._v_number_stream()


    def z_matrix(self, pf):
        """

        ArrayAd< double, 2 > LRadDriver::z_matrix(const ArrayAd< double, 3 > &pf) const
        Calculate the z matrix.

        This does the full calculation, and is somewhat expensive to call.

        Must be called after setup_geometry since this calculation uses
        geometry values setup from it. 
        """
        return _l_rad_driver.LRadDriver_z_matrix(self, pf)


    def setup_geometry(self, alt, sza, zen, azm):
        """

        void LRadDriver::setup_geometry(blitz::Array< double, 1 > alt, double sza, double zen, double azm)
        const
        Setup viewing geometry, should only be called once per instance or if
        the viewing geometry changes.

        Update the altitude information.

        This can change the number of layers if desired. 
        """
        return _l_rad_driver.LRadDriver_setup_geometry(self, alt, sza, zen, azm)


    def setup_surface_params(self, surface_param):
        """

        void LRadDriver::setup_surface_params(const blitz::Array< double, 1 > &surface_param)
        Set up surface parameters for spectral point. 
        """
        return _l_rad_driver.LRadDriver_setup_surface_params(self, surface_param)


    def setup_optical_inputs(self, od, ssa, pf, zmat):
        """

        void LRadDriver::setup_optical_inputs(const blitz::Array< double, 1 > &od, const blitz::Array< double, 1 >
        &ssa, const blitz::Array< double, 3 > &pf, const blitz::Array< double,
        2 > &zmat)
        Set up optical depth, single scattering albedo and scattering matrix
        Should be called per spectral point. 
        """
        return _l_rad_driver.LRadDriver_setup_optical_inputs(self, od, ssa, pf, zmat)


    def clear_linear_inputs(self):
        """

        void LRadDriver::clear_linear_inputs()
        Mark that we are not retrieving weighting functions. 
        """
        return _l_rad_driver.LRadDriver_clear_linear_inputs(self)


    def setup_linear_inputs(self, od, ssa, pf, zmat):
        """

        void LRadDriver::setup_linear_inputs(const ArrayAd< double, 1 > &od, const ArrayAd< double, 1 > &ssa,
        const ArrayAd< double, 3 > &pf, const ArrayAd< double, 2 > &zmat)
        Set up linearization, weighting functions. 
        """
        return _l_rad_driver.LRadDriver_setup_linear_inputs(self, od, ssa, pf, zmat)


    def calculate_first_order(self):
        """

        void LRadDriver::calculate_first_order()
        Perform radiative transfer calculation with the values setup by
        setup_optical_inputs and setup_linear_inputs for the first order of
        scattering. 
        """
        return _l_rad_driver.LRadDriver_calculate_first_order(self)


    def calculate_second_order(self):
        """

        void LRadDriver::calculate_second_order()
        Perform radiative transfer calculation with the values setup by
        setup_optical_inputs and setup_linear_inputs for the second order of
        scattering. 
        """
        return _l_rad_driver.LRadDriver_calculate_second_order(self)


    def stokes(self):
        """

        Array< double, 1 > LRadDriver::stokes() const
        Retrieve the stokes values calculated. 
        """
        return _l_rad_driver.LRadDriver_stokes(self)


    def atmospheric_jacobian(self):
        """

        Array< double, 3 > LRadDriver::atmospheric_jacobian() const
        Atmospheric jacobian from last calculation. 
        """
        return _l_rad_driver.LRadDriver_atmospheric_jacobian(self)


    def surface_jacobian(self):
        """

        Array< double, 2 > LRadDriver::surface_jacobian() const
        Surface jacobian. 
        """
        return _l_rad_driver.LRadDriver_surface_jacobian(self)


    def print_desc(self, Os, Short_form=False):
        """

        void LRadDriver::print(std::ostream &Os, bool Short_form=false) const

        """
        return _l_rad_driver.LRadDriver_print_desc(self, Os, Short_form)

    __swig_destroy__ = _l_rad_driver.delete_LRadDriver
LRadDriver._v_number_stokes = new_instancemethod(_l_rad_driver.LRadDriver__v_number_stokes, None, LRadDriver)
LRadDriver._v_number_stream = new_instancemethod(_l_rad_driver.LRadDriver__v_number_stream, None, LRadDriver)
LRadDriver.z_matrix = new_instancemethod(_l_rad_driver.LRadDriver_z_matrix, None, LRadDriver)
LRadDriver.setup_geometry = new_instancemethod(_l_rad_driver.LRadDriver_setup_geometry, None, LRadDriver)
LRadDriver.setup_surface_params = new_instancemethod(_l_rad_driver.LRadDriver_setup_surface_params, None, LRadDriver)
LRadDriver.setup_optical_inputs = new_instancemethod(_l_rad_driver.LRadDriver_setup_optical_inputs, None, LRadDriver)
LRadDriver.clear_linear_inputs = new_instancemethod(_l_rad_driver.LRadDriver_clear_linear_inputs, None, LRadDriver)
LRadDriver.setup_linear_inputs = new_instancemethod(_l_rad_driver.LRadDriver_setup_linear_inputs, None, LRadDriver)
LRadDriver.calculate_first_order = new_instancemethod(_l_rad_driver.LRadDriver_calculate_first_order, None, LRadDriver)
LRadDriver.calculate_second_order = new_instancemethod(_l_rad_driver.LRadDriver_calculate_second_order, None, LRadDriver)
LRadDriver.stokes = new_instancemethod(_l_rad_driver.LRadDriver_stokes, None, LRadDriver)
LRadDriver.atmospheric_jacobian = new_instancemethod(_l_rad_driver.LRadDriver_atmospheric_jacobian, None, LRadDriver)
LRadDriver.surface_jacobian = new_instancemethod(_l_rad_driver.LRadDriver_surface_jacobian, None, LRadDriver)
LRadDriver.print_desc = new_instancemethod(_l_rad_driver.LRadDriver_print_desc, None, LRadDriver)
LRadDriver_swigregister = _l_rad_driver.LRadDriver_swigregister
LRadDriver_swigregister(LRadDriver)



