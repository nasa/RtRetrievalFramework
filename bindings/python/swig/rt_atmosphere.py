# This file was automatically generated by SWIG (http://www.swig.org).
# Version 3.0.12
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.

from sys import version_info as _swig_python_version_info
if _swig_python_version_info >= (3, 0, 0):
    new_instancemethod = lambda func, inst, cls: _rt_atmosphere.SWIG_PyInstanceMethod_New(func)
else:
    from new import instancemethod as new_instancemethod
if _swig_python_version_info >= (2, 7, 0):
    def swig_import_helper():
        import importlib
        pkg = __name__.rpartition('.')[0]
        mname = '.'.join((pkg, '_rt_atmosphere')).lstrip('.')
        try:
            return importlib.import_module(mname)
        except ImportError:
            return importlib.import_module('_rt_atmosphere')
    _rt_atmosphere = swig_import_helper()
    del swig_import_helper
elif _swig_python_version_info >= (2, 6, 0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_rt_atmosphere', [dirname(__file__)])
        except ImportError:
            import _rt_atmosphere
            return _rt_atmosphere
        try:
            _mod = imp.load_module('_rt_atmosphere', fp, pathname, description)
        finally:
            if fp is not None:
                fp.close()
        return _mod
    _rt_atmosphere = swig_import_helper()
    del swig_import_helper
else:
    import _rt_atmosphere
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


SHARED_PTR_DISOWN = _rt_atmosphere.SHARED_PTR_DISOWN

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
class ObservableRtAtmosphere(full_physics_swig.generic_object.GenericObject):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')

    def __init__(self, *args, **kwargs):
        raise AttributeError("No constructor defined - class is abstract")
    __repr__ = _swig_repr
    __swig_destroy__ = _rt_atmosphere.delete_ObservableRtAtmosphere
ObservableRtAtmosphere.add_observer_and_keep_reference = new_instancemethod(_rt_atmosphere.ObservableRtAtmosphere_add_observer_and_keep_reference, None, ObservableRtAtmosphere)
ObservableRtAtmosphere.add_observer = new_instancemethod(_rt_atmosphere.ObservableRtAtmosphere_add_observer, None, ObservableRtAtmosphere)
ObservableRtAtmosphere.remove_observer = new_instancemethod(_rt_atmosphere.ObservableRtAtmosphere_remove_observer, None, ObservableRtAtmosphere)
ObservableRtAtmosphere_swigregister = _rt_atmosphere.ObservableRtAtmosphere_swigregister
ObservableRtAtmosphere_swigregister(ObservableRtAtmosphere)

class ObserverRtAtmosphere(full_physics_swig.generic_object.GenericObject):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr

    def __init__(self):
        _rt_atmosphere.ObserverRtAtmosphere_swiginit(self, _rt_atmosphere.new_ObserverRtAtmosphere())
    __swig_destroy__ = _rt_atmosphere.delete_ObserverRtAtmosphere
ObserverRtAtmosphere.notify_update = new_instancemethod(_rt_atmosphere.ObserverRtAtmosphere_notify_update, None, ObserverRtAtmosphere)
ObserverRtAtmosphere.notify_add = new_instancemethod(_rt_atmosphere.ObserverRtAtmosphere_notify_add, None, ObserverRtAtmosphere)
ObserverRtAtmosphere.notify_remove = new_instancemethod(_rt_atmosphere.ObserverRtAtmosphere_notify_remove, None, ObserverRtAtmosphere)
ObserverRtAtmosphere_swigregister = _rt_atmosphere.ObserverRtAtmosphere_swigregister
ObserverRtAtmosphere_swigregister(ObserverRtAtmosphere)

class RtAtmosphere(full_physics_swig.state_vector.StateVectorObserver, ObservableRtAtmosphere):
    """

    This class is responsible for setting up the atmosphere and ground
    information needed to run the Radiative transfer code.

    There are many, many properties associated with the atmosphere. This
    class is not meant to model these properties, it is really the very
    limited information needed to run the Radiative transfer code.

    Note that this includes both the atmosphere and surface parameters
    needed by the RT code.

    The calculation of the Jacobians in LIDORT takes a time directly
    proportional to the number of variables we are taking the Jacobian
    with respect to, we use an "intermediate" set of variables for some
    of the reported gradients (e.g., AtmosphereOco uses taur, taug, and
    tau for each of the aerosol). To support future Atmosphere classes, we
    are purposely vague on exactly what these intermediate variables are,
    at least through the RtAtmosphere interface. The
    "intermediate_variable" function can be used to get the value of
    these intermediate variables and Jacobian with the state vector
    variables.

    A description of the intermediate variables can be found in
    doc/LIDORT_Jacobian.pdf.

    Note that it is assumed by the LSI that averaging these intermediate
    variables to get average optical properties makes sense. This is true
    if the variables are taur etc., but might not be true in general. If
    we add a class derived from RtAtmosphere where this doesn't make
    sense, we will need to rework this interface.

    Other objects may depend on the RtAtmosphere, and should be updated
    when the RtAtmosphere is updated. To facilitate that, this class is an
    Oberverable, and objects can add themselves as Observers to be
    notified when the RtAtmosphere is updated.

    Because the absorber calculation tends to be a bottle neck, we keep a
    timer in this class. This class keeps track of the time used in the
    atmosphere calculations. Other classes can make use of this
    information for logging if desired.

    C++ includes: rt_atmosphere.h 
    """

    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')

    def __init__(self, *args, **kwargs):
        raise AttributeError("No constructor defined - class is abstract")
    __repr__ = _swig_repr
    __swig_destroy__ = _rt_atmosphere.delete_RtAtmosphere

    def _v_timer_info(self):
        """

        std::string RtAtmosphere::timer_info() const
        Return timer information. 
        """
        return _rt_atmosphere.RtAtmosphere__v_timer_info(self)


    @property
    def timer_info(self):
        return self._v_timer_info()


    def _v_number_layer(self):
        """

        virtual int FullPhysics::RtAtmosphere::number_layer() const =0
        Number of layers we currently have. 
        """
        return _rt_atmosphere.RtAtmosphere__v_number_layer(self)


    @property
    def number_layer(self):
        return self._v_number_layer()


    def _v_number_spectrometer(self):
        """

        virtual int FullPhysics::RtAtmosphere::number_spectrometer() const =0
        Number of spectrometers we have. 
        """
        return _rt_atmosphere.RtAtmosphere__v_number_spectrometer(self)


    @property
    def number_spectrometer(self):
        return self._v_number_spectrometer()


    def altitude(self, spec_index):
        """

        virtual ArrayAdWithUnit<double, 1> FullPhysics::RtAtmosphere::altitude(int spec_index) const =0
        Altitude grid for current pressure grid. 
        """
        return _rt_atmosphere.RtAtmosphere_altitude(self, spec_index)


    def column_optical_depth(self, wn, spec_index, Gas_name):
        """

        virtual AutoDerivative<double> FullPhysics::RtAtmosphere::column_optical_depth(double wn, int spec_index, const std::string &Gas_name) const =0
        Total column optical depth for the given gas.

        This is 0 if the band isn't one that sees that gas. 
        """
        return _rt_atmosphere.RtAtmosphere_column_optical_depth(self, wn, spec_index, Gas_name)


    def optical_depth_wrt_iv(self, *args):
        """

        virtual ArrayAd<double, 1> FullPhysics::RtAtmosphere::optical_depth_wrt_iv(double wn, int spec_index) const =0
        The optical depth for each layer, for the given wave number.

        The derivatives of the optical depth are with respect to the
        intermediate variables, rather than the state vector variables (see
        description of RtAtmosphere class for discussion of this).

        Parameters:
        -----------

        wn:  The wave number to calculate parameters for.

        spec_index:  The spectrometer index

        Optical depth for each layer. This is number_layer() in size 
        """
        return _rt_atmosphere.RtAtmosphere_optical_depth_wrt_iv(self, *args)


    def single_scattering_albedo_wrt_iv(self, *args):
        """

        virtual ArrayAd<double, 1> FullPhysics::RtAtmosphere::single_scattering_albedo_wrt_iv(double wn, int spec_index, const ArrayAd< double, 2 > &iv) const =0
        This is a variation of single_scattering_albedo that takes the
        supplied value for the intermediate variables rather than calculating
        it own value.

        This is used by the LSI to get "average optical properties".

        Parameters:
        -----------

        wn:  The wave number to calculate parameters for.

        spec_index:  The spectrometer index

        iv:  Intermediate variable values to use.

        Single scattering albedo for each layer. This is number_layer() in
        size 
        """
        return _rt_atmosphere.RtAtmosphere_single_scattering_albedo_wrt_iv(self, *args)


    def scattering_moment_wrt_iv(self, *args):
        """

        virtual ArrayAd<double, 3> FullPhysics::RtAtmosphere::scattering_moment_wrt_iv(double wn, int spec_index, const ArrayAd< double, 2 > &iv, int
        nummom=-1, int numscat=-1) const =0
        This is a variation of scattering_moment that takes the supplied value
        for the intermediate variables rather than calculating it own value.

        This is used by the LSI to get "average optical properties".

        Parameters:
        -----------

        wn:  The wave number to calculate parameters for.

        spec_index:  The spectrometer index

        iv:  Intermediate variable values to use.

        nummom:  Number of moments to include in scatt_mom_each_layer, the
        default it to include all of them.

        numscat:  Number of scattering matrix elements to include in
        scatt_mom_each_layer, the default it to include all of them.

        Scattering moments for each layer. This is number_moment + 1 x
        number_layer() x number scattering matrix elements 
        """
        return _rt_atmosphere.RtAtmosphere_scattering_moment_wrt_iv(self, *args)


    def optical_depth_wrt_state_vector(self, wn, spec_index):
        """

        ArrayAd< double, 1 > RtAtmosphere::optical_depth_wrt_state_vector(double wn, int spec_index) const
        The optical depth for each layer, for the given wave number.

        This variation gives the derivatives with respect to the state vector,
        this just combines optical_depth with the Jacobian of the intermediate
        variables given by intermediate_variable.

        Parameters:
        -----------

        wn:  The wave number to calculate parameters for.

        spec_index:  The spectrometer index

        Optical depth for each layer. This is number_layer() in size 
        """
        return _rt_atmosphere.RtAtmosphere_optical_depth_wrt_state_vector(self, wn, spec_index)


    def single_scattering_albedo_wrt_state_vector(self, wn, spec_index):
        """

        ArrayAd< double, 1 > RtAtmosphere::single_scattering_albedo_wrt_state_vector(double wn, int spec_index) const
        The single scattering albedo for each layer, for the given wave
        number.

        This variation gives the derivatives with respect to the state vector,
        this just combines single_scattering_albedo with the Jacobian of the
        intermediate variables given by intermediate_variable.

        Parameters:
        -----------

        wn:  The wave number to calculate parameters for.

        spec_index:  The spectrometer index

        Single scattering albedo for each layer. This is number_layer() in
        size 
        """
        return _rt_atmosphere.RtAtmosphere_single_scattering_albedo_wrt_state_vector(self, wn, spec_index)


    def scattering_moment_wrt_state_vector(self, wn, spec_index, nummom=-1, numscat=-1):
        """

        ArrayAd< double, 3 > RtAtmosphere::scattering_moment_wrt_state_vector(double wn, int spec_index, int nummom=-1, int numscat=-1) const
        The scattering moments for for each layer, for the given wave number.

        The scattering moments use the de Rooij convention for the 6
        scattering matrix element.

        This variation gives the derivatives with respect to the state vector,
        this just combines single_scattering_albedo with the Jacobian of the
        intermediate variables given by intermediate_variable.

        Parameters:
        -----------

        wn:  The wave number to calculate parameters for.

        spec_index:  The spectrometer index

        nummom:  Number of moments to include in scatt_mom_each_layer, the
        default it to include all of them.

        numscat:  Number of scattering matrix elements to include in
        scatt_mom_each_layer, the default it to include all of them.

        Scattering moments for each layer. This is number_moment + 1 x
        number_layer() x number scattering matrix elements 
        """
        return _rt_atmosphere.RtAtmosphere_scattering_moment_wrt_state_vector(self, wn, spec_index, nummom, numscat)


    def intermediate_variable(self, wn, spec_index):
        """

        virtual ArrayAd<double, 2> FullPhysics::RtAtmosphere::intermediate_variable(double wn, int spec_index) const =0
        This gives the values of the intermediate variables and the Jacobian
        with respect to the state vector.

        This is number_layer() x number variables 
        """
        return _rt_atmosphere.RtAtmosphere_intermediate_variable(self, wn, spec_index)


    def _v_ground(self):
        """

        virtual const boost::shared_ptr<Ground> FullPhysics::RtAtmosphere::ground() const =0
        Object that represents the ground surface.

        If null then there is no surface for this atmosphere. 
        """
        return _rt_atmosphere.RtAtmosphere__v_ground(self)


    @property
    def ground(self):
        return self._v_ground()


    def _v_uplooking(self):
        """

        virtual bool FullPhysics::RtAtmosphere::uplooking() const =0
        Return true if we have an atmosphere for uplooking mode, i.e., we
        don't have a ground defined. 
        """
        return _rt_atmosphere.RtAtmosphere__v_uplooking(self)


    @property
    def uplooking(self):
        return self._v_uplooking()


    def reset_timer(self):
        """

        virtual void FullPhysics::RtAtmosphere::reset_timer()
        Reset timer. 
        """
        return _rt_atmosphere.RtAtmosphere_reset_timer(self)

RtAtmosphere._v_timer_info = new_instancemethod(_rt_atmosphere.RtAtmosphere__v_timer_info, None, RtAtmosphere)
RtAtmosphere._v_number_layer = new_instancemethod(_rt_atmosphere.RtAtmosphere__v_number_layer, None, RtAtmosphere)
RtAtmosphere._v_number_spectrometer = new_instancemethod(_rt_atmosphere.RtAtmosphere__v_number_spectrometer, None, RtAtmosphere)
RtAtmosphere.altitude = new_instancemethod(_rt_atmosphere.RtAtmosphere_altitude, None, RtAtmosphere)
RtAtmosphere.column_optical_depth = new_instancemethod(_rt_atmosphere.RtAtmosphere_column_optical_depth, None, RtAtmosphere)
RtAtmosphere.optical_depth_wrt_iv = new_instancemethod(_rt_atmosphere.RtAtmosphere_optical_depth_wrt_iv, None, RtAtmosphere)
RtAtmosphere.single_scattering_albedo_wrt_iv = new_instancemethod(_rt_atmosphere.RtAtmosphere_single_scattering_albedo_wrt_iv, None, RtAtmosphere)
RtAtmosphere.scattering_moment_wrt_iv = new_instancemethod(_rt_atmosphere.RtAtmosphere_scattering_moment_wrt_iv, None, RtAtmosphere)
RtAtmosphere.optical_depth_wrt_state_vector = new_instancemethod(_rt_atmosphere.RtAtmosphere_optical_depth_wrt_state_vector, None, RtAtmosphere)
RtAtmosphere.single_scattering_albedo_wrt_state_vector = new_instancemethod(_rt_atmosphere.RtAtmosphere_single_scattering_albedo_wrt_state_vector, None, RtAtmosphere)
RtAtmosphere.scattering_moment_wrt_state_vector = new_instancemethod(_rt_atmosphere.RtAtmosphere_scattering_moment_wrt_state_vector, None, RtAtmosphere)
RtAtmosphere.intermediate_variable = new_instancemethod(_rt_atmosphere.RtAtmosphere_intermediate_variable, None, RtAtmosphere)
RtAtmosphere._v_ground = new_instancemethod(_rt_atmosphere.RtAtmosphere__v_ground, None, RtAtmosphere)
RtAtmosphere._v_uplooking = new_instancemethod(_rt_atmosphere.RtAtmosphere__v_uplooking, None, RtAtmosphere)
RtAtmosphere.reset_timer = new_instancemethod(_rt_atmosphere.RtAtmosphere_reset_timer, None, RtAtmosphere)
RtAtmosphere.__str__ = new_instancemethod(_rt_atmosphere.RtAtmosphere___str__, None, RtAtmosphere)
RtAtmosphere_swigregister = _rt_atmosphere.RtAtmosphere_swigregister
RtAtmosphere_swigregister(RtAtmosphere)



