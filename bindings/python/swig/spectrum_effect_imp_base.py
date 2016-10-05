# This file was automatically generated by SWIG (http://www.swig.org).
# Version 2.0.9
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.



from sys import version_info
if version_info >= (3,0,0):
    new_instancemethod = lambda func, inst, cls: _spectrum_effect_imp_base.SWIG_PyInstanceMethod_New(func)
else:
    from new import instancemethod as new_instancemethod
if version_info >= (2,6,0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_spectrum_effect_imp_base', [dirname(__file__)])
        except ImportError:
            import _spectrum_effect_imp_base
            return _spectrum_effect_imp_base
        if fp is not None:
            try:
                _mod = imp.load_module('_spectrum_effect_imp_base', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _spectrum_effect_imp_base = swig_import_helper()
    del swig_import_helper
else:
    import _spectrum_effect_imp_base
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


SHARED_PTR_DISOWN = _spectrum_effect_imp_base.SHARED_PTR_DISOWN
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

import full_physics_swig.sub_state_vector_array
import full_physics_swig.generic_object
import full_physics_swig.spectrum_effect
class SpectrumEffectImpBase(full_physics_swig.spectrum_effect.SubStateVectorArraySpectrumEffect):
    """
    As a design principle, we have each base class with the absolutely
    minimum interface needed for use from the rest of the system.

    This allows us to support any future code that supports this minimum
    interface.

    However, almost always you will want to derive from this class
    instead. See PressureImpBase for a more complete discussion of this.

    C++ includes: spectrum_effect_imp_base.h 
    """
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    __swig_destroy__ = _spectrum_effect_imp_base.delete_SpectrumEffectImpBase
    def clone(self):
        """
        virtual boost::shared_ptr<SpectrumEffect> FullPhysics::SpectrumEffectImpBase::clone() const =0

        """
        return _spectrum_effect_imp_base.SpectrumEffectImpBase_clone(self)

    def print_desc(self, *args):
        """
        virtual void FullPhysics::SpectrumEffectImpBase::print(std::ostream &Os, bool Short_form=false) const
        Print to stream.

        The default calls the function "desc" that returns a string. This
        gives cleaner interface for deriving from this class in python, but
        most C++ classes will want to override this function rather than using
        desc. 
        """
        return _spectrum_effect_imp_base.SpectrumEffectImpBase_print_desc(self, *args)

    def _v_desc(self):
        """
        virtual std::string FullPhysics::SpectrumEffectImpBase::desc() const
        Description of object, to be printed to stream.

        This gives a cleaner interface for deriving from python. 
        """
        return _spectrum_effect_imp_base.SpectrumEffectImpBase__v_desc(self)

    @property
    def desc(self):
        return self._v_desc()

    def __init__(self, *args): 
        if self.__class__ == SpectrumEffectImpBase:
            _self = None
        else:
            _self = self
        _spectrum_effect_imp_base.SpectrumEffectImpBase_swiginit(self,_spectrum_effect_imp_base.new_SpectrumEffectImpBase(_self, *args))
    def __disown__(self):
        self.this.disown()
        _spectrum_effect_imp_base.disown_SpectrumEffectImpBase(self)
        return weakref_proxy(self)
SpectrumEffectImpBase.clone = new_instancemethod(_spectrum_effect_imp_base.SpectrumEffectImpBase_clone,None,SpectrumEffectImpBase)
SpectrumEffectImpBase.apply_effect = new_instancemethod(_spectrum_effect_imp_base.SpectrumEffectImpBase_apply_effect,None,SpectrumEffectImpBase)
SpectrumEffectImpBase.add_observer = new_instancemethod(_spectrum_effect_imp_base.SpectrumEffectImpBase_add_observer,None,SpectrumEffectImpBase)
SpectrumEffectImpBase.remove_observer = new_instancemethod(_spectrum_effect_imp_base.SpectrumEffectImpBase_remove_observer,None,SpectrumEffectImpBase)
SpectrumEffectImpBase.update_sub_state_hook = new_instancemethod(_spectrum_effect_imp_base.SpectrumEffectImpBase_update_sub_state_hook,None,SpectrumEffectImpBase)
SpectrumEffectImpBase.print_desc = new_instancemethod(_spectrum_effect_imp_base.SpectrumEffectImpBase_print_desc,None,SpectrumEffectImpBase)
SpectrumEffectImpBase._v_desc = new_instancemethod(_spectrum_effect_imp_base.SpectrumEffectImpBase__v_desc,None,SpectrumEffectImpBase)
SpectrumEffectImpBase.mark_used = new_instancemethod(_spectrum_effect_imp_base.SpectrumEffectImpBase_mark_used,None,SpectrumEffectImpBase)
SpectrumEffectImpBase.state_vector_name = new_instancemethod(_spectrum_effect_imp_base.SpectrumEffectImpBase_state_vector_name,None,SpectrumEffectImpBase)
SpectrumEffectImpBase.notify_update = new_instancemethod(_spectrum_effect_imp_base.SpectrumEffectImpBase_notify_update,None,SpectrumEffectImpBase)
SpectrumEffectImpBase.notify_add = new_instancemethod(_spectrum_effect_imp_base.SpectrumEffectImpBase_notify_add,None,SpectrumEffectImpBase)
SpectrumEffectImpBase.notify_remove = new_instancemethod(_spectrum_effect_imp_base.SpectrumEffectImpBase_notify_remove,None,SpectrumEffectImpBase)
SpectrumEffectImpBase.update_sub_state = new_instancemethod(_spectrum_effect_imp_base.SpectrumEffectImpBase_update_sub_state,None,SpectrumEffectImpBase)
SpectrumEffectImpBase.state_vector_name_i = new_instancemethod(_spectrum_effect_imp_base.SpectrumEffectImpBase_state_vector_name_i,None,SpectrumEffectImpBase)
SpectrumEffectImpBase.state_vector_name_sub = new_instancemethod(_spectrum_effect_imp_base.SpectrumEffectImpBase_state_vector_name_sub,None,SpectrumEffectImpBase)
SpectrumEffectImpBase_swigregister = _spectrum_effect_imp_base.SpectrumEffectImpBase_swigregister
SpectrumEffectImpBase_swigregister(SpectrumEffectImpBase)



