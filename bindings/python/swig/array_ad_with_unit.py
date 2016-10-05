# This file was automatically generated by SWIG (http://www.swig.org).
# Version 2.0.9
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.



from sys import version_info
if version_info >= (3,0,0):
    new_instancemethod = lambda func, inst, cls: _array_ad_with_unit.SWIG_PyInstanceMethod_New(func)
else:
    from new import instancemethod as new_instancemethod
if version_info >= (2,6,0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_array_ad_with_unit', [dirname(__file__)])
        except ImportError:
            import _array_ad_with_unit
            return _array_ad_with_unit
        if fp is not None:
            try:
                _mod = imp.load_module('_array_ad_with_unit', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _array_ad_with_unit = swig_import_helper()
    del swig_import_helper
else:
    import _array_ad_with_unit
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


SHARED_PTR_DISOWN = _array_ad_with_unit.SHARED_PTR_DISOWN
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
class ArrayAdWithUnitDouble_1(object):
    """
    This is a ArrayAd that also has units associated with it.

    This is a simple structure that just keeps these two things together.

    C++ includes: array_ad_with_unit.h 
    """
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        """
        FullPhysics::ArrayAdWithUnit< T, D >::ArrayAdWithUnit(const ArrayAd< T, D > &V)

        """
        _array_ad_with_unit.ArrayAdWithUnitDouble_1_swiginit(self,_array_ad_with_unit.new_ArrayAdWithUnitDouble_1(*args))
    def convert(self, *args):
        """
        ArrayAdWithUnit<T, D> FullPhysics::ArrayAdWithUnit< T, D >::convert(const Unit &R) const
        Convert to the given units. 
        """
        return _array_ad_with_unit.ArrayAdWithUnitDouble_1_convert(self, *args)

    def _v_rows(self):
        """
        int FullPhysics::ArrayAdWithUnit< T, D >::rows() const

        """
        return _array_ad_with_unit.ArrayAdWithUnitDouble_1__v_rows(self)

    @property
    def rows(self):
        return self._v_rows()

    def _v_cols(self):
        """
        int FullPhysics::ArrayAdWithUnit< T, D >::cols() const

        """
        return _array_ad_with_unit.ArrayAdWithUnitDouble_1__v_cols(self)

    @property
    def cols(self):
        return self._v_cols()

    def _v_depth(self):
        """
        int FullPhysics::ArrayAdWithUnit< T, D >::depth() const

        """
        return _array_ad_with_unit.ArrayAdWithUnitDouble_1__v_depth(self)

    @property
    def depth(self):
        return self._v_depth()

    def _v_is_constant(self):
        """
        bool FullPhysics::ArrayAdWithUnit< T, D >::is_constant() const

        """
        return _array_ad_with_unit.ArrayAdWithUnitDouble_1__v_is_constant(self)

    @property
    def is_constant(self):
        return self._v_is_constant()

    def _v_number_variable(self):
        """
        int FullPhysics::ArrayAdWithUnit< T, D >::number_variable() const

        """
        return _array_ad_with_unit.ArrayAdWithUnitDouble_1__v_number_variable(self)

    @property
    def number_variable(self):
        return self._v_number_variable()

    def reference(self, *args):
        """
        void FullPhysics::ArrayAdWithUnit< T, D >::reference(const ArrayAdWithUnit< T, D > &V)

        """
        return _array_ad_with_unit.ArrayAdWithUnitDouble_1_reference(self, *args)

    @property
    def value(self):
      return self._value()

    @value.setter
    def value(self, val):
      self._value_set(val)

    @property
    def units(self):
      return self._units()

    @units.setter
    def units(self,val):
        self._units_set(val)
      
    __swig_destroy__ = _array_ad_with_unit.delete_ArrayAdWithUnitDouble_1
ArrayAdWithUnitDouble_1.convert = new_instancemethod(_array_ad_with_unit.ArrayAdWithUnitDouble_1_convert,None,ArrayAdWithUnitDouble_1)
ArrayAdWithUnitDouble_1._v_rows = new_instancemethod(_array_ad_with_unit.ArrayAdWithUnitDouble_1__v_rows,None,ArrayAdWithUnitDouble_1)
ArrayAdWithUnitDouble_1._v_cols = new_instancemethod(_array_ad_with_unit.ArrayAdWithUnitDouble_1__v_cols,None,ArrayAdWithUnitDouble_1)
ArrayAdWithUnitDouble_1._v_depth = new_instancemethod(_array_ad_with_unit.ArrayAdWithUnitDouble_1__v_depth,None,ArrayAdWithUnitDouble_1)
ArrayAdWithUnitDouble_1._v_is_constant = new_instancemethod(_array_ad_with_unit.ArrayAdWithUnitDouble_1__v_is_constant,None,ArrayAdWithUnitDouble_1)
ArrayAdWithUnitDouble_1._v_number_variable = new_instancemethod(_array_ad_with_unit.ArrayAdWithUnitDouble_1__v_number_variable,None,ArrayAdWithUnitDouble_1)
ArrayAdWithUnitDouble_1.reference = new_instancemethod(_array_ad_with_unit.ArrayAdWithUnitDouble_1_reference,None,ArrayAdWithUnitDouble_1)
ArrayAdWithUnitDouble_1._value = new_instancemethod(_array_ad_with_unit.ArrayAdWithUnitDouble_1__value,None,ArrayAdWithUnitDouble_1)
ArrayAdWithUnitDouble_1._value_set = new_instancemethod(_array_ad_with_unit.ArrayAdWithUnitDouble_1__value_set,None,ArrayAdWithUnitDouble_1)
ArrayAdWithUnitDouble_1._units = new_instancemethod(_array_ad_with_unit.ArrayAdWithUnitDouble_1__units,None,ArrayAdWithUnitDouble_1)
ArrayAdWithUnitDouble_1._units_set = new_instancemethod(_array_ad_with_unit.ArrayAdWithUnitDouble_1__units_set,None,ArrayAdWithUnitDouble_1)
ArrayAdWithUnitDouble_1_swigregister = _array_ad_with_unit.ArrayAdWithUnitDouble_1_swigregister
ArrayAdWithUnitDouble_1_swigregister(ArrayAdWithUnitDouble_1)

class ArrayAdWithUnitDouble_2(object):
    """
    This is a ArrayAd that also has units associated with it.

    This is a simple structure that just keeps these two things together.

    C++ includes: array_ad_with_unit.h 
    """
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        """
        FullPhysics::ArrayAdWithUnit< T, D >::ArrayAdWithUnit(const ArrayAd< T, D > &V)

        """
        _array_ad_with_unit.ArrayAdWithUnitDouble_2_swiginit(self,_array_ad_with_unit.new_ArrayAdWithUnitDouble_2(*args))
    def convert(self, *args):
        """
        ArrayAdWithUnit<T, D> FullPhysics::ArrayAdWithUnit< T, D >::convert(const Unit &R) const
        Convert to the given units. 
        """
        return _array_ad_with_unit.ArrayAdWithUnitDouble_2_convert(self, *args)

    def _v_rows(self):
        """
        int FullPhysics::ArrayAdWithUnit< T, D >::rows() const

        """
        return _array_ad_with_unit.ArrayAdWithUnitDouble_2__v_rows(self)

    @property
    def rows(self):
        return self._v_rows()

    def _v_cols(self):
        """
        int FullPhysics::ArrayAdWithUnit< T, D >::cols() const

        """
        return _array_ad_with_unit.ArrayAdWithUnitDouble_2__v_cols(self)

    @property
    def cols(self):
        return self._v_cols()

    def _v_depth(self):
        """
        int FullPhysics::ArrayAdWithUnit< T, D >::depth() const

        """
        return _array_ad_with_unit.ArrayAdWithUnitDouble_2__v_depth(self)

    @property
    def depth(self):
        return self._v_depth()

    def _v_is_constant(self):
        """
        bool FullPhysics::ArrayAdWithUnit< T, D >::is_constant() const

        """
        return _array_ad_with_unit.ArrayAdWithUnitDouble_2__v_is_constant(self)

    @property
    def is_constant(self):
        return self._v_is_constant()

    def _v_number_variable(self):
        """
        int FullPhysics::ArrayAdWithUnit< T, D >::number_variable() const

        """
        return _array_ad_with_unit.ArrayAdWithUnitDouble_2__v_number_variable(self)

    @property
    def number_variable(self):
        return self._v_number_variable()

    def reference(self, *args):
        """
        void FullPhysics::ArrayAdWithUnit< T, D >::reference(const ArrayAdWithUnit< T, D > &V)

        """
        return _array_ad_with_unit.ArrayAdWithUnitDouble_2_reference(self, *args)

    @property
    def value(self):
      return self._value()

    @value.setter
    def value(self, val):
      self._value_set(val)

    @property
    def units(self):
      return self._units()

    @units.setter
    def units(self,val):
        self._units_set(val)
      
    __swig_destroy__ = _array_ad_with_unit.delete_ArrayAdWithUnitDouble_2
ArrayAdWithUnitDouble_2.convert = new_instancemethod(_array_ad_with_unit.ArrayAdWithUnitDouble_2_convert,None,ArrayAdWithUnitDouble_2)
ArrayAdWithUnitDouble_2._v_rows = new_instancemethod(_array_ad_with_unit.ArrayAdWithUnitDouble_2__v_rows,None,ArrayAdWithUnitDouble_2)
ArrayAdWithUnitDouble_2._v_cols = new_instancemethod(_array_ad_with_unit.ArrayAdWithUnitDouble_2__v_cols,None,ArrayAdWithUnitDouble_2)
ArrayAdWithUnitDouble_2._v_depth = new_instancemethod(_array_ad_with_unit.ArrayAdWithUnitDouble_2__v_depth,None,ArrayAdWithUnitDouble_2)
ArrayAdWithUnitDouble_2._v_is_constant = new_instancemethod(_array_ad_with_unit.ArrayAdWithUnitDouble_2__v_is_constant,None,ArrayAdWithUnitDouble_2)
ArrayAdWithUnitDouble_2._v_number_variable = new_instancemethod(_array_ad_with_unit.ArrayAdWithUnitDouble_2__v_number_variable,None,ArrayAdWithUnitDouble_2)
ArrayAdWithUnitDouble_2.reference = new_instancemethod(_array_ad_with_unit.ArrayAdWithUnitDouble_2_reference,None,ArrayAdWithUnitDouble_2)
ArrayAdWithUnitDouble_2._value = new_instancemethod(_array_ad_with_unit.ArrayAdWithUnitDouble_2__value,None,ArrayAdWithUnitDouble_2)
ArrayAdWithUnitDouble_2._value_set = new_instancemethod(_array_ad_with_unit.ArrayAdWithUnitDouble_2__value_set,None,ArrayAdWithUnitDouble_2)
ArrayAdWithUnitDouble_2._units = new_instancemethod(_array_ad_with_unit.ArrayAdWithUnitDouble_2__units,None,ArrayAdWithUnitDouble_2)
ArrayAdWithUnitDouble_2._units_set = new_instancemethod(_array_ad_with_unit.ArrayAdWithUnitDouble_2__units_set,None,ArrayAdWithUnitDouble_2)
ArrayAdWithUnitDouble_2_swigregister = _array_ad_with_unit.ArrayAdWithUnitDouble_2_swigregister
ArrayAdWithUnitDouble_2_swigregister(ArrayAdWithUnitDouble_2)

class ArrayAdWithUnitDouble_3(object):
    """
    This is a ArrayAd that also has units associated with it.

    This is a simple structure that just keeps these two things together.

    C++ includes: array_ad_with_unit.h 
    """
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        """
        FullPhysics::ArrayAdWithUnit< T, D >::ArrayAdWithUnit(const ArrayAd< T, D > &V)

        """
        _array_ad_with_unit.ArrayAdWithUnitDouble_3_swiginit(self,_array_ad_with_unit.new_ArrayAdWithUnitDouble_3(*args))
    def convert(self, *args):
        """
        ArrayAdWithUnit<T, D> FullPhysics::ArrayAdWithUnit< T, D >::convert(const Unit &R) const
        Convert to the given units. 
        """
        return _array_ad_with_unit.ArrayAdWithUnitDouble_3_convert(self, *args)

    def _v_rows(self):
        """
        int FullPhysics::ArrayAdWithUnit< T, D >::rows() const

        """
        return _array_ad_with_unit.ArrayAdWithUnitDouble_3__v_rows(self)

    @property
    def rows(self):
        return self._v_rows()

    def _v_cols(self):
        """
        int FullPhysics::ArrayAdWithUnit< T, D >::cols() const

        """
        return _array_ad_with_unit.ArrayAdWithUnitDouble_3__v_cols(self)

    @property
    def cols(self):
        return self._v_cols()

    def _v_depth(self):
        """
        int FullPhysics::ArrayAdWithUnit< T, D >::depth() const

        """
        return _array_ad_with_unit.ArrayAdWithUnitDouble_3__v_depth(self)

    @property
    def depth(self):
        return self._v_depth()

    def _v_is_constant(self):
        """
        bool FullPhysics::ArrayAdWithUnit< T, D >::is_constant() const

        """
        return _array_ad_with_unit.ArrayAdWithUnitDouble_3__v_is_constant(self)

    @property
    def is_constant(self):
        return self._v_is_constant()

    def _v_number_variable(self):
        """
        int FullPhysics::ArrayAdWithUnit< T, D >::number_variable() const

        """
        return _array_ad_with_unit.ArrayAdWithUnitDouble_3__v_number_variable(self)

    @property
    def number_variable(self):
        return self._v_number_variable()

    def reference(self, *args):
        """
        void FullPhysics::ArrayAdWithUnit< T, D >::reference(const ArrayAdWithUnit< T, D > &V)

        """
        return _array_ad_with_unit.ArrayAdWithUnitDouble_3_reference(self, *args)

    @property
    def value(self):
      return self._value()

    @value.setter
    def value(self, val):
      self._value_set(val)

    @property
    def units(self):
      return self._units()

    @units.setter
    def units(self,val):
        self._units_set(val)
      
    __swig_destroy__ = _array_ad_with_unit.delete_ArrayAdWithUnitDouble_3
ArrayAdWithUnitDouble_3.convert = new_instancemethod(_array_ad_with_unit.ArrayAdWithUnitDouble_3_convert,None,ArrayAdWithUnitDouble_3)
ArrayAdWithUnitDouble_3._v_rows = new_instancemethod(_array_ad_with_unit.ArrayAdWithUnitDouble_3__v_rows,None,ArrayAdWithUnitDouble_3)
ArrayAdWithUnitDouble_3._v_cols = new_instancemethod(_array_ad_with_unit.ArrayAdWithUnitDouble_3__v_cols,None,ArrayAdWithUnitDouble_3)
ArrayAdWithUnitDouble_3._v_depth = new_instancemethod(_array_ad_with_unit.ArrayAdWithUnitDouble_3__v_depth,None,ArrayAdWithUnitDouble_3)
ArrayAdWithUnitDouble_3._v_is_constant = new_instancemethod(_array_ad_with_unit.ArrayAdWithUnitDouble_3__v_is_constant,None,ArrayAdWithUnitDouble_3)
ArrayAdWithUnitDouble_3._v_number_variable = new_instancemethod(_array_ad_with_unit.ArrayAdWithUnitDouble_3__v_number_variable,None,ArrayAdWithUnitDouble_3)
ArrayAdWithUnitDouble_3.reference = new_instancemethod(_array_ad_with_unit.ArrayAdWithUnitDouble_3_reference,None,ArrayAdWithUnitDouble_3)
ArrayAdWithUnitDouble_3._value = new_instancemethod(_array_ad_with_unit.ArrayAdWithUnitDouble_3__value,None,ArrayAdWithUnitDouble_3)
ArrayAdWithUnitDouble_3._value_set = new_instancemethod(_array_ad_with_unit.ArrayAdWithUnitDouble_3__value_set,None,ArrayAdWithUnitDouble_3)
ArrayAdWithUnitDouble_3._units = new_instancemethod(_array_ad_with_unit.ArrayAdWithUnitDouble_3__units,None,ArrayAdWithUnitDouble_3)
ArrayAdWithUnitDouble_3._units_set = new_instancemethod(_array_ad_with_unit.ArrayAdWithUnitDouble_3__units_set,None,ArrayAdWithUnitDouble_3)
ArrayAdWithUnitDouble_3_swigregister = _array_ad_with_unit.ArrayAdWithUnitDouble_3_swigregister
ArrayAdWithUnitDouble_3_swigregister(ArrayAdWithUnitDouble_3)

class ArrayAdWithUnitDouble_4(object):
    """
    This is a ArrayAd that also has units associated with it.

    This is a simple structure that just keeps these two things together.

    C++ includes: array_ad_with_unit.h 
    """
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        """
        FullPhysics::ArrayAdWithUnit< T, D >::ArrayAdWithUnit(const ArrayAd< T, D > &V)

        """
        _array_ad_with_unit.ArrayAdWithUnitDouble_4_swiginit(self,_array_ad_with_unit.new_ArrayAdWithUnitDouble_4(*args))
    def convert(self, *args):
        """
        ArrayAdWithUnit<T, D> FullPhysics::ArrayAdWithUnit< T, D >::convert(const Unit &R) const
        Convert to the given units. 
        """
        return _array_ad_with_unit.ArrayAdWithUnitDouble_4_convert(self, *args)

    def _v_rows(self):
        """
        int FullPhysics::ArrayAdWithUnit< T, D >::rows() const

        """
        return _array_ad_with_unit.ArrayAdWithUnitDouble_4__v_rows(self)

    @property
    def rows(self):
        return self._v_rows()

    def _v_cols(self):
        """
        int FullPhysics::ArrayAdWithUnit< T, D >::cols() const

        """
        return _array_ad_with_unit.ArrayAdWithUnitDouble_4__v_cols(self)

    @property
    def cols(self):
        return self._v_cols()

    def _v_depth(self):
        """
        int FullPhysics::ArrayAdWithUnit< T, D >::depth() const

        """
        return _array_ad_with_unit.ArrayAdWithUnitDouble_4__v_depth(self)

    @property
    def depth(self):
        return self._v_depth()

    def _v_is_constant(self):
        """
        bool FullPhysics::ArrayAdWithUnit< T, D >::is_constant() const

        """
        return _array_ad_with_unit.ArrayAdWithUnitDouble_4__v_is_constant(self)

    @property
    def is_constant(self):
        return self._v_is_constant()

    def _v_number_variable(self):
        """
        int FullPhysics::ArrayAdWithUnit< T, D >::number_variable() const

        """
        return _array_ad_with_unit.ArrayAdWithUnitDouble_4__v_number_variable(self)

    @property
    def number_variable(self):
        return self._v_number_variable()

    def reference(self, *args):
        """
        void FullPhysics::ArrayAdWithUnit< T, D >::reference(const ArrayAdWithUnit< T, D > &V)

        """
        return _array_ad_with_unit.ArrayAdWithUnitDouble_4_reference(self, *args)

    @property
    def value(self):
      return self._value()

    @value.setter
    def value(self, val):
      self._value_set(val)

    @property
    def units(self):
      return self._units()

    @units.setter
    def units(self,val):
        self._units_set(val)
      
    __swig_destroy__ = _array_ad_with_unit.delete_ArrayAdWithUnitDouble_4
ArrayAdWithUnitDouble_4.convert = new_instancemethod(_array_ad_with_unit.ArrayAdWithUnitDouble_4_convert,None,ArrayAdWithUnitDouble_4)
ArrayAdWithUnitDouble_4._v_rows = new_instancemethod(_array_ad_with_unit.ArrayAdWithUnitDouble_4__v_rows,None,ArrayAdWithUnitDouble_4)
ArrayAdWithUnitDouble_4._v_cols = new_instancemethod(_array_ad_with_unit.ArrayAdWithUnitDouble_4__v_cols,None,ArrayAdWithUnitDouble_4)
ArrayAdWithUnitDouble_4._v_depth = new_instancemethod(_array_ad_with_unit.ArrayAdWithUnitDouble_4__v_depth,None,ArrayAdWithUnitDouble_4)
ArrayAdWithUnitDouble_4._v_is_constant = new_instancemethod(_array_ad_with_unit.ArrayAdWithUnitDouble_4__v_is_constant,None,ArrayAdWithUnitDouble_4)
ArrayAdWithUnitDouble_4._v_number_variable = new_instancemethod(_array_ad_with_unit.ArrayAdWithUnitDouble_4__v_number_variable,None,ArrayAdWithUnitDouble_4)
ArrayAdWithUnitDouble_4.reference = new_instancemethod(_array_ad_with_unit.ArrayAdWithUnitDouble_4_reference,None,ArrayAdWithUnitDouble_4)
ArrayAdWithUnitDouble_4._value = new_instancemethod(_array_ad_with_unit.ArrayAdWithUnitDouble_4__value,None,ArrayAdWithUnitDouble_4)
ArrayAdWithUnitDouble_4._value_set = new_instancemethod(_array_ad_with_unit.ArrayAdWithUnitDouble_4__value_set,None,ArrayAdWithUnitDouble_4)
ArrayAdWithUnitDouble_4._units = new_instancemethod(_array_ad_with_unit.ArrayAdWithUnitDouble_4__units,None,ArrayAdWithUnitDouble_4)
ArrayAdWithUnitDouble_4._units_set = new_instancemethod(_array_ad_with_unit.ArrayAdWithUnitDouble_4__units_set,None,ArrayAdWithUnitDouble_4)
ArrayAdWithUnitDouble_4_swigregister = _array_ad_with_unit.ArrayAdWithUnitDouble_4_swigregister
ArrayAdWithUnitDouble_4_swigregister(ArrayAdWithUnitDouble_4)



