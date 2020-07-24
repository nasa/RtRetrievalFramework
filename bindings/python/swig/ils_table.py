# This file was automatically generated by SWIG (http://www.swig.org).
# Version 3.0.12
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.

from sys import version_info as _swig_python_version_info
if _swig_python_version_info >= (3, 0, 0):
    new_instancemethod = lambda func, inst, cls: _ils_table.SWIG_PyInstanceMethod_New(func)
else:
    from new import instancemethod as new_instancemethod
if _swig_python_version_info >= (2, 7, 0):
    def swig_import_helper():
        import importlib
        pkg = __name__.rpartition('.')[0]
        mname = '.'.join((pkg, '_ils_table')).lstrip('.')
        try:
            return importlib.import_module(mname)
        except ImportError:
            return importlib.import_module('_ils_table')
    _ils_table = swig_import_helper()
    del swig_import_helper
elif _swig_python_version_info >= (2, 6, 0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_ils_table', [dirname(__file__)])
        except ImportError:
            import _ils_table
            return _ils_table
        try:
            _mod = imp.load_module('_ils_table', fp, pathname, description)
        finally:
            if fp is not None:
                fp.close()
        return _mod
    _ils_table = swig_import_helper()
    del swig_import_helper
else:
    import _ils_table
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


SHARED_PTR_DISOWN = _ils_table.SHARED_PTR_DISOWN

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
import full_physics_swig.ils_function
class IlsTableLinear(full_physics_swig.ils_function.IlsFunction):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr

    def __init__(self, *args):
        _ils_table.IlsTableLinear_swiginit(self, _ils_table.new_IlsTableLinear(*args))
    __swig_destroy__ = _ils_table.delete_IlsTableLinear

    @property
    def wavenumber(self):
        return self._v_wavenumber()


    @property
    def delta_lambda(self):
        return self._v_delta_lambda()


    @property
    def response(self):
        return self._v_response()

IlsTableLinear._v_wavenumber = new_instancemethod(_ils_table.IlsTableLinear__v_wavenumber, None, IlsTableLinear)
IlsTableLinear._v_delta_lambda = new_instancemethod(_ils_table.IlsTableLinear__v_delta_lambda, None, IlsTableLinear)
IlsTableLinear._v_response = new_instancemethod(_ils_table.IlsTableLinear__v_response, None, IlsTableLinear)
IlsTableLinear.create_delta_lambda_to_response = new_instancemethod(_ils_table.IlsTableLinear_create_delta_lambda_to_response, None, IlsTableLinear)
IlsTableLinear.clone = new_instancemethod(_ils_table.IlsTableLinear_clone, None, IlsTableLinear)
IlsTableLinear_swigregister = _ils_table.IlsTableLinear_swigregister
IlsTableLinear_swigregister(IlsTableLinear)

class IlsTableLog(full_physics_swig.ils_function.IlsFunction):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr

    def __init__(self, *args):
        _ils_table.IlsTableLog_swiginit(self, _ils_table.new_IlsTableLog(*args))
    __swig_destroy__ = _ils_table.delete_IlsTableLog

    @property
    def wavenumber(self):
        return self._v_wavenumber()


    @property
    def delta_lambda(self):
        return self._v_delta_lambda()


    @property
    def response(self):
        return self._v_response()

IlsTableLog._v_wavenumber = new_instancemethod(_ils_table.IlsTableLog__v_wavenumber, None, IlsTableLog)
IlsTableLog._v_delta_lambda = new_instancemethod(_ils_table.IlsTableLog__v_delta_lambda, None, IlsTableLog)
IlsTableLog._v_response = new_instancemethod(_ils_table.IlsTableLog__v_response, None, IlsTableLog)
IlsTableLog.create_delta_lambda_to_response = new_instancemethod(_ils_table.IlsTableLog_create_delta_lambda_to_response, None, IlsTableLog)
IlsTableLog.clone = new_instancemethod(_ils_table.IlsTableLog_clone, None, IlsTableLog)
IlsTableLog_swigregister = _ils_table.IlsTableLog_swigregister
IlsTableLog_swigregister(IlsTableLog)



