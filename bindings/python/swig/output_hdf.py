# This file was automatically generated by SWIG (http://www.swig.org).
# Version 2.0.9
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.



from sys import version_info
if version_info >= (3,0,0):
    new_instancemethod = lambda func, inst, cls: _output_hdf.SWIG_PyInstanceMethod_New(func)
else:
    from new import instancemethod as new_instancemethod
if version_info >= (2,6,0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_output_hdf', [dirname(__file__)])
        except ImportError:
            import _output_hdf
            return _output_hdf
        if fp is not None:
            try:
                _mod = imp.load_module('_output_hdf', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _output_hdf = swig_import_helper()
    del swig_import_helper
else:
    import _output_hdf
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


SHARED_PTR_DISOWN = _output_hdf.SHARED_PTR_DISOWN
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

import full_physics_swig.output
import full_physics_swig.generic_object
class OutputHdf(full_physics_swig.output.Output):
    """
    This write the output of the Level 2 Full physics.

    This particular implementation writes out a HDF5 file.

    Note that we make a few assumptions to simplify the implementation. We
    could relax any of these constraints, we would just need to modify the
    implementation.

    We assume that there are only a few hardcoded Shapes and Dimensions.
    The current implementation just uses a fixed hardcoded set. We could
    create a more general, flexible (but also more complicated)
    implementation if needed in the future.

    We determine the shape metadata of any particular field by looking at
    the size. If we don't recognize the size, we simply silently leave off
    the shape metadata information. This allows new fields to be added
    without needing to necessarily update the shape information (e.g., a
    new diagnostic field). A consequence of this design decision is that
    actual mistakes (e.g., wrong number of pressures level reported) won't
    be caught, and new field may be missing shape information until we
    change this class. I think this is the right trade, but we may need to
    reevaluate this at some point in the future.

    We assume that the different dimensions that can make the same rank
    shape are different sizes (e.g., Retrieval_Level_Array and
    Retrieval_StateVectorElement_Array). It is hard to think of a case
    where this wouldn't be true, but this assumption is entirely to make
    the implementation easier. There is no intrinsic reason why we
    couldn't support these being the same. We check and trigger an error
    if these are the same, but we can change the code to support this if
    needed.

    An alternative implementation would be to have a table mapping each
    field to the shape it is. This isn't hugely complicated, but would
    require us to keep a table of all possible fields. For now, we will do
    the simpler implementation.

    C++ includes: output_hdf.h 
    """
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        """
        OutputHdf::OutputHdf(const boost::shared_ptr< HdfFileGenerating > &H, int Num_level, int
        Statevector_size, int Num_aerosol, int Number_band)
        Constructor. This takes the file to write to. 
        """
        _output_hdf.OutputHdf_swiginit(self,_output_hdf.new_OutputHdf(*args))
    __swig_destroy__ = _output_hdf.delete_OutputHdf
OutputHdf_swigregister = _output_hdf.OutputHdf_swigregister
OutputHdf_swigregister(OutputHdf)



