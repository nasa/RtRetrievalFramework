# This file was automatically generated by SWIG (http://www.swig.org).
# Version 2.0.9
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.



from sys import version_info
if version_info >= (3,0,0):
    new_instancemethod = lambda func, inst, cls: _lidort_interface_masters.SWIG_PyInstanceMethod_New(func)
else:
    from new import instancemethod as new_instancemethod
if version_info >= (2,6,0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_lidort_interface_masters', [dirname(__file__)])
        except ImportError:
            import _lidort_interface_masters
            return _lidort_interface_masters
        if fp is not None:
            try:
                _mod = imp.load_module('_lidort_interface_masters', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _lidort_interface_masters = swig_import_helper()
    del swig_import_helper
else:
    import _lidort_interface_masters
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


SHARED_PTR_DISOWN = _lidort_interface_masters.SHARED_PTR_DISOWN
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

class Brdf_Linsup_Masters(object):
    """
    C++ includes: lidort_interface_masters.h

    """
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        """
        FullPhysics::Brdf_Linsup_Masters::Brdf_Linsup_Masters(const int &thread_in)

        """
        _lidort_interface_masters.Brdf_Linsup_Masters_swiginit(self,_lidort_interface_masters.new_Brdf_Linsup_Masters(*args))
    def _v_brdf_sup_in(self):
        """
        const Brdf_Sup_Inputs& FullPhysics::Brdf_Linsup_Masters::brdf_sup_in() const

        """
        return _lidort_interface_masters.Brdf_Linsup_Masters__v_brdf_sup_in(self)

    @property
    def brdf_sup_in(self):
        return self._v_brdf_sup_in()

    def _v_brdf_sup_inputstatus(self):
        """
        const Brdf_Input_Exception_Handling& FullPhysics::Brdf_Linsup_Masters::brdf_sup_inputstatus() const

        """
        return _lidort_interface_masters.Brdf_Linsup_Masters__v_brdf_sup_inputstatus(self)

    @property
    def brdf_sup_inputstatus(self):
        return self._v_brdf_sup_inputstatus()

    def _v_thread(self):
        """
        const int& FullPhysics::Brdf_Linsup_Masters::thread() const

        """
        return _lidort_interface_masters.Brdf_Linsup_Masters__v_thread(self)

    @property
    def thread(self):
        return self._v_thread()

    def _v_brdf_sup_out(self):
        """
        const Brdf_Sup_Outputs& FullPhysics::Brdf_Linsup_Masters::brdf_sup_out() const

        """
        return _lidort_interface_masters.Brdf_Linsup_Masters__v_brdf_sup_out(self)

    @property
    def brdf_sup_out(self):
        return self._v_brdf_sup_out()

    def _v_brdf_linsup_out(self):
        """
        const Brdf_Linsup_Outputs& FullPhysics::Brdf_Linsup_Masters::brdf_linsup_out() const

        """
        return _lidort_interface_masters.Brdf_Linsup_Masters__v_brdf_linsup_out(self)

    @property
    def brdf_linsup_out(self):
        return self._v_brdf_linsup_out()

    def read_config(self, *args):
        """
        void FullPhysics::Brdf_Linsup_Masters::read_config(const std::string &filnam_in)

        """
        return _lidort_interface_masters.Brdf_Linsup_Masters_read_config(self, *args)

    def run(self, *args):
        """
        void FullPhysics::Brdf_Linsup_Masters::run(const bool &do_debug_restoration_in, const int &nmoments_input_in)

        """
        return _lidort_interface_masters.Brdf_Linsup_Masters_run(self, *args)

    __swig_destroy__ = _lidort_interface_masters.delete_Brdf_Linsup_Masters
Brdf_Linsup_Masters._v_brdf_sup_in = new_instancemethod(_lidort_interface_masters.Brdf_Linsup_Masters__v_brdf_sup_in,None,Brdf_Linsup_Masters)
Brdf_Linsup_Masters._v_brdf_sup_inputstatus = new_instancemethod(_lidort_interface_masters.Brdf_Linsup_Masters__v_brdf_sup_inputstatus,None,Brdf_Linsup_Masters)
Brdf_Linsup_Masters._v_thread = new_instancemethod(_lidort_interface_masters.Brdf_Linsup_Masters__v_thread,None,Brdf_Linsup_Masters)
Brdf_Linsup_Masters._v_brdf_sup_out = new_instancemethod(_lidort_interface_masters.Brdf_Linsup_Masters__v_brdf_sup_out,None,Brdf_Linsup_Masters)
Brdf_Linsup_Masters._v_brdf_linsup_out = new_instancemethod(_lidort_interface_masters.Brdf_Linsup_Masters__v_brdf_linsup_out,None,Brdf_Linsup_Masters)
Brdf_Linsup_Masters.read_config = new_instancemethod(_lidort_interface_masters.Brdf_Linsup_Masters_read_config,None,Brdf_Linsup_Masters)
Brdf_Linsup_Masters.run = new_instancemethod(_lidort_interface_masters.Brdf_Linsup_Masters_run,None,Brdf_Linsup_Masters)
Brdf_Linsup_Masters_swigregister = _lidort_interface_masters.Brdf_Linsup_Masters_swigregister
Brdf_Linsup_Masters_swigregister(Brdf_Linsup_Masters)

class Brdf_Sup_Masters(object):
    """
    C++ includes: lidort_interface_masters.h

    """
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        """
        FullPhysics::Brdf_Sup_Masters::Brdf_Sup_Masters(const int &thread_in)

        """
        _lidort_interface_masters.Brdf_Sup_Masters_swiginit(self,_lidort_interface_masters.new_Brdf_Sup_Masters(*args))
    def _v_brdf_sup_in(self):
        """
        const Brdf_Sup_Inputs& FullPhysics::Brdf_Sup_Masters::brdf_sup_in() const

        """
        return _lidort_interface_masters.Brdf_Sup_Masters__v_brdf_sup_in(self)

    @property
    def brdf_sup_in(self):
        return self._v_brdf_sup_in()

    def _v_brdf_sup_inputstatus(self):
        """
        const Brdf_Input_Exception_Handling& FullPhysics::Brdf_Sup_Masters::brdf_sup_inputstatus() const

        """
        return _lidort_interface_masters.Brdf_Sup_Masters__v_brdf_sup_inputstatus(self)

    @property
    def brdf_sup_inputstatus(self):
        return self._v_brdf_sup_inputstatus()

    def _v_thread(self):
        """
        const int& FullPhysics::Brdf_Sup_Masters::thread() const

        """
        return _lidort_interface_masters.Brdf_Sup_Masters__v_thread(self)

    @property
    def thread(self):
        return self._v_thread()

    def _v_brdf_sup_out(self):
        """
        const Brdf_Sup_Outputs& FullPhysics::Brdf_Sup_Masters::brdf_sup_out() const

        """
        return _lidort_interface_masters.Brdf_Sup_Masters__v_brdf_sup_out(self)

    @property
    def brdf_sup_out(self):
        return self._v_brdf_sup_out()

    def read_config(self, *args):
        """
        void FullPhysics::Brdf_Sup_Masters::read_config(const std::string &filnam_in)

        """
        return _lidort_interface_masters.Brdf_Sup_Masters_read_config(self, *args)

    def run(self, *args):
        """
        void FullPhysics::Brdf_Sup_Masters::run(const bool &do_debug_restoration_in, const int &nmoments_input_in)

        """
        return _lidort_interface_masters.Brdf_Sup_Masters_run(self, *args)

    __swig_destroy__ = _lidort_interface_masters.delete_Brdf_Sup_Masters
Brdf_Sup_Masters._v_brdf_sup_in = new_instancemethod(_lidort_interface_masters.Brdf_Sup_Masters__v_brdf_sup_in,None,Brdf_Sup_Masters)
Brdf_Sup_Masters._v_brdf_sup_inputstatus = new_instancemethod(_lidort_interface_masters.Brdf_Sup_Masters__v_brdf_sup_inputstatus,None,Brdf_Sup_Masters)
Brdf_Sup_Masters._v_thread = new_instancemethod(_lidort_interface_masters.Brdf_Sup_Masters__v_thread,None,Brdf_Sup_Masters)
Brdf_Sup_Masters._v_brdf_sup_out = new_instancemethod(_lidort_interface_masters.Brdf_Sup_Masters__v_brdf_sup_out,None,Brdf_Sup_Masters)
Brdf_Sup_Masters.read_config = new_instancemethod(_lidort_interface_masters.Brdf_Sup_Masters_read_config,None,Brdf_Sup_Masters)
Brdf_Sup_Masters.run = new_instancemethod(_lidort_interface_masters.Brdf_Sup_Masters_run,None,Brdf_Sup_Masters)
Brdf_Sup_Masters_swigregister = _lidort_interface_masters.Brdf_Sup_Masters_swigregister
Brdf_Sup_Masters_swigregister(Brdf_Sup_Masters)

class Lidort_Lcs_Masters(object):
    """
    C++ includes: lidort_interface_masters.h

    """
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        """
        FullPhysics::Lidort_Lcs_Masters::Lidort_Lcs_Masters(const int &thread_in)

        """
        _lidort_interface_masters.Lidort_Lcs_Masters_swiginit(self,_lidort_interface_masters.new_Lidort_Lcs_Masters(*args))
    def _v_thread(self):
        """
        const int& FullPhysics::Lidort_Lcs_Masters::thread() const

        """
        return _lidort_interface_masters.Lidort_Lcs_Masters__v_thread(self)

    @property
    def thread(self):
        return self._v_thread()

    def _v_lidort_fixin(self):
        """
        void FullPhysics::Lidort_Lcs_Masters::lidort_fixin(Lidort_Fixed_Inputs &lidort_fixin_in)

        """
        return _lidort_interface_masters.Lidort_Lcs_Masters__v_lidort_fixin(self)

    @property
    def lidort_fixin(self):
        return self._v_lidort_fixin()

    def _v_lidort_modin(self):
        """
        void FullPhysics::Lidort_Lcs_Masters::lidort_modin(Lidort_Modified_Inputs &lidort_modin_in)

        """
        return _lidort_interface_masters.Lidort_Lcs_Masters__v_lidort_modin(self)

    @property
    def lidort_modin(self):
        return self._v_lidort_modin()

    def _v_lidort_sup(self):
        """
        void FullPhysics::Lidort_Lcs_Masters::lidort_sup(Lidort_Sup_Inout &lidort_sup_in)

        """
        return _lidort_interface_masters.Lidort_Lcs_Masters__v_lidort_sup(self)

    @property
    def lidort_sup(self):
        return self._v_lidort_sup()

    def _v_lidort_out(self):
        """
        void FullPhysics::Lidort_Lcs_Masters::lidort_out(Lidort_Outputs &lidort_out_in)

        """
        return _lidort_interface_masters.Lidort_Lcs_Masters__v_lidort_out(self)

    @property
    def lidort_out(self):
        return self._v_lidort_out()

    def _v_lidort_linfixin(self):
        """
        void FullPhysics::Lidort_Lcs_Masters::lidort_linfixin(Lidort_Fixed_Lininputs &lidort_linfixin_in)

        """
        return _lidort_interface_masters.Lidort_Lcs_Masters__v_lidort_linfixin(self)

    @property
    def lidort_linfixin(self):
        return self._v_lidort_linfixin()

    def _v_lidort_linmodin(self):
        """
        void FullPhysics::Lidort_Lcs_Masters::lidort_linmodin(Lidort_Modified_Lininputs &lidort_linmodin_in)

        """
        return _lidort_interface_masters.Lidort_Lcs_Masters__v_lidort_linmodin(self)

    @property
    def lidort_linmodin(self):
        return self._v_lidort_linmodin()

    def _v_lidort_linsup(self):
        """
        void FullPhysics::Lidort_Lcs_Masters::lidort_linsup(Lidort_Linsup_Inout &lidort_linsup_in)

        """
        return _lidort_interface_masters.Lidort_Lcs_Masters__v_lidort_linsup(self)

    @property
    def lidort_linsup(self):
        return self._v_lidort_linsup()

    def _v_lidort_linout(self):
        """
        void FullPhysics::Lidort_Lcs_Masters::lidort_linout(Lidort_Linoutputs &lidort_linout_in)

        """
        return _lidort_interface_masters.Lidort_Lcs_Masters__v_lidort_linout(self)

    @property
    def lidort_linout(self):
        return self._v_lidort_linout()

    def run(self):
        """
        void FullPhysics::Lidort_Lcs_Masters::run()

        """
        return _lidort_interface_masters.Lidort_Lcs_Masters_run(self)

    __swig_destroy__ = _lidort_interface_masters.delete_Lidort_Lcs_Masters
Lidort_Lcs_Masters._v_thread = new_instancemethod(_lidort_interface_masters.Lidort_Lcs_Masters__v_thread,None,Lidort_Lcs_Masters)
Lidort_Lcs_Masters._v_lidort_fixin = new_instancemethod(_lidort_interface_masters.Lidort_Lcs_Masters__v_lidort_fixin,None,Lidort_Lcs_Masters)
Lidort_Lcs_Masters._v_lidort_modin = new_instancemethod(_lidort_interface_masters.Lidort_Lcs_Masters__v_lidort_modin,None,Lidort_Lcs_Masters)
Lidort_Lcs_Masters._v_lidort_sup = new_instancemethod(_lidort_interface_masters.Lidort_Lcs_Masters__v_lidort_sup,None,Lidort_Lcs_Masters)
Lidort_Lcs_Masters._v_lidort_out = new_instancemethod(_lidort_interface_masters.Lidort_Lcs_Masters__v_lidort_out,None,Lidort_Lcs_Masters)
Lidort_Lcs_Masters._v_lidort_linfixin = new_instancemethod(_lidort_interface_masters.Lidort_Lcs_Masters__v_lidort_linfixin,None,Lidort_Lcs_Masters)
Lidort_Lcs_Masters._v_lidort_linmodin = new_instancemethod(_lidort_interface_masters.Lidort_Lcs_Masters__v_lidort_linmodin,None,Lidort_Lcs_Masters)
Lidort_Lcs_Masters._v_lidort_linsup = new_instancemethod(_lidort_interface_masters.Lidort_Lcs_Masters__v_lidort_linsup,None,Lidort_Lcs_Masters)
Lidort_Lcs_Masters._v_lidort_linout = new_instancemethod(_lidort_interface_masters.Lidort_Lcs_Masters__v_lidort_linout,None,Lidort_Lcs_Masters)
Lidort_Lcs_Masters.run = new_instancemethod(_lidort_interface_masters.Lidort_Lcs_Masters_run,None,Lidort_Lcs_Masters)
Lidort_Lcs_Masters_swigregister = _lidort_interface_masters.Lidort_Lcs_Masters_swigregister
Lidort_Lcs_Masters_swigregister(Lidort_Lcs_Masters)

class Lidort_Lps_Masters(object):
    """
    C++ includes: lidort_interface_masters.h

    """
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        """
        FullPhysics::Lidort_Lps_Masters::Lidort_Lps_Masters(const int &thread_in)

        """
        _lidort_interface_masters.Lidort_Lps_Masters_swiginit(self,_lidort_interface_masters.new_Lidort_Lps_Masters(*args))
    def _v_thread(self):
        """
        const int& FullPhysics::Lidort_Lps_Masters::thread() const

        """
        return _lidort_interface_masters.Lidort_Lps_Masters__v_thread(self)

    @property
    def thread(self):
        return self._v_thread()

    def _v_lidort_fixin(self):
        """
        void FullPhysics::Lidort_Lps_Masters::lidort_fixin(Lidort_Fixed_Inputs &lidort_fixin_in)

        """
        return _lidort_interface_masters.Lidort_Lps_Masters__v_lidort_fixin(self)

    @property
    def lidort_fixin(self):
        return self._v_lidort_fixin()

    def _v_lidort_modin(self):
        """
        void FullPhysics::Lidort_Lps_Masters::lidort_modin(Lidort_Modified_Inputs &lidort_modin_in)

        """
        return _lidort_interface_masters.Lidort_Lps_Masters__v_lidort_modin(self)

    @property
    def lidort_modin(self):
        return self._v_lidort_modin()

    def _v_lidort_sup(self):
        """
        void FullPhysics::Lidort_Lps_Masters::lidort_sup(Lidort_Sup_Inout &lidort_sup_in)

        """
        return _lidort_interface_masters.Lidort_Lps_Masters__v_lidort_sup(self)

    @property
    def lidort_sup(self):
        return self._v_lidort_sup()

    def _v_lidort_out(self):
        """
        void FullPhysics::Lidort_Lps_Masters::lidort_out(Lidort_Outputs &lidort_out_in)

        """
        return _lidort_interface_masters.Lidort_Lps_Masters__v_lidort_out(self)

    @property
    def lidort_out(self):
        return self._v_lidort_out()

    def _v_lidort_linfixin(self):
        """
        void FullPhysics::Lidort_Lps_Masters::lidort_linfixin(Lidort_Fixed_Lininputs &lidort_linfixin_in)

        """
        return _lidort_interface_masters.Lidort_Lps_Masters__v_lidort_linfixin(self)

    @property
    def lidort_linfixin(self):
        return self._v_lidort_linfixin()

    def _v_lidort_linmodin(self):
        """
        void FullPhysics::Lidort_Lps_Masters::lidort_linmodin(Lidort_Modified_Lininputs &lidort_linmodin_in)

        """
        return _lidort_interface_masters.Lidort_Lps_Masters__v_lidort_linmodin(self)

    @property
    def lidort_linmodin(self):
        return self._v_lidort_linmodin()

    def _v_lidort_linsup(self):
        """
        void FullPhysics::Lidort_Lps_Masters::lidort_linsup(Lidort_Linsup_Inout &lidort_linsup_in)

        """
        return _lidort_interface_masters.Lidort_Lps_Masters__v_lidort_linsup(self)

    @property
    def lidort_linsup(self):
        return self._v_lidort_linsup()

    def _v_lidort_linout(self):
        """
        void FullPhysics::Lidort_Lps_Masters::lidort_linout(Lidort_Linoutputs &lidort_linout_in)

        """
        return _lidort_interface_masters.Lidort_Lps_Masters__v_lidort_linout(self)

    @property
    def lidort_linout(self):
        return self._v_lidort_linout()

    def run(self):
        """
        void FullPhysics::Lidort_Lps_Masters::run()

        """
        return _lidort_interface_masters.Lidort_Lps_Masters_run(self)

    __swig_destroy__ = _lidort_interface_masters.delete_Lidort_Lps_Masters
Lidort_Lps_Masters._v_thread = new_instancemethod(_lidort_interface_masters.Lidort_Lps_Masters__v_thread,None,Lidort_Lps_Masters)
Lidort_Lps_Masters._v_lidort_fixin = new_instancemethod(_lidort_interface_masters.Lidort_Lps_Masters__v_lidort_fixin,None,Lidort_Lps_Masters)
Lidort_Lps_Masters._v_lidort_modin = new_instancemethod(_lidort_interface_masters.Lidort_Lps_Masters__v_lidort_modin,None,Lidort_Lps_Masters)
Lidort_Lps_Masters._v_lidort_sup = new_instancemethod(_lidort_interface_masters.Lidort_Lps_Masters__v_lidort_sup,None,Lidort_Lps_Masters)
Lidort_Lps_Masters._v_lidort_out = new_instancemethod(_lidort_interface_masters.Lidort_Lps_Masters__v_lidort_out,None,Lidort_Lps_Masters)
Lidort_Lps_Masters._v_lidort_linfixin = new_instancemethod(_lidort_interface_masters.Lidort_Lps_Masters__v_lidort_linfixin,None,Lidort_Lps_Masters)
Lidort_Lps_Masters._v_lidort_linmodin = new_instancemethod(_lidort_interface_masters.Lidort_Lps_Masters__v_lidort_linmodin,None,Lidort_Lps_Masters)
Lidort_Lps_Masters._v_lidort_linsup = new_instancemethod(_lidort_interface_masters.Lidort_Lps_Masters__v_lidort_linsup,None,Lidort_Lps_Masters)
Lidort_Lps_Masters._v_lidort_linout = new_instancemethod(_lidort_interface_masters.Lidort_Lps_Masters__v_lidort_linout,None,Lidort_Lps_Masters)
Lidort_Lps_Masters.run = new_instancemethod(_lidort_interface_masters.Lidort_Lps_Masters_run,None,Lidort_Lps_Masters)
Lidort_Lps_Masters_swigregister = _lidort_interface_masters.Lidort_Lps_Masters_swigregister
Lidort_Lps_Masters_swigregister(Lidort_Lps_Masters)

class Lidort_Inputs(object):
    """
    C++ includes: lidort_interface_masters.h

    """
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self): 
        """
        FullPhysics::Lidort_Inputs::Lidort_Inputs()

        """
        _lidort_interface_masters.Lidort_Inputs_swiginit(self,_lidort_interface_masters.new_Lidort_Inputs())
    def _v_lidort_fixin(self):
        """
        const Lidort_Fixed_Inputs& FullPhysics::Lidort_Inputs::lidort_fixin() const

        """
        return _lidort_interface_masters.Lidort_Inputs__v_lidort_fixin(self)

    @property
    def lidort_fixin(self):
        return self._v_lidort_fixin()

    def _v_lidort_modin(self):
        """
        const Lidort_Modified_Inputs& FullPhysics::Lidort_Inputs::lidort_modin() const

        """
        return _lidort_interface_masters.Lidort_Inputs__v_lidort_modin(self)

    @property
    def lidort_modin(self):
        return self._v_lidort_modin()

    def _v_lidort_inputstatus(self):
        """
        const Lidort_Input_Exception_Handling& FullPhysics::Lidort_Inputs::lidort_inputstatus() const

        """
        return _lidort_interface_masters.Lidort_Inputs__v_lidort_inputstatus(self)

    @property
    def lidort_inputstatus(self):
        return self._v_lidort_inputstatus()

    def read_config(self, *args):
        """
        void FullPhysics::Lidort_Inputs::read_config(const std::string &filnam_in)

        """
        return _lidort_interface_masters.Lidort_Inputs_read_config(self, *args)

    __swig_destroy__ = _lidort_interface_masters.delete_Lidort_Inputs
Lidort_Inputs._v_lidort_fixin = new_instancemethod(_lidort_interface_masters.Lidort_Inputs__v_lidort_fixin,None,Lidort_Inputs)
Lidort_Inputs._v_lidort_modin = new_instancemethod(_lidort_interface_masters.Lidort_Inputs__v_lidort_modin,None,Lidort_Inputs)
Lidort_Inputs._v_lidort_inputstatus = new_instancemethod(_lidort_interface_masters.Lidort_Inputs__v_lidort_inputstatus,None,Lidort_Inputs)
Lidort_Inputs.read_config = new_instancemethod(_lidort_interface_masters.Lidort_Inputs_read_config,None,Lidort_Inputs)
Lidort_Inputs_swigregister = _lidort_interface_masters.Lidort_Inputs_swigregister
Lidort_Inputs_swigregister(Lidort_Inputs)

class Lidort_Masters(object):
    """
    C++ includes: lidort_interface_masters.h

    """
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        """
        FullPhysics::Lidort_Masters::Lidort_Masters(const int &thread_in)

        """
        _lidort_interface_masters.Lidort_Masters_swiginit(self,_lidort_interface_masters.new_Lidort_Masters(*args))
    def _v_thread(self):
        """
        const int& FullPhysics::Lidort_Masters::thread() const

        """
        return _lidort_interface_masters.Lidort_Masters__v_thread(self)

    @property
    def thread(self):
        return self._v_thread()

    def _v_lidort_fixin(self):
        """
        void FullPhysics::Lidort_Masters::lidort_fixin(Lidort_Fixed_Inputs &lidort_fixin_in)

        """
        return _lidort_interface_masters.Lidort_Masters__v_lidort_fixin(self)

    @property
    def lidort_fixin(self):
        return self._v_lidort_fixin()

    def _v_lidort_modin(self):
        """
        void FullPhysics::Lidort_Masters::lidort_modin(Lidort_Modified_Inputs &lidort_modin_in)

        """
        return _lidort_interface_masters.Lidort_Masters__v_lidort_modin(self)

    @property
    def lidort_modin(self):
        return self._v_lidort_modin()

    def _v_lidort_sup(self):
        """
        void FullPhysics::Lidort_Masters::lidort_sup(Lidort_Sup_Inout &lidort_sup_in)

        """
        return _lidort_interface_masters.Lidort_Masters__v_lidort_sup(self)

    @property
    def lidort_sup(self):
        return self._v_lidort_sup()

    def _v_lidort_out(self):
        """
        void FullPhysics::Lidort_Masters::lidort_out(Lidort_Outputs &lidort_out_in)

        """
        return _lidort_interface_masters.Lidort_Masters__v_lidort_out(self)

    @property
    def lidort_out(self):
        return self._v_lidort_out()

    def run(self):
        """
        void FullPhysics::Lidort_Masters::run()

        """
        return _lidort_interface_masters.Lidort_Masters_run(self)

    __swig_destroy__ = _lidort_interface_masters.delete_Lidort_Masters
Lidort_Masters._v_thread = new_instancemethod(_lidort_interface_masters.Lidort_Masters__v_thread,None,Lidort_Masters)
Lidort_Masters._v_lidort_fixin = new_instancemethod(_lidort_interface_masters.Lidort_Masters__v_lidort_fixin,None,Lidort_Masters)
Lidort_Masters._v_lidort_modin = new_instancemethod(_lidort_interface_masters.Lidort_Masters__v_lidort_modin,None,Lidort_Masters)
Lidort_Masters._v_lidort_sup = new_instancemethod(_lidort_interface_masters.Lidort_Masters__v_lidort_sup,None,Lidort_Masters)
Lidort_Masters._v_lidort_out = new_instancemethod(_lidort_interface_masters.Lidort_Masters__v_lidort_out,None,Lidort_Masters)
Lidort_Masters.run = new_instancemethod(_lidort_interface_masters.Lidort_Masters_run,None,Lidort_Masters)
Lidort_Masters_swigregister = _lidort_interface_masters.Lidort_Masters_swigregister
Lidort_Masters_swigregister(Lidort_Masters)

class Lidort_Sup_Accessories(object):
    """
    C++ includes: lidort_interface_masters.h

    """
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        """
        FullPhysics::Lidort_Sup_Accessories::Lidort_Sup_Accessories(boost::shared_ptr< Brdf_Sup_Inputs > &brdf_sup_in_in,
        boost::shared_ptr< Lidort_Fixed_Inputs > &lidort_fixin_in,
        boost::shared_ptr< Lidort_Modified_Inputs > &lidort_modin_in)

        """
        _lidort_interface_masters.Lidort_Sup_Accessories_swiginit(self,_lidort_interface_masters.new_Lidort_Sup_Accessories(*args))
    def _v_brdf_sup_in(self):
        """
        const Brdf_Sup_Inputs& FullPhysics::Lidort_Sup_Accessories::brdf_sup_in() const

        """
        return _lidort_interface_masters.Lidort_Sup_Accessories__v_brdf_sup_in(self)

    @property
    def brdf_sup_in(self):
        return self._v_brdf_sup_in()

    def _v_lidort_fixin(self):
        """
        const Lidort_Fixed_Inputs& FullPhysics::Lidort_Sup_Accessories::lidort_fixin() const

        """
        return _lidort_interface_masters.Lidort_Sup_Accessories__v_lidort_fixin(self)

    @property
    def lidort_fixin(self):
        return self._v_lidort_fixin()

    def _v_lidort_modin(self):
        """
        const Lidort_Modified_Inputs& FullPhysics::Lidort_Sup_Accessories::lidort_modin() const

        """
        return _lidort_interface_masters.Lidort_Sup_Accessories__v_lidort_modin(self)

    @property
    def lidort_modin(self):
        return self._v_lidort_modin()

    def _v_lidort_brdfcheck_status(self):
        """
        const Lidort_Exception_Handling& FullPhysics::Lidort_Sup_Accessories::lidort_brdfcheck_status() const

        """
        return _lidort_interface_masters.Lidort_Sup_Accessories__v_lidort_brdfcheck_status(self)

    @property
    def lidort_brdfcheck_status(self):
        return self._v_lidort_brdfcheck_status()

    def brdf_input_checker(self):
        """
        void FullPhysics::Lidort_Sup_Accessories::brdf_input_checker()

        """
        return _lidort_interface_masters.Lidort_Sup_Accessories_brdf_input_checker(self)

    __swig_destroy__ = _lidort_interface_masters.delete_Lidort_Sup_Accessories
Lidort_Sup_Accessories._v_brdf_sup_in = new_instancemethod(_lidort_interface_masters.Lidort_Sup_Accessories__v_brdf_sup_in,None,Lidort_Sup_Accessories)
Lidort_Sup_Accessories._v_lidort_fixin = new_instancemethod(_lidort_interface_masters.Lidort_Sup_Accessories__v_lidort_fixin,None,Lidort_Sup_Accessories)
Lidort_Sup_Accessories._v_lidort_modin = new_instancemethod(_lidort_interface_masters.Lidort_Sup_Accessories__v_lidort_modin,None,Lidort_Sup_Accessories)
Lidort_Sup_Accessories._v_lidort_brdfcheck_status = new_instancemethod(_lidort_interface_masters.Lidort_Sup_Accessories__v_lidort_brdfcheck_status,None,Lidort_Sup_Accessories)
Lidort_Sup_Accessories.brdf_input_checker = new_instancemethod(_lidort_interface_masters.Lidort_Sup_Accessories_brdf_input_checker,None,Lidort_Sup_Accessories)
Lidort_Sup_Accessories_swigregister = _lidort_interface_masters.Lidort_Sup_Accessories_swigregister
Lidort_Sup_Accessories_swigregister(Lidort_Sup_Accessories)


