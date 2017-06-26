from builtins import object
# Test that we do casting correctly, where we take a function returning a
# boost::shared_ptr<Base> and return the exact derived type to python
# (e.g., a function returning a boost::shared_ptr<Pressure>, which is actually
# a PressureSigma, will return the PressureSigma type to python.

from nose.tools import *
from full_physics import *
from numpy.testing import *
from nose.plugins.skip import Skip, SkipTest

if(not have_full_physics_swig):
    class PressureImpBase(object):
        pass

# A simple python based class, to 
class PythonPressureSigma(PressureImpBase):
    '''This class uses pressure sigma levels to determine the pressure
    levels we do the RT on.'''
    def __init__(self, a, b, surface_pressure, pressure_flag = True):
        coef = np.array([surface_pressure])
        flag = np.array([pressure_flag])
        PressureImpBase.__init__(self,coef, flag)
        self.a = a
        self.b = b

    def clone(self):
        res = PressureSigma(self.a, self.b,
                            self.coefficient().value()[0],
                            self.used_flag_value()[0])
        return res
    
    def calc_pressure_grid(self):
        t = self.b * self.coefficient()[0] + self.a
        self.pgrid.reference(np_to_array_ad(t))

    def state_vector_name_i(self, i):
        return "Surface Pressure (Pascals)"

    def my_func(self):
        '''Function just to make sure we have an object this type.'''
        return 102

    def desc(self):
        return '''
Pressure Sigma:
  Surface pressure: %f
  Retrieval flag: %s
  a: %s
  b: %s
''' %(self.coefficient().value()[0], self.used_flag_value()[0].__str__(),
      self.a.__str__(), self.b.__str__())

def test_cast_cpp():
    '''Test returning a C++ class'''
    if(not have_full_physics_swig):
        raise SkipTest
    psigma = PressureSigma([0,0,0], [0.3, 0.6, 1.0], 10, True)
    pwrap = PressureHolder(psigma)
    # Test functions only in PressureSigma
    assert_almost_equal(pwrap.p.b, [0.3, 0.6, 1.0])
    pinp = PressureLevelInput([1, 2, 3])
    plevel = PressureFixedLevel(False, pinp, 2.5)
    pwrap.p = plevel
    # Function only in PressureFixedLevel
    assert pwrap.p.number_active_level == 3

@raises(AttributeError)
def test_cast_cpp_excep():
    '''Test returning a C++ class'''
    if(not have_full_physics_swig):
        raise SkipTest
    psigma = PressureSigma([0,0,0], [0.3, 0.6, 1.0], 10, True)
    pwrap = PressureHolder(psigma)
    pinp = PressureLevelInput([1, 2, 3])
    plevel = PressureFixedLevel(False, pinp, 2.5)
    pwrap.p = plevel
    # This will cause an exception
    pwrap.p.b

def test_cast_python():
    '''Make sure we handle classes that are actually python correctly.'''
    if(not have_full_physics_swig):
        raise SkipTest
    psigma = PythonPressureSigma([0,0,0], [0.3, 0.6, 1.0], 10, True)
    pwrap = PressureHolder(psigma)
    assert pwrap.p.my_func() == 102
