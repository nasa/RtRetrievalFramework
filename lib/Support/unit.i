// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "unit.h"
%}

%base_import(generic_object)
%fp_shared_ptr(FullPhysics::Unit)
namespace FullPhysics {
class Unit : public GenericObject {
public:
  std::string print_to_string() const;
  const int number_base_unit = 11;
  Unit(const std::string& Name, double Conversion_to_si, 
	      const boost::rational<int>& Length_power,
	      const boost::rational<int>& Mass_power,
	      const boost::rational<int>& Time_power,
	      const boost::rational<int>& Current_power,
	      const boost::rational<int>& Temperature_power,
	      const boost::rational<int>& Amount_power,
	      const boost::rational<int>& Luminous_intensity_power,
	      const boost::rational<int>& Solid_angle_power,
	      const boost::rational<int>& Angle_power,
              const boost::rational<int>& Photon_count_power,
              const boost::rational<int>& Sample_index);
  Unit(const std::string& Name, const Unit& Dunit);
  Unit(const std::string& Name_to_parse);
  %python_attribute(base_unit_powers, const boost::array<boost::rational<int>, number_base_unit>&);
  %python_attribute(conversion_to_si, double)
  %python_attribute(name, std::string)
  bool is_commensurate(const Unit& Units) const;
  Unit operator*=(const Unit& Dunit);
  Unit operator*=(double Scale_factor);
  Unit operator/=(const Unit& Dunit);
  %extend {
    Unit __mul__(const Unit& Y) 
    { return *$self * Y; }
    Unit __mul__(double Y) 
    { return *$self * Y; }
    Unit __rmul__(double X) 
    { return X * *$self;}
    // Python 2 division operator name 
    Unit __div__(const Unit& Y) 
    { return *$self / Y; }
    // Python 3 division operator name 
    Unit __truediv__(const Unit& Y) 
    { return *$self / Y; }
    Unit __rdiv__(double X) 
    { return X / *$self;}
    Unit __pow__(int X) 
    { return FullPhysics::pow(*$self,  X);}
  }

  // 1 here is the pickle format version, so we can tell if we try to
  // read data with a different format version than the code here.
  %pickle_init(1, self.name)
};

double conversion(const Unit& Dunit_from, const Unit& Dunit_to);
}

