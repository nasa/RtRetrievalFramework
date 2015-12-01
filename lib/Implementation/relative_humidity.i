// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include <std_vector.i>
%include "common.i"
%{
#include "relative_humidity.h"
%}
%base_import(generic_object)
%import "absorber.i"
%import "temperature.i"
%import "pressure.i"
%fp_shared_ptr(FullPhysics::RelativeHumidity);
namespace FullPhysics {
class RelativeHumidity {
public:
  RelativeHumidity(const boost::shared_ptr<Absorber>& Abs, 
		   const boost::shared_ptr<Temperature>& Temp,
		   const boost::shared_ptr<Pressure>& Press);
  ArrayAd<double, 1> relative_humidity_grid() const;
  std::string print_to_string() const;
};
}
