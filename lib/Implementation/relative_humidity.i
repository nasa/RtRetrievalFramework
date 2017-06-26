// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include <std_vector.i>
%include "common.i"
%{
#include "relative_humidity.h"
#include "sub_state_vector_array.h"
#include "altitude.h"
%}
%base_import(generic_object)
%import "absorber.i"
%import "temperature.i"
%import "pressure.i"
%import "array_ad.i"
%fp_shared_ptr(FullPhysics::RelativeHumidity);
namespace FullPhysics {
class RelativeHumidity {
public:
  RelativeHumidity(const boost::shared_ptr<Absorber>& Abs, 
		   const boost::shared_ptr<Temperature>& Temp,
		   const boost::shared_ptr<Pressure>& Press);
  virtual boost::shared_ptr<RelativeHumidity> clone() const;
  virtual boost::shared_ptr<RelativeHumidity> 
  clone(const boost::shared_ptr<Absorber>& Abs, 
	const boost::shared_ptr<Temperature>& Temp,
	const boost::shared_ptr<Pressure>& Press) const;
  ArrayAd<double, 1> relative_humidity_grid() const;
  ArrayAd<double, 1> relative_humidity_layer() const;
  ArrayAd<double, 1> specific_humidity_grid() const;
  std::string print_to_string() const;
};
}
