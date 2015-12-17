// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "ecmwf.h"
%}
%base_import(generic_object)
%import "array_ad.i"

%fp_shared_ptr(FullPhysics::Ecmwf);
%nodefaultctor FullPhysics::Ecmwf;

namespace FullPhysics {
class Ecmwf : public GenericObject {
public:
  virtual ~Ecmwf();
  blitz::Array<double, 1> specific_humidity(const 
      blitz::Array<double, 1>& Pressure_level) const;
  ArrayAd<double, 1> specific_humidity(const 
       ArrayAd<double, 1>& Pressure_level) const;
  blitz::Array<double, 1> h2o_vmr(const
	  blitz::Array<double, 1>& Pressure_level) const;
  ArrayAd<double, 1> h2o_vmr(const
	  ArrayAd<double, 1>& Pressure_level) const;
  void temperature_grid(blitz::Array<double, 1>& OUTPUT,
			blitz::Array<double, 1>& OUTPUT) const;
  void specific_humidity_grid(blitz::Array<double, 1>& OUTPUT,
			     blitz::Array<double, 1>& OUTPUT) const;
  blitz::Array<double, 1> temperature(const 
      blitz::Array<double, 1>& Pressure_level) const;
  ArrayAd<double, 1> temperature(const 
      ArrayAd<double, 1>& Pressure_level) const;
  %python_attribute(surface_pressure, double);
  %python_attribute(windspeed, double);
  std::string print_to_string();
};
}

