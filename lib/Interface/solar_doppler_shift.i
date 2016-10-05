// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "solar_doppler_shift.h"
%}

%base_import(generic_object)
%import "spectral_domain.i"
%import "double_with_unit.i"

%fp_shared_ptr(FullPhysics::SolarDopplerShift);

namespace FullPhysics {
class SolarDopplerShift  : public GenericObject {
public:
  virtual ~SolarDopplerShift();
  std::string print_to_string() const;
  %python_attribute_abstract(solar_distance, DoubleWithUnit);
  virtual SpectralDomain doppler_stretch(
     const SpectralDomain& Spec_domain) const = 0;
};
}
