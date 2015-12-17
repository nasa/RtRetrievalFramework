// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include <std_vector.i>
%include "common.i"

%{
#include "ils.h"
%}

%base_import(state_vector)
%import "array_ad.i"
%import "spectral_domain.i"
%import "double_with_unit.i"

%fp_shared_ptr(FullPhysics::Ils);

namespace FullPhysics {
class Ils : public StateVectorObserver {
public:
  virtual ~Ils();
  std::string print_to_string() const;
  virtual blitz::Array<double, 1> apply_ils
  (const blitz::Array<double, 1>& High_resolution_wave_number,
   const blitz::Array<double, 1>& High_resolution_radiance,
   const std::vector<int>& Pixel_list) const = 0;
  virtual ArrayAd<double, 1> apply_ils
  (const blitz::Array<double, 1>& High_resolution_wave_number,
   const ArrayAd<double, 1>& High_resolution_radiance,
   const std::vector<int>& Pixel_list) const = 0;
  virtual boost::shared_ptr<Ils> clone() const = 0;
  %python_attribute(band_name, virtual std::string);
  %python_attribute(hdf_band_name, virtual std::string);
  %python_attribute(pixel_grid, virtual SpectralDomain);
  %python_attribute_with_set(ils_half_width, DoubleWithUnit);
};
}
%template(vector_ils) std::vector<boost::shared_ptr<FullPhysics::Ils> >;
