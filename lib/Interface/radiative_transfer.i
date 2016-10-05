// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "radiative_transfer.h"
%}
%base_import(generic_object)
%fp_shared_ptr(FullPhysics::RadiativeTransfer);
%import "spectrum.i"
%import "spectral_domain.i"

namespace FullPhysics {
class RadiativeTransfer : public GenericObject {
public:
  virtual ~RadiativeTransfer();
  std::string print_to_string() const;
  %python_attribute_abstract(number_stokes, int)
  %python_attribute_abstract(number_spectrometer, int)
  virtual Spectrum reflectance
  (const SpectralDomain& Spec_domain, int Spec_index, 
   bool Skip_jacobian = false) const = 0;
  virtual blitz::Array<double, 2> stokes(const FullPhysics::SpectralDomain& Spec_domain,
					 int Spec_index) const = 0;
  virtual ArrayAd<double, 2> stokes_and_jacobian
  (const FullPhysics::SpectralDomain& Spec_domain, int Spec_index) const = 0;
};
}
