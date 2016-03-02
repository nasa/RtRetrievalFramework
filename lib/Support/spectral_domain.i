// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "spectral_domain.h"
%}
%base_import(generic_object)
%import "array_ad.i"
%import "unit.i"
%import "array_with_unit.i"

%fp_shared_ptr(FullPhysics::SpectralDomain)
namespace FullPhysics {
class SpectralDomain : public GenericObject {
public:
  enum TypePreference {PREFER_WAVENUMBER, PREFER_WAVELENGTH};
  SpectralDomain(const ArrayAd<double, 1>& Data, 
		 const Unit& U);
  SpectralDomain(const ArrayAd<double, 1>& Data, 
		 const std::string& U);
  SpectralDomain(const blitz::Array<double, 1>& FORCE_COPY,
		 const Unit& Units = units::inv_cm);
  SpectralDomain(const ArrayWithUnit<double, 1>& Data);
  SpectralDomain(const ArrayAd<double, 1>& Data, 
		 const blitz::Array<int, 1>& FORCE_COPY,
		 const Unit& U);
  SpectralDomain(const ArrayAd<double, 1>& Data, 
		 const blitz::Array<int, 1>& FORCE_COPY,
		 const std::string& U);
  SpectralDomain(const blitz::Array<double, 1>& FORCE_COPY,
		 const blitz::Array<int, 1>& FORCE_COPY,
		 const Unit& Units);
  SpectralDomain(const ArrayWithUnit<double, 1>& Data,
		 const blitz::Array<int, 1>& FORCE_COPY);
  %python_attribute(data, blitz::Array<double, 1>)
  %python_attribute(sample_index, blitz::Array<int, 1>)
  %python_attribute(units, Unit);
  %python_attribute(type_preference, TypePreference);
  blitz::Array<double, 1> convert_wave(const Unit& Units) const;
  blitz::Array<double, 1> convert_wave(const std::string& Units) const;
  blitz::Array<double, 1> wavenumber(const Unit& Units = units::inv_cm) const;
  blitz::Array<double, 1> wavenumber(const std::string&) const;
  blitz::Array<double, 1> wavelength(const Unit& Units = units::micron) const;
  blitz::Array<double, 1> wavelength(const std::string& Units) const;
  ArrayWithUnit<double, 1> photon_to_radiance_factor() const;
  std::string print_to_string() const;
};
}
