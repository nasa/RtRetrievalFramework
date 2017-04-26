// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"
%include "swig_optional.i"

%{
#include "forward_model.h"
%}

%base_import(generic_object)
%import "state_vector.i"
%import "spectrum.i"

%fp_shared_ptr(FullPhysics::ForwardModel)

namespace FullPhysics {
class StateVector;

class ForwardModel : public GenericObject {
public:
  virtual ~ForwardModel();
  std::string print_to_string() const;
  %python_attribute_abstract(state_vector, boost::shared_ptr<StateVector>)
  %python_attribute(number_spectrometer, virtual int)
  virtual std::string hdf_band_name(int Spec_index) const;
  virtual SpectralDomain spectral_domain(int Spec_index) const;
  virtual Spectrum radiance(int Spec_index, bool Skip_jacobian = false) 
    const = 0;
  virtual Spectrum measured_radiance(int Spec_index) 
    const = 0;
  virtual void setup_grid();
  Spectrum radiance_all(bool Skip_jacobian = false) const;
  %python_attribute(measured_radiance_all, Spectrum);
  boost::optional<blitz::Range> pixel_range(int Spec_index) const;
  %python_attribute_with_set(input_file_description, std::string);
};
}
