// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include <std_vector.i>
%include "common.i"
%{
#include "aerosol_property.h"
%}
%base_import(generic_object)

%fp_shared_ptr(FullPhysics::AerosolProperty)
namespace FullPhysics {
class AerosolProperty : public GenericObject {
public:
  virtual ~AerosolProperty();
  virtual double extinction_coefficient(double wn) const = 0;
  virtual double scattering_coefficient(double wn) const = 0;
  virtual blitz::Array<double, 2> phase_function_moment(double wn, 
			int nmom = -1, int nscatt = -1) const = 0;
  std::string print_to_string() const;
};
}
%template(vector_aerosol_property) std::vector<boost::shared_ptr<FullPhysics::AerosolProperty> >;
