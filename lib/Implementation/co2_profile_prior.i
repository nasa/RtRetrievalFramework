// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "co2_profile_prior.h"
%}
%base_import(generic_object)
%import "hdf_file.i"
%import "pressure.i"
%import "oco_met_file.i"
%fp_shared_ptr(FullPhysics::CO2ProfilePrior);

namespace FullPhysics {
class CO2ProfilePrior: public GenericObject {
public:
  CO2ProfilePrior(const OcoMetFile& Met_file,
		  const HdfFile& Profile_file) {}
  blitz::Array<double, 1> apriori_vmr(const Pressure& pressure) const;
};
}
