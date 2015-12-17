// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "empirical_orthogonal_function_output.h"
%}
%base_import(register_output_base)
%import "output.i"
%import "empirical_orthogonal_function.i"

%fp_shared_ptr(FullPhysics::EmpiricalOrthogonalFunctionOutput);

namespace FullPhysics {
class EmpiricalOrthogonalFunctionOutput : public RegisterOutputBase {
public:
  EmpiricalOrthogonalFunctionOutput(const boost::shared_ptr<EmpiricalOrthogonalFunction>& Eof,
			     const std::string& Hdf_band_name);
  virtual void register_output(const boost::shared_ptr<Output>& out) const;
  virtual void register_output_apriori(const boost::shared_ptr<Output>& out) const;
};
}


