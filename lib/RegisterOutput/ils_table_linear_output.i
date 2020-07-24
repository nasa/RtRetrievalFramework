// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "ils_table_linear_output.h"
%}
%base_import(register_output_base)
%import "output.i"
%import "ils_table_linear.i"
%fp_shared_ptr(FullPhysics::IlsTableLinearOutput);

namespace FullPhysics {
class IlsTableLinearOutput : public RegisterOutputBase {
public:
  IlsTableLinearOutput(const boost::shared_ptr<IlsTableLinear>& D,
			     const std::string& Hdf_band_name);
  virtual void register_output(const boost::shared_ptr<Output>& out) const;
  virtual void register_output_apriori(const boost::shared_ptr<Output>& out) const;
};
}


