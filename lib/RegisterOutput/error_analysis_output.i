// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "error_analysis_output.h"
#include "sub_state_vector_array.h"
%}
%base_import(register_output_base)
%import "output.i"
%import "error_analysis.i"

%fp_shared_ptr(FullPhysics::ErrorAnalysisOutput);

namespace FullPhysics {
class ErrorAnalysisOutput : public RegisterOutputBase {
public:
  ErrorAnalysisOutput(const boost::shared_ptr<ErrorAnalysis>& E, 
		      const blitz::Array<bool, 1>& Spec_flag,
		      bool Have_co2 = false);

  virtual void register_output(const boost::shared_ptr<Output>& out) const;
  virtual void register_output_apriori(const boost::shared_ptr<Output>& out) const;
};
}


