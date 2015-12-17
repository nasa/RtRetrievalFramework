// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "forward_model_cost_function.h"
%}
%base_import(cost_function)
%import "forward_model.i"
%fp_shared_ptr(FullPhysics::ForwardModelCostFunction);

namespace FullPhysics {
class ForwardModelCostFunction : public CostFunction {
public:
  ForwardModelCostFunction(const boost::shared_ptr<ForwardModel>& fm);
  virtual void cost_function(const blitz::Array<double, 1>& X,
			     blitz::Array<double, 1>& OUTPUT,
			     blitz::Array<double, 1>& OUTPUT,
			     blitz::Array<double, 2>& OUTPUT) const;
};
}
