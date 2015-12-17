// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "max_a_posteriori_oco.h"
%}
%base_import(max_a_posteriori)
%base_import(model_measure_oco)
%import "forward_model.i"
%fp_shared_ptr(FullPhysics::MaxAPosterioriOCO);

namespace FullPhysics {
class MaxAPosterioriOCO : public MaxAPosteriori, public ModelMeasureOCO {
public:
  MaxAPosterioriOCO(const boost::shared_ptr<ForwardModel>& fm,
                    const blitz::Array<double, 1> a_priori_params,
                    const blitz::Array<double, 2> a_priori_cov);
  virtual ~MaxAPosterioriOCO();

};
}
