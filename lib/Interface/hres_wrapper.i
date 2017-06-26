// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "hres_wrapper.h"
#include "sub_state_vector_array.h"
#include "pressure.h"
%}

%base_import(radiative_transfer_single_wn)
%fp_shared_ptr(FullPhysics::HresWrapper);

namespace FullPhysics {
class HresWrapper : public RadiativeTransferSingleWn {
public:
  HresWrapper(const boost::shared_ptr<RadiativeTransferSingleWn>& Rt);
  %python_attribute(number_stokes, virtual int)
  %python_attribute(number_stream, virtual int)
  virtual blitz::Array<double, 1> stokes_single_wn
  (double Wn, int Spec_index, const ArrayAd<double, 2>& Iv) const;
  virtual ArrayAd<double, 1> stokes_and_jacobian_single_wn
  (double Wn, int Spec_index, const ArrayAd<double, 2>& Iv) const;
  %python_attribute(rt, boost::shared_ptr<RadiativeTransfer>)
};
}
