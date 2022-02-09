// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "twostream_rt.h"
#include "sub_state_vector_array.h"
#include "pressure.h"
%}
%base_import(spurr_rt)
%import "twostream_driver.i"

%fp_shared_ptr(FullPhysics::TwostreamRt);

namespace FullPhysics {

// Force to be not abstract
%feature("notabstract") TwostreamRt;

class TwostreamRt : public SpurrRt {
public:
  TwostreamRt(const boost::shared_ptr<RtAtmosphere>& Atm,
              const boost::shared_ptr<StokesCoefficient>& Stokes_coef,
              const blitz::Array<double, 1>& Sza, 
              const blitz::Array<double, 1>& Zen, 
              const blitz::Array<double, 1>& Azm,
              bool do_fullquadrature = true);
  %python_attribute(number_stream, int)
  %python_attribute(number_moment, int)
  %python_attribute(brdf_driver, boost::shared_ptr<TwostreamBrdfDriver>)
  %python_attribute(rt_driver, boost::shared_ptr<TwostreamRtDriver>)
  virtual void print(std::ostream& Os, bool Short_form = false) const;
};
}
