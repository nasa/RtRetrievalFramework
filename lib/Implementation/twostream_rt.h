#ifndef TWOSTREAM_RT_H
#define TWOSTREAM_RT_H

#include "twostream_driver.h"
#include "spurr_rt.h"

namespace FullPhysics {

/****************************************************************//**
  Uses the Spurr interfaces to construct a radiative transfer
  class connecting L2 FP and TwoStream
 *******************************************************************/
class TwostreamRt : public SpurrRt {
public:
  TwostreamRt(const boost::shared_ptr<RtAtmosphere>& Atm,
              const boost::shared_ptr<StokesCoefficient>& Stokes_coef,
              const blitz::Array<double, 1>& Sza, 
              const blitz::Array<double, 1>& Zen, 
              const blitz::Array<double, 1>& Azm,
              bool do_fullquadrature = true);

  /// Number of quadtature streams in the cosine half space
  virtual int number_stream() const { return 1; }
  
  /// Number of moments for scattering matrix.
  /// 2stream natuarally uses up to 3 moments 
  int number_moment() const { return 3; };
  
  /// Convenience routine to get brdf driver object
  const boost::shared_ptr<TwostreamBrdfDriver> brdf_driver() const { return rt_driver()->twostream_brdf_driver(); }

  /// Convenience routine to get rt driver object
  const boost::shared_ptr<TwostreamRtDriver> rt_driver() const 
  { return boost::shared_ptr<TwostreamRtDriver>(boost::dynamic_pointer_cast<TwostreamRtDriver>(rt_driver_)); }

  virtual void print(std::ostream& Os, bool Short_form = false) const;
};
}
#endif
