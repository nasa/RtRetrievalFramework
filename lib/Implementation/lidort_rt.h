#ifndef LIDORT_RT_H
#define LIDORT_RT_H

#include "lidort_driver.h"
#include "spurr_rt.h"

namespace FullPhysics {

/****************************************************************//**
  Uses the Spurr interfaces to construct a radiative transfer
  class connecting L2 FP and LIDORT 3.5
 *******************************************************************/
class LidortRt : public SpurrRt {
public:
  LidortRt(const boost::shared_ptr<RtAtmosphere>& Atm,
	   const boost::shared_ptr<StokesCoefficient>& Stokes_coef,
	   const blitz::Array<double, 1>& Sza, 
	   const blitz::Array<double, 1>& Zen, 
	   const blitz::Array<double, 1>& Azm,
	   bool Pure_nadir,
	   int Number_streams, 
	   int Number_moments, 
	   bool Do_multi_scatt_only
	   );

  /// Number of quadtature streams in the cosine half space
  virtual int number_stream() const { return rt_driver()->number_stream(); }
  
  /// Number of moments for scattering matrix.
  int number_moment() const { return rt_driver()->number_moment(); };
  
  /// Convenience routine to get brdf driver object
  const boost::shared_ptr<LidortBrdfDriver> brdf_driver() const { return rt_driver()->lidort_brdf_driver(); }

  /// Convenience routine to get rt  driver object
  const boost::shared_ptr<LidortRtDriver> rt_driver() const 
  { return boost::shared_ptr<LidortRtDriver>(boost::dynamic_pointer_cast<LidortRtDriver>(rt_driver_)); }
  
  virtual void print(std::ostream& Os, bool Short_form = false) const;
};
}
#endif
