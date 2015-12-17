// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "lidort_rt.h"
#include "sub_state_vector_array.h"
#include "pressure.h"
%}
%base_import(observer)
%base_import(rt_atmosphere)
%base_import(spurr_rt)
%import "lidort_driver.i"
%fp_shared_ptr(FullPhysics::LidortRt);

namespace FullPhysics {
class LidortRt : public SpurrRt  {
public:
  LidortRt(const boost::shared_ptr<RtAtmosphere>& Atm,
	   const blitz::Array<double, 2>& Stokes_coef,
	   const blitz::Array<double, 1>& Sza, 
	   const blitz::Array<double, 1>& Zen, 
	   const blitz::Array<double, 1>& Azm,
       bool Pure_nadir,
	   int Number_streams, 
	   int Number_moments, 
	   bool Do_multi_scatt_only
	   );

  %python_attribute(number_stream, virtual int)
  %python_attribute(number_moment, int)
  
  %python_attribute(brdf_driver, boost::shared_ptr<LidortBrdfDriver>)
  %python_attribute(rt_driver, boost::shared_ptr<LidortRtDriver>)
  
  virtual void print(std::ostream& Os, bool Short_form = false) const;
};
}
